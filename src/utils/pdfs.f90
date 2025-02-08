module PDFS
    use TLab_Constants, only: wp, wi, big_wp
#ifdef USE_MPI
    use mpi_f08
#endif
    implicit none
    private

    public :: PDF1V2D, PDF1V2D1G, PDF2V2D, PDF_ANALIZE

    integer(wi) i, k, up, vp, ip, offset
    real(wp) umin, umax, ustep

#ifdef USE_MPI
    integer ims_err, impi
    real(wp) umin_p, umax_p
#endif

contains
    !########################################################################
    !#
    !# Calculate the PDF over plane of an array u using nbins bins.
    !#
    !# ilim     In    0, externally forced through umin_ext/umax_ext
    !#                otherwise, calculate locally the min/max
    !#
    !########################################################################
    subroutine PDF1V2D(ilim, nx, ny, nz, j, umin_ext, umax_ext, u, nbins, pdf, wrk1d, a, avg)
        implicit none

        integer(wi), intent(IN) :: ilim, nx, ny, nz, j, nbins
        real(wp), intent(IN) :: umin_ext, umax_ext
        real(wp), intent(IN) :: u(nx, ny, nz)
        real(wp), intent(OUT) :: pdf(nbins + 2)            ! Space at the end for min/max bins of u
        real(wp), intent(INOUT) :: wrk1d(nbins)
        real(wp), optional :: a(nx, ny, nz), avg(nbins) ! For conditional average, if needed

        ! ###################################################################
        pdf = 0.0_wp
        if (present(avg)) avg = 0.0_wp

        ! -------------------------------------------------------------------
        ! Calculate Minimum and Maximum
        ! -------------------------------------------------------------------
        if (ilim == 0) then
            umin = umin_ext
            umax = umax_ext

        else
            umin = u(1, j, 1)
            umax = u(1, j, 1)
            do k = 1, nz
                do i = 1, nx
                    umin = min(umin, u(i, j, k))
                    umax = max(umax, u(i, j, k))
                end do
            end do

#ifdef USE_MPI
            call MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
            umin = umin_p
            call MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
            umax = umax_p
#endif

        end if

        ustep = (umax - umin)/real(nbins, wp)    ! Calculate step in histogram
        pdf(nbins + 1) = umin + 0.5_wp*ustep    ! Calculate coordinate of histogram
        pdf(nbins + 2) = umax - 0.5_wp*ustep
        if (ustep == 0.0_wp) ustep = 1.0_wp   ! Just 1 point, prevent division by zero and force all in first bin

        ! -------------------------------------------------------------------
        ! Calculate Histogram
        ! -------------------------------------------------------------------
        do k = 1, nz
            do i = 1, nx
                up = int((u(i, j, k) - umin)/ustep) + 1
                if (ilim == 0) then
                    if (up <= nbins .and. up >= 1) then
                        pdf(up) = pdf(up) + 1.0_wp
                        if (present(a) .and. present(avg)) avg(up) = avg(up) + a(i, j, k)
                    end if
                else ! put last point in the last bin
                    up = min(up, nbins)
                    pdf(up) = pdf(up) + 1.0_wp
                    if (present(a) .and. present(avg)) avg(up) = avg(up) + a(i, j, k)
                end if
            end do
        end do

#ifdef USE_MPI
        impi = nbins
        call MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        pdf(1:nbins) = wrk1d(1:nbins)
        if (present(avg)) then
            call MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
            avg(1:nbins) = wrk1d(1:nbins)
        end if
#endif

        if (present(avg)) then              ! Save avg data in pdf array
            do up = 1, nbins
                if (pdf(up) > 0.0_wp) then       ! Avg remains zero if there is no point in this interval
                    avg(up) = avg(up)/pdf(up)
                end if
            end do
        end if

        return
    end subroutine PDF1V2D

    !########################################################################
    !#
    !# Conditioned on the intermittency field gate. Same as before, but with a
    !# conditional inside the loop.
    !#
    !# igate  In   Level of the gate signal to use as intermittency function
    !#
    !########################################################################
    subroutine PDF1V2D1G(ilim, nx, ny, nz, j, igate, gate, umin_ext, umax_ext, u, nbins, pdf, wrk1d, a, avg)
        implicit none

        integer(wi), intent(IN) :: ilim, nx, ny, nz, j, nbins
        real(wp), intent(IN) :: umin_ext, umax_ext
        real(wp), intent(IN) :: u(nx, ny, nz)
        real(wp), intent(OUT) :: pdf(nbins + 2)            ! Space at the end for min/max bins of u
        real(wp), intent(INOUT) :: wrk1d(nbins)
        integer(1), intent(IN) :: gate(nx, ny, nz), igate
        real(wp), optional :: a(nx, ny, nz), avg(nbins) ! For conditional average, if needed

        ! ###################################################################
        pdf = 0.0_wp
        if (present(avg)) avg = 0.0_wp

        ! -------------------------------------------------------------------
        ! Calculate Minimum and Maximum
        ! -------------------------------------------------------------------
        if (ilim == 0) then
            umin = umin_ext
            umax = umax_ext

        else
            umin = big_wp
            umax = -big_wp
            do k = 1, nz
                do i = 1, nx
                    if (gate(i, j, k) == igate) then
                        umin = min(umin, u(i, j, k))
                        umax = max(umax, u(i, j, k))
                    end if
                end do
            end do

#ifdef USE_MPI
            call MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
            umin = umin_p
            call MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
            umax = umax_p
#endif

        end if

        ustep = (umax - umin)/real(nbins, wp)    ! Calculate step in histogram
        pdf(nbins + 1) = umin + 0.5_wp*ustep    ! Calculate coordinate of histogram
        pdf(nbins + 2) = umax - 0.5_wp*ustep
        if (ustep == 0.0_wp) ustep = 1.0_wp   ! Just 1 point, prevent division by zero and force all in first bin

        ! -------------------------------------------------------------------
        ! Calculate Histogram
        ! -------------------------------------------------------------------
        do k = 1, nz
            do i = 1, nx
                if (gate(i, j, k) == igate) then
                    up = int((u(i, j, k) - umin)/ustep) + 1
                    if (ilim == 0) then
                        if (up <= nbins .and. up >= 1) then
                            pdf(up) = pdf(up) + 1.0_wp
                            if (present(a) .and. present(avg)) avg(up) = avg(up) + a(i, j, k)
                        end if
                    else ! put last point in the last bin
                        up = min(up, nbins)
                        pdf(up) = pdf(up) + 1.0_wp
                        if (present(a) .and. present(avg)) avg(up) = avg(up) + a(i, j, k)
                    end if
                end if

            end do
        end do

#ifdef USE_MPI
        impi = nbins
        call MPI_ALLREDUCE(pdf, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        pdf(1:nbins) = wrk1d(1:nbins)
        if (present(avg)) then
            call MPI_ALLREDUCE(avg, wrk1d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
            avg(1:nbins) = wrk1d(1:nbins)
        end if
#endif

        if (present(avg)) then              ! Save avg data in pdf array
            do up = 1, nbins
                if (pdf(up) > 0.0_wp) then       ! Avg remains zero if there is no point in this interval
                    avg(up) = avg(up)/pdf(up)
                end if
            end do
        end if

        return
    end subroutine PDF1V2D1G

    !########################################################################
    ! Joint PDFs
    !########################################################################
    subroutine PDF2V2D(nx, ny, nz, j, u, v, nbins, pdf, wrk2d, a, avg)
        implicit none

        integer(wi), intent(IN) :: nx, ny, nz, j, nbins(2)
        real(wp), intent(IN) :: u(nx, ny, nz), v(nx, ny, nz)
        real(wp), intent(OUT) :: pdf(nbins(1)*nbins(2) + 2 + 2*nbins(1)) ! Space at the end for min/max bins of u,v
        real(wp), intent(INOUT) :: wrk2d(nbins(1), nbins(2))              ! nbins(2) should be greater than 2 for enough memory space
        real(wp), optional :: a(nx, ny, nz), avg(nbins(1)*nbins(2))   ! For conditional average, if needed

        ! ###################################################################
        pdf = 0.0_wp
        if (present(avg)) avg = 0.0_wp

        offset = nbins(1)*nbins(2) + 2

        ! -------------------------------------------------------------------
        ! Calculate Minimum and Maximum
        ! -------------------------------------------------------------------
        ! First variable
        umin = u(1, j, 1)
        umax = u(1, j, 1)
        do k = 1, nz
            do i = 1, nx
                umin = min(umin, u(i, j, k))
                umax = max(umax, u(i, j, k))
            end do
        end do

#ifdef USE_MPI
        call MPI_ALLREDUCE(umin, umin_p, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
        umin = umin_p
        call MPI_ALLREDUCE(umax, umax_p, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        umax = umax_p
#endif

        ustep = (umax - umin)/real(nbins(1), wp)         ! Calculate step in histogram
        pdf(nbins(1)*nbins(2) + 1) = umin + 0.5_wp*ustep ! Calculate coordinate of histogram
        pdf(nbins(1)*nbins(2) + 2) = umax - 0.5_wp*ustep
        if (ustep == 0.0_wp) ustep = 1.0_wp            ! Just 1 point, prevent division by zero and force all in first bin

        ! Second variable
#define vmin(j)   wrk2d(j,1)
#define vmax(j)   wrk2d(j,2)
#define vstep(j)  wrk2d(j,3)
        vmin(1:nbins(1)) = big_wp  ! To calculate min
        vmax(1:nbins(1)) = -big_wp  ! To calculate max
        do k = 1, nz
            do i = 1, nx
                up = int((u(i, j, k) - umin)/ustep) + 1
                up = min(up, nbins(1))
                vmin(up) = min(vmin(up), v(i, j, k))
                vmax(up) = max(vmax(up), v(i, j, k))
            end do
        end do
#ifdef USE_MPI
        impi = nbins(1)
        call MPI_ALLREDUCE(vmin(1), vstep(1), impi, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
        vmin(1:nbins(1)) = vstep(1:nbins(1))
        call MPI_ALLREDUCE(vmax(1), vstep(1), impi, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        vmax(1:nbins(1)) = vstep(1:nbins(1))
#endif

        do up = 1, nbins(1)                                          ! Calculate step in histogram
            vstep(up) = (vmax(up) - vmin(up))/real(nbins(2), wp)
            ip = offset + up; pdf(ip) = vmin(up) + 0.5_wp*vstep(up)   ! Calculate coordinate of histogram
            ip = ip + nbins(1); pdf(ip) = vmax(up) - 0.5_wp*vstep(up)
            if (vstep(up) == 0.0_wp) vstep(up) = 1.0_wp               ! Just 1 point, prevent division by zero and force all in first bin
        end do

        ! -------------------------------------------------------------------
        ! Calculate Histogram
        ! -------------------------------------------------------------------
        do k = 1, nz
            do i = 1, nx
                up = int((u(i, j, k) - umin)/ustep) + 1
                up = min(up, nbins(1))
                vp = int((v(i, j, k) - vmin(up))/vstep(up)) + 1
                vp = min(vp, nbins(2))
                ip = (vp - 1)*nbins(1) + up
                pdf(ip) = pdf(ip) + 1.0_wp
                if (present(a) .and. present(avg)) avg(ip) = avg(ip) + a(i, j, k)
            end do
        end do

#undef vmin
#undef vmax
#undef vstep

#ifdef USE_MPI
        impi = nbins(1)*nbins(2)
        call MPI_ALLREDUCE(pdf, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        pdf(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2), 1)
        if (present(avg)) then
            call MPI_ALLREDUCE(avg, wrk2d, impi, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
            avg(1:nbins(1)*nbins(2)) = wrk2d(1:nbins(1)*nbins(2), 1)
        end if
#endif

        if (present(avg)) then              ! Save avg data in pdf array
            do ip = 1, nbins(1)*nbins(2)
                if (pdf(ip) > 0.0_wp) then       ! Avg remains zero if there is no point in this interval
                    avg(ip) = avg(ip)/pdf(ip)
                end if
            end do
        end if

        return
    end subroutine PDF2V2D

    !########################################################################
    !# Recalculating max/min for a given relative threshold plim.
    !# Adding nplim, the count of points above threshold.
    !# BCs flag ibc to drop extreme points.
    !########################################################################
    subroutine PDF_ANALIZE(ibc, nbins, pdf, umin_ext, umax_ext, plim, nplim)
        implicit none

        integer(wi), intent(IN) :: ibc, nbins
        real(wp), intent(IN) :: pdf(nbins + 2), plim
        real(wp), intent(INOUT) :: umin_ext, umax_ext
        integer(wi), intent(OUT) :: nplim

        ! -------------------------------------------------------------------
        integer(wi) upmin, upmax
        real(wp) pdf_threshold

        ! ###################################################################
        ustep = (pdf(nbins + 2) - pdf(nbins + 1))/real(nbins - 1, wp)

        upmin = 1                                       ! eliminate the outer bins according to BCs
        upmax = nbins
        if (ibc == 1 .or. ibc == 3) then
            upmin = upmin + 1
        end if
        if (ibc == 2 .or. ibc == 3) then
            upmax = upmax - 1
        end if
        umin_ext = pdf(nbins + 1) - 0.5_wp*ustep + ustep*real(upmin - 1, wp)
        umax_ext = pdf(nbins + 1) - 0.5_wp*ustep + ustep*real(upmax, wp)

        pdf_threshold = plim*maxval(pdf(upmin:upmax))  ! get absolute threshold

        nplim = 0                                       ! count the number of points above threshold
        do up = upmin, upmax
            if (pdf(up) > pdf_threshold) then
                nplim = nplim + 1
            end if
        end do

        do up = upmin, upmax                             ! elmininate smallest u-values if their probability is below threshold
            if (pdf(up) > pdf_threshold) then
                umin_ext = pdf(nbins + 1) - 0.5_wp*ustep + ustep*real(up - 1, wp)
                exit
            end if
        end do

        do up = upmax, upmin, -1                          ! elmininate largest u-values if their probability is below threshold
            if (pdf(up) > pdf_threshold) then
                umax_ext = pdf(nbins + 1) - 0.5_wp*ustep + ustep*real(up, wp)
                exit
            end if
        end do

        return
    end subroutine PDF_ANALIZE

end module PDFS
