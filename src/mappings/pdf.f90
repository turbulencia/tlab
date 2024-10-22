!########################################################################
!#
!# Calcualte {pdf(u_i), i=1,...,nv] in ny planes. (Histograms, not normalized.)
!# A last j-plane is added with the PDF in all the volume.
!#
!# ibc       In    BCs: 0 homogeneous interval
!#                      1 local interval
!#                      2 local interval, analysis and drop no point
!#                      3 local interval, analysis and drop left point
!#                      4 local interval, analysis and drop right point
!#                      5 local interval, analysis and drop both points
!#
!########################################################################
subroutine PDF1V_N(fname, time, nx, ny, nz, nv, nbins, ibc, umin, umax, u, igate, gate, y, pdf)
    use TLAB_CONSTANTS, only: lfile, wp, wi
    use TLAB_TYPES, only: pointers_dt
    use TLAB_ARRAYS, only: wrk1d
    use TLab_WorkFlow
    use PDFS
#ifdef USE_MPI
    use MPI
#endif

    implicit none

    character*(*), intent(IN) :: fname
    real(wp), intent(IN) :: time
    integer(wi), intent(IN) :: nx, ny, nz, nv, nbins, ibc(nv)
    real(wp), intent(IN) :: umin(nv), umax(nv)            ! Random variables
    type(pointers_dt), intent(IN) :: u(nv)
    integer(1), intent(IN) :: gate(*), igate               ! discrete conditioning criteria
    real(wp), intent(IN) :: y(ny)                        ! heights of each plane
    real(wp), intent(OUT) :: pdf(nbins + 2, ny + 1, nv)         ! last 2 bins contain the interval bounds

    ! -------------------------------------------------------------------
    integer(wi) iv, j, nplim, ibc_loc
    real(wp) plim, umin_loc, umax_loc

    character*64 name

#ifdef USE_MPI
    integer ims_pro, ims_err
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)
#endif

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    plim = 1.0e-4_wp                 ! relative threshold in PDF analysis; adapt to sample sizeo

    do iv = 1, nv

        do j = 1, ny                   ! calculation in planes
            if (igate == 0) then
                call PDF1V2D(ibc(iv), nx, ny, nz, j, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
            else
                call PDF1V2D1G(ibc(iv), nx, ny, nz, j, igate, gate, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
            end if

            if (ibc(iv) > 1) then
                ibc_loc = ibc(iv) - 2; umin_loc = umin(iv); umax_loc = umax(iv)
                call PDF_ANALIZE(ibc_loc, nbins, pdf(1, j, iv), umin_loc, umax_loc, plim, nplim)
                if (igate == 0) then
                    call PDF1V2D(0, nx, ny, nz, j, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
                else
                    call PDF1V2D1G(0, nx, ny, nz, j, igate, gate, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
                end if
            end if

        end do

        if (ny > 1) then            ! calculation in whole volume, saved as plane j=ny+1
            if (igate == 0) then
                call PDF1V2D(ibc(iv), nx*ny, 1, nz, 1, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
            else
                call PDF1V2D1G(ibc(iv), nx*ny, 1, nz, 1, igate, gate, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
            end if

            if (ibc(iv) > 1) then
                ibc_loc = ibc(iv) - 2; umin_loc = umin(iv); umax_loc = umax(iv)
                call PDF_ANALIZE(ibc_loc, nbins, pdf(1, j, iv), umin_loc, umax_loc, plim, nplim)
                if (igate == 0) then
                    call PDF1V2D(0, nx*ny, 1, nz, 1, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
                else
                    call PDF1V2D1G(0, nx*ny, 1, nz, 1, igate, gate, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, j, iv), wrk1d)
                end if
            end if

        end if

    end do

    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
        do iv = 1, nv
            name = trim(adjustl(fname))
            if (u(iv)%tag /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(u(iv)%tag))
            call TLAB_WRITE_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')
#include "dns_open_file.h"
            if (ny > 1) then
                write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:, :, iv))
            else
                write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:, 1, iv))
            end if
            close (LOC_UNIT_ID)
        end do

#ifdef USE_MPI
    end if
#endif

    return

end subroutine PDF1V_N

!########################################################################
!########################################################################
subroutine PDF2V(fname, time, nx, ny, nz, nbins, u, v, y, pdf)
    use TLAB_CONSTANTS, only: lfile, wp, wi
    use TLAB_ARRAYS, only: wrk2d
    use TLab_WorkFlow
    use PDFS
#ifdef USE_MPI
    use MPI
#endif

    implicit none

    character*(*), intent(IN) :: fname
    real(wp), intent(IN) :: time
    integer(wi), intent(IN) :: nx, ny, nz, nbins(2)
    real(wp), intent(IN) :: u(nx*ny*nz), v(nx*ny*nz)
    real(wp), intent(IN) :: y(ny)
    real(wp), intent(OUT) :: pdf(nbins(1)*nbins(2) + 2 + 2*nbins(1), ny + 1)

    ! -------------------------------------------------------------------
    integer(wi) j
    character*64 name

#ifdef USE_MPI
    integer ims_pro, ims_err
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)
#endif

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    do j = 1, ny               ! calculation in planes
        call PDF2V2D(nx, ny, nz, j, u, v, nbins, pdf(1, j), wrk2d)
    end do

    if (ny > 1) then        ! calculation in whole volume, saved as plane ny+1
        call PDF2V2D(nx*ny, 1, nz, 1, u, v, nbins, pdf(1, j), wrk2d)
    end if

    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
        name = trim(adjustl(fname))
        call TLAB_WRITE_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')
#include "dns_open_file.h"
        if (ny > 1) then
            write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:, :))
        else
            write (LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:, 1))
        end if
        close (LOC_UNIT_ID)
#ifdef USE_MPI
    end if
#endif

    return

end subroutine PDF2V
