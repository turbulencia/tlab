#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!########################################################################
subroutine INTEGRATE_SPECTRUM(nx, ny, nz, kr_total, isize_aux, &
                              spec_2d, data_x, data_z, spec_r, tmp_x, tmp_z, wrk2d)
    use TLab_Constants, only: wp, wi
    use FDM, only: g
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_err
    use TLabMPI_VARS, only: ims_npro_k
    use TLabMPI_VARS, only: ims_comm_x, ims_comm_z
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
    use TLabMPI_PROCS
#endif

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz, kr_total, isize_aux
    real(wp), dimension(nx, ny, nz), intent(IN) :: spec_2d ! power spectral density

    real(wp), dimension(kr_total, ny), intent(OUT) :: spec_r

    real(wp), dimension(nx, ny), intent(OUT) :: data_x
    real(wp), dimension(nx, ny, 2), intent(INOUT) :: tmp_x
    real(wp), dimension(nz/2, ny), intent(OUT) :: data_z
    real(wp), dimension(isize_aux, nz, 2), intent(INOUT) :: tmp_z ! need space for transpostion

#ifdef USE_MPI
    real(wp), dimension(isize_aux/ims_npro_k, g(3)%size, 2), intent(INOUT) :: wrk2d
#else
    real(wp), dimension(ny, g(3)%size, 2), intent(INOUT) :: wrk2d
#endif

! -----------------------------------------------------------------------
    integer(wi) :: i, k, kx_global, kz_global, kr_global, flag, ny_local
#ifdef USE_MPI
    integer(wi) count, id
#endif

! #######################################################################
    tmp_x = 0; tmp_z = 0; wrk2d = 0

    do k = 1, nz
#ifdef USE_MPI
        kz_global = k + ims_offset_k
#else
        kz_global = k
#endif

        tmp_x(:, :, 1) = tmp_x(:, :, 1) + spec_2d(:, :, k)

        if (kz_global <= g(3)%size/2) then; kz_global = g(3)%size/2 - kz_global + 2; flag = 1
        else; kz_global = kz_global - g(3)%size/2; flag = 0
        end if
! drop the Nyquist frequency; we add it to the previous mode to keep structure below
        if (kz_global == g(3)%size/2 + 1) kz_global = max(g(3)%size/2, 1)

        do i = 1, nx
#ifdef USE_MPI
            kx_global = i + ims_offset_i/2
#else
            kx_global = i
#endif

            tmp_z(1:ny, k, 1) = tmp_z(1:ny, k, 1) + spec_2d(i, 1:ny, k)

            kr_global = int(sqrt(real((kx_global - 1)**2 + (kz_global - 1)**2, wp))) + 1
            if (kr_global <= kr_total) &
                spec_r(kr_global, :) = spec_r(kr_global, :) + spec_2d(i, :, k)

! correction; half of the modes kx_global=1 have been counted twice
            if (kx_global == 1 .and. flag == 1) then
                tmp_z(1:ny, k, 1) = tmp_z(1:ny, k, 1) - spec_2d(i, 1:ny, k)

                if (kr_global <= kr_total) &
                    spec_r(kr_global, :) = spec_r(kr_global, :) - spec_2d(i, :, k)
            end if

! correction; half of the modes kz_global=1 need to be counted twice
            if (kz_global == 1 .and. kx_global > 1) then
                tmp_z(1:ny, k, 1) = tmp_z(1:ny, k, 1) + spec_2d(i, 1:ny, k)
            end if

        end do
    end do

! Finalize Ox spectrum
#ifdef USE_MPI
    count = nx*ny
    call MPI_ALLREDUCE(tmp_x(:, :, 1), tmp_x(:, :, 2), count, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
    data_x(:, :) = data_x(:, :) + tmp_x(:, :, 2)

#else
    data_x(:, :) = data_x(:, :) + tmp_x(:, :, 1)

#endif

! Finalize Oz spectrum
#ifdef USE_MPI
    count = isize_aux*nz
    call MPI_ALLREDUCE(tmp_z(:, :, 1), tmp_z(:, :, 2), count, MPI_REAL8, MPI_SUM, ims_comm_x, ims_err)

    if (ims_npro_k > 1) then
        id = TLabMPI_K_AUX2
        call TLabMPI_TRPF_K(tmp_z(:, :, 2), wrk2d(:, :, 1), id)

    else
        wrk2d(1:ny*nz, 1, 1) = tmp_z(1:ny*nz, 1, 2)

    end if

    ny_local = isize_aux/ims_npro_k

#else
    wrk2d(1:ny*nz, 1, 1) = tmp_z(1:ny*nz, 1, 1)

    ny_local = ny

#endif

    do k = 1, g(3)%size
        kz_global = k
        if (kz_global <= g(3)%size/2) then; kz_global = g(3)%size/2 - kz_global + 2
        else; kz_global = kz_global - g(3)%size/2
        end if
! drop the Nyquist frequency; we add it to the previous mode to keep structure below
        if (kz_global == g(3)%size/2 + 1) kz_global = max(g(3)%size/2, 1)

        wrk2d(1:ny_local, kz_global, 2) = wrk2d(1:ny_local, kz_global, 2) + wrk2d(1:ny_local, k, 1)

    end do

#ifdef USE_MPI
    if (ims_npro_k > 1) then
        count = g(3)%size/2/ims_npro_k ! add strides for the transposition
        do k = 1, g(3)%size/2, count
            wrk2d(1:ny_local*count, (k - 1)*2 + 1, 1) = wrk2d(1:ny_local*count, k, 2)
        end do

        call TLabMPI_TRPB_K(wrk2d(:, :, 1), tmp_z(:, :, 1), id)

    else
#endif

        tmp_z(1:isize_aux*nz, 1, 1) = wrk2d(1:ny_local*g(3)%size, 1, 2)

#ifdef USE_MPI
    end if
#endif

    do k = 1, nz/2
        data_z(k, 1:ny) = data_z(k, 1:ny) + tmp_z(1:ny, k, 1)
    end do

    return
end subroutine INTEGRATE_SPECTRUM

!########################################################################
!########################################################################
subroutine REDUCE_SPECTRUM(nx, ny, nz, nblock, in, out, tmp1, variance)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: isize_txc_dimz

! need to know about domain decomposition in x b/o
! nyquist frequency and zero frequency account different for the variance
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_pro_i, ims_npro_i, ims_err
#endif

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz, nblock
    complex(wp), dimension(isize_txc_dimz/2, nz), intent(IN) :: in, tmp1
    real(wp), dimension(nx/2, ny/nblock, 2*nz), intent(OUT) :: out ! Amplitude (1:nz) and Phase (nz+1:2*nz)
    real(wp), dimension(ny, 2) :: variance

! -----------------------------------------------------------------------
    integer(wi) :: kx, y, kz, ipy, ip, kx_global, ny_loc
    complex(wp) :: cdummy
    real(wp) :: power

! #######################################################################
! Calculate PSD and phase
! #######################################################################
    variance = 0.0_wp ! use variance to control result

! Drop the uppermost ny%nblock lines as there would not be
! nblock levels contributing to the output
    ny_loc = ny - mod(ny, nblock)

    do kz = 1, nz

        do y = 1, ny_loc
            ipy = (y - 1)/nblock + 1

! Drop the Nyquist frequency nx/2+1;
! keeping it would make writing inhomogeneous across processors
            do kx = 1, nx/2
#ifdef USE_MPI
                kx_global = kx + ims_offset_i
#else
                kx_global = kx
#endif
                ip = (nx/2 + 1)*(y - 1) + kx

                power = real(in(ip, kz))
                out(kx, ipy, kz) = out(kx, ipy, kz) + power

! phase; calculate phase only if there is some power
                if (power >= 1.e-16) then
                    cdummy = tmp1(ip, kz)
                    out(kx, ipy, kz + nz) = out(kx, ipy, kz + nz) + atan2(aimag(cdummy), real(cdummy))
                end if

! use variance to control result
                if (kx_global == 1) then; variance(y, 1) = variance(y, 1) + power
                else; variance(y, 1) = variance(y, 1) + power*2.0_wp
                end if

            end do

! add variance from nyquist frequency for the Parseval identity to hold
#ifdef USE_MPI
            if (ims_pro_i == ims_npro_i - 1) then
#endif
                ip = (nx/2 + 1)*(y - 1) + nx/2 + 1
                variance(y, 1) = variance(y, 1) + real(in(ip, kz))
#ifdef USE_MPI
            end if
#endif

        end do
    end do

! use variance to control result
#ifdef USE_MPI
    variance(:, 2) = variance(:, 1)
    call MPI_AllReduce(variance(1, 2), variance(1, 1), ny_loc, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif

    return

end subroutine REDUCE_SPECTRUM

!########################################################################
!########################################################################
subroutine REDUCE_CORRELATION(nx, ny, nz, nblock, nr_total, &
                              in, data_2d, data_x, data_z, data_r, variance1, variance2, icalc_radial)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: isize_wrk1d
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz, nblock, nr_total
    integer(wi), optional, intent(IN) :: icalc_radial ! whether to reduce radial correlations or not
    real(wp), dimension(nx, ny, nz), intent(IN) :: in
    real(wp), dimension(nx, ny/nblock, nz), intent(OUT) :: data_2d
    real(wp), dimension(nx, ny/nblock), intent(OUT) :: data_x
    real(wp), dimension(nz, ny/nblock), intent(OUT) :: data_z
    real(wp), dimension(nr_total, ny/nblock), intent(OUT) :: data_r
    real(wp), dimension(isize_wrk1d, 2), intent(IN) :: variance1 ! to normalize
    real(wp), dimension(isize_wrk1d), intent(OUT) :: variance2 ! to validate

! -----------------------------------------------------------------------
    integer(wi) i, j, k, ipy, ny_loc, i_global, k_global, r_global
    real(wp) norm

! #######################################################################
    variance2(:) = 0.0_wp ! use variance to control result

! Drop the uppermost ny%nblock lines as there would not be
! nblock levels contributing to the output
    ny_loc = ny - mod(ny, nblock)

    do j = 1, ny_loc
        ipy = (j - 1)/nblock + 1
        norm = sqrt(variance1(j, 1))*sqrt(variance1(j, 2))
        if (norm > 0.0_wp) then; norm = 1.0_wp/norm
        else; norm = 1.0_wp
        end if

        do k = 1, nz
#ifdef USE_MPI
            k_global = k + ims_offset_k
#else
            k_global = k
#endif

            do i = 1, nx
#ifdef USE_MPI
                i_global = i + ims_offset_i
#else
                i_global = i
#endif

! Reduce 2D data
                data_2d(i, ipy, k) = data_2d(i, ipy, k) + in(i, j, k)*norm

! Reduce 1D data
                if (k_global == 1) data_x(i, ipy) = data_x(i, ipy) + in(i, j, k)*norm

                if (i_global == 1) data_z(k, ipy) = data_z(k, ipy) + in(i, j, k)*norm

                if (icalc_radial == 1) then
                    r_global = int(sqrt(real((k_global - 1)**2 + (i_global - 1)**2, wp))) + 1
                    if (r_global <= nr_total) data_r(r_global, ipy) = data_r(r_global, ipy) + in(i, j, k)*norm

                end if
! use variance to control result
                if (i_global == 1 .and. k_global == 1) variance2(j) = in(i, j, k)

            end do
        end do
    end do

    return
end subroutine REDUCE_CORRELATION

!########################################################################
!########################################################################
subroutine RADIAL_SAMPLESIZE(nx, nz, nr_total, samplesize)
    use TLab_Constants, only: wp, wi

#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    integer(wi), intent(IN) :: nx, nz, nr_total
    real(wp), dimension(nr_total), intent(OUT) :: samplesize

! -----------------------------------------------------------------------
    integer(wi) i, k, i_global, k_global, r_global

! #######################################################################
    do k = 1, nz
#ifdef USE_MPI
        k_global = k + ims_offset_k
#else
        k_global = k
#endif

        do i = 1, nx
#ifdef USE_MPI
            i_global = i + ims_offset_i
#else
            i_global = i
#endif
            r_global = int(sqrt(real((k_global - 1)**2 + (i_global - 1)**2, wp))) + 1
            if (r_global <= nr_total) samplesize(r_global) = samplesize(r_global) + 1.0_wp

        end do
    end do

    return
end subroutine RADIAL_SAMPLESIZE

!########################################################################
!########################################################################
subroutine WRITE_SPECTRUM1D(fname, varname, nxy, nvar, pow)
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif

    implicit none

    character*(*), intent(IN) :: fname
    character*32, dimension(nvar), intent(IN) :: varname(nvar)
    integer(wi), intent(IN) :: nxy, nvar
    real(wp), dimension(nxy, nvar), intent(IN) :: pow

! -----------------------------------------------------------------------
    integer(wi) iv
    character*64 name

! #######################################################################
    if (nxy == 1) return

#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

        do iv = 1, nvar
            name = trim(adjustl(fname))
            if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))

            call TLab_Write_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')

#include "dns_open_file.h"
            write (LOC_UNIT_ID) SNGL(pow(1:nxy, iv))
            close (LOC_UNIT_ID)

        end do

#ifdef USE_MPI
    end if
#endif

    return

end subroutine WRITE_SPECTRUM1D

!########################################################################
!########################################################################
#ifdef USE_MPI

subroutine SPECTRA_MPIO_AUX(opt_main, nblock)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax
    use IO_FIELDS, only: io_aux
    use MPI
    use TLabMPI_VARS

    implicit none

    integer(wi), intent(IN) :: opt_main, nblock

! -----------------------------------------------------------------------
    integer(wi) :: ndims
    integer(wi), dimension(3) :: sizes, locsize, offset

    integer ims_color

! #######################################################################
    io_aux(:)%active = .false.
    io_aux(:)%offset = 0
    io_aux(:)%precision = IO_TYPE_SINGLE

! #######################################################################
    if (opt_main == 1 .or. opt_main == 2) then

! 1. Ox spectrum
        if (ims_pro_k == 0) io_aux(1)%active = .true.
        io_aux(1)%communicator = ims_comm_x

        ndims = 2 ! Subarray for the output of the Ox spectrum
        sizes(1) = imax*ims_npro_i/2; sizes(2) = jmax/nblock
        locsize(1) = imax/2; locsize(2) = jmax/nblock
        offset(1) = ims_offset_i/2; offset(2) = ims_offset_j/nblock

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(1)%subarray, ims_err)
        call MPI_Type_commit(io_aux(1)%subarray, ims_err)

!     subarray(1) = io_aux(1)%subarray ! to be removed

! 2. Oz spectrum
        if (ims_pro_i == 0) io_aux(2)%active = .true.
        io_aux(2)%communicator = ims_comm_z

        ndims = 2 ! Subarray for the output of the Oz spectrum
        sizes(1) = kmax*ims_npro_k/2; sizes(2) = jmax/nblock
        locsize(1) = kmax/2; locsize(2) = jmax/nblock
        offset(1) = ims_offset_k/2; offset(2) = ims_offset_j/nblock

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(2)%subarray, ims_err)
        call MPI_Type_commit(io_aux(2)%subarray, ims_err)

!     subarray(2) = io_aux(2)%subarray ! to be removed

! 3. Full 2D spectrum
        io_aux(3)%active = .true.
        io_aux(3)%communicator = MPI_COMM_WORLD

        ndims = 3 ! Subarray for the output of the 2D data
        sizes(1) = imax*ims_npro_i/2; sizes(2) = jmax/nblock; sizes(3) = kmax*ims_npro_k
        locsize(1) = imax/2; locsize(2) = jmax/nblock; locsize(3) = kmax
        offset(1) = ims_offset_i/2; offset(2) = ims_offset_j/nblock; offset(3) = ims_offset_k

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(3)%subarray, ims_err)
        call MPI_Type_commit(io_aux(3)%subarray, ims_err)

!     subarray(3) = io_aux(3)%subarray ! to be removed

! #######################################################################
    else if (opt_main == 3) then

! 1. Ox auto-correlation

! Create new communicator with 1/2 of the Ox domain
! Assuming even number of processes in each direction
        ims_color = 1
        if (ims_pro_i <= (ims_npro_i - 1)/2) &
            ims_color = 0
        call MPI_COMM_SPLIT(ims_comm_x, ims_color, ims_pro, ims_comm_x_aux, ims_err)

        if (ims_color == 0) then
            if (ims_pro_k == 0) io_aux(1)%active = .true.
            io_aux(1)%communicator = ims_comm_x_aux

            ndims = 2 ! Subarray for the output of the Ox cross-correlation
            sizes(1) = imax*ims_npro_i/2; sizes(2) = jmax/nblock
            locsize(1) = imax; locsize(2) = jmax/nblock
            offset(1) = ims_offset_i; offset(2) = ims_offset_j/nblock

            call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(1)%subarray, ims_err)
            call MPI_Type_commit(io_aux(1)%subarray, ims_err)

!        subarray(1) = io_aux(1)%subarray ! to be removed

        end if

! 2. Oz auto-correlation

! Create new communicator with 1/2 of the Oz domain
! Assuming even number of processes in each direction
        ims_color = 1
        if (ims_pro_k <= (ims_npro_k - 1)/2) &
            ims_color = 0
        call MPI_COMM_SPLIT(ims_comm_z, ims_color, ims_pro, ims_comm_z_aux, ims_err)

        if (ims_color == 0) then
            if (ims_pro_i == 0) io_aux(2)%active = .true.
            io_aux(2)%communicator = ims_comm_z_aux

            ndims = 2 ! Subarray for the output of the Oz cross-correlation
            sizes(1) = kmax*ims_npro_k/2; sizes(2) = jmax/nblock
            locsize(1) = kmax; locsize(2) = jmax/nblock
            offset(1) = ims_offset_k; offset(2) = ims_offset_j/nblock

            call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(2)%subarray, ims_err)
            call MPI_Type_commit(io_aux(2)%subarray, ims_err)

!        subarray(2) = io_aux(2)%subarray ! to be removed

        end if

! 2. Full 2D auto-correlation

! Create new communicator with the first 1/4 of the xOz domain
! Assuming even number of processes in each direction
        ims_color = 1
        if (ims_pro_i <= (ims_npro_i - 1)/2 .and. ims_pro_k <= (ims_npro_k - 1)/2) &
            ims_color = 0
        call MPI_COMM_SPLIT(MPI_COMM_WORLD, ims_color, ims_pro, ims_comm_xz_aux, ims_err)

        if (ims_color == 0) then
            io_aux(3)%active = .true.
            io_aux(3)%communicator = ims_comm_xz_aux

            ndims = 3 ! Subarray for the output of the 2D data
            sizes(1) = imax*ims_npro_i/2; sizes(2) = jmax/nblock; sizes(3) = kmax*ims_npro_k/2
            locsize(1) = imax; locsize(2) = jmax/nblock; locsize(3) = kmax
            offset(1) = ims_offset_i; offset(2) = ims_offset_j/nblock; offset(3) = ims_offset_k

            call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                          MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(3)%subarray, ims_err)
            call MPI_Type_commit(io_aux(3)%subarray, ims_err)

!        subarray(3) = io_aux(3)%subarray ! to be removed

        end if

! #######################################################################
    else if (opt_main == 4) then

! 1. Ox cross-correlation
        if (ims_pro_k == 0) io_aux(1)%active = .true.
        io_aux(1)%communicator = ims_comm_x

        ndims = 2 ! Subarray for the output of the Ox cross-correlation
        sizes(1) = imax*ims_npro_i; sizes(2) = jmax/nblock
        locsize(1) = imax; locsize(2) = jmax/nblock
        offset(1) = ims_offset_i; offset(2) = ims_offset_j/nblock

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(1)%subarray, ims_err)
        call MPI_Type_commit(io_aux(1)%subarray, ims_err)

!     subarray(1) = io_aux(1)%subarray ! to be removed

! 2. Oz cross-correlation
        if (ims_pro_i == 0) io_aux(2)%active = .true.
        io_aux(2)%communicator = ims_comm_z

        ndims = 2 ! Subarray for the output of the Oz cross-correlation
        sizes(1) = kmax*ims_npro_k; sizes(2) = jmax/nblock
        locsize(1) = kmax; locsize(2) = jmax/nblock
        offset(1) = ims_offset_k; offset(2) = ims_offset_j/nblock

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(2)%subarray, ims_err)
        call MPI_Type_commit(io_aux(2)%subarray, ims_err)

!     subarray(2) = io_aux(2)%subarray ! to be removed

! 3. Full 2D cross-correlation
        io_aux(3)%active = .true.
        io_aux(3)%communicator = MPI_COMM_WORLD

        ndims = 3 ! Subarray for the output of the 2D data
        sizes(1) = imax*ims_npro_i; sizes(2) = jmax/nblock; sizes(3) = kmax*ims_npro_k
        locsize(1) = imax; locsize(2) = jmax/nblock; locsize(3) = kmax
        offset(1) = ims_offset_i; offset(2) = ims_offset_j/nblock; offset(3) = ims_offset_k

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(3)%subarray, ims_err)
        call MPI_Type_commit(io_aux(3)%subarray, ims_err)

!     subarray(3) = io_aux(3)%subarray ! to be removed

    end if

    return
end subroutine SPECTRA_MPIO_AUX

#endif
