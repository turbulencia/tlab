#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

module OPR_FOURIER
    use TLAB_CONSTANTS, only: wp, wi, efile
    use TLAB_VARS, only: isize_txc_field, isize_txc_dimz, isize_wrk2d
    use TLAB_VARS, only: imax, jmax
    use TLAB_VARS, only: g
    use TLAB_VARS, only: ivfilter
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only: ims_npro_i, ims_npro_k
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k, ims_pro, ims_err
    use TLAB_MPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLAB_MPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLAB_MPI_PROCS
#endif
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

    implicit none
    private

    integer(8) :: fft_plan_fx, fft_plan_bx, fft_plan_fx_bcs
    integer(8) :: fft_plan_fy, fft_plan_by!, fft_plan_fy1d, fft_plan_by1d
    integer(8) :: fft_plan_fz, fft_plan_bz

    logical :: fft_reordering

    integer(wi) k
    complex(wp), pointer :: c_in(:,:) => null(), c_out(:,:) => null(), c_tmp1(:,:) => null(), c_tmp2(:,:) => null(), c_wrk3d(:,:) => null(), c_in1(:,:) => null()

    ! public :: fft_plan_fy1d, fft_plan_by1d ! vertical spectral pressure filter
    public :: OPR_FOURIER_INITIALIZE
    public :: OPR_FOURIER_F_X_EXEC, OPR_FOURIER_B_X_EXEC    ! Main routines, used in Poisson solver, frequently used and hence optimized
    public :: OPR_FOURIER_F_Z_EXEC, OPR_FOURIER_B_Z_EXEC
    public :: OPR_FOURIER_F                                 ! Minor routines, less frequently used
    public :: OPR_FOURIER_B
    public :: OPR_FOURIER_CONVOLUTION_FXZ
    public :: OPR_FOURIER_SPECTRA_3D

contains
! #######################################################################
! #######################################################################
    subroutine OPR_FOURIER_INITIALIZE()
        use TLAB_ARRAYS, only: wrk1d, wrk3d, txc
#ifdef USE_FFTW
#include "fftw3.f"
#endif

        ! -----------------------------------------------------------------------
        integer(wi) isize_stride, isize_disp, isize_fft_z, isize_fft_y, isize_fft_x

        ! #######################################################################
        ! FFTW library
        ! #######################################################################
#ifdef USE_FFTW

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            if (ims_size_i(TLAB_MPI_I_POISSON1) /= ims_size_i(TLAB_MPI_I_POISSON2)) then
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Error in the size in the transposition arrays.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if
        end if
#endif

        ! -----------------------------------------------------------------------
        ! Oz direction
        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            isize_fft_z = ims_size_k(TLAB_MPI_K_POISSON)/2 ! divide by 2 bcs. we work w complex #
        else
#endif
            isize_fft_z = (imax/2 + 1)*(jmax + 2)
#ifdef USE_MPI
        end if
#endif

        isize_stride = isize_fft_z

        if (g(3)%size > 1) then
#ifdef _DEBUG
            call dfftw_plan_many_dft(fft_plan_fz, 1, g(3)%size, isize_fft_z, &
                                     txc(:,1), g(3)%size, isize_stride, 1, &
                                     wrk3d, g(3)%size, isize_stride, 1, FFTW_FORWARD, FFTW_ESTIMATE)

            call dfftw_plan_many_dft(fft_plan_bz, 1, g(3)%size, isize_fft_z, &
                                     txc(:,1), g(3)%size, isize_stride, 1, &
                                     wrk3d, g(3)%size, isize_stride, 1, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
            call dfftw_plan_many_dft(fft_plan_fz, 1, g(3)%size, isize_fft_z, &
                                     txc(:,1), g(3)%size, isize_stride, 1, &
                                     wrk3d, g(3)%size, isize_stride, 1, FFTW_FORWARD, FFTW_MEASURE)

            call dfftw_plan_many_dft(fft_plan_bz, 1, g(3)%size, isize_fft_z, &
                                     txc(:,1), g(3)%size, isize_stride, 1, &
                                     wrk3d, g(3)%size, isize_stride, 1, FFTW_BACKWARD, FFTW_MEASURE)
#endif
        end if

        ! -----------------------------------------------------------------------
        ! Ox direction
        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            isize_fft_x = ims_size_i(TLAB_MPI_I_POISSON1)
        else
#endif
            isize_fft_x = jmax
#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            isize_disp = (imax/2 + 1)*ims_npro_i
        else
#endif
            isize_disp = g(1)%size/2 + 1
#ifdef _DEBUG
            ! wrk1d(:,2) should be complex with size n/2+1, i.e., n+2 real, but there is space at the end of wrk1d(:,2)
            call dfftw_plan_dft_r2c_1d(fft_plan_fx_bcs, g(1)%size, wrk1d(:,1), wrk1d(:,2), FFTW_ESTIMATE)
#else
            call dfftw_plan_dft_r2c_1d(fft_plan_fx_bcs, g(1)%size, wrk1d(:,1), wrk1d(:,2), FFTW_MEASURE)
#endif
#ifdef USE_MPI
        end if
#endif

#ifdef _DEBUG
        call dfftw_plan_many_dft_r2c(fft_plan_fx, 1, g(1)%size, isize_fft_x, &
                                     txc(:,1), g(1)%size, 1, g(1)%size, &
                                     wrk3d, g(1)%size/2 + 1, 1, isize_disp, FFTW_ESTIMATE)
        call dfftw_plan_many_dft_c2r(fft_plan_bx, 1, g(1)%size, isize_fft_x, &
                                     txc(:,1), g(1)%size/2 + 1, 1, isize_disp, &
                                     wrk3d, g(1)%size, 1, g(1)%size, FFTW_ESTIMATE)
#else
        call dfftw_plan_many_dft_r2c(fft_plan_fx, 1, g(1)%size, isize_fft_x, &
                                     txc(:,1), g(1)%size, 1, g(1)%size, &
                                     wrk3d, g(1)%size/2 + 1, 1, isize_disp, FFTW_MEASURE)
        call dfftw_plan_many_dft_c2r(fft_plan_bx, 1, g(1)%size, isize_fft_x, &
                                     txc(:,1), g(1)%size/2 + 1, 1, isize_disp, &
                                     wrk3d, g(1)%size, 1, g(1)%size, FFTW_MEASURE)
#endif

        ! -----------------------------------------------------------------------
        ! Oy direction
        ! -----------------------------------------------------------------------
        isize_fft_y = imax/2 + 1

        isize_stride = isize_fft_y

        if (g(2)%size > 1) then
#ifdef _DEBUG
            call dfftw_plan_many_dft(fft_plan_fy, 1, g(2)%size, isize_fft_y, &
                                     txc(:,1), g(2)%size, isize_stride, 1, &
                                     wrk3d, g(2)%size, isize_stride, 1, FFTW_FORWARD, FFTW_ESTIMATE)

            call dfftw_plan_many_dft(fft_plan_by, 1, g(2)%size, isize_fft_y, &
                                     txc(:,1), g(2)%size, isize_stride, 1, &
                                     wrk3d, g(2)%size, isize_stride, 1, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
            call dfftw_plan_many_dft(fft_plan_fy, 1, g(2)%size, isize_fft_y, &
                                     txc(:,1), g(2)%size, isize_stride, 1, &
                                     wrk3d, g(2)%size, isize_stride, 1, FFTW_FORWARD, FFTW_MEASURE)

            call dfftw_plan_many_dft(fft_plan_by, 1, g(2)%size, isize_fft_y, &
                                     txc(:,1), g(2)%size, isize_stride, 1, &
                                     wrk3d, g(2)%size, isize_stride, 1, FFTW_BACKWARD, FFTW_MEASURE)
#endif
        end if

        ! -----------------------------------------------------------------------
        ! Oy direction - spectral pressure filter (with pre. grid staggering)
        ! (only with quasi-periodic Bcs in the vertical, e.g. closed channel flow)
        ! -----------------------------------------------------------------------
        ! if (ivfilter == 1) then
        !     call OPR_STAGGERING_INITIALIZE(wrk3d)
        ! end if

#else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. FFTW needed for POISSON solver.')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)

#endif

        return
    end subroutine OPR_FOURIER_INITIALIZE

! #######################################################################
! #######################################################################
    subroutine OPR_FOURIER_F(flag_mode, nx, ny, nz, in, out, tmp1)
        use TLAB_ARRAYS, only: wrk2d
        use TLAB_POINTERS_C, only: c_wrk3d
        integer(wi) :: flag_mode                                ! 1D, 2D or 3D
        integer(wi) :: nx, ny, nz
        real(wp), intent(INOUT) :: in(isize_txc_field)          ! extended w/ BCs below
        real(wp), intent(OUT) :: out(isize_txc_dimz, nz)
        real(wp), intent(INOUT) :: tmp1(isize_txc_dimz, nz)

        target out, tmp1

        ! #######################################################################
        ! Pass memory address from real array to complex array
        call c_f_pointer(c_loc(out), c_out, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])

        wrk2d(:,1:2) = 0.0_wp          ! BCs

        if (flag_mode == 3 .and. g(2)%size > 1) then ! 3D FFT (unless 2D sim)
            if (g(3)%size > 1) then
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in, wrk2d(1, 1), wrk2d(1, 2), c_out, c_tmp1, c_wrk3d)
                call OPR_FOURIER_F_Z_EXEC(c_out, c_tmp1)
            else
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in, wrk2d(1, 1), wrk2d(1, 2), c_tmp1, c_out, c_wrk3d)
            end if

            do k = 1, nz
                call dfftw_execute_dft(fft_plan_fy, c_tmp1(:, k), c_out(:, k))
            end do

        else
            if (flag_mode == 2 .and. g(3)%size > 1) then ! 2D FFT (unless 1D sim)
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in, wrk2d(1, 1), wrk2d(1, 2), c_tmp1, c_out, c_wrk3d)
                call OPR_FOURIER_F_Z_EXEC(c_tmp1, c_out)
            else
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in, wrk2d(1, 1), wrk2d(1, 2), c_out, c_tmp1, c_wrk3d)
            end if
        end if

        return
    end subroutine OPR_FOURIER_F

! #######################################################################
! #######################################################################
    subroutine OPR_FOURIER_B(flag_mode, nx, ny, nz, in, out)
        use TLAB_POINTERS_C, only: c_wrk3d
        integer(wi) :: flag_mode  ! 1D, 2D or 3D
        integer(wi) :: nx, ny, nz
        real(wp), intent(INOUT) :: in(isize_txc_dimz, nz)
        real(wp), intent(OUT) :: out(isize_txc_field)

        target in

        ! #######################################################################
        ! Pass memory address from real array to complex array
        call c_f_pointer(c_loc(in), c_in, shape=[isize_txc_dimz/2, nz])

        if (flag_mode == 3 .and. g(2)%size > 1) then ! 3D FFT (unless 2D sim)
            do k = 1, nz
                call dfftw_execute_dft(fft_plan_by, c_in(:, k), c_wrk3d(:, k))
            end do

            if (g(3)%size > 1) then
                call OPR_FOURIER_B_Z_EXEC(c_wrk3d, c_in)
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_in, out, c_wrk3d)
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, out, c_in)
            end if

        else
            if (flag_mode == 2 .and. g(3)%size > 1) then ! 2D FFT (unless 1D sim)
                call OPR_FOURIER_B_Z_EXEC(c_in, c_wrk3d)
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, out, c_in)
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_in, out, c_wrk3d)
            end if

        end if

        return
    end subroutine OPR_FOURIER_B

! #######################################################################
! #######################################################################
! Calculate spectrum in array in1
! Calculate correlation in array in2
    subroutine OPR_FOURIER_CONVOLUTION_FXZ(flag1, flag2, nx, ny, nz, in1, in2, tmp1, tmp2, wrk2d, wrk3d)

        character(len=*), intent(in) :: flag1
        integer(wi), intent(in) :: flag2, nx, ny, nz
        real(wp), dimension(isize_txc_field), intent(INOUT) :: in1, in2, tmp1, tmp2, wrk3d
        real(wp), dimension(isize_wrk2d, 2), intent(INOUT) :: wrk2d ! BCs padding

        target in1, tmp1, tmp2, wrk3d

! #######################################################################
        ! Pass memory address from real array to complex array
        call c_f_pointer(c_loc(in1), c_in1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(wrk3d), c_wrk3d, shape=[isize_txc_dimz/2, nz])

        wrk2d = 0.0_wp
        fft_reordering = .true.

        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in1, wrk2d(1, 1), wrk2d(1, 2), c_wrk3d, c_tmp1, c_tmp2)
            call OPR_FOURIER_F_Z_EXEC(c_wrk3d, c_tmp1)
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in1, wrk2d(1, 1), wrk2d(1, 2), c_tmp1, c_tmp2, c_wrk3d)
        end if

        select case (trim(adjustl(flag1)))

        case ('auto')        ! Auto-spectra
            c_in1 = c_tmp1*conjg(c_tmp1)

        case ('cross')       ! Cross-spectra
            if (g(3)%size > 1) then
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in2, wrk2d(1, 1), wrk2d(1, 2), c_wrk3d, c_in1, c_tmp2)
                call OPR_FOURIER_F_Z_EXEC(c_wrk3d, c_in1)
            else
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, in2, wrk2d(1, 1), wrk2d(1, 2), c_in1, c_tmp2, c_wrk3d)
            end if
            c_in1 = c_in1*conjg(c_tmp1)

        end select

! -----------------------------------------------------------------------
        if (flag2 == 2) then         ! Calculate correlation in array in2
            tmp2 = in1                    ! the routines below can overwrite the entry array
            if (g(3)%size > 1) then
                call OPR_FOURIER_B_Z_EXEC(c_tmp2, c_wrk3d)
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, in2, c_tmp1)
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp2, in2, c_wrk3d)
            end if
        end if

! -----------------------------------------------------------------------
        fft_reordering = .false.

        return
    end subroutine OPR_FOURIER_CONVOLUTION_FXZ
!########################################################################
!# Forward/Backward Fourier transform of the array a extended by two planes in Oy direction (that is, jmax+2).
!#
!# In the case of the pressure solver these planes represent the boundary conditions.
!#
!# The transformed complex field is saved to an array with shape
!# (isize_txc_dimz/2,nz), that is, a stride isize_txc_dimz/2 between z-planes
!# (z-planes need not be contiguous).
!#
!# Each z-plane contains jmax+2 lines of nx/2+1 complex numbers. The first
!# nx/2 Fourier coefficients are contiguous and the element nx/2+1 is the
!# Nyquist frequency; if OX MPI decomposition, this element is just a padding used
!# to homogeneize arrays across PEs with ims_npro_i-PE, which contains then the
!# Nyquist frequency of the corresponding line.
!#
!########################################################################
    subroutine OPR_FOURIER_F_X_EXEC(nx, ny, nz, in, in_bcs_hb, in_bcs_ht, out, wrk1, wrk2)

        integer(wi) nx, ny, nz
        real(wp), dimension(nx*ny, *) :: in
        real(wp), dimension(nx, nz) :: in_bcs_hb, in_bcs_ht
        complex(wp), dimension(isize_txc_dimz/2, nz) :: out
        complex(wp), dimension(nx/2 + 1, *) :: wrk2
#ifdef USE_MPI
        complex(wp), dimension((nx/2 + 1)*ims_npro_i, *) :: wrk1
#else
        complex(wp), dimension(g(1)%size/2 + 1, *) :: wrk1
#endif

        target wrk1, wrk2

        ! -----------------------------------------------------------------------
        integer(wi) j, ip, isize_page, isize_line

#ifdef USE_MPI
        integer(wi) i, id, iold, inew
        real(wp), pointer :: r_wrk1(:) => null(), r_wrk2(:) => null()
#endif

        ! #######################################################################
        isize_line = nx/2 + 1; isize_page = isize_line*ny

        ! #######################################################################
#ifdef USE_MPI
        if (ims_npro_i > 1) then

            ! Pass memory address from complex array to real array
            call c_f_pointer(c_loc(wrk1), r_wrk1, shape=[isize_txc_field])
            call c_f_pointer(c_loc(wrk2), r_wrk2, shape=[isize_txc_field])

            ! Add bcs into array a; there must be space !
            ip = 1
            ip = ip + nx*ny*nz; in(ip:ip + nx*nz - 1, 1) = in_bcs_hb(1:nx*nz, 1)
            ip = ip + nx*nz; in(ip:ip + nx*nz - 1, 1) = in_bcs_ht(1:nx*nz, 1)

            ! Transpose array a into b
            id = TLAB_MPI_I_POISSON1
            call TLAB_MPI_TRPF_I(in, r_wrk2, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

            ! ims_size_k(id) FFTWs
            call dfftw_execute_dft_r2c(fft_plan_fx, r_wrk2, wrk1)

            ! reorganize wrk1 (FFTW make a stride in wrk1 already before)
            id = TLAB_MPI_I_POISSON1
            do k = 1, ims_size_i(id)
                inew = (nx/2 + 1)*ims_npro_i
                iold = g(1)%size/2 + 1
                wrk1(inew, k) = wrk1(iold, k)
                do ip = ims_npro_i, 2, -1
                    do i = nx/2, 1, -1
                        inew = (ip - 1)*(nx/2 + 1) + i
                        iold = (ip - 1)*nx/2 + i
                        wrk1(inew, k) = wrk1(iold, k)
                    end do
                end do
            end do

            ! Transpose array back
            id = TLAB_MPI_I_POISSON2
            call TLAB_MPI_TRPB_I(r_wrk1, r_wrk2, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

            ! reorganize wrk2 into b
            do k = 1, nz
                ip = 1 + ny*(k - 1); out(1:isize_page, k) = wrk2(1:isize_page, ip)
                ip = k + ny*nz; out(isize_page + 1:, k) = wrk2(1:isize_line, ip)
                ip = ip + nz; out(isize_page + isize_line + 1:, k) = wrk2(1:isize_line, ip)
            end do

            nullify (r_wrk1, r_wrk2)

        else
#endif

            ! #######################################################################
            do k = 1, nz
                call dfftw_execute_dft_r2c(fft_plan_fx, in(:, k), out(:, k))
                j = 1 + ny; ip = (j - 1)*isize_line + 1
                call dfftw_execute_dft_r2c(fft_plan_fx_bcs, in_bcs_hb(:, k), out(ip:, k))
                j = 2 + ny; ip = (j - 1)*isize_line + 1
                call dfftw_execute_dft_r2c(fft_plan_fx_bcs, in_bcs_ht(:, k), out(ip:, k))
            end do

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_FOURIER_F_X_EXEC

    !########################################################################
    !########################################################################
    subroutine OPR_FOURIER_B_X_EXEC(nx, ny, nz, in, out, wrk)
        integer(wi) nx, ny, nz
        complex(wp), dimension(isize_txc_dimz/2, nz), intent(IN) :: in
        real(wp), dimension(nx*ny, nz), intent(OUT) :: out
        complex(wp), dimension(nx/2 + 1, *), intent(inout) :: wrk

        target out, wrk

        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        integer(wi) i, ip, id, iold, inew, isize_page
        real(wp), pointer :: r_wrk(:) => null()
#endif

        !########################################################################
        ! Ox Parallel Decomposition
        !########################################################################
#ifdef USE_MPI
        if (ims_npro_i > 1) then

            ! Pass memory address from complex array to real array
            call c_f_pointer(c_loc(wrk), r_wrk, shape=[isize_txc_field])
            call c_f_pointer(c_loc(out), c_out, shape=[(nx/2 + 1)*ims_npro_i, nz])

            ! reorganize in into wrk
            isize_page = (nx/2 + 1)*ny
            do k = 1, nz
                ip = 1 + ny*(k - 1); wrk(1:isize_page, ip) = in(1:isize_page, k)
                !     ip = 1  +ny*(nx/2+1)* nz;   wrk(1: nx/2+1,    ip) = in(   ny   *(nx/2+1)+1:,k) !Idonotneed
                !     ip = ip +   (nx/2+1)* nz;   wrk(1: nx/2+1,    ip) = in(  (ny+1)*(nx/2+1)+1:,k) !Idonotneed
            end do

            ! Transpose array
            id = TLAB_MPI_I_POISSON2
            call TLAB_MPI_TRPF_I(r_wrk, out, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

            ! reorganize a (FFTW make a stride in a already before)
            id = TLAB_MPI_I_POISSON1
            do k = 1, ims_size_i(id)
                do ip = 2, ims_npro_i
                    do i = 1, nx/2
                        iold = (ip - 1)*(nx/2 + 1) + i
                        inew = (ip - 1)*nx/2 + i
                        c_out(inew, k) = c_out(iold, k)
                    end do
                end do
                iold = ims_npro_i*(nx/2 + 1)
                inew = g(1)%size/2 + 1
                c_out(inew, k) = c_out(iold, k)
            end do

            ! ims_size_i(id) FFTWs
            call dfftw_execute_dft_c2r(fft_plan_bx, c_out, r_wrk)

            ! Transpose array wrk into out
            id = TLAB_MPI_I_POISSON1
            call TLAB_MPI_TRPB_I(r_wrk, out, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

            nullify (r_wrk, c_out)

        else
#endif

            !########################################################################
            ! No Ox Parallel Decomposition
            !########################################################################
            do k = 1, nz
                call dfftw_execute_dft_c2r(fft_plan_bx, in(:,k), out(:, k))
            end do

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_FOURIER_B_X_EXEC

    !########################################################################
    !########################################################################
    subroutine OPR_FOURIER_F_Z_EXEC(in, out)

#ifdef USE_MPI
        complex(wp), dimension(ims_size_k(TLAB_MPI_K_POISSON)/2, g(3)%size), target :: in, out
#else
        complex(wp), dimension(isize_txc_dimz/2, g(3)%size), target :: in, out
#endif

        ! -----------------------------------------------------------------------
        complex(wp), dimension(:, :), pointer :: p_org, p_dst

        integer(wi) k_old1, k_old2, k_new1, k_new2
#ifdef USE_MPI
        integer(wi) id
        real(wp), pointer :: r_in(:) => null(), r_out(:) => null()
#endif

        ! #######################################################################
        ! Forward complex FFT in z
#ifdef USE_MPI
        id = TLAB_MPI_K_POISSON

        if (ims_npro_k > 1) then
            ! Pass memory address from complex array to real array
            call c_f_pointer(c_loc(in), r_in, shape=[isize_txc_field])
            call c_f_pointer(c_loc(out), r_out, shape=[isize_txc_field])

            call TLAB_MPI_TRPF_K(r_in, r_out, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            p_org => out
            p_dst => in
        else
#endif
            p_org => in
            p_dst => out
#ifdef USE_MPI
        end if
#endif

        call dfftw_execute_dft(fft_plan_fz, p_org, p_dst)

        if (fft_reordering) then ! re-shuffle spectra in z
            do k = 1, g(3)%size/2
                k_old1 = k + g(3)%size/2
                k_new1 = k
                k_old2 = k
                k_new2 = k + g(3)%size/2

                p_org(:, k_new1) = p_dst(:, k_old1)
                p_org(:, k_new2) = p_dst(:, k_old2)
            end do
            do k = 1, g(3)%size
                p_dst(:, k) = p_org(:, k)
            end do
        end if

#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPB_K(r_in, r_out, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            nullify (r_in, r_out)
        end if
#endif

        nullify (p_org, p_dst)

        return
    end subroutine OPR_FOURIER_F_Z_EXEC

    !########################################################################
    !########################################################################
    subroutine OPR_FOURIER_B_Z_EXEC(in, out)

#ifdef USE_MPI
        complex(wp), dimension(ims_size_k(TLAB_MPI_K_POISSON)/2, g(3)%size), target :: in, out
#else
        complex(wp), dimension(isize_txc_dimz/2, g(3)%size), target :: in, out
#endif

        ! -----------------------------------------------------------------------
        complex(wp), dimension(:, :), pointer :: p_org, p_dst

        integer(wi) k_old1, k_old2, k_new1, k_new2
#ifdef USE_MPI
        integer(wi) id
        real(wp), pointer :: r_in(:) => null(), r_out(:) => null()
#endif

        ! #######################################################################
        ! Forward complex FFT in z
#ifdef USE_MPI
        id = TLAB_MPI_K_POISSON

        if (ims_npro_k > 1) then
            ! Pass memory address from complex array to real array
            call c_f_pointer(c_loc(in), r_in, shape=[isize_txc_field])
            call c_f_pointer(c_loc(out), r_out, shape=[isize_txc_field])

            call TLAB_MPI_TRPF_K(r_in, r_out, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            p_org => out
            p_dst => in
        else
#endif
            p_org => in
            p_dst => out
#ifdef USE_MPI
        end if
#endif

        if (fft_reordering) then ! re-shuffle spectra in z
            do k = 1, g(3)%size/2
                k_new1 = k + g(3)%size/2
                k_old1 = k
                k_new2 = k
                k_old2 = k + g(3)%size/2

                p_dst(:, k_new1) = p_org(:, k_old1)
                p_dst(:, k_new2) = p_org(:, k_old2)
            end do
            do k = 1, g(3)%size
                p_org(:, k) = p_dst(:, k)
            end do
        end if

        call dfftw_execute_dft(fft_plan_bz, p_org, p_dst)

#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLAB_MPI_TRPB_K(r_in, r_out, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            nullify (r_in, r_out)
        end if
#endif

        nullify (p_org, p_dst)

        return
    end subroutine OPR_FOURIER_B_Z_EXEC

    ! #######################################################################
    ! #######################################################################
    subroutine OPR_FOURIER_SPECTRA_3D(nx, ny, nz, isize_psd, u, psd, wrk1d)

        integer(wi), intent(IN) :: nx, ny, nz, isize_psd
        real(wp), dimension(isize_txc_dimz, nz), intent(IN) :: u
        real(wp), dimension(isize_psd), intent(OUT) :: psd
        real(wp), dimension(isize_psd), intent(INOUT) :: wrk1d

        ! -----------------------------------------------------------------------
        integer(wi) i, j, r, iglobal, kglobal, ip
        real(wp) fr, fi, fj, fk

        ! #######################################################################
        psd = 0.0_wp

        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif
            if (kglobal <= g(3)%size/2 + 1) then; fk = real(kglobal - 1, wp)
            else; fk = -real(g(3)%size + 1 - kglobal, wp)
            end if

            do j = 1, ny
                if (j <= g(2)%size/2 + 1) then; fj = real(j - 1, wp)
                else; fj = -real(g(2)%size + 1 - j, wp)
                end if

                do i = 1, nx/2 + 1
#ifdef USE_MPI
                    iglobal = i + ims_offset_i/2
#else
                    iglobal = i
#endif
                    if (iglobal <= g(1)%size/2 + 1) then; fi = real(iglobal - 1, wp)
                    else; fi = -real(g(1)%size + 1 - iglobal, wp)
                    end if

                    fr = ceiling(sqrt(real(fi**2 + fj**2 + fk**2, wp)))

                    ! -----------------------------------------------------------------------
                    ! psd
                    ! -----------------------------------------------------------------------
                    ip = (nx + 2)*(j - 1) + 2*i

                    r = int(fr)
                    if (r == 0) then ! zero is not written
                        !              print*, i,j,k, r
                    else
                        psd(r) = psd(r) + u(ip - 1, k)**2 + u(ip, k)**2
                    end if

                end do
            end do
        end do

#ifdef USE_MPI
        wrk1d = psd
        call MPI_Reduce(wrk1d, psd, isize_psd, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
        if (ims_pro == 0) then
            psd = wrk1d
        end if
#endif

        return
    end subroutine OPR_FOURIER_SPECTRA_3D

end module OPR_FOURIER
