#include "dns_error.h"

module OPR_Fourier
    use TLab_Constants, only: wp, wi, pi_wp, efile
    use TLab_Memory, only: isize_txc_field, isize_txc_dimz
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_Arrays, only: wrk3d
    use TLab_Pointers_C, only: c_wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Grid
    use, intrinsic :: iso_c_binding
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
    use TLabMPI_Transpose
#endif
    implicit none
    private

    public :: OPR_Fourier_Initialize
    public :: OPR_Fourier_X_Forward, OPR_Fourier_X_Backward     ! Main routines, used in Poisson solver; optimize
    public :: OPR_Fourier_Z_Forward, OPR_Fourier_Z_Backward
    public :: OPR_Fourier_F                                     ! Minor routines, less frequently used
    public :: OPR_Fourier_B
    public :: OPR_Fourier_CONVOLUTION_FXZ
    public :: OPR_Fourier_ComputePSD
    public :: OPR_Fourier_SetPSD, OPR_Fourier_SetPSD_2d

    logical, public :: fft_x_on = .false.                       ! Switch on and off the FFT in the corresponding directions.
    logical, public :: fft_y_on = .false.
    logical, public :: fft_z_on = .false.

    ! -----------------------------------------------------------------------
    integer(wi) size_fft_x, size_fft_y, size_fft_z

    type(c_ptr) :: fft_plan_fx, fft_plan_bx
    type(c_ptr) :: fft_plan_fy, fft_plan_by!, fft_plan_fy1d, fft_plan_by1d
    type(c_ptr) :: fft_plan_fz, fft_plan_bz
#ifdef USE_MPI
    type(tmpi_transpose_dt), public :: tmpi_plan_fftx
    type(tmpi_transpose_dt), public :: tmpi_plan_fftz
#endif

    logical :: fft_reordering_k = .false.
    logical :: fft_reordering_i = .false.

    integer(wi) k
    complex(wp), pointer :: c_tmp1(:, :) => null(), c_in(:, :) => null(), c_out(:, :) => null(), c_in1(:, :) => null()
    complex(wp), pointer :: c_tmp2(:) => null()

contains
    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_Initialize()
#ifdef USE_MPI
        use mpi_f08, only: MPI_DOUBLE_COMPLEX
#endif
        use TLab_Arrays, only: txc

#ifdef USE_FFTW
#include "fftw3.f03"
#endif
        ! -----------------------------------------------------------------------
        integer(wi) stride, isize_disp, nlines

        ! #######################################################################
#ifndef USE_FFTW
        call TLab_Write_ASCII(efile, __FILE__//'. FFTW needed for POISSON solver.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
#endif

        if (mod(imax, 2) /= 0) then
            call TLab_Write_ASCII(efile, __FILE__//'. Imax must be a multiple of 2 for the FFT operations.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        ! -----------------------------------------------------------------------
        ! Oz direction
        size_fft_z = z%size

        if (size_fft_z > 1) then
            fft_z_on = .true.

#ifdef USE_MPI
            if (ims_npro_k > 1) then
                tmpi_plan_fftz = TLabMPI_Trp_PlanK(kmax, isize_txc_dimz/2, &
                                                   locType=MPI_DOUBLE_COMPLEX, &
                                                   message='Oz FFTW in Poisson solver.')
                nlines = tmpi_plan_fftz%nlines

            else
#endif
                nlines = (imax/2 + 1)*jmax
#ifdef USE_MPI
            end if
#endif

            stride = nlines

#ifdef _DEBUG
            call dfftw_plan_many_dft(fft_plan_fz, 1, size_fft_z, nlines, &
                                     txc(:, 1), size_fft_z, stride, 1, &
                                     wrk3d, size_fft_z, stride, 1, &
                                     FFTW_FORWARD, FFTW_ESTIMATE)

            call dfftw_plan_many_dft(fft_plan_bz, 1, size_fft_z, nlines, &
                                     txc(:, 1), size_fft_z, stride, 1, &
                                     wrk3d, size_fft_z, stride, 1, &
                                     FFTW_BACKWARD, FFTW_ESTIMATE)
#else
            call dfftw_plan_many_dft(fft_plan_fz, 1, size_fft_z, nlines, &
                                     txc(:, 1), size_fft_z, stride, 1, &
                                     wrk3d, size_fft_z, stride, 1, &
                                     FFTW_FORWARD, FFTW_MEASURE)

            call dfftw_plan_many_dft(fft_plan_bz, 1, size_fft_z, nlines, &
                                     txc(:, 1), size_fft_z, stride, 1, &
                                     wrk3d, size_fft_z, stride, 1, &
                                     FFTW_BACKWARD, FFTW_MEASURE)
#endif
        end if

        ! -----------------------------------------------------------------------
        ! Ox direction
        size_fft_x = x%size
        if (size_fft_x > 1) then
            fft_x_on = .true.

#ifdef USE_MPI
            if (ims_npro_i > 1) then
                ! Extended with the Nyquist frequency
                tmpi_plan_fftx = TLabMPI_Trp_PlanI(imax/2 + 1, jmax*kmax, &
                                                   locType=MPI_DOUBLE_COMPLEX, &
                                                   message='extended Ox FFTW in Poisson solver.')

                if (tmpi_plan_fftx%nlines /= tmpi_plan_dx%nlines) then
                    call TLab_Write_ASCII(efile, __FILE__//'. Error in the size in the transposition arrays.')
                    call TLab_Stop(DNS_ERROR_UNDEVELOP)
                end if

                nlines = tmpi_plan_dx%nlines
                isize_disp = (imax/2 + 1)*ims_npro_i

            else
#endif
                nlines = jmax*kmax
                isize_disp = size_fft_x/2 + 1

#ifdef USE_MPI
            end if
#endif

#ifdef _DEBUG
            call dfftw_plan_many_dft_r2c(fft_plan_fx, 1, size_fft_x, nlines, &
                                         txc(:, 1), size_fft_x, 1, size_fft_x, &
                                         wrk3d, size_fft_x/2 + 1, 1, isize_disp, &
                                         FFTW_ESTIMATE)
            call dfftw_plan_many_dft_c2r(fft_plan_bx, 1, size_fft_x, nlines, &
                                         txc(:, 1), size_fft_x/2 + 1, 1, isize_disp, &
                                         wrk3d, size_fft_x, 1, size_fft_x, &
                                         FFTW_ESTIMATE)
#else
            call dfftw_plan_many_dft_r2c(fft_plan_fx, 1, size_fft_x, nlines, &
                                         txc(:, 1), size_fft_x, 1, size_fft_x, &
                                         wrk3d, size_fft_x/2 + 1, 1, isize_disp, &
                                         FFTW_MEASURE)
            call dfftw_plan_many_dft_c2r(fft_plan_bx, 1, size_fft_x, nlines, &
                                         txc(:, 1), size_fft_x/2 + 1, 1, isize_disp, &
                                         wrk3d, size_fft_x, 1, size_fft_x, &
                                         FFTW_MEASURE)
#endif

        end if

        ! -----------------------------------------------------------------------
        ! Oy direction
        size_fft_y = y%size
        if (size_fft_y > 1) then
            fft_y_on = .true.

            nlines = imax/2 + 1    !(imax/2 + 1)*kmax
            stride = imax/2 + 1

#ifdef _DEBUG
            call dfftw_plan_many_dft(fft_plan_fy, 1, size_fft_y, nlines, &
                                     txc(:, 1), size_fft_y, stride, 1, &
                                     wrk3d, size_fft_y, stride, 1, &
                                     FFTW_FORWARD, FFTW_ESTIMATE)

            call dfftw_plan_many_dft(fft_plan_by, 1, size_fft_y, nlines, &
                                     txc(:, 1), size_fft_y, stride, 1, &
                                     wrk3d, size_fft_y, stride, 1, &
                                     FFTW_BACKWARD, FFTW_ESTIMATE)
#else
            call dfftw_plan_many_dft(fft_plan_fy, 1, size_fft_y, nlines, &
                                     txc(:, 1), size_fft_y, stride, 1, &
                                     wrk3d, size_fft_y, stride, 1, &
                                     FFTW_FORWARD, FFTW_MEASURE)

            call dfftw_plan_many_dft(fft_plan_by, 1, size_fft_y, nlines, &
                                     txc(:, 1), size_fft_y, stride, 1, &
                                     wrk3d, size_fft_y, stride, 1, &
                                     FFTW_BACKWARD, FFTW_MEASURE)
#endif
        end if

        return
    end subroutine OPR_Fourier_Initialize

    !########################################################################
    !# Forward/Backward Fourier transform.
    !#
    !# Each z-plane contains jmax lines of nx/2+1 complex numbers. The element nx/2+1 is the
    !# Nyquist frequency.
    !# If OX MPI decomposition, one can reorganize modes
    !# to homogeneize arrays across PEs.
    !#
    !########################################################################
    subroutine OPR_Fourier_X_Forward(nx, ny, nz, in, out)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: in(nx*ny*nz)
        complex(wp), intent(out) :: out(isize_txc_field)

        target out

        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        integer(wi) i, iold, inew, ip, isize_line
        complex(wp), pointer :: wrk1(:, :) => null()
        real(wp), pointer :: r_out(:) => null()
#endif

        ! #######################################################################

#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call c_f_pointer(c_loc(wrk3d), wrk1, shape=[(nx/2 + 1)*ims_npro_i, tmpi_plan_dx%nlines])
            call c_f_pointer(c_loc(out), r_out, shape=[isize_txc_field])

            call TLabMPI_Trp_ExecI_Forward(in(:), r_out(:), tmpi_plan_dx)

            call dfftw_execute_dft_r2c(fft_plan_fx, r_out, wrk1)

            if (fft_reordering_i) then      ! reorganize a (FFTW make a stride in a already before)
                isize_line = nx/2 + 1
                do k = 1, tmpi_plan_dx%nlines
                    inew = isize_line*ims_npro_i
                    iold = size_fft_x/2 + 1
                    wrk1(inew, k) = wrk1(iold, k)
                    do ip = ims_npro_i, 2, -1
                        do i = nx/2, 1, -1
                            inew = (ip - 1)*isize_line + i
                            iold = (ip - 1)*nx/2 + i
                            wrk1(inew, k) = wrk1(iold, k)
                        end do
                    end do
                end do
            end if

            call TLabMPI_Trp_ExecI_Backward(wrk1(:, 1), out(:), tmpi_plan_fftx)

            nullify (wrk1, r_out)

        else
#endif
            call dfftw_execute_dft_r2c(fft_plan_fx, in, out)

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_Fourier_X_Forward

    !########################################################################
    !########################################################################
    subroutine OPR_Fourier_X_Backward(nx, ny, nz, in, out)
        integer(wi) nx, ny, nz
        complex(wp), intent(in) :: in(isize_txc_field)
        real(wp), intent(out) :: out(nx*ny*nz)

        target in, out

        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        integer(wi) i, ip, iold, inew, isize_line
        real(wp), pointer :: r_in(:) => null()
#endif

        !########################################################################
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            call c_f_pointer(c_loc(in), r_in, shape=[isize_txc_field])
            call c_f_pointer(c_loc(out), c_out, shape=[(nx/2 + 1)*ims_npro_i, nz])

            call TLabMPI_Trp_ExecI_Forward(in(:), c_out(:, 1), tmpi_plan_fftx)

            if (fft_reordering_i) then      ! reorganize a (FFTW make a stride in a already before)
                isize_line = nx/2 + 1
                do k = 1, tmpi_plan_dx%nlines !tmpi_plan_fftx1%nlines
                    do ip = 2, ims_npro_i
                        do i = 1, nx/2
                            iold = (ip - 1)*isize_line + i
                            inew = (ip - 1)*nx/2 + i
                            c_out(inew, k) = c_out(iold, k)
                        end do
                    end do
                    iold = ims_npro_i*isize_line
                    inew = size_fft_x/2 + 1
                    c_out(inew, k) = c_out(iold, k)
                end do
            end if

            call dfftw_execute_dft_c2r(fft_plan_bx, c_out, r_in)

            call TLabMPI_Trp_ExecI_Backward(r_in(:), out(:), tmpi_plan_dx) !tmpi_plan_fftx1)

            nullify (r_in, c_out)

        else
#endif
            call dfftw_execute_dft_c2r(fft_plan_bx, in, out)

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine OPR_Fourier_X_Backward

    !########################################################################
    !########################################################################
    subroutine OPR_Fourier_Z_Forward(in, out)
        complex(wp), intent(inout), target :: in(*)     ! in might be overwritten
        complex(wp), intent(out), target :: out(*)

        ! -----------------------------------------------------------------------
        complex(wp), pointer :: p_org(:, :), p_dst(:, :)
        integer(wi) k_old1, k_old2, k_new1, k_new2

        ! #######################################################################
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_Trp_ExecK_Forward(in, out, tmpi_plan_fftz)
            p_org(1:tmpi_plan_fftz%nlines, 1:size_fft_z) => out(1:isize_txc_field)
            p_dst(1:tmpi_plan_fftz%nlines, 1:size_fft_z) => in(1:isize_txc_field)
        else
#endif
            p_org(1:isize_txc_dimz/2, 1:size_fft_z) => in(1:isize_txc_field)
            p_dst(1:isize_txc_dimz/2, 1:size_fft_z) => out(1:isize_txc_field)
#ifdef USE_MPI
        end if
#endif

        call dfftw_execute_dft(fft_plan_fz, p_org, p_dst)

        if (fft_reordering_k) then                    ! re-shuffle spectra in z
            do k = 1, size_fft_z/2
                k_old1 = k + size_fft_z/2
                k_new1 = k
                k_old2 = k
                k_new2 = k + size_fft_z/2

                p_org(:, k_new1) = p_dst(:, k_old1)
                p_org(:, k_new2) = p_dst(:, k_old2)
            end do
            do k = 1, size_fft_z
                p_dst(:, k) = p_org(:, k)
            end do
        end if

#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_Trp_ExecK_Backward(in, out, tmpi_plan_fftz)
        end if
#endif

        nullify (p_org, p_dst)

        return
    end subroutine OPR_Fourier_Z_Forward

    !########################################################################
    !########################################################################
    subroutine OPR_Fourier_Z_Backward(in, out)
        complex(wp), intent(inout), target :: in(*)     ! in might be overwritten
        complex(wp), intent(out), target :: out(*)

        ! -----------------------------------------------------------------------
        complex(wp), pointer :: p_org(:, :), p_dst(:, :)
        integer(wi) k_old1, k_old2, k_new1, k_new2

        ! #######################################################################
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_Trp_ExecK_Forward(in, out, tmpi_plan_fftz)
            p_org(1:tmpi_plan_fftz%nlines, 1:size_fft_z) => out(1:isize_txc_field)
            p_dst(1:tmpi_plan_fftz%nlines, 1:size_fft_z) => in(1:isize_txc_field)
        else
#endif
            p_org(1:isize_txc_dimz/2, 1:size_fft_z) => in(1:isize_txc_field)
            p_dst(1:isize_txc_dimz/2, 1:size_fft_z) => out(1:isize_txc_field)
#ifdef USE_MPI
        end if
#endif

        if (fft_reordering_k) then                    ! re-shuffle spectra in z
            do k = 1, size_fft_z/2
                k_new1 = k + size_fft_z/2
                k_old1 = k
                k_new2 = k
                k_old2 = k + size_fft_z/2

                p_dst(:, k_new1) = p_org(:, k_old1)
                p_dst(:, k_new2) = p_org(:, k_old2)
            end do
            do k = 1, size_fft_z
                p_org(:, k) = p_dst(:, k)
            end do
        end if

        call dfftw_execute_dft(fft_plan_bz, p_org, p_dst)

#ifdef USE_MPI
        if (ims_npro_k > 1) then
            call TLabMPI_Trp_ExecK_Backward(in, out, tmpi_plan_fftz)
        end if
#endif

        nullify (p_org, p_dst)

        return
    end subroutine OPR_Fourier_Z_Backward

    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_F(flag_mode, nx, ny, nz, in, out, tmp1)
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

        fft_reordering_i = .true.

        if (flag_mode == 3 .and. size_fft_y > 1) then ! 3D FFT (unless 2D sim)
            if (size_fft_z > 1) then
                call OPR_Fourier_X_Forward(nx, ny, nz, in, c_out)
                call OPR_Fourier_Z_Forward(c_out, c_tmp1) ! out might be overwritten
            else
                call OPR_Fourier_X_Forward(nx, ny, nz, in, c_tmp1)
            end if

            do k = 1, nz
                call dfftw_execute_dft(fft_plan_fy, c_tmp1(:, k), c_out(:, k))
            end do
            ! call dfftw_execute_dft(fft_plan_fy, c_tmp1, c_out)

        else
            if (flag_mode == 2 .and. size_fft_z > 1) then ! 2D FFT (unless 1D sim)
                call OPR_Fourier_X_Forward(nx, ny, nz, in, c_tmp1)
                call OPR_Fourier_Z_Forward(c_tmp1, c_out) ! tmp1 might be overwritten
            else
                call OPR_Fourier_X_Forward(nx, ny, nz, in, c_out)
            end if
        end if

        fft_reordering_i = .false.

        return
    end subroutine OPR_Fourier_F

    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_B(flag_mode, nx, ny, nz, in, out)
        integer(wi) :: flag_mode  ! 1D, 2D or 3D
        integer(wi) :: nx, ny, nz
        real(wp), intent(INOUT) :: in(isize_txc_dimz, nz)
        real(wp), intent(OUT) :: out(isize_txc_field)

        target in

        ! #######################################################################
        ! Pass memory address from real array to complex array
        call c_f_pointer(c_loc(in), c_in, shape=[isize_txc_dimz/2, nz])

        fft_reordering_i = .true.

        if (flag_mode == 3 .and. size_fft_y > 1) then ! 3D FFT (unless 2D sim)
            do k = 1, nz
                call dfftw_execute_dft(fft_plan_by, c_in(:, k), c_wrk3d(:, k))
            end do
            ! call dfftw_execute_dft(fft_plan_by, c_in, c_wrk3d)

            if (size_fft_z > 1) then
                call OPR_Fourier_Z_Backward(c_wrk3d, c_in)            ! wrk3d might be overwritten
                call OPR_Fourier_X_Backward(nx, ny, nz, c_in, out)    ! c_in might be overwritten
            else
                call OPR_Fourier_X_Backward(nx, ny, nz, c_wrk3d, out) ! wrk3d might be overwritten
            end if

        else
            if (flag_mode == 2 .and. size_fft_z > 1) then ! 2D FFT (unless 1D sim)
                call OPR_Fourier_Z_Backward(c_in, c_wrk3d)            ! in might be overwritten
                call OPR_Fourier_X_Backward(nx, ny, nz, c_wrk3d, out) ! wrk3d might be overwritten
            else
                call OPR_Fourier_X_Backward(nx, ny, nz, c_in, out)    ! c_in might be overwritten
            end if

        end if

        fft_reordering_i = .false.

        return
    end subroutine OPR_Fourier_B

    ! #######################################################################
    ! #######################################################################
    ! Calculate spectrum in array in1
    ! Calculate correlation in array in2
    subroutine OPR_Fourier_CONVOLUTION_FXZ(flag1, flag2, nx, ny, nz, in1, in2, tmp1, tmp2)
        character(len=*), intent(in) :: flag1
        integer(wi), intent(in) :: flag2, nx, ny, nz
        real(wp), dimension(isize_txc_field), intent(INOUT) :: in1, in2, tmp1, tmp2

        target in1, tmp1, tmp2

        ! #######################################################################
        ! Pass memory address from real array to complex array
        call c_f_pointer(c_loc(in1), c_in1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field])

        fft_reordering_k = .true.
        fft_reordering_i = .true.

        if (size_fft_z > 1) then
            call OPR_Fourier_X_Forward(nx, ny, nz, in1, c_tmp2)
            call OPR_Fourier_Z_Forward(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_Fourier_X_Forward(nx, ny, nz, in1, c_tmp1)
        end if

        select case (trim(adjustl(flag1)))

        case ('auto')        ! Auto-spectra
            c_in1 = c_tmp1*conjg(c_tmp1)

        case ('cross')       ! Cross-spectra
            if (size_fft_z > 1) then
                call OPR_Fourier_X_Forward(nx, ny, nz, in2, c_tmp2)
                call OPR_Fourier_Z_Forward(c_tmp2, c_in1) ! tmp2 might be overwritten
            else
                call OPR_Fourier_X_Forward(nx, ny, nz, in2, c_in1)
            end if
            c_in1 = c_in1*conjg(c_tmp1)

        end select

        ! -----------------------------------------------------------------------
        if (flag2 == 2) then         ! Calculate correlation in array in2
            if (size_fft_z > 1) then
                call OPR_Fourier_Z_Backward(c_in1, c_tmp1)            ! in1 might be overwritten
                call OPR_Fourier_X_Backward(nx, ny, nz, c_tmp1, in2)  ! tmp1 might be overwritten
            else
                call OPR_Fourier_X_Backward(nx, ny, nz, c_in1, in2)   ! c_in1 might be overwritten
            end if
        end if

        ! -----------------------------------------------------------------------
        fft_reordering_k = .false.
        fft_reordering_i = .false.

        return
    end subroutine OPR_Fourier_CONVOLUTION_FXZ

    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_ComputePSD(nx, ny, nz, u, psd)
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS, only: ims_pro, ims_err
#endif
        integer(wi), intent(in) :: nx, ny, nz
        complex(wp), intent(in) :: u(nx/2 + 1, ny, nz)
        real(wp), intent(out) :: psd(:)

        ! -----------------------------------------------------------------------
        integer(wi) i, j, r, iglobal, kglobal
        real(wp) fr, fi, fj, fk
#ifdef USE_MPI
        integer(wi) isize_psd
#endif

        ! #######################################################################
        psd(:) = 0.0_wp

        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif
            if (kglobal <= size_fft_z/2 + 1) then
                fk = real(kglobal - 1, wp)
            else
                fk = -real(size_fft_z + 1 - kglobal, wp)
            end if

            do j = 1, ny
                if (j <= size_fft_y/2 + 1) then
                    fj = real(j - 1, wp)
                else
                    fj = -real(size_fft_y + 1 - j, wp)
                end if

                do i = 1, nx/2 + 1
#ifdef USE_MPI
                    iglobal = i + ims_offset_i/2
#else
                    iglobal = i
#endif
                    fi = real(iglobal - 1, wp)

                    fr = ceiling(sqrt(real(fi**2 + fj**2 + fk**2, wp)))

                    ! -----------------------------------------------------------------------
                    ! psd
                    r = int(fr)
                    if (r == 0) then ! zero is not written
                        ! print *, i, j, k, r
                    else
                        psd(r) = psd(r) + abs(u(i, j, k))**2
                    end if

                end do
            end do
        end do

#ifdef USE_MPI
        isize_psd = size(psd)
        call MPI_Reduce(psd, wrk3d, isize_psd, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
        if (ims_pro == 0) then
            psd(1:isize_psd) = wrk3d(1:isize_psd)
        end if
#endif

        return
    end subroutine OPR_Fourier_ComputePSD

    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_SetPSD(nx, ny, nz, u, psd, locPhase)
        use Distributions
        use FDM, only: g

        integer(wi) nx, ny, nz
        complex(wp), intent(inout) :: u(nx/2 + 1, ny, nz)
        type(distributions_dt), intent(in) :: psd
        real(wp), intent(in), optional :: locPhase(nx/2 + 1, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) i, j, iglobal, jglobal, kglobal
        real(wp) pow_dst, pow_org, phase
        real(wp) f, fi, fj, fk

        ! #######################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif
            if (kglobal <= size_fft_z/2 + 1) then
                fk = real(kglobal - 1, wp)/g(3)%scale
            else
                fk = -real(size_fft_z + 1 - kglobal, wp)/g(3)%scale
            end if

            do j = 1, ny
                jglobal = j ! No MPI decomposition along Oy
                if (jglobal <= size_fft_y/2 + 1) then
                    fj = real(jglobal - 1, wp)/g(2)%scale
                else
                    fj = -real(size_fft_y + 1 - jglobal, wp)/g(2)%scale
                end if

                do i = 1, nx/2 + 1
#ifdef USE_MPI
                    iglobal = i + ims_offset_i/2
#else
                    iglobal = i
#endif
                    fi = real(iglobal - 1, wp)/g(1)%scale

                    f = sqrt(fi**2 + fj**2 + fk**2)

                    ! -----------------------------------------------------------------------
                    ! target psd
                    pow_dst = Distributions_Compute(psd, f)

                    if (f == 0.0_wp) then
                        pow_dst = 0.0_wp

                    else
                        if (size_fft_y == 1 .or. size_fft_z == 1) then ! 2D spectrum
                            pow_dst = pow_dst/(pi_wp*f)
                        else
                            pow_dst = pow_dst/(2*pi_wp*f**2)
                        end if

                    end if

                    ! -----------------------------------------------------------------------
                    ! phase and scaling of complex data
                    if (pow_dst > 0.0_wp) pow_dst = sqrt(pow_dst)

                    if (present(locPhase)) then
                        if (iglobal == 1 .or. iglobal == size_fft_x/2 + 1) then
                            phase = 0.0_wp
                        else
                            phase = (locPhase(i, j, k) - 0.5_wp)*2.0_wp*pi_wp
                        end if
                        u(i, j, k) = cmplx(pow_dst*cos(phase), pow_dst*sin(phase), wp)

                    else
                        pow_org = abs(u(i, j, k))

                        if (pow_org > 0.0_wp) pow_dst = pow_dst/pow_org

                        u(i, j, k) = u(i, j, k)*pow_dst

                    end if

                end do
            end do
        end do

        return
    end subroutine OPR_Fourier_SetPSD

    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Fourier_SetPSD_2d(nx, ny, nz, u, psd)
        use Distributions
        use FDM, only: g

        integer(wi) nx, ny, nz
        complex(wp), intent(inout) :: u(nx/2 + 1, ny, nz)
        type(distributions_dt), intent(in) :: psd

        ! -----------------------------------------------------------------------
        integer(wi) i, iglobal, kglobal
        real(wp) pow_dst
        real(wp) f, fi, fk

        ! #######################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif
            if (kglobal <= size_fft_z/2 + 1) then
                fk = real(kglobal - 1, wp)/g(3)%scale
            else
                fk = -real(size_fft_z + 1 - kglobal, wp)/g(3)%scale
            end if

            do i = 1, nx/2 + 1
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif
                fi = real(iglobal - 1, wp)/g(1)%scale

                f = sqrt(fi**2 + fk**2)

                ! -----------------------------------------------------------------------
                ! target psd
                pow_dst = Distributions_Compute(psd, f)

                ! -----------------------------------------------------------------------
                ! scaling of complex data
                u(i, :, k) = u(i, :, k)*pow_dst

            end do
        end do

        return
    end subroutine OPR_Fourier_SetPSD_2d

end module OPR_Fourier
