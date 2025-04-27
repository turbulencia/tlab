#include "dns_error.h"

! Split the routines into the ones that are initialized and the ones that not?
! If not initialized, you can enter with any jmax, but the periodic directions need to be the global ones because of OPR_Fourier.
module OPR_Elliptic
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: lfile
    use TLab_Memory, only: TLab_Allocate_Real
    use TLab_Memory, only: imax, jmax, kmax, isize_txc_field
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLab_Arrays, only: wrk1d, wrk2d, wrk3d
    use TLab_Pointers_C, only: c_wrk3d
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k, ims_pro_i
#endif
    use FDM, only: fdm_dt
    use FDM_Integral
    use OPR_Fourier
    use OPR_ODES
    use OPR_Partial, only: OPR_Partial_Y, OPR_P1
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    public :: OPR_Elliptic_Initialize
    public :: OPR_Poisson
    public :: OPR_Helmholtz

    ! -----------------------------------------------------------------------
    procedure(OPR_Poisson_interface) :: OPR_Poisson_dt      ! Implicit pointer (Procedure type)
    abstract interface
        subroutine OPR_Poisson_interface(nx, ny, nz, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
            use TLab_Constants, only: wi, wp
            use FDM, only: fdm_dt
            integer(wi), intent(in) :: nx, ny, nz
            integer, intent(in) :: ibc                                      ! Dirichlet/Neumman BCs at jmin/jmax: BCS_DD, BCS_ND, BCS_DN, BCS_NN
            real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
            real(wp), intent(inout), target :: tmp1(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
            real(wp), intent(inout), target :: tmp2(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
            real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
            real(wp), intent(out), optional :: dpdy(nx, ny, nz)             ! Vertical derivative of solution
        end subroutine
    end interface
    procedure(OPR_Poisson_dt), pointer :: OPR_Poisson

    procedure(OPR_Helmholtz_interface) :: OPR_Helmholtz_dt  ! Implicit pointer (Procedure type)
    abstract interface
        subroutine OPR_Helmholtz_interface(nx, ny, nz, ibc, alpha, p, tmp1, tmp2, bcs_hb, bcs_ht)
            use TLab_Constants, only: wi, wp
            use FDM, only: fdm_dt
            integer(wi), intent(in) :: nx, ny, nz
            integer, intent(in) :: ibc                                      ! Dirichlet/Neumman BCs at jmin/jmax: BCS_DD, BCS_ND, BCS_DN, BCS_NN
            real(wp), intent(in) :: alpha
            real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
            real(wp), intent(inout), target :: tmp1(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
            real(wp), intent(inout), target :: tmp2(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
            real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
        end subroutine
    end interface
    procedure(OPR_Helmholtz_dt), pointer :: OPR_Helmholtz

    real(wp) norm
    integer(wi) i_sing(2), k_sing(2)                                ! singular modes
    integer(wi) i, k, i_max, isize_line

    type(fdm_dt) fdm_loc                                            ! scheme used for the elliptic solvers

    type(fdm_integral_dt), allocatable :: fdm_int1(:, :, :)         ! factorized method
    real(wp), allocatable, target :: rhs_b(:, :), rhs_t(:, :)       ! rhs to free memory space
    type(fdm_integral_dt) :: fdm_int1_loc(2)

    type(fdm_integral_dt), allocatable :: fdm_int2(:, :)            ! direct method
    real(wp), allocatable, target :: rhs_d(:, :)                    ! rhs to free memory space
    type(fdm_integral_dt) :: fdm_int2_loc

    real(wp), allocatable :: lambda(:, :)

    complex(wp), pointer :: c_tmp1(:) => null(), c_tmp2(:) => null()
    real(wp), pointer :: p_wrk3d(:, :, :) => null()

contains
    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Elliptic_Initialize(inifile)
        use FDM, only: g, FDM_CreatePlan
        use FDM_Derivative, only: FDM_COM4_DIRECT, FDM_COM6_DIRECT

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer imode_elliptic
        integer, parameter :: TYPE_FACTORIZE = 1
        integer, parameter :: TYPE_DIRECT = 2

        integer(wi) :: ndl, ndr, nd
        character*512 sRes
        character*32 bakfile

        integer(wi) iglobal, kglobal
        integer(wi) fft_offset_i, fft_offset_k

        ! ###################################################################
        ! Reading
        bakfile = trim(adjustl(inifile))//'.bak'

        imode_elliptic = TYPE_FACTORIZE         ! default is the finite-difference method used for the derivatives
        fdm_loc%der1%mode_fdm = g(2)%der1%mode_fdm      ! to impose zero divergence down to round-off error in the interior points
        fdm_loc%der2%mode_fdm = g(2)%der2%mode_fdm

        call ScanFile_Char(bakfile, inifile, 'Main', 'EllipticOrder', 'void', sRes)
        if (trim(adjustl(sRes)) == 'compactdirect4') then
            imode_elliptic = TYPE_DIRECT
            fdm_loc%der1%mode_fdm = FDM_COM4_DIRECT
            fdm_loc%der2%mode_fdm = FDM_COM4_DIRECT
        else if (trim(adjustl(sRes)) == 'compactdirect6') then
            imode_elliptic = TYPE_DIRECT
            fdm_loc%der1%mode_fdm = FDM_COM6_DIRECT
            fdm_loc%der2%mode_fdm = FDM_COM6_DIRECT
        end if

        ! ###################################################################
        ! Initializing
        call FDM_CreatePlan(g(2)%nodes, fdm_loc)

        isize_line = imax/2 + 1

        allocate (lambda(kmax, isize_line))
        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        select case (imode_elliptic)
        case (TYPE_FACTORIZE)
            OPR_Poisson => OPR_Poisson_FourierXZ_Factorize
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Factorize

            ndl = fdm_loc%der1%nb_diag(1)
            ndr = fdm_loc%der1%nb_diag(2)
            nd = ndl
            allocate (fdm_int1(2, kmax, isize_line))
            call TLab_Allocate_Real(__FILE__, rhs_b, [g(2)%size, nd], 'rhs_b')
            call TLab_Allocate_Real(__FILE__, rhs_t, [g(2)%size, nd], 'rhs_t')

            if (stagger_on) then
                i_sing = [1, 1]                 ! only one singular mode + different modified wavenumbers
                k_sing = [1, 1]
            else
                i_sing = [1, g(1)%size/2 + 1]   ! global indexes, transformed below to task-local indexes.
                k_sing = [1, g(3)%size/2 + 1]
            end if

        case (TYPE_DIRECT)
            OPR_Poisson => OPR_Poisson_FourierXZ_Direct
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Direct

            ndl = fdm_loc%der2%nb_diag(1)
            ndr = fdm_loc%der2%nb_diag(2)
            nd = ndl
            allocate (fdm_int2(kmax, isize_line))
            call TLab_Allocate_Real(__FILE__, rhs_d, [g(2)%size, nd], 'rhs_d')

            i_sing = [1, 1]                     ! 2nd order FDMs are non-zero at Nyquist
            k_sing = [1, 1]

        end select

#ifdef USE_MPI
        ! fft_offset_i = ims_offset_i/2
        fft_offset_i = ims_pro_i*isize_line
        fft_offset_k = ims_offset_k

#else
        fft_offset_i = 0
        fft_offset_k = 0
#endif

        i_sing = i_sing - [fft_offset_i, fft_offset_i]              ! Singular modes in task-local variables
        k_sing = k_sing - [fft_offset_k, fft_offset_k]
        i_max = min(g(1)%size/2 + 1 - fft_offset_i, isize_line)     ! Maximum mode is x direction

        do i = 1, i_max
#ifdef USE_MPI

            iglobal = i + fft_offset_i
#else
            iglobal = i
#endif

            do k = 1, kmax
#ifdef USE_MPI
                kglobal = k + fft_offset_k
#else
                kglobal = k
#endif

                select case (imode_elliptic)
                case (TYPE_FACTORIZE)
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda(k, i) = g(1)%der1%mwn(iglobal)**2.0 + g(3)%der1%mwn(kglobal)**2.0
                    else
                        lambda(k, i) = g(1)%der1%mwn(iglobal)**2.0
                    end if

                    call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%der1, &
                                             sqrt(lambda(k, i)), BCS_MIN, fdm_int1(BCS_MIN, k, i))

                    call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%der1, &
                                             -sqrt(lambda(k, i)), BCS_MAX, fdm_int1(BCS_MAX, k, i))

                    if (any(i_sing == i) .and. any(k_sing == k)) then
                    else                                        ! free memory that is independent of lambda
                        rhs_b(:, :) = fdm_int1(BCS_MIN, k, i)%rhs(:, :)
                        if (allocated(fdm_int1(BCS_MIN, k, i)%rhs)) deallocate (fdm_int1(BCS_MIN, k, i)%rhs)

                        rhs_t(:, :) = fdm_int1(BCS_MAX, k, i)%rhs(:, :)
                        if (allocated(fdm_int1(BCS_MAX, k, i)%rhs)) deallocate (fdm_int1(BCS_MAX, k, i)%rhs)

                    end if

                    ! idr = ndr/2 + 1
                    ! fdm_int1(:, k, i)%lhs(:, idr) = fdm_int1(:, k, i)%lhs(:, idr)*norm
                    ! fdm_int1(:, k, i)%lhs(:, idr + 1:ndr) = fdm_int1(:, k, i)%lhs(:, idr + 1:ndr)/norm

                case (TYPE_DIRECT)     ! only for case BCS_NN
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda(k, i) = g(1)%der2%mwn(iglobal) + g(3)%der2%mwn(kglobal)
                    else
                        lambda(k, i) = g(1)%der2%mwn(iglobal)
                    end if

                    ! Compatibility constraint. The reference value of p at the lower boundary will be set to zero
                    if (any(i_sing == i) .and. any(k_sing == k)) then
                        call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2, lambda(k, i), BCS_DN, fdm_int2(k, i))
                    else
                        call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2, lambda(k, i), BCS_NN, fdm_int2(k, i))
                    end if

                    ! free memory that is independent of lambda
                    rhs_d(:, :) = fdm_int2(k, i)%rhs(:, :)
                    if (allocated(fdm_int2(k, i)%rhs)) deallocate (fdm_int2(k, i)%rhs)

                end select

            end do
        end do

        return
    end subroutine OPR_Elliptic_Initialize

    !########################################################################
    !#
    !# Solve Lap p = f using Fourier in xOz planes, to rewrite the problem as
    !#
    !#     \hat{p}''-\lambda \hat{p} = \hat{f}
    !#
    !# where \lambda = kx^2+kz^2
    !#
    !# The reference value of p at the lower boundary is set to zero
    !#
    !########################################################################
    subroutine OPR_Poisson_FourierXZ_Factorize(nx, ny, nz, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(inout) :: tmp2(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
        real(wp), intent(out), optional :: dpdy(nx, ny, nz)             ! Vertical derivative of solution

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        real(wp) bcs(2, 2)

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field])
        p_wrk3d(1:2*ny, 1:nz, 1:nx/2 + 1) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        p(1:nx, 1, 1:nz) = bcs_hb(1:nx, 1:nz)           ! Passing boundary conditions in forcing array
        p(1:nx, ny, 1:nz) = bcs_ht(1:nx, 1:nz)

        if (fft_z_on) then
            call OPR_Fourier_X_Forward(nx, ny, nz, p, c_tmp2)
            call OPR_Fourier_Z_Forward(c_tmp2, c_tmp1)   ! tmp2 might be overwritten; cannot use wrk3d
        else
            call OPR_Fourier_X_Forward(nx, ny, nz, p, c_tmp1)
        end if

        tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
        ! ###################################################################
        ! Make x direction last one and leave y direction first
        call TLab_Transpose_COMPLEX(c_tmp1, isize_line, ny*nz, isize_line, c_tmp2, ny*nz)

#define f(j,k,i) tmp2(j,k,i)
#define v(j,k,i) tmp1(j,k,i)
#define u(j,k,i) p_wrk3d(j,k,i)

        ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
        do i = 1, i_max
            do k = 1, nz
                bcs(1:2, 1) = f(1:2, k, i)                  ! bottom boundary conditions
                bcs(1:2, 2) = f(2*ny - 1:2*ny, k, i)        ! top boundary conditions

                select case (ibc)
                case (BCS_NN)       ! Neumann & Neumann boundary conditions
                    if (any(i_sing == i) .and. any(k_sing == k)) then
                        call OPR_ODE2_Factorize_NN_Sing(2, fdm_int1(:, k, i), u(:, k, i), f(:, k, i), bcs, v(:, k, i), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_Factorize_NN(2, fdm_int1(:, k, i), rhs_b, rhs_t, &
                                                   u(:, k, i), f(:, k, i), bcs, v(1:, k, i), wrk1d, wrk2d)
                    end if

                case (BCS_DD)       ! Dirichlet & Dirichlet boundary conditions
                    if (any(i_sing == i) .and. any(k_sing == k)) then
                        call OPR_ODE2_Factorize_DD_Sing(2, fdm_int1(:, k, i), u(:, k, i), f(:, k, i), bcs, v(:, k, i), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_Factorize_DD(2, fdm_int1(:, k, i), rhs_b, rhs_t, &
                                                   u(:, k, i), f(:, k, i), bcs, v(:, k, i), wrk1d, wrk2d)
                    end if

                end select

            end do
        end do

        if (present(dpdy)) call TLab_Transpose_COMPLEX(c_tmp1, ny*nz, isize_line, ny*nz, c_tmp2, isize_line)
        call TLab_Transpose_COMPLEX(c_wrk3d, ny*nz, isize_line, ny*nz, c_tmp1, isize_line)

        ! ###################################################################
        ! Fourier field p (based on array tmp1)
        ! ###################################################################
        if (fft_z_on) then
            call OPR_Fourier_Z_Backward(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
            call OPR_Fourier_X_Backward(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
        else
            call OPR_Fourier_X_Backward(nx, ny, nz, c_tmp1, p)    ! tmp2 might be overwritten
        end if

        ! Fourier derivatives (based on array tmp2)
        if (present(dpdy)) then
            if (fft_z_on) then
                call OPR_Fourier_Z_Backward(c_tmp2, c_wrk3d)              ! tmp2 might be overwritten
                call OPR_Fourier_X_Backward(nx, ny, nz, c_wrk3d, dpdy)    ! wrk3d might be overwritten
            else
                call OPR_Fourier_X_Backward(nx, ny, nz, c_tmp2, dpdy)     ! tmp2 might be overwritten
            end if
        end if

        nullify (c_tmp1, c_tmp2, p_wrk3d)
#undef f
#undef v
#undef u

        return
    end subroutine OPR_Poisson_FourierXZ_Factorize

    !########################################################################
    !########################################################################
    subroutine OPR_Poisson_FourierXZ_Direct(nx, ny, nz, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        use FDM, only: g
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(inout) :: tmp2(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
        real(wp), intent(out), optional :: dpdy(nx, ny, nz)             ! Vertical derivative of solution

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        integer(wi), parameter :: bcs_p(2, 2) = 0                       ! For partial_y at the end

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field])
        p_wrk3d(1:2*ny, 1:nz, 1:nx/2 + 1) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        p(1:nx, 1, 1:nz) = bcs_hb(1:nx, 1:nz)       ! Passing boundary conditions in forcing array
        p(1:nx, ny, 1:nz) = bcs_ht(1:nx, 1:nz)

        if (fft_z_on) then
            call OPR_Fourier_X_Forward(nx, ny, nz, p, c_tmp2)
            call OPR_Fourier_Z_Forward(c_tmp2, c_tmp1) ! tmp2 might be overwritten; cannot use wrk3d
        else
            call OPR_Fourier_X_Forward(nx, ny, nz, p, c_tmp1)
        end if

        tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
        ! ###################################################################
        ! Make x direction last one and leave y direction first
        call TLab_Transpose_COMPLEX(c_tmp1, isize_line, ny*nz, isize_line, c_tmp2, ny*nz)

#define f(j,k,i) tmp2(j,k,i)
#define u(j,k,i) p_wrk3d(j,k,i)

        ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
        do i = 1, i_max
            do k = 1, nz
                u(1:2, k, i) = f(1:2, k, i)                         ! bottom boundary conditions
                u(2*ny - 1:2*ny, k, i) = f(2*ny - 1:2*ny, k, i)     ! top boundary conditions

                select case (ibc)
                case (BCS_NN)           ! use precalculated LU factorization
                    ! Compatibility constraint for singular modes. The reference value of p at bottom is set to zero
                    if (any(i_sing == i) .and. any(k_sing == k)) u(1:2, k, i) = 0.0_wp

                    call FDM_Int2_Solve(2, fdm_int2(k, i), rhs_d, f(:, k, i), u(:, k, i), wrk2d)

                case default            ! Need to calculate and factorize LHS
                    call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2, lambda(k, i), ibc, fdm_int2_loc)
                    call FDM_Int2_Solve(2, fdm_int2_loc, fdm_int2_loc%rhs, f(:, k, i), u(:, k, i), wrk2d)

                end select

            end do
        end do

        call TLab_Transpose_COMPLEX(c_wrk3d, ny*nz, isize_line, ny*nz, c_tmp1, isize_line)

        ! ###################################################################
        ! Fourier field p (based on array tmp1)
        ! ###################################################################
        if (fft_z_on) then
            call OPR_Fourier_Z_Backward(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
            call OPR_Fourier_X_Backward(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
        else
            call OPR_Fourier_X_Backward(nx, ny, nz, c_tmp1, p)    ! tmp1 might be overwritten
        end if

        if (present(dpdy)) then
            call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs_p, g(2), p, dpdy)
        end if

        nullify (c_tmp1, c_tmp2, p_wrk3d)
#undef f
#undef u

        return
    end subroutine OPR_Poisson_FourierXZ_Direct

    !########################################################################
    !#
    !# Solve Lap a + \alpha a = f using Fourier in xOz planes, to rewrite the problem as
    !#
    !#      \hat{a}''-(\lambda-\alpha) \hat{a} = \hat{f}
    !#
    !# where \lambda = kx^2+kz^2
    !#
    !########################################################################
    subroutine OPR_Helmholtz_FourierXZ_Factorize(nx, ny, nz, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                        ! Forcing term, and solution field
        real(wp), intent(inout) :: tmp1(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(inout) :: tmp2(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        real(wp) bcs(2, 2)

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field])
        p_wrk3d(1:2*ny, 1:nz, 1:nx/2 + 1) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        a(1:nx, 1, 1:nz) = bcs_hb(1:nx, 1:nz)       ! Passing boundary conditions in forcing array
        a(1:nx, ny, 1:nz) = bcs_ht(1:nx, 1:nz)

        if (fft_z_on) then
            call OPR_Fourier_X_Forward(nx, ny, nz, a, c_tmp2)
            call OPR_Fourier_Z_Forward(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_Fourier_X_Forward(nx, ny, nz, a, c_tmp1)
        end if

        tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE (\hat{p}')'-(\lambda+\alpha) \hat{p} = \hat{f}
        ! ###################################################################
        ! Make x direction last one and leave y direction first
        call TLab_Transpose_COMPLEX(c_tmp1, isize_line, ny*nz, isize_line, c_tmp2, ny*nz)

#define f(j,k,i) tmp2(j,k,i)
#define v(j,k,i) tmp1(j,k,i)
#define u(j,k,i) p_wrk3d(j,k,i)

        ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
        do i = 1, i_max
            do k = 1, nz
                bcs(1:2, 1) = f(1:2, k, i)                  ! bottom boundary conditions
                bcs(1:2, 2) = f(2*ny - 1:2*ny, k, i)        ! top boundary conditions

                call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%der1, &
                                         sqrt(lambda(k, i) - alpha), BCS_MIN, fdm_int1_loc(BCS_MIN))

                call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%der1, &
                                         -sqrt(lambda(k, i) - alpha), BCS_MAX, fdm_int1_loc(BCS_MAX))

                select case (ibc)
                case (BCS_NN)
                    call OPR_ODE2_Factorize_NN(2, fdm_int1_loc, fdm_int1_loc(BCS_MIN)%rhs, fdm_int1_loc(BCS_MAX)%rhs, &
                                               u(:, k, i), f(:, k, i), bcs, v(:, k, i), wrk1d, wrk2d)

                case (BCS_DD)
                    call OPR_ODE2_Factorize_DD(2, fdm_int1_loc, fdm_int1_loc(BCS_MIN)%rhs, fdm_int1_loc(BCS_MAX)%rhs, &
                                               u(:, k, i), f(:, k, i), bcs, v(:, k, i), wrk1d, wrk2d)
                end select

            end do
        end do

        call TLab_Transpose_COMPLEX(c_wrk3d, ny*nz, isize_line, ny*nz, c_tmp1, isize_line)

        ! ###################################################################
        ! Fourier field a (based on array tmp1)
        ! ###################################################################
        if (fft_z_on) then
            call OPR_Fourier_Z_Backward(c_tmp1, c_wrk3d) ! tmp1 might be overwritten
            call OPR_Fourier_X_Backward(nx, ny, nz, c_wrk3d, a)
        else
            call OPR_Fourier_X_Backward(nx, ny, nz, c_tmp1, a)
        end if

        nullify (c_tmp1, c_tmp2, p_wrk3d)
#undef f
#undef v
#undef u

        return
    end subroutine OPR_Helmholtz_FourierXZ_Factorize

    !########################################################################
    !########################################################################
    subroutine OPR_Helmholtz_FourierXZ_Direct(nx, ny, nz, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                        ! Forcing term, and solution field
        real(wp), intent(inout) :: tmp1(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(inout) :: tmp2(2*ny, nz, nx/2 + 1)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields

        target tmp1, tmp2

        ! -----------------------------------------------------------------------

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_field])
        p_wrk3d(1:2*ny, 1:nz, 1:nx/2 + 1) => wrk3d(1:isize_txc_field)

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        a(1:nx, 1, 1:nz) = bcs_hb(1:nx, 1:nz)       ! Passing boundary conditions in forcing array
        a(1:nx, ny, 1:nz) = bcs_ht(1:nx, 1:nz)

        if (fft_z_on) then
            call OPR_Fourier_X_Forward(nx, ny, nz, a, c_tmp2)
            call OPR_Fourier_Z_Forward(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_Fourier_X_Forward(nx, ny, nz, a, c_tmp1)
        end if

        tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE \hat{p}''-(\lambda+\alpha) \hat{p} = \hat{f}
        ! ###################################################################
        ! Make x direction last one and leave y direction first
        call TLab_Transpose_COMPLEX(c_tmp1, isize_line, ny*nz, isize_line, c_tmp2, ny*nz)

#define f(j,k,i) tmp2(j,k,i)
#define u(j,k,i) p_wrk3d(j,k,i)

        ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
        do i = 1, i_max
            do k = 1, nz
                u(1:2, k, i) = f(1:2, k, i)
                u(2*ny - 1:2*ny, k, i) = f(2*ny - 1:2*ny, k, i)

                call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2, lambda(k, i) - alpha, ibc, fdm_int2_loc)
                call FDM_Int2_Solve(2, fdm_int2_loc, fdm_int2_loc%rhs, f(:, k, i), u(:, k, i), wrk2d)

            end do
        end do

        call TLab_Transpose_COMPLEX(c_wrk3d, ny*nz, isize_line, ny*nz, c_tmp1, isize_line)

        ! ###################################################################
        ! Fourier field a (based on array tmp1)
        ! ###################################################################
        if (fft_z_on) then
            call OPR_Fourier_Z_Backward(c_tmp1, c_wrk3d) ! tmp1 might be overwritten
            call OPR_Fourier_X_Backward(nx, ny, nz, c_wrk3d, a)
        else
            call OPR_Fourier_X_Backward(nx, ny, nz, c_tmp1, a)
        end if

        nullify (c_tmp1, c_tmp2, p_wrk3d)
#undef f
#undef u

        return
    end subroutine OPR_Helmholtz_FourierXZ_Direct

end module OPR_Elliptic
