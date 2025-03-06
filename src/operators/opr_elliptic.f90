#include "dns_const.h"
#include "dns_error.h"
! You need to split the routines into the ones that are initialized and the ones that not.
module OPR_ELLIPTIC
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: efile
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use TLab_Memory, only: TLab_Allocate_Real
    use TLab_Memory, only: isize_txc_dimz, imax, jmax, kmax
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLab_Arrays, only: wrk1d, wrk2d, wrk3d
    use TLab_Pointers_3D, only: p_wrk1d, p_wrk2d
    use TLab_Pointers_C, only: c_wrk1d, c_wrk3d
    use FDM, only: fdm_dt, FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA, FDM_COM4_DIRECT, FDM_COM6_DIRECT
    use FDM_Integral
    use FDM_MatMul
    use OPR_FOURIER
    use OPR_ODES
    use OPR_PARTIAL
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    procedure(OPR_Poisson_interface) :: OPR_Poisson_dt ! Implicit pointer (Procedure type)
    abstract interface
        subroutine OPR_Poisson_interface(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
            use TLab_Constants, only: wi, wp
            use FDM, only: fdm_dt
            use TLab_Memory, only: isize_txc_dimz
            integer(wi), intent(in) :: nx, ny, nz
            integer, intent(in) :: ibc                                      ! Dirichlet/Neumman BCs at jmin/jmax: BCS_DD, BCS_ND, BCS_DN, BCS_NN
            type(fdm_dt), intent(in) :: g(3)
            real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
            real(wp), intent(inout), target :: tmp1(isize_txc_dimz, nz)     ! FFT of forcing term
            real(wp), intent(inout), target :: tmp2(isize_txc_dimz, nz)     ! Aux array for FFT
            real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
            real(wp), intent(out), optional :: dpdy(nx, ny, nz)             ! Vertical derivative of solution
        end subroutine
    end interface

    procedure(OPR_Helmholtz_interface) :: OPR_Helmholtz_dt ! Implicit pointer (Procedure type)
    abstract interface
        subroutine OPR_Helmholtz_interface(nx, ny, nz, g, ibc, alpha, p, tmp1, tmp2, bcs_hb, bcs_ht)
            use TLab_Constants, only: wi, wp
            use FDM, only: fdm_dt
            use TLab_Memory, only: isize_txc_dimz
            integer(wi), intent(in) :: nx, ny, nz
            integer, intent(in) :: ibc                                      ! Dirichlet/Neumman BCs at jmin/jmax: BCS_DD, BCS_ND, BCS_DN, BCS_NN
            type(fdm_dt), intent(in) :: g(3)
            real(wp), intent(in) :: alpha
            real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
            real(wp), intent(inout), target :: tmp1(isize_txc_dimz, nz)     ! FFT of forcing term
            real(wp), intent(inout), target :: tmp2(isize_txc_dimz, nz)     ! Aux array for FFT
            real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
        end subroutine
    end interface

    procedure(OPR_Poisson_dt), pointer :: OPR_Poisson
    procedure(OPR_Helmholtz_dt), pointer :: OPR_Helmholtz

    public :: OPR_Elliptic_Initialize
    public :: OPR_Poisson
    public :: OPR_Helmholtz

    ! -----------------------------------------------------------------------
    complex(wp), target :: bcs(3)
    real(wp), pointer :: r_bcs(:) => null()
    complex(wp), pointer :: c_tmp1(:, :) => null(), c_tmp2(:, :) => null()
    integer(wi) i, j, k, iglobal, kglobal, ip, isize_line
    real(wp) norm
    integer(wi) i_sing(2), k_sing(2)    ! singular global modes
    type(fdm_dt) fdm_loc
    type(fdm_integral_dt), allocatable :: fdm_int_loc(:, :, :)      ! factorized method
    real(wp), allocatable, target :: rhs_b(:, :), rhs_t(:, :)
    real(wp), allocatable, target :: lu_d(:, :, :, :)               ! direct method
    type(fdm_integral_dt) :: fdm_int_helmholtz(2)
    real(wp), allocatable :: lambda(:, :)

#define p_a(icpp,jcpp,kcpp)   lu_d(icpp,1,jcpp,kcpp)
#define p_b(icpp,jcpp,kcpp)   lu_d(icpp,2,jcpp,kcpp)
#define p_c(icpp,jcpp,kcpp)   lu_d(icpp,3,jcpp,kcpp)
#define p_d(icpp,jcpp,kcpp)   lu_d(icpp,4,jcpp,kcpp)
#define p_e(icpp,jcpp,kcpp)   lu_d(icpp,5,jcpp,kcpp)
#define p_f1(icpp,jcpp,kcpp)  lu_d(icpp,6,jcpp,kcpp)
#define p_f2(icpp,jcpp,kcpp)  lu_d(icpp,7,jcpp,kcpp)
#define p_rhs1(icpp,jcpp,kcpp)  lu_d(icpp,8,jcpp,kcpp)
#define p_rhs2(icpp,jcpp,kcpp)  lu_d(icpp,9,jcpp,kcpp)

contains
! #######################################################################
! #######################################################################
! We precalculate the LU factorization for the case BCS_NN, which is the one used in the pressure-Poisson equation
    subroutine OPR_Elliptic_Initialize(inifile)
        use FDM, only: g, FDM_Initialize

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer imode_elliptic           ! finite-difference method for pressure-Poisson and Helmholtz equations
        integer ibc_loc
        integer(wi) :: ndl, ndr, nd
        integer, parameter :: i1 = 1
        character*512 sRes
        character*32 bakfile

! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        call ScanFile_Char(bakfile, inifile, 'Main', 'EllipticOrder', 'compactjacobian6', sRes)
        if (trim(adjustl(sRes)) == 'compactjacobian4') then; imode_elliptic = FDM_COM4_JACOBIAN
        else if (trim(adjustl(sRes)) == 'compactjacobian6') then; imode_elliptic = FDM_COM6_JACOBIAN
        else if (trim(adjustl(sRes)) == 'compactjacobian6penta') then; imode_elliptic = FDM_COM6_JACOBIAN_PENTA
        else if (trim(adjustl(sRes)) == 'compactdirect4') then; imode_elliptic = FDM_COM4_DIRECT
        else if (trim(adjustl(sRes)) == 'compactdirect6') then; imode_elliptic = FDM_COM6_DIRECT
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong Main.EllipticOrder option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! It should probably be better to say factorize or not factorize and use the global fdm plan
        fdm_loc%mode_fdm1 = imode_elliptic
        fdm_loc%mode_fdm2 = imode_elliptic
        call FDM_Initialize(g(2)%nodes, fdm_loc)

        isize_line = imax/2 + 1

        allocate (lambda(isize_line, kmax))
        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        select case (imode_elliptic)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA)
            OPR_Poisson => OPR_Poisson_FourierXZ_Factorize
            ! OPR_Poisson => OPR_Poisson_FourierXZ_Factorize_Old
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Factorize

            ndl = fdm_loc%nb_diag_1(1)
            ndr = fdm_loc%nb_diag_1(2)
            nd = ndl
            allocate (fdm_int_loc(2, isize_line, kmax))
            call TLab_Allocate_Real(__FILE__, rhs_b, [g(2)%size, nd, isize_line, kmax], 'rhs_b')
            call TLab_Allocate_Real(__FILE__, rhs_t, [g(2)%size, nd, isize_line, kmax], 'rhs_t')

            if (.not. stagger_on) then
                i_sing = [1, g(1)%size/2 + 1]
                k_sing = [1, g(3)%size/2 + 1]
            else                    ! In case of staggering only one singular mode + different modified wavenumbers
                i_sing = [1, 1]
                k_sing = [1, 1]
            end if

        case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
            OPR_Poisson => OPR_Poisson_FourierXZ_Direct
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Direct

            ndl = fdm_loc%nb_diag_2(1)
            ndr = fdm_loc%nb_diag_2(2)
            nd = ndr + ndl - 1 + 2     ! The rhs diagonal is 1 and not need to store; we add 2 independent terms
            call TLab_Allocate_Real(__FILE__, lu_d, [g(2)%size, nd, isize_line, kmax], 'lu_d')

        end select

        do k = 1, kmax
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif

            do i = 1, isize_line
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif

                select case (imode_elliptic)
                case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA)
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda(i, k) = g(1)%mwn1(iglobal)**2.0 + g(3)%mwn1(kglobal)**2.0
                    else
                        lambda(i, k) = g(1)%mwn1(iglobal)**2.0
                    end if

                    fdm_int_loc(BCS_MIN, i, k)%bc = BCS_MIN
                    fdm_int_loc(BCS_MIN, i, k)%mode_fdm1 = fdm_loc%mode_fdm1
                    call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                             sqrt(lambda(i, k)), fdm_int_loc(BCS_MIN, i, k))

                    rhs_b(:, :) = fdm_int_loc(BCS_MIN, i, k)%rhs(:, :)          ! free memory that is independent of lambda
                    if (allocated(fdm_int_loc(BCS_MIN, i, k)%rhs)) deallocate (fdm_int_loc(BCS_MIN, i, k)%rhs)

                    fdm_int_loc(BCS_MAX, i, k)%bc = BCS_MAX
                    fdm_int_loc(BCS_MAX, i, k)%mode_fdm1 = fdm_loc%mode_fdm1
                    call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                             -sqrt(lambda(i, k)), fdm_int_loc(BCS_MAX, i, k))

                    rhs_t(:, :) = fdm_int_loc(BCS_MAX, i, k)%rhs(:, :)          ! free memory that is independent of lambda
                    if (allocated(fdm_int_loc(BCS_MAX, i, k)%rhs)) deallocate (fdm_int_loc(BCS_MAX, i, k)%rhs)

                case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)     ! only for case BCS_NN
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda(i, k) = g(1)%mwn2(iglobal) + g(3)%mwn2(kglobal)
                    else
                        lambda(i, k) = g(1)%mwn2(iglobal)
                    end if

                    ! Compatibility constraint. The reference value of p at the lower boundary is set to zero
                    if (iglobal == 1 .and. kglobal == 1) then
                        ibc_loc = BCS_DN
                    else
                        ibc_loc = BCS_NN
                    end if

                    ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                    call FDM_Int2_CreateSystem(g(2)%size, g(2)%nodes, ibc_loc, fdm_loc%lhs2, fdm_loc%rhs2, lambda(i, k), &
                                               lu_d(:, 1:5, i, k), lu_d(:, 6:7, i, k), lu_d(:, 8:9, i, k))

                    ! LU decomposizion
                    ! We rely on this routines not changing a(2:3), b(2), e(ny-2:ny-1), d(ny-1)
                    call PENTADFS(g(2)%size - 2, p_a(2, i, k), p_b(2, i, k), p_c(2, i, k), p_d(2, i, k), p_e(2, i, k))

                    ! Particular solutions
                    call PENTADSS(g(2)%size - 2, i1, p_a(2, i, k), p_b(2, i, k), p_c(2, i, k), p_d(2, i, k), p_e(2, i, k), p_f1(2, i, k))
                    call PENTADSS(g(2)%size - 2, i1, p_a(2, i, k), p_b(2, i, k), p_c(2, i, k), p_d(2, i, k), p_e(2, i, k), p_f2(2, i, k))

                end select

            end do
        end do

        return
    end subroutine OPR_Elliptic_Initialize

!########################################################################
!#
!# Solve Lap p = f using Fourier in xOz planes, to rewritte the problem as
!#     \hat{p}''-\lambda \hat{p} = \hat{f}
!#
!# where \lambda = kx^2+kz^2
!#
!# The reference value of p at the lower boundary is set to zero
!#
!# The global variable isize_txc_field defines the size of array txc equal
!# to (imax+2)*(jmax+2)*kmax, or larger if PARALLEL mode
!#
!########################################################################
    subroutine OPR_Poisson_FourierXZ_Factorize(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        type(fdm_dt), intent(in) :: g(3)
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
        real(wp), intent(out), optional :: dpdy(nx, ny, nz)             ! Vertical derivative of solution

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        real(wp), pointer :: u(:, :) => null(), v(:, :) => null(), f(:, :) => null()

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp1)
        end if

        tmp1 = tmp1*norm

        ! ###################################################################
        ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
        ! ###################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif

            ! Make x direction last one and leave y direction first
            call TLab_Transpose_COMPLEX(c_tmp1(:, k), isize_line, ny + 2, isize_line, c_tmp2(:, k), ny + 2)

            f(1:2*(ny + 2), 1:isize_line) => tmp2(1:2*(ny + 2)*isize_line, k)
            v(1:2*ny, 1:isize_line) => tmp1(1:2*ny*isize_line, k)
            u(1:2*ny, 1:isize_line) => wrk3d(1:2*ny*isize_line)

            do i = 1, isize_line
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                select case (ibc)
                case (BCS_NN) ! Neumann   & Neumann   BCs
                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                        call OPR_ODE2_SINGULAR_NN(2, fdm_Int0, u(:, i), f(:, i), f(2*ny + 1:, i), v(:, i), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_NN(2, fdm_int_loc(:, i, k), rhs_b, rhs_t, &
                                         u(:, i), f(:, i), f(2*ny + 1:, i), v(:, i), wrk1d, wrk2d)
                    end if

                case (BCS_DD) ! Dirichlet & Dirichlet BCs
                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                        call OPR_ODE2_SINGULAR_DD(2, fdm_Int0, u(:, i), f(:, i), f(2*ny + 1:, i), v(:, i), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_DD(2, fdm_int_loc(:, i, k), rhs_b, rhs_t, &
                                         u(:, i), f(:, i), f(2*ny + 1:, i), v(:, i), wrk1d, wrk2d)
                    end if

                end select

            end do

            if (present(dpdy)) call TLab_Transpose_COMPLEX(c_tmp1(:, k), ny, isize_line, ny, c_tmp2(:, k), isize_line)
            call TLab_Transpose_COMPLEX(c_wrk3d, ny, isize_line, ny, c_tmp1(:, k), isize_line)

        end do

        ! ###################################################################
        ! Fourier field p (based on array tmp1)
        ! ###################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, p)    ! tmp2 might be overwritten
        end if

        ! Fourier derivatives (based on array tmp2)
        if (present(dpdy)) then
            if (g(3)%size > 1) then
                call OPR_FOURIER_B_Z_EXEC(c_tmp2, c_wrk3d)              ! tmp2 might be overwritten
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, dpdy)    ! wrk3d might be overwritten
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp2, dpdy)     ! tmp2 might be overwritten
            end if
        end if

        nullify (u, v, f)
        nullify (c_tmp1, c_tmp2)

        return
    end subroutine OPR_Poisson_FourierXZ_Factorize

!########################################################################
!########################################################################
    subroutine OPR_Poisson_FourierXZ_Direct(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc
        type(fdm_dt), intent(in) :: g(3)
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
        real(wp), intent(out), optional :: dpdy(nx, ny, nz)             ! Vertical derivative of solution

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        integer(wi) ibc_loc, bcs_p(2, 2)
        integer, parameter :: i1 = 1, i2 = 2

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp1)
        end if

        ! ###################################################################
        ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
        ! ###################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif

            do i = 1, isize_line
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif

                ! forcing term in c_wrk1d(:,5), i.e. p_wrk1d(:,9), solution will be in c_wrk1d(:,6), i.e., p_wrk1d(:,11)
                do j = 1, ny
                    ip = (j - 1)*isize_line + i; c_wrk1d(j, 5) = c_tmp1(ip, k)
                end do

                ! BCs
                j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
                j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

                ! Compatibility constraint for singular modes. 2nd order FDMs are non-zero at Nyquist
                ! The reference value of p at the lower boundary is set to zero
                if (iglobal == 1 .and. kglobal == 1 .and. ibc == BCS_NN) then
                    ibc_loc = BCS_DN
                    bcs(1) = (0.0_wp, 0.0_wp)
                else
                    ibc_loc = ibc
                end if

                ! -----------------------------------------------------------------------
                if (ibc /= BCS_NN) then     ! Need to calculate and factorize LHS
                    ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                    call FDM_Int2_CreateSystem(ny, g(2)%nodes, ibc_loc, fdm_loc%lhs2, fdm_loc%rhs2, lambda(i, k), &
                                               p_wrk1d(:, 1:5), p_wrk1d(:, 6:7), p_wrk1d(:, 13:14))

                    ! LU factorization
                    call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

                    ! Particular solutions
                    call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
                    call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

                    ! Construct rhs
                    p_wrk1d(1:2, 11) = 0.0_wp       ! This element is simply the solution at imin of p(0)
                    p_wrk1d(ny - 1:ny, 12) = 0.0_wp ! This element is simply the solution at imax of p(0)
                    call MatMul_3d(ny - 2, 2, p_wrk1d(2:, 13), p_wrk1d(2:, 14), p_wrk1d(3:, 9), p_wrk1d(3:, 11))

                    ! Solve pentadiagonal linear system
                    call PENTADSS(ny - 2, i2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(3, 11))

                    c_wrk1d(:, 6) = (c_wrk1d(:, 6) + bcs(1)*p_wrk1d(:, 6) + bcs(2)*p_wrk1d(:, 7))*norm

                    !   Corrections to the BCS_DD to account for Neumann
                    if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                        c_wrk1d(1, 6) = c_wrk1d(1, 6) + p_wrk1d(1, 3)*c_wrk1d(2, 6) &
                                        + p_wrk1d(1, 4)*c_wrk1d(3, 6) + p_wrk1d(1, 5)*c_wrk1d(4, 6) &
                                        + p_wrk1d(1, 2)*c_wrk1d(2, 5)*norm
                    end if

                    if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                        c_wrk1d(ny, 6) = c_wrk1d(ny, 6) + p_wrk1d(ny, 3)*c_wrk1d(ny - 1, 6) &
                                         + p_wrk1d(ny, 2)*c_wrk1d(ny - 2, 6) + p_wrk1d(ny, 1)*c_wrk1d(ny - 3, 6) &
                                         + p_wrk1d(ny, 5)*c_wrk1d(ny - 1, 5)*norm
                    end if

                else                        ! use precalculated LU factorization
                    ! Construct rhs
                    p_wrk1d(1:2, 11) = 0.0_wp       ! This element is simply the solution at imin of p(0)
                    p_wrk1d(ny - 1:ny, 12) = 0.0_wp ! This element is simply the solution at imax of p(0)
                    call MatMul_3d(ny - 2, 2, p_rhs1(2, i, k), p_rhs2(2, i, k), p_wrk1d(3:, 9), p_wrk1d(3:, 11))

                    ! Solve pentadiagonal linear system
                    call PENTADSS(ny - 2, i2, p_a(2, i, k), p_b(2, i, k), p_c(2, i, k), p_d(2, i, k), p_e(2, i, k), p_wrk1d(3, 11))

                    c_wrk1d(:, 6) = (c_wrk1d(:, 6) + bcs(1)*p_f1(:, i, k) + bcs(2)*p_f2(:, i, k))*norm

                    !   Corrections to the BCS_DD to account for Neumann
                    if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                        c_wrk1d(1, 6) = c_wrk1d(1, 6) + p_c(1, i, k)*c_wrk1d(2, 6) &
                                        + p_d(1, i, k)*c_wrk1d(3, 6) + p_e(1, i, k)*c_wrk1d(4, 6) &
                                        + p_b(1, i, k)*c_wrk1d(2, 5)*norm
                    end if

                    if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                        c_wrk1d(ny, 6) = c_wrk1d(ny, 6) + p_c(ny, i, k)*c_wrk1d(ny - 1, 6) &
                                         + p_b(ny, i, k)*c_wrk1d(ny - 2, 6) + p_a(ny, i, k)*c_wrk1d(ny - 3, 6) &
                                         + p_e(ny, i, k)*c_wrk1d(ny - 1, 5)*norm
                    end if

                end if

                ! Rearrange in memory and normalize
                do j = 1, ny
                    ip = (j - 1)*isize_line + i
                    c_tmp1(ip, k) = c_wrk1d(j, 6)
                end do

            end do
        end do

        ! ###################################################################
        ! Fourier field p (based on array tmp1)
        ! ###################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, p)    ! tmp2 might be overwritten
        end if

        ! Fourier derivatives (based on array tmp2)
        if (present(dpdy)) then
            bcs_p = 0
            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs_p, g(2), p, dpdy)
        end if

        nullify (c_tmp1, c_tmp2)

        return
    end subroutine OPR_Poisson_FourierXZ_Direct

!     subroutine OPR_Poisson_FourierXZ_Factorize_Old(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
!         integer(wi), intent(in) :: nx, ny, nz
!         integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
!         !                                                   1, for Neumann   & Dirichlet
!         !                                                   2, for Dirichlet & Neumann
!         !                                                   3, for Neumann   & Neumann
!         type(fdm_dt), intent(in) :: g(3)
!         real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
!         real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
!         real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
!         real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields
!         real(wp), intent(out), optional :: dpdy(nx, ny, nz)           ! Vertical derivative of solution

!         target tmp1, tmp2

!         ! -----------------------------------------------------------------------
!         integer(wi) i_sing(2), k_sing(2)    ! singular global modes

!         ! #######################################################################
!         call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
!         call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
!         call c_f_pointer(c_loc(bcs), r_bcs, shape=[3*2])

!         norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

!         isize_line = nx/2 + 1

!         if (.not. stagger_on) then
!             i_sing = [1, g(1)%size/2 + 1]
!             k_sing = [1, g(3)%size/2 + 1]
!         else                    ! In case of staggering only one singular mode + different modified wavenumbers
!             i_sing = [1, 1]
!             k_sing = [1, 1]
!         end if

!         ! #######################################################################
!         ! Fourier transform of forcing term; output of this section in array tmp1
!         ! #######################################################################
!         if (g(3)%size > 1) then
!             call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2)
!             call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
!         else
!             call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp1)
!         end if

!         ! ###################################################################
!         ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
!         ! ###################################################################
!         do k = 1, nz
! #ifdef USE_MPI
!             kglobal = k + ims_offset_k
! #else
!             kglobal = k
! #endif

!             do i = 1, isize_line
! #ifdef USE_MPI
!                 iglobal = i + ims_offset_i/2
! #else
!                 iglobal = i
! #endif

!                 ! forcing term
!                 do j = 1, ny
!                     ip = (j - 1)*isize_line + i; c_wrk1d(j, 1) = c_tmp1(ip, k)
!                 end do

!                 ! BCs
!                 j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
!                 j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

!                 ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
!                 select case (ibc)
!                 case (BCS_NN) ! Neumann   & Neumann   BCs
!                     if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
!                         ! call OPR_ODE2_1_SINGULAR_NN_OLD(g(2)%mode_fdm1, ny, 2, &
!                         !                             g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_SINGULAR_NN(2, fdm_Int0, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), p_wrk2d)
!                     else
!                         ! call OPR_ODE2_1_REGULAR_NN_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k), &
!                         !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_NN(2, fdm_int_loc(:, i, k), rhs_b, rhs_t, &
!                                          p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), p_wrk2d)
!                     end if

!                 case (BCS_DD) ! Dirichlet & Dirichlet BCs
!                     if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
!                         ! call OPR_ODE2_1_SINGULAR_DD_OLD(g(2)%mode_fdm1, ny, 2, &
!                         !                                 g(2)%nodes, g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_SINGULAR_DD(2, fdm_Int0, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), p_wrk2d)
!                     else
!                         ! call OPR_ODE2_1_REGULAR_DD_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k), &
!                         !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_DD(2, fdm_int_loc(:, i, k), rhs_b, rhs_t, &
!                                          p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), p_wrk2d)
!                     end if

!                 end select

!                 ! Rearrange in memory and normalize
!                 do j = 1, ny
!                     ip = (j - 1)*isize_line + i
!                     c_tmp1(ip, k) = c_wrk1d(j, 2)*norm ! solution
!                     c_tmp2(ip, k) = c_wrk1d(j, 3)*norm ! Oy derivative
!                 end do

!             end do
!         end do

    !     ! ###################################################################
    !     ! Fourier field p (based on array tmp1)
    !     ! ###################################################################
    !     if (g(3)%size > 1) then
    !         call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
    !         call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
    !     else
    !         call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, p)    ! tmp2 might be overwritten
    !     end if

    !     ! Fourier derivatives (based on array tmp2)
    !     if (present(dpdy)) then
    !         if (g(3)%size > 1) then
    !             call OPR_FOURIER_B_Z_EXEC(c_tmp2, c_wrk3d)              ! tmp2 might be overwritten
    !             call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, dpdy)    ! wrk3d might be overwritten
    !         else
    !             call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp2, dpdy)     ! tmp2 might be overwritten
    !         end if
    !     end if

    !     nullify (c_tmp1, c_tmp2, r_bcs)

    !     return
    ! end subroutine OPR_Poisson_FourierXZ_Factorize_Old

!########################################################################
!#
!# Solve Lap a + lpha a = f using Fourier in xOz planes, to rewritte
!# the problem as
!#
!#      \hat{a}''-(\lambda-lpha) \hat{a} = \hat{f}
!#
!# where \lambda = kx^2+kz^2
!#
!# The global variable isize_txc_field defines the size of array txc equal
!# to (imax+2)*(jmax+2)*kmax, or larger if PARALLEL mode
!#
!########################################################################
    subroutine OPR_Helmholtz_FourierXZ_Factorize(nx, ny, nz, g, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(fdm_dt), intent(in) :: g(3)
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                       ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        integer(wi) ndl, ndr

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(bcs), r_bcs, shape=[3*2])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp2)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp1)
        end if

        ! ###################################################################
        ! Solve FDE (\hat{p}')'-(\lambda+lpha) \hat{p} = \hat{f}
        ! ###################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif

            do i = 1, isize_line
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif

                ! forcing term
                do j = 1, ny
                    ip = (j - 1)*isize_line + i; c_wrk1d(j, 1) = c_tmp1(ip, k)
                end do

                ! BCs
                j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
                j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                ndl = fdm_loc%nb_diag_1(1)
                ndr = fdm_loc%nb_diag_2(2)

                fdm_int_helmholtz(BCS_MIN)%bc = BCS_MIN
                fdm_int_helmholtz(BCS_MIN)%mode_fdm1 = fdm_loc%mode_fdm1
                call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                         sqrt(lambda(i, k) - alpha), fdm_int_helmholtz(BCS_MIN))

                fdm_int_helmholtz(BCS_MAX)%bc = BCS_MAX
                fdm_int_helmholtz(BCS_MAX)%mode_fdm1 = fdm_loc%mode_fdm1
                call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                         -sqrt(lambda(i, k) - alpha), fdm_int_helmholtz(BCS_MAX))

                select case (ibc)
                case (3) ! Neumann   & Neumann   BCs
                    ! call OPR_ODE2_1_REGULAR_NN_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k)-alpha, &
                    !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    call OPR_ODE2_NN(2, fdm_int_helmholtz, fdm_int_helmholtz(BCS_MIN)%rhs, fdm_int_helmholtz(BCS_MAX)%rhs, &
                                     p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), p_wrk2d)

                case (0) ! Dirichlet & Dirichlet BCs
                    ! call OPR_ODE2_1_REGULAR_DD_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k)-alpha, &
                    !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))

                    call OPR_ODE2_DD(2, fdm_int_helmholtz, fdm_int_helmholtz(BCS_MIN)%rhs, fdm_int_helmholtz(BCS_MAX)%rhs, &
                                     p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), p_wrk2d)
                end select

                ! Rearrange in memory and normalize
                do j = 1, ny
                    ip = (j - 1)*isize_line + i
                    c_tmp1(ip, k) = c_wrk1d(j, 2)*norm ! solution
                end do

            end do
        end do

        ! ###################################################################
        ! Fourier field a (based on array tmp1)
        ! ###################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d) ! tmp1 might be overwritten
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, a)
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, a)
        end if

        nullify (c_tmp1, c_tmp2, r_bcs)

        return
    end subroutine OPR_Helmholtz_FourierXZ_Factorize

!########################################################################
!########################################################################
! Same, but using the direct mode of FDM
! Opposite to previous routine, here we use the first 8 wrk1d arrays for the diagonals of the LHS,
! and the last ones for the forcing and solution. The reason is the routine after this one.
    subroutine OPR_Helmholtz_FourierXZ_Direct(nx, ny, nz, g, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(fdm_dt), intent(in) :: g(3)
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                       ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields

        target tmp1, tmp2

        ! -----------------------------------------------------------------------

        integer, parameter :: i1 = 1, i2 = 2

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp2)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp1)
        end if

        ! ###################################################################
        ! Solve FDE \hat{p}''-(\lambda+lpha) \hat{p} = \hat{f}
        ! ###################################################################
        do k = 1, nz
#ifdef USE_MPI
            kglobal = k + ims_offset_k
#else
            kglobal = k
#endif

            do i = 1, isize_line
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif
                ! forcing term in c_wrk1d(:,5), i.e. p_wrk1d(:,9), solution will be in c_wrk1d(:,6), i.e., p_wrk1d(:,11)
                do j = 1, ny
                    ip = (j - 1)*isize_line + i; c_wrk1d(j, 5) = c_tmp1(ip, k)
                end do

                ! BCs
                j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
                j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                p_wrk1d(:, 1:7) = 0.0_wp
                call FDM_Int2_CreateSystem(ny, g(2)%nodes, ibc, fdm_loc%lhs2, fdm_loc%rhs2, lambda(i, k) - alpha, &
                                           p_wrk1d(:, 1:5), p_wrk1d(:, 6:7), p_wrk1d(:, 13:14))

                ! LU factorization
                call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

                ! Particular solutions
                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

                ! Construct rhs
                p_wrk1d(1:2, 11) = 0.0_wp       ! This element is simply the solution at imin of p(0)
                p_wrk1d(ny - 1:ny, 12) = 0.0_wp ! This element is simply the solution at imax of p(0)
                call MatMul_3d(ny - 2, 2, p_wrk1d(2:, 13), p_wrk1d(2:, 14), p_wrk1d(3:, 9), p_wrk1d(3:, 11))

                ! Solve pentadiagonal linear system
                call PENTADSS(ny - 2, i2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(3, 11))

                c_wrk1d(:, 6) = c_wrk1d(:, 6) + bcs(1)*p_wrk1d(:, 6) + bcs(2)*p_wrk1d(:, 7)

                !   Corrections to the BCS_DD to account for Neumann
                if (any([BCS_ND, BCS_NN] == ibc)) then
                    c_wrk1d(1, 6) = c_wrk1d(1, 6) + p_wrk1d(1, 3)*c_wrk1d(2, 6) &
                                    + p_wrk1d(1, 4)*c_wrk1d(3, 6) + p_wrk1d(1, 5)*c_wrk1d(4, 6) &
                                    + p_wrk1d(1, 2)*c_wrk1d(2, 5)
                end if

                if (any([BCS_DN, BCS_NN] == ibc)) then
                    c_wrk1d(ny, 6) = c_wrk1d(ny, 6) + p_wrk1d(ny, 3)*c_wrk1d(ny - 1, 6) &
                                     + p_wrk1d(ny, 2)*c_wrk1d(ny - 2, 6) + p_wrk1d(ny, 1)*c_wrk1d(ny - 3, 6) &
                                     + p_wrk1d(ny, 5)*c_wrk1d(ny - 1, 5)
                end if

                ! Rearrange in memory and normalize
                do j = 1, ny
                    ip = (j - 1)*isize_line + i
                    c_tmp1(ip, k) = c_wrk1d(j, 6)*norm ! solution
                end do

            end do
        end do

        ! ###################################################################
        ! Fourier field a (based on array tmp1)
        ! ###################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d) ! tmp1 might be overwritten
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, a)
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, a)
        end if

        nullify (c_tmp1, c_tmp2)

        return
    end subroutine OPR_Helmholtz_FourierXZ_Direct

! !########################################################################
! !########################################################################
! ! Same, but for n fields
! ! I THINK THIS VERSION FIXES A PREVIOUS BUG BUT NEEDS TO BE TESTED
!     subroutine OPR_HELMHOLTZ_FXZ_D_N(nx, ny, nz, nfield, g, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
!         use TLab_Pointers, only: pointers_dt

!         integer(wi), intent(in) :: nx, ny, nz, nfield
!         integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
!         !                                                   1, for Neumann   & Dirichlet
!         !                                                   2, for Dirichlet & Neumann
!         !                                                   3, for Neumann   & Neumann
!         type(fdm_dt), intent(in) :: g(3)
!         real(wp), intent(in) :: alpha
!         type(pointers_dt), intent(in) :: a(nfield)                      ! Forcing term, and solution field p
!         real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz, nfield)     ! FFT of forcing term
!         real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
!         real(wp), intent(in) :: bcs_hb(nx, nz, nfield), bcs_ht(nx, nz, nfield)      ! Boundary-condition fields

!         target tmp1, tmp2

!         ! -----------------------------------------------------------------------
!         integer ifield, ip_sol
!         complex(wp) :: bcs_n(nfield, 2)
!         complex(wp), pointer :: aux_n(:, :, :) => null()
!         complex(wp), pointer :: c_tmp1_n(:, :, :) => null()

!         integer, parameter :: i1 = 1, i2 = 2

!         ! #######################################################################
!         if (ibc /= 0) then ! So far only implemented for Dirichlet BCs
!             call TLab_Write_ASCII(efile, 'OPR_HELMHOLT_FXZ_D. Undeveloped BCs.')
!             call TLab_Stop(DNS_ERROR_UNDEVELOP)
!         end if

!         call c_f_pointer(c_loc(tmp1), c_tmp1_n, shape=[isize_txc_dimz/2, nz, nfield])
!         call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
!         call c_f_pointer(c_loc(p_wrk1d(1, 9)), aux_n, shape=[nfield, ny, 2]) ! lines of forcing and solution

!         norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

!         isize_line = nx/2 + 1

!         ! #######################################################################
!         ! Fourier transform of forcing term; output of this section in array tmp1
!         ! #######################################################################
!         do ifield = 1, nfield
!             if (g(3)%size > 1) then
!                 call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, &
!                                           bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), c_tmp2)
!                 call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1_n(:, :, ifield)) ! tmp2 might be overwritten
!             else
!                 call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, &
!                                           bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), c_tmp1_n(:, :, ifield))
!             end if
!         end do

!         ! ###################################################################
!         ! Solve FDE \hat{p}''-(\lambda+lpha) \hat{p} = \hat{f}
!         ! ###################################################################
!         do k = 1, nz
! #ifdef USE_MPI
!             kglobal = k + ims_offset_k
! #else
!             kglobal = k
! #endif

!             do i = 1, isize_line
! #ifdef USE_MPI
!                 iglobal = i + ims_offset_i/2
! #else
!                 iglobal = i
! #endif

!                 ! Define \lambda based on modified wavenumbers (real)
!                 if (g(3)%size > 1) then
!                     lambda = g(1)%mwn2(iglobal) + g(3)%mwn2(kglobal)
!                 else
!                     lambda = g(1)%mwn2(iglobal)
!                 end if

!                 lambda = lambda - alpha

!                 ! forcing term in aux_n(:,:,1), i.e. p_wrk1d(:,9), solution will be in aux_n(:,:,2)
!                 do ifield = 1, nfield
!                     do j = 1, ny
!                         ip = (j - 1)*isize_line + i; aux_n(ifield, j, 1) = c_tmp1_n(ip, k, ifield)
!                     end do

!                     ! BCs
!                     j = ny + 1; ip = (j - 1)*isize_line + i; bcs_n(ifield, 1) = c_tmp1_n(ip, k, ifield) ! Dirichlet or Neumann
!                     j = ny + 2; ip = (j - 1)*isize_line + i; bcs_n(ifield, 2) = c_tmp1_n(ip, k, ifield) ! Dirichlet or Neumann

!                 end do
!                 ip_sol = 9 + nfield*2

!                 ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
!                 ! if (ibc == 0) then ! Dirichlet BCs
!                 p_wrk1d(:, 1:7) = 0.0_wp
!                 call FDM_Int2_CreateSystem(ny, g(2)%nodes, ibc, lhs, rhs, lambda, &
!                                     p_wrk1d(:, 1:5), p_wrk1d(:, 6:7), p_wrk1d(:, 13:14))
!                 call INT_C2NX_RHS(ny, i2, lhs, p_wrk1d(1, 9), p_wrk1d(1, 11))

!                 call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

!                 call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
!                 call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

!      call PENTADSS(ny - 2, i2*nfield, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(3, ip_sol))

!                 do ifield = 1, nfield
!                     ! BCa
!                     aux_n(ifield, :, 2) = aux_n(ifield, :, 2) &
!                                           + bcs_n(ifield, 1)*p_wrk1d(:, 6) + bcs_n(ifield, 2)*p_wrk1d(:, 7)

!                     ! Rearrange in memory and normalize
!                     do j = 1, ny
!                         ip = (j - 1)*isize_line + i
!                         c_tmp1_n(ip, k, ifield) = aux_n(ifield, j, 2)*norm ! solution
!                     end do

!                 end do

!             end do
!         end do

!         ! ###################################################################
!         ! Fourier field a (based on array tmp1)
!         ! ###################################################################
!         do ifield = 1, nfield
!             if (g(3)%size > 1) then
!                 call OPR_FOURIER_B_Z_EXEC(c_tmp1_n(:, :, ifield), c_wrk3d) ! tmp1 might be overwritten
!                 call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, a(ifield)%field) ! wrk3d might be overwritten
!             else
!                 call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1_n(:, :, ifield), a(ifield)%field)
!             end if
!         end do

!         nullify (c_tmp1_n, c_tmp2, aux_n)

!         return
!     end subroutine OPR_HELMHOLTZ_FXZ_D_N

end module OPR_ELLIPTIC
