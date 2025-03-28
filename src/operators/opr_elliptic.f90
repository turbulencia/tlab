#include "dns_const.h"
#include "dns_error.h"
! You need to split the routines into the ones that are initialized and the ones that not.
! If not initialized, you can enter with any jmax, but the periodic directions need to be the global ones because of OPR_Fourier.
module OPR_ELLIPTIC
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: lfile
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use TLab_Memory, only: TLab_Allocate_Real
    use TLab_Memory, only: isize_txc_dimz, imax, jmax, kmax
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLab_Arrays, only: wrk1d, wrk2d, wrk3d
    use TLab_Pointers_3D, only: p_wrk1d
    use TLab_Pointers_C, only: c_wrk1d, c_wrk3d
    use FDM, only: fdm_dt, FDM_COM4_DIRECT, FDM_COM6_DIRECT
    use FDM_Integral
    use OPR_FOURIER
    use OPR_ODES
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
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
    real(wp) norm
    integer(wi) i_sing(2), k_sing(2)                                ! singular global modes

    type(fdm_dt) fdm_loc                                            ! scheme used for the elliptic solvers

    type(fdm_integral_dt), allocatable :: fdm_int1(:, :, :)         ! factorized method
    real(wp), allocatable, target :: rhs_b(:, :), rhs_t(:, :)       ! rhs to free memory space
    type(fdm_integral_dt) :: fdm_int1_loc(2)

    type(fdm_integral_dt), allocatable :: fdm_int2(:, :)            ! direct method
    real(wp), allocatable, target :: rhs_d(:, :)                    ! rhs to free memory space
    type(fdm_integral_dt) :: fdm_int2_loc

    real(wp), allocatable :: lambda(:, :)

    complex(wp), target :: bcs(3)
    real(wp), pointer :: r_bcs(:) => null()
    complex(wp), pointer :: c_tmp1(:, :) => null(), c_tmp2(:, :) => null()
    integer(wi) i, j, k, iglobal, kglobal, ip, isize_line

contains
    ! #######################################################################
    ! #######################################################################
    subroutine OPR_Elliptic_Initialize(inifile)
        use FDM, only: g, FDM_Initialize

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer imode_elliptic
        integer, parameter :: TYPE_FACTORIZE = 1
        integer, parameter :: TYPE_DIRECT = 2

        integer(wi) :: ndl, ndr, nd
        character*512 sRes
        character*32 bakfile

        ! ###################################################################
        ! Reading
        bakfile = trim(adjustl(inifile))//'.bak'

        imode_elliptic = TYPE_FACTORIZE         ! default is the finite-difference method used for the derivatives
        fdm_loc%mode_fdm1 = g(2)%mode_fdm1      ! to impose zero divergence down to round-off error in the interior points
        fdm_loc%der2%mode_fdm = g(2)%der2%mode_fdm

        call ScanFile_Char(bakfile, inifile, 'Main', 'EllipticOrder', 'void', sRes)
        if (trim(adjustl(sRes)) == 'compactdirect4') then
            imode_elliptic = TYPE_DIRECT
            fdm_loc%mode_fdm1 = FDM_COM4_DIRECT
            fdm_loc%der2%mode_fdm = FDM_COM4_DIRECT
        else if (trim(adjustl(sRes)) == 'compactdirect6') then
            imode_elliptic = TYPE_DIRECT
            fdm_loc%mode_fdm1 = FDM_COM6_DIRECT
            fdm_loc%der2%mode_fdm = FDM_COM6_DIRECT
        end if

        ! ###################################################################
        ! Initializing
        call FDM_Initialize(g(2)%nodes, fdm_loc)

        isize_line = imax/2 + 1

        allocate (lambda(isize_line, kmax))
        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        select case (imode_elliptic)
        case (TYPE_FACTORIZE)
            OPR_Poisson => OPR_Poisson_FourierXZ_Factorize
            ! OPR_Poisson => OPR_Poisson_FourierXZ_Factorize_Old
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Factorize

            ndl = fdm_loc%nb_diag_1(1)
            ndr = fdm_loc%nb_diag_1(2)
            nd = ndl
            allocate (fdm_int1(2, isize_line, kmax))
            call TLab_Allocate_Real(__FILE__, rhs_b, [g(2)%size, nd], 'rhs_b')
            call TLab_Allocate_Real(__FILE__, rhs_t, [g(2)%size, nd], 'rhs_t')

            if (.not. stagger_on) then
                i_sing = [1, g(1)%size/2 + 1]
                k_sing = [1, g(3)%size/2 + 1]
            else                    ! In case of staggering only one singular mode + different modified wavenumbers
                i_sing = [1, 1]
                k_sing = [1, 1]
            end if

        case (TYPE_DIRECT)
            OPR_Poisson => OPR_Poisson_FourierXZ_Direct
            ! OPR_Poisson => OPR_Poisson_FourierXZ_Direct_Old
            OPR_Helmholtz => OPR_Helmholtz_FourierXZ_Direct

            ndl = fdm_loc%der2%nb_diag(1)
            ndr = fdm_loc%der2%nb_diag(2)
            nd = ndl
            allocate (fdm_int2(isize_line, kmax))
            call TLab_Allocate_Real(__FILE__, rhs_d, [g(2)%size, nd], 'rhs_d')

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
                case (TYPE_FACTORIZE)
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda(i, k) = g(1)%mwn1(iglobal)**2.0 + g(3)%mwn1(kglobal)**2.0
                    else
                        lambda(i, k) = g(1)%mwn1(iglobal)**2.0
                    end if

                    fdm_int1(BCS_MIN, i, k)%bc = BCS_MIN
                    fdm_int1(BCS_MIN, i, k)%mode_fdm = fdm_loc%mode_fdm1
                    call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                             sqrt(lambda(i, k)), fdm_int1(BCS_MIN, i, k))

                    fdm_int1(BCS_MAX, i, k)%bc = BCS_MAX
                    fdm_int1(BCS_MAX, i, k)%mode_fdm = fdm_loc%mode_fdm1
                    call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                             -sqrt(lambda(i, k)), fdm_int1(BCS_MAX, i, k))

                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                    else                                        ! free memory that is independent of lambda
                        rhs_b(:, :) = fdm_int1(BCS_MIN, i, k)%rhs(:, :)
                        if (allocated(fdm_int1(BCS_MIN, i, k)%rhs)) deallocate (fdm_int1(BCS_MIN, i, k)%rhs)

                        rhs_t(:, :) = fdm_int1(BCS_MAX, i, k)%rhs(:, :)
                        if (allocated(fdm_int1(BCS_MAX, i, k)%rhs)) deallocate (fdm_int1(BCS_MAX, i, k)%rhs)

                    end if

                    ! idr = ndr/2 + 1
                    ! fdm_int1(:, i, k)%lhs(:, idr) = fdm_int1(:, i, k)%lhs(:, idr)*norm
                    ! fdm_int1(:, i, k)%lhs(:, idr + 1:ndr) = fdm_int1(:, i, k)%lhs(:, idr + 1:ndr)/norm

                case (TYPE_DIRECT)     ! only for case BCS_NN
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda(i, k) = g(1)%der2%mwn(iglobal) + g(3)%der2%mwn(kglobal)
                    else
                        lambda(i, k) = g(1)%der2%mwn(iglobal)
                    end if

                    ! Compatibility constraint. The reference value of p at the lower boundary will be set to zero
                    if (iglobal == 1 .and. kglobal == 1) then
                        fdm_int2(i, k)%bc = BCS_DN
                    else
                        fdm_int2(i, k)%bc = BCS_NN
                    end if

                    ! free memory that is independent of lambda
                    call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2%lhs(:, 1:ndl), fdm_loc%der2%rhs(:, 1:ndr), lambda(i, k), fdm_int2(i, k))
                    rhs_d(:, :) = fdm_int2(i, k)%rhs(:, :)
                    if (allocated(fdm_int2(i, k)%rhs)) deallocate (fdm_int2(i, k)%rhs)

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
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten; cannot use wrk3d
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
                        call OPR_ODE2_Factorize_NN_Sing(2, fdm_int1(:, i, k), u(:, i), f(:, i), f(2*ny + 1:, i), v(:, i), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_Factorize_NN(2, fdm_int1(:, i, k), rhs_b, rhs_t, &
                                                   u(:, i), f(:, i), f(2*ny + 1:, i), v(:, i), wrk1d, wrk2d)
                    end if

                case (BCS_DD) ! Dirichlet & Dirichlet BCs
                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                        call OPR_ODE2_Factorize_DD_Sing(2, fdm_int1(:, i, k), u(:, i), f(:, i), f(2*ny + 1:, i), v(:, i), wrk1d, wrk2d)
                    else
                        call OPR_ODE2_Factorize_DD(2, fdm_int1(:, i, k), rhs_b, rhs_t, &
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
        integer(wi) ibc_loc, bcs_p(2, 2), ndl, ndr
        real(wp), pointer :: u(:, :) => null(), f(:, :) => null()

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])

        ndl = fdm_loc%der2%nb_diag(1)
        ndr = fdm_loc%der2%nb_diag(2)

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten; cannot use wrk3d
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
            u(1:2*ny, 1:isize_line) => wrk3d(1:2*ny*isize_line)

            do i = 1, isize_line
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif

                ! BCs
                ip = 2*ny
                u(1:2, i) = f(ip + 1:ip + 2, i)
                u(ip - 1:ip, i) = f(ip + 3:ip + 4, i)

                ! Compatibility constraint for singular modes. 2nd order FDMs are non-zero at Nyquist
                ! The reference value of p at the lower boundary is set to zero
                if (iglobal == 1 .and. kglobal == 1 .and. ibc == BCS_NN) then
                    ibc_loc = BCS_DN
                    u(1:2, i) = 0.0_wp
                else
                    ibc_loc = ibc
                end if

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                if (ibc /= BCS_NN) then     ! Need to calculate and factorize LHS
                    fdm_int2_loc%bc = ibc_loc
                    call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2%lhs(:, 1:ndl), fdm_loc%der2%rhs(:, 1:ndr), lambda(i, k), fdm_int2_loc)
                    call FDM_Int2_Solve(2, fdm_int2_loc, fdm_int2_loc%rhs, f(:, i), u(:, i), wrk2d)

                else                        ! use precalculated LU factorization
                    call FDM_Int2_Solve(2, fdm_int2(i, k), rhs_d, f(:, i), u(:, i), wrk2d)

                end if

            end do

            call TLab_Transpose_COMPLEX(c_wrk3d, ny, isize_line, ny, c_tmp1(:, k), isize_line)

        end do

        ! ###################################################################
        ! Fourier field p (based on array tmp1)
        ! ###################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, p)    ! tmp1 might be overwritten
        end if

        if (present(dpdy)) then
            bcs_p = 0
            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs_p, g(2), p, dpdy)
        end if

        nullify (u, f)
        nullify (c_tmp1, c_tmp2)

        return
    end subroutine OPR_Poisson_FourierXZ_Direct

!     ! ###################################################################
!     ! ###################################################################
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
!                         ! call OPR_ODE2_Factorize_1_SINGULAR_NN_OLD(g(2)%mode_fdm1, ny, 2, &
!                         !                             g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_Factorize_NN_Sing(2, fdm_Int0, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), wrk2d)
!                     else
!                         ! call OPR_ODE2_Factorize_1_REGULAR_NN_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k), &
!                         !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_Factorize_NN(2, fdm_int1(:, i, k), rhs_b, rhs_t, &
!                                                    p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), wrk2d)
!                     end if

!                 case (BCS_DD) ! Dirichlet & Dirichlet BCs
!                     if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
!                         ! call OPR_ODE2_Factorize_1_SINGULAR_DD_OLD(g(2)%mode_fdm1, ny, 2, &
!                         !                                 g(2)%nodes, g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_Factorize_DD_Sing(2, fdm_Int0, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), wrk2d)
!                     else
!                         ! call OPR_ODE2_Factorize_1_REGULAR_DD_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k), &
!                         !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
!                         call OPR_ODE2_Factorize_DD(2, fdm_int1(:, i, k), rhs_b, rhs_t, &
!                                                    p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), wrk2d)
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

!         ! ###################################################################
!         ! Fourier field p (based on array tmp1)
!         ! ###################################################################
!         if (g(3)%size > 1) then
!             call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
!             call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
!         else
!             call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, p)    ! tmp2 might be overwritten
!         end if

!         ! Fourier derivatives (based on array tmp2)
!         if (present(dpdy)) then
!             if (g(3)%size > 1) then
!                 call OPR_FOURIER_B_Z_EXEC(c_tmp2, c_wrk3d)              ! tmp2 might be overwritten
!                 call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, dpdy)    ! wrk3d might be overwritten
!             else
!                 call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp2, dpdy)     ! tmp2 might be overwritten
!             end if
!         end if

!         nullify (c_tmp1, c_tmp2, r_bcs)

!         return
!     end subroutine OPR_Poisson_FourierXZ_Factorize_Old

!     !########################################################################
!     !########################################################################
!     subroutine OPR_Poisson_FourierXZ_Direct_Old(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
!         integer(wi), intent(in) :: nx, ny, nz
!         integer, intent(in) :: ibc
!         type(fdm_dt), intent(in) :: g(3)
!         real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
!         real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
!         real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
!         real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)          ! Boundary-condition fields
!         real(wp), intent(out), optional :: dpdy(nx, ny, nz)             ! Vertical derivative of solution

!         target tmp1, tmp2

!         ! -----------------------------------------------------------------------
!         integer(wi) ibc_loc, bcs_p(2, 2), ndl, ndr

!         ! #######################################################################
!         call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
!         call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
!         call c_f_pointer(c_loc(bcs), r_bcs, shape=[3*2])

!         norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

!         isize_line = nx/2 + 1
!         ndl = fdm_loc%der2%nb_diag(1)
!         ndr = fdm_loc%der2%nb_diag(2)

!         ! #######################################################################
!         ! Fourier transform of forcing term; output of this section in array tmp1
!         ! #######################################################################
!         if (g(3)%size > 1) then
!             call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2)
!             call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten; cannot use wrk3d
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

!                 ! forcing term in c_wrk1d(:,5), i.e. p_wrk1d(:,9), solution will be in c_wrk1d(:,6), i.e., p_wrk1d(:,11)
!                 do j = 1, ny
!                     ip = (j - 1)*isize_line + i; c_wrk1d(j, 5) = c_tmp1(ip, k)
!                 end do

!                 ! BCs
!                 j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
!                 j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

!                 p_wrk1d(1:2, 11) = r_bcs(1:2)
!                 p_wrk1d(ny - 1:ny, 12) = r_bcs(3:4)

!                 ! Compatibility constraint for singular modes. 2nd order FDMs are non-zero at Nyquist
!                 ! The reference value of p at the lower boundary is set to zero
!                 if (iglobal == 1 .and. kglobal == 1 .and. ibc == BCS_NN) then
!                     ibc_loc = BCS_DN
!                     ! bcs(1) = (0.0_wp, 0.0_wp)
!                     p_wrk1d(1:2, 11) = 0.0_wp
!                 else
!                     ibc_loc = ibc
!                 end if

!                 ! -----------------------------------------------------------------------
!                 if (ibc /= BCS_NN) then     ! Need to calculate and factorize LHS
!                     ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
!                     fdm_int2_loc%bc = ibc_loc
!                     call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2%lhs(:, 1:ndl), fdm_loc%der2%rhs(:, 1:ndr), lambda(i, k), fdm_int2_loc)
!                     call FDM_Int2_Solve(2, fdm_int2_loc, fdm_int2_loc%rhs, p_wrk1d(:, 9), p_wrk1d(:, 11), wrk2d)

!                 else                        ! use precalculated LU factorization
!                     call FDM_Int2_Solve(2, fdm_int2(i, k), rhs_d, p_wrk1d(:, 9), p_wrk1d(:, 11), wrk2d)

!                 end if

!                 ! Rearrange in memory and normalize
!                 do j = 1, ny
!                     ip = (j - 1)*isize_line + i
!                     c_tmp1(ip, k) = c_wrk1d(j, 6)*norm
!                 end do

!             end do
!         end do

!         ! ###################################################################
!         ! Fourier field p (based on array tmp1)
!         ! ###################################################################
!         if (g(3)%size > 1) then
!             call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)          ! tmp1 might be overwritten
!             call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, p)   ! wrk3d might be overwritten
!         else
!             call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, p)    ! tmp2 might be overwritten
!         end if

!         ! Fourier derivatives (based on array tmp2)
!         if (present(dpdy)) then
!             bcs_p = 0
!             call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs_p, g(2), p, dpdy)
!         end if

!         nullify (c_tmp1, c_tmp2)

!         return
!     end subroutine OPR_Poisson_FourierXZ_Direct_Old

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
        integer, intent(in) :: ibc
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
                ndr = fdm_loc%der2%nb_diag(2)

                fdm_int1_loc(BCS_MIN)%bc = BCS_MIN
                fdm_int1_loc(BCS_MIN)%mode_fdm = fdm_loc%mode_fdm1
                call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                         sqrt(lambda(i, k) - alpha), fdm_int1_loc(BCS_MIN))

                fdm_int1_loc(BCS_MAX)%bc = BCS_MAX
                fdm_int1_loc(BCS_MAX)%mode_fdm = fdm_loc%mode_fdm1
                call FDM_Int1_Initialize(fdm_loc%nodes(:), fdm_loc%lhs1(:, 1:ndl), fdm_loc%rhs1(:, 1:ndr), &
                                         -sqrt(lambda(i, k) - alpha), fdm_int1_loc(BCS_MAX))

                select case (ibc)
                case (3) ! Neumann   & Neumann   BCs
                    ! call OPR_ODE2_Factorize_1_REGULAR_NN_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k)-alpha, &
                    !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    call OPR_ODE2_Factorize_NN(2, fdm_int1_loc, fdm_int1_loc(BCS_MIN)%rhs, fdm_int1_loc(BCS_MAX)%rhs, &
                                               p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), wrk2d)

                case (0) ! Dirichlet & Dirichlet BCs
                    ! call OPR_ODE2_Factorize_1_REGULAR_DD_OLD(g(2)%mode_fdm1, ny, 2, lambda(i,k)-alpha, &
                    !                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))

                    call OPR_ODE2_Factorize_DD(2, fdm_int1_loc, fdm_int1_loc(BCS_MIN)%rhs, fdm_int1_loc(BCS_MAX)%rhs, &
                                               p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7), wrk2d)
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
        integer, intent(in) :: ibc
        type(fdm_dt), intent(in) :: g(3)
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                       ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        integer ndl, ndr
        integer, parameter :: i1 = 1, i2 = 2

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(bcs), r_bcs, shape=[3*2])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1
        ndl = fdm_loc%der2%nb_diag(1)
        ndr = fdm_loc%der2%nb_diag(2)

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
                fdm_int2_loc%bc = ibc
                call FDM_Int2_Initialize(fdm_loc%nodes(:), fdm_loc%der2%lhs(:, 1:ndl), fdm_loc%der2%rhs(:, 1:ndr), lambda(i, k) - alpha, fdm_int2_loc)
                p_wrk1d(1:2, 11) = r_bcs(1:2)
                p_wrk1d(ny - 1:ny, 12) = r_bcs(3:4)
                call FDM_Int2_Solve(2, fdm_int2_loc, fdm_int2_loc%rhs, p_wrk1d(:, 9), p_wrk1d(:, 11), wrk2d)

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

end module OPR_ELLIPTIC
