#include "dns_const.h"
#include "dns_error.h"

module OPR_ELLIPTIC
    use TLAB_CONSTANTS
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: isize_txc_dimz, imax, jmax, kmax
    use TLAB_VARS, only: stagger_on
    use TLAB_POINTERS_3D, only: p_wrk1d
    use TLAB_POINTERS_C, only: c_wrk1d, c_wrk3d
    use TLAB_PROCS
    use OPR_FOURIER
    use OPR_FDE
    use OPR_PARTIAL
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    complex(wp), target :: bcs(3)
    real(wp), pointer :: r_bcs(:) => null()
    complex(wp), pointer :: c_tmp1(:, :) => null(), c_tmp2(:, :) => null()
    integer(wi) i, j, k, iglobal, kglobal, ip, isize_line
    real(wp) lambda, norm
    real(wp), allocatable :: lhs(:, :), rhs(:, :)
    real(wp), allocatable, target :: lu_poisson(:, :)       ! 3D array; here or in TLAB_ARRAYS?
    real(wp), pointer :: a(:, :, :), b(:, :, :), c(:, :, :), d(:, :, :), e(:, :, :), f1(:, :, :), f2(:, :, :)

    public :: OPR_ELLIPTIC_INITIALIZE
    public :: OPR_POISSON_FXZ
    public :: OPR_POISSON_FXZ_D         ! Using direct formulation of FDM schemes
    ! public :: OPR_POISSON_FXZ_D_TRANSPOSE   ! Not yet working; txc shape (isize_txc_dimz,nz) cannot always be reshape to contiguous (nx+2)*(ny+2)*nz
    public :: OPR_HELMHOLTZ_FXZ
    public :: OPR_HELMHOLTZ_FXZ_D       ! Using direct formulation of FDM schemes
    public :: OPR_HELMHOLTZ_FXZ_D_N     ! For N fields

contains
! #######################################################################
! #######################################################################
! We precalculate the LU factorization for the case BCS_NN, which is the one used in the pressure-Poisson equation
    subroutine OPR_ELLIPTIC_INITIALIZE()
        use TLAB_VARS, only: g, imode_elliptic
        use FDM_COM_DIRECT
        
        integer ibc_loc
        integer, parameter :: i1 = 1, i2 = 2

        ! -----------------------------------------------------------------------
        select case (imode_elliptic)
        case (FDM_COM4_DIRECT)
            allocate (lhs(g(2)%size, 3), rhs(g(2)%size, 4))
            call FDM_C2N4ND_INITIALIZE(g(2)%size, g(2)%nodes, lhs, rhs)

        case default !(FDM_COM6_DIRECT) ! I need it for helmholtz
            allocate (lhs(g(2)%size, 3), rhs(g(2)%size, 4))
            call FDM_C2N6ND_INITIALIZE(g(2)%size, g(2)%nodes, lhs, rhs)

        end select

        ! -----------------------------------------------------------------------
        ! LU factorization for direct cases in case BCS_NN, the one for the pressure equation; needs 5 3D arrays
        select case (imode_elliptic)
        case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
            isize_line = imax/2 + 1

            call TLAB_ALLOCATE_ARRAY_DOUBLE(__FILE__, lu_poisson, [g(2)%size*isize_line*kmax, 7], 'lu_poisson')

            a(1:g(2)%size, 1:isize_line, 1:kmax) => lu_poisson(1:g(2)%size*isize_line*kmax, 1)
            b(1:g(2)%size, 1:isize_line, 1:kmax) => lu_poisson(1:g(2)%size*isize_line*kmax, 2)
            c(1:g(2)%size, 1:isize_line, 1:kmax) => lu_poisson(1:g(2)%size*isize_line*kmax, 3)
            d(1:g(2)%size, 1:isize_line, 1:kmax) => lu_poisson(1:g(2)%size*isize_line*kmax, 4)
            e(1:g(2)%size, 1:isize_line, 1:kmax) => lu_poisson(1:g(2)%size*isize_line*kmax, 5)
            f1(1:g(2)%size, 1:isize_line, 1:kmax) => lu_poisson(1:g(2)%size*isize_line*kmax, 6)
            f2(1:g(2)%size, 1:isize_line, 1:kmax) => lu_poisson(1:g(2)%size*isize_line*kmax, 7)

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

                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda = g(1)%mwn(iglobal, 2) + g(3)%mwn(kglobal, 2)
                    else
                        lambda = g(1)%mwn(iglobal, 2)
                    end if

                    ! Compatibility constraint. The reference value of p at the lower boundary is set to zero
                    if (iglobal == 1 .and. kglobal == 1) then
                        ibc_loc = BCS_DN
                    else
                        ibc_loc = BCS_NN
                    end if

                    ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                    if (g(2)%uniform) then          ! FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA
                        call INT_C2N6_LHS_E(g(2)%size, g(2)%jac, ibc_loc, lambda, &
                                            a(1, i, k), b(1, i, k), c(1, i, k), d(1, i, k), e(1, i, k), f1(1, i, k), f2(1, i, k))

                    else                            ! FDM_COM6_DIRECT, although this is = to FDM_COM6_JACOBIAN if uniform
                        call INT_C2NXND_LHS_E(g(2)%size, g(2)%nodes, ibc_loc, lhs, rhs, lambda, &
                                              a(1, i, k), b(1, i, k), c(1, i, k), d(1, i, k), e(1, i, k), f1(1, i, k), f2(1, i, k))

                    end if

                    ! LU decomposizion
                    call PENTADFS(g(2)%size - 2, a(2, i, k), b(2, i, k), c(2, i, k), d(2, i, k), e(2, i, k))

                    ! Particular solutions
                    call PENTADSS(g(2)%size - 2, i1, a(2, i, k), b(2, i, k), c(2, i, k), d(2, i, k), e(2, i, k), f1(2, i, k))
                    call PENTADSS(g(2)%size - 2, i1, a(2, i, k), b(2, i, k), c(2, i, k), d(2, i, k), e(2, i, k), f2(2, i, k))

                end do
            end do

        end select

        return
    end subroutine OPR_ELLIPTIC_INITIALIZE

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
!# We use c_wrk1d and p_wrk1d for complex and real reference to same data (see tlab_procs%define_pointers_c)
!#
!########################################################################
    subroutine OPR_POISSON_FXZ(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(grid_dt), intent(in) :: g(3)
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields
        real(wp), intent(out), optional :: dpdy(nx, ny, nz)           ! Vertical derivative of solution

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        integer(wi) i_sing(2), k_sing(2)    ! singular global modes

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(bcs), r_bcs, shape=[3*2])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        if (.not. stagger_on) then
            i_sing = [1, g(1)%size/2 + 1]
            k_sing = [1, g(3)%size/2 + 1]
        else                    ! In case of staggering only one singular mode + different modified wavenumbers
            i_sing = [1, 1]
            k_sing = [1, 1]
        end if

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

                ! Define \lambda based on modified wavenumbers (real)
                if (g(3)%size > 1) then
                    lambda = g(1)%mwn(iglobal, 1) + g(3)%mwn(kglobal, 1)
                else
                    lambda = g(1)%mwn(iglobal, 1)
                end if

                ! forcing term
                do j = 1, ny
                    ip = (j - 1)*isize_line + i; c_wrk1d(j, 1) = c_tmp1(ip, k)
                end do

                ! BCs
                j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
                j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                select case (ibc)
                case (BCS_NN) ! Neumann   & Neumann   BCs
                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                        call FDE_BVP_SINGULAR_NN(g(2)%mode_fdm, ny, 2, &
                                                 g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    else
                        call FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    end if

                case (BCS_DD) ! Dirichlet & Dirichlet BCs
                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                        call FDE_BVP_SINGULAR_DD(g(2)%mode_fdm, ny, 2, &
                                            g(2)%nodes, g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    else
                        call FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    end if

                end select

                ! Rearrange in memory and normalize
                do j = 1, ny
                    ip = (j - 1)*isize_line + i
                    c_tmp1(ip, k) = c_wrk1d(j, 2)*norm ! solution
                    c_tmp2(ip, k) = c_wrk1d(j, 3)*norm ! Oy derivative
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
            if (g(3)%size > 1) then
                call OPR_FOURIER_B_Z_EXEC(c_tmp2, c_wrk3d)              ! tmp2 might be overwritten
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, dpdy)    ! wrk3d might be overwritten
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp2, dpdy)     ! tmp2 might be overwritten
            end if
        end if

        nullify (c_tmp1, c_tmp2, r_bcs)

        return
    end subroutine OPR_POISSON_FXZ

!########################################################################
!########################################################################
! The tranpose step does not work because the txc arrays are not necessarily contiguous...
! Same, but using the direct mode of FDM
! Opposite to previous routine, here we use the first 8 wrk1d arrays for the diagonals of the LHS,
! and the last ones for the forcing and solution. The reason is the routine after this one.
    subroutine OPR_POISSON_FXZ_D_TRANSPOSE(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(grid_dt), intent(in) :: g(3)
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields
        real(wp), intent(out), optional :: dpdy(nx, ny, nz)           ! Vertical derivative of solution

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        integer(wi) ibc_loc, bcs_p(2, 2)
        integer, parameter :: i1 = 1, i2 = 2

        integer(wi) ny_txc, nyz_txc
        complex(wp), pointer :: c_tmp1_t(:, :, :) => null(), c_tmp2_t(:, :, :) => null()
        real(wp), pointer :: r_tmp1_t(:, :, :) => null(), r_tmp2_t(:, :, :) => null()

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ny_txc = isize_txc_dimz/2/isize_line
        print *, isize_txc_dimz/2, isize_line, isize_line*(ny+2)
        nyz_txc = ny_txc*nz
        call c_f_pointer(c_loc(tmp2), c_tmp2_t, shape=[ny_txc, nz, isize_line])
        call c_f_pointer(c_loc(tmp2), r_tmp2_t, shape=[ny_txc*2, nz, isize_line])
        call c_f_pointer(c_loc(tmp1), c_tmp1_t, shape=[ny_txc, nz, isize_line])
        call c_f_pointer(c_loc(tmp1), r_tmp1_t, shape=[ny_txc*2, nz, isize_line])

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp1)
            call OPR_FOURIER_F_Z_EXEC(c_tmp1, c_tmp2) ! tmp1 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2)
        end if

        ! ###################################################################
        ! Local transposition: make x-direction the last one
        ! ###################################################################
#ifdef USE_ESSL
        call ZGECMO(c_tmp2, isize_line, isize_line, nyz_txc, c_tmp1_t, nyz_txc)
#else
        call DNS_TRANSPOSE_COMPLEX(c_tmp2, isize_line, nyz_txc, isize_line, c_tmp1_t, nyz_txc)
#endif

        ! ###################################################################
        ! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{f}
        ! ###################################################################
        do i = 1, isize_line
#ifdef USE_MPI
            iglobal = i + ims_offset_i/2
#else
            iglobal = i
#endif

            do k = 1, nz
#ifdef USE_MPI
                kglobal = k + ims_offset_k
#else
                kglobal = k
#endif

                bcs(1) = c_tmp1_t(ny + 1, k, i)
                bcs(2) = c_tmp1_t(ny + 2, k, i)

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
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda = g(1)%mwn(iglobal, 2) + g(3)%mwn(kglobal, 2)
                    else
                        lambda = g(1)%mwn(iglobal, 2)
                    end if

                    ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                    if (g(2)%uniform) then         ! FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA
                        call INT_C2N6_LHS_E(ny, g(2)%jac, ibc_loc, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))

                    else                            ! FDM_COM6_DIRECT, although this is = to FDM_COM6_JACOBIAN if uniform
                        call INT_C2NXND_LHS_E(ny, g(2)%nodes, ibc_loc, lhs, rhs, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))

                    end if

                    ! LU factorization
                    call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

                    ! Parciular solutions
                 call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
                 call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

                end if

                ! -----------------------------------------------------------------------
                if (g(2)%uniform) then         ! FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA
                    call INT_C2N6_RHS(ny, i2, g(2)%jac, r_tmp1_t(1, k, i), r_tmp2_t(1, k, i))

                else                            ! FDM_COM6_DIRECT, although this is = to FDM_COM6_JACOBIAN if uniform
                    call INT_C2NXND_RHS(ny, i2, lhs, r_tmp1_t(1, k, i), r_tmp2_t(1, k, i))

                end if

                if (ibc /= BCS_NN) then     ! use local LU factorization
                    call PENTADSS(ny - 2, i2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), r_tmp2_t(3, k, i))

                    c_tmp2_t(1:ny, k, i) = (c_tmp2_t(1:ny, k, i) + bcs(1)*p_wrk1d(:, 6) + bcs(2)*p_wrk1d(:, 7))*norm

                    !   Corrections to the BCS_DD to account for Neumann
                    if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                        c_tmp2_t(1, k, i) = c_tmp2_t(1, k, i) + p_wrk1d(1, 3)*c_tmp2_t(2, k, i) &
                                            + p_wrk1d(1, 4)*c_tmp2_t(3, k, i) + p_wrk1d(1, 5)*c_tmp2_t(4, k, i)
                    end if

                    if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                        c_tmp2_t(ny, k, i) = c_tmp2_t(ny, k, i) + p_wrk1d(ny, 3)*c_tmp2_t(ny - 1, k, i) &
                                             + p_wrk1d(ny, 2)*c_tmp2_t(ny - 2, k, i) + p_wrk1d(ny, 1)*c_tmp2_t(ny - 3, k, i)
                    end if

                else                        ! use precalculated LU factorization
                    call PENTADSS(ny - 2, i2, a(2, i, k), b(2, i, k), c(2, i, k), d(2, i, k), e(2, i, k), r_tmp2_t(3, k, i))

                    c_tmp2_t(1:ny, k, i) = (c_tmp2_t(1:ny, k, i) + bcs(1)*f1(1:ny, i, k) + bcs(2)*f2(1:ny, i, k))*norm

                    !   Corrections to the BCS_DD to account for Neumann
                    if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                        c_tmp2_t(1, k, i) = c_tmp2_t(1, k, i) + c(1, i, k)*c_tmp2_t(2, k, i) &
                                            + d(1, i, k)*c_tmp2_t(3, k, i) + e(1, i, k)*c_tmp2_t(4, k, i)
                    end if

                    if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                        c_tmp2_t(ny, k, i) = c_tmp2_t(ny, k, i) + c(ny, i, k)*c_tmp2_t(ny - 1, k, i) &
                                             + b(ny, i, k)*c_tmp2_t(ny - 2, k, i) + a(ny, i, k)*c_tmp2_t(ny - 3, k, i)
                    end if

                end if

            end do
        end do

        ! ###################################################################
        ! Put arrays back in the order in which they came in
        ! ###################################################################
#ifdef USE_ESSL
        call ZGECMO(c_tmp2_t, nyz_txc, nyz_txc, isize_line, c_tmp1, isize_line)
#else
        call DNS_TRANSPOSE_COMPLEX(c_tmp2_t, nyz_txc, isize_line, nyz_txc, c_tmp1, isize_line)
#endif

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

        nullify (c_tmp1, c_tmp2, c_tmp1_t, c_tmp2_t, r_tmp1_t, r_tmp2_t)

        return
    end subroutine OPR_POISSON_FXZ_D_TRANSPOSE

!########################################################################
!########################################################################
! Same, but using the direct mode of FDM
! Opposite to previous routine, here we use the first 8 wrk1d arrays for the diagonals of the LHS,
! and the last ones for the forcing and solution. The reason is the routine after this one.
    subroutine OPR_POISSON_FXZ_D(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(grid_dt), intent(in) :: g(3)
        real(wp), intent(inout) :: p(nx, ny, nz)                        ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields
        real(wp), intent(out), optional :: dpdy(nx, ny, nz)           ! Vertical derivative of solution

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
                    ! Define \lambda based on modified wavenumbers (real)
                    if (g(3)%size > 1) then
                        lambda = g(1)%mwn(iglobal, 2) + g(3)%mwn(kglobal, 2)
                    else
                        lambda = g(1)%mwn(iglobal, 2)
                    end if

                    ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                    if (g(2)%uniform) then         ! FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA
                        call INT_C2N6_LHS_E(ny, g(2)%jac, ibc_loc, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))

                    else                            ! FDM_COM6_DIRECT, although this is = to FDM_COM6_JACOBIAN if uniform
                        call INT_C2NXND_LHS_E(ny, g(2)%nodes, ibc_loc, lhs, rhs, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))

                    end if

                    ! LU factorization
                    call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

                    ! Parciular solutions
                    call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
                    call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

                end if

                ! -----------------------------------------------------------------------
                if (g(2)%uniform) then         ! FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA
                    call INT_C2N6_RHS(ny, i2, g(2)%jac, p_wrk1d(1, 9), p_wrk1d(1, 11))

                else                            ! FDM_COM6_DIRECT, although this is = to FDM_COM6_JACOBIAN if uniform
                    call INT_C2NXND_RHS(ny, i2, lhs, p_wrk1d(1, 9), p_wrk1d(1, 11))

                end if

                if (ibc /= BCS_NN) then     ! use local LU factorization
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
                    call PENTADSS(ny - 2, i2, a(2, i, k), b(2, i, k), c(2, i, k), d(2, i, k), e(2, i, k), p_wrk1d(3, 11))

                    c_wrk1d(:, 6) = (c_wrk1d(:, 6) + bcs(1)*f1(:, i, k) + bcs(2)*f2(:, i, k))*norm

                    !   Corrections to the BCS_DD to account for Neumann
                    if (any([BCS_ND, BCS_NN] == ibc_loc)) then
                        c_wrk1d(1, 6) = c_wrk1d(1, 6) + c(1, i, k)*c_wrk1d(2, 6) &
                                        + d(1, i, k)*c_wrk1d(3, 6) + e(1, i, k)*c_wrk1d(4, 6) &
                                        + b(1, i, k)*c_wrk1d(2, 5)*norm
                    end if

                    if (any([BCS_DN, BCS_NN] == ibc_loc)) then
                        c_wrk1d(ny, 6) = c_wrk1d(ny, 6) + c(ny, i, k)*c_wrk1d(ny - 1, 6) &
                                         + b(ny, i, k)*c_wrk1d(ny - 2, 6) + a(ny, i, k)*c_wrk1d(ny - 3, 6) &
                                         + e(ny, i, k)*c_wrk1d(ny - 1, 5)*norm
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
    end subroutine OPR_POISSON_FXZ_D

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
    subroutine OPR_HELMHOLTZ_FXZ(nx, ny, nz, g, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(grid_dt), intent(in) :: g(3)
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                       ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields

        target tmp1, tmp2

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

                ! Define \lambda based on modified wavenumbers (real)
                if (g(3)%size > 1) then
                    lambda = g(1)%mwn(iglobal, 1) + g(3)%mwn(kglobal, 1)
                else
                    lambda = g(1)%mwn(iglobal, 1)
                end if

                lambda = lambda - alpha

                ! forcing term
                do j = 1, ny
                    ip = (j - 1)*isize_line + i; c_wrk1d(j, 1) = c_tmp1(ip, k)
                end do

                ! BCs
                j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
                j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                select case (ibc)
                case (3) ! Neumann   & Neumann   BCs
                    call FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny, 2, lambda, &
                                            g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))

                case (0) ! Dirichlet & Dirichlet BCs
                    call FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny, 2, lambda, &
                                            g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))

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
    end subroutine OPR_HELMHOLTZ_FXZ

!########################################################################
!########################################################################
! Same, but using the direct mode of FDM
! Opposite to previous routine, here we use the first 8 wrk1d arrays for the diagonals of the LHS,
! and the last ones for the forcing and solution. The reason is the routine after this one.
    subroutine OPR_HELMHOLTZ_FXZ_D(nx, ny, nz, g, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(grid_dt), intent(in) :: g(3)
        real(wp), intent(in) :: alpha
        real(wp), intent(inout) :: a(nx, ny, nz)                       ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz)             ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(inout) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields

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

                ! Define \lambda based on modified wavenumbers (real)
                if (g(3)%size > 1) then
                    lambda = g(1)%mwn(iglobal, 2) + g(3)%mwn(kglobal, 2)
                else
                    lambda = g(1)%mwn(iglobal, 2)
                end if

                lambda = lambda - alpha

                ! forcing term in c_wrk1d(:,5), i.e. p_wrk1d(:,9), solution will be in c_wrk1d(:,6), i.e., p_wrk1d(:,11)
                do j = 1, ny
                    ip = (j - 1)*isize_line + i; c_wrk1d(j, 5) = c_tmp1(ip, k)
                end do

                ! BCs
                j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = c_tmp1(ip, k) ! Dirichlet or Neumann
                j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = c_tmp1(ip, k) ! Dirichlet or Neumann

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                if (g(2)%uniform) then          ! FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA
                    call INT_C2N6_LHS_E(ny, g(2)%jac, ibc, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2N6_RHS(ny, i2, g(2)%jac, p_wrk1d(1, 9), p_wrk1d(1, 11))

                else                            ! FDM_COM6_DIRECT, although this is = to FDM_COM6_JACOBIAN if uniform
                    p_wrk1d(:, 1:7) = 0.0_wp
                    call INT_C2NXND_LHS_E(ny, g(2)%nodes, ibc, lhs, rhs, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2NXND_RHS(ny, i2, lhs, p_wrk1d(1, 9), p_wrk1d(1, 11))

                end if

                call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

                call PENTADSS(ny - 2, i2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(3, 11))

                c_wrk1d(:, 6) = c_wrk1d(:, 6) + bcs(1)*p_wrk1d(:, 6) + bcs(2)*p_wrk1d(:, 7)

                !   Corrections to the BCS_DD to account for Neumann
                if (any([BCS_ND, BCS_NN] == ibc)) then
             c_wrk1d(1, 6) = c_wrk1d(1, 6) + p_wrk1d(1, 3)*c_wrk1d(2, 6) + p_wrk1d(1, 4)*c_wrk1d(3, 6) + p_wrk1d(1, 5)*c_wrk1d(4, 6)
                end if

                if (any([BCS_DN, BCS_NN] == ibc)) then
                    c_wrk1d(ny, 6) = c_wrk1d(ny, 6) + p_wrk1d(ny, 3)*c_wrk1d(ny - 1, 6) + p_wrk1d(ny, 2)*c_wrk1d(ny - 2, 6) + p_wrk1d(ny, 1)*c_wrk1d(ny - 3, 6)
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
    end subroutine OPR_HELMHOLTZ_FXZ_D

!########################################################################
!########################################################################
! Same, but for n fields
! I THINK THIS VERSION FIXES A PREVIOUS BUG BUT NEEDS TO BE TESTED
    subroutine OPR_HELMHOLTZ_FXZ_D_N(nx, ny, nz, nfield, g, ibc, alpha, a, tmp1, tmp2, bcs_hb, bcs_ht)
        use TLAB_TYPES, only: pointers_dt

        integer(wi), intent(in) :: nx, ny, nz, nfield
        integer, intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(grid_dt), intent(in) :: g(3)
        real(wp), intent(in) :: alpha
        type(pointers_dt), intent(in) :: a(nfield)                      ! Forcing term, and solution field p
        real(wp), intent(inout) :: tmp1(isize_txc_dimz, nz, nfield)     ! FFT of forcing term
        real(wp), intent(inout) :: tmp2(isize_txc_dimz, nz)             ! Aux array for FFT
        real(wp), intent(in) :: bcs_hb(nx, nz, nfield), bcs_ht(nx, nz, nfield)      ! Boundary-condition fields

        target tmp1, tmp2

        ! -----------------------------------------------------------------------
        integer ifield, ip_sol
        complex(wp) :: bcs_n(nfield, 2)
        complex(wp), pointer :: aux_n(:, :, :) => null()
        complex(wp), pointer :: c_tmp1_n(:, :, :) => null()

        integer, parameter :: i1 = 1, i2 = 2

        ! #######################################################################
        if (ibc /= 0) then ! So far only implemented for Dirichlet BCs
            call TLAB_WRITE_ASCII(efile, 'OPR_HELMHOLT_FXZ_D. Undeveloped BCs.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        call c_f_pointer(c_loc(tmp1), c_tmp1_n, shape=[isize_txc_dimz/2, nz, nfield])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(p_wrk1d(1, 9)), aux_n, shape=[nfield, ny, 2]) ! lines of forcing and solution

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        do ifield = 1, nfield
            if (g(3)%size > 1) then
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, &
                                          bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), c_tmp2)
                call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1_n(:, :, ifield)) ! tmp2 might be overwritten
            else
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, &
                                          bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), c_tmp1_n(:, :, ifield))
            end if
        end do

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

                ! Define \lambda based on modified wavenumbers (real)
                if (g(3)%size > 1) then
                    lambda = g(1)%mwn(iglobal, 2) + g(3)%mwn(kglobal, 2)
                else
                    lambda = g(1)%mwn(iglobal, 2)
                end if

                lambda = lambda - alpha

                ! forcing term in aux_n(:,:,1), i.e. p_wrk1d(:,9), solution will be in aux_n(:,:,2)
                do ifield = 1, nfield
                    do j = 1, ny
                        ip = (j - 1)*isize_line + i; aux_n(ifield, j, 1) = c_tmp1_n(ip, k, ifield)
                    end do

                    ! BCs
                    j = ny + 1; ip = (j - 1)*isize_line + i; bcs_n(ifield, 1) = c_tmp1_n(ip, k, ifield) ! Dirichlet or Neumann
                    j = ny + 2; ip = (j - 1)*isize_line + i; bcs_n(ifield, 2) = c_tmp1_n(ip, k, ifield) ! Dirichlet or Neumann

                end do
                ip_sol = 9 + nfield*2

                ! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
                ! if (ibc == 0) then ! Dirichlet BCs
                if (g(2)%uniform) then         ! FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA
                    call INT_C2N6_LHS_E(ny, g(2)%jac, ibc, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2N6_RHS(ny, i2*nfield, g(2)%jac, p_wrk1d(1, 9), p_wrk1d(1, ip_sol))

                else                        ! FDM_COM6_DIRECT, although this is = to FDM_COM6_JACOBIAN if uniform
                    p_wrk1d(:, 1:7) = 0.0_wp
                    call INT_C2NXND_LHS_E(ny, g(2)%nodes, ibc, lhs, rhs, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2NXND_RHS(ny, i2, lhs, p_wrk1d(1, 9), p_wrk1d(1, 11))

                end if

                call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

     call PENTADSS(ny - 2, i2*nfield, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(3, ip_sol))

                do ifield = 1, nfield
                    ! BCa
                    aux_n(ifield, :, 2) = aux_n(ifield, :, 2) &
                                          + bcs_n(ifield, 1)*p_wrk1d(:, 6) + bcs_n(ifield, 2)*p_wrk1d(:, 7)

                    ! Rearrange in memory and normalize
                    do j = 1, ny
                        ip = (j - 1)*isize_line + i
                        c_tmp1_n(ip, k, ifield) = aux_n(ifield, j, 2)*norm ! solution
                    end do

                end do

            end do
        end do

        ! ###################################################################
        ! Fourier field a (based on array tmp1)
        ! ###################################################################
        do ifield = 1, nfield
            if (g(3)%size > 1) then
                call OPR_FOURIER_B_Z_EXEC(c_tmp1_n(:, :, ifield), c_wrk3d) ! tmp1 might be overwritten
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, a(ifield)%field) ! wrk3d might be overwritten
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1_n(:, :, ifield), a(ifield)%field)
            end if
        end do

        nullify (c_tmp1_n, c_tmp2, aux_n)

        return
    end subroutine OPR_HELMHOLTZ_FXZ_D_N

end module OPR_ELLIPTIC
