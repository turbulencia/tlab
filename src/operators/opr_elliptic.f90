#include "dns_const.h"
#include "dns_error.h"

module OPR_ELLIPTIC
    use TLAB_CONSTANTS
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: stagger_on
    use TLAB_POINTERS_3D, only: p_wrk1d
    use TLAB_POINTERS_C, only: c_wrk1d, c_wrk3d
    use TLAB_PROCS
    use OPR_FOURIER
    use OPR_FDE
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

    public :: OPR_POISSON_FXZ
    public :: OPR_HELMHOLTZ_FXZ
    public :: OPR_HELMHOLTZ_FXZ_D       ! Using direct formulation of FDM schemes
    public :: OPR_HELMHOLTZ_FXZ_D_N     ! For N fields

contains
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
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2) !, c_tmp1, c_wrk3d)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp1) !, c_tmp2, c_wrk3d)
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
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, p, c_tmp1)
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, p, c_wrk3d)
        end if

        ! Fourier derivatives (based on array tmp2)
        if (present(dpdy)) then
            if (g(3)%size > 1) then
                call OPR_FOURIER_B_Z_EXEC(c_tmp2, c_wrk3d)
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, dpdy, c_tmp2)
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp2, dpdy, c_wrk3d)
            end if
        end if

        nullify (c_tmp1, c_tmp2, r_bcs)

        return
    end subroutine OPR_POISSON_FXZ

!########################################################################
!#
!# Solve Lap a + \alpha a = f using Fourier in xOz planes, to rewritte
!# the problem as
!#
!#      \hat{a}''-(\lambda-\alpha) \hat{a} = \hat{f}
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

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp2) !, c_tmp1, c_wrk3d)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp1) !, c_tmp2, c_wrk3d)
        end if

        ! ###################################################################
        ! Solve FDE (\hat{p}')'-(\lambda+\alpha) \hat{p} = \hat{f}
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
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, a, c_tmp1)
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, a, c_wrk3d)
        end if

        nullify (c_tmp1, c_tmp2)

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
        real(wp), intent(in) :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields

        target tmp1, tmp2

        ! -----------------------------------------------------------------------

        integer, parameter :: i1 = 1, i2 = 2

        ! #######################################################################
        if (ibc /= BCS_DD) then ! So far only implemented for Dirichlet BCs
            call TLAB_WRITE_ASCII(efile, 'OPR_HELMHOLT_FXZ_D. Undeveloped BCs.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp2) !, c_tmp1, c_wrk3d)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, c_tmp1) !, c_tmp2, c_wrk3d)
        end if

        ! ###################################################################
        ! Solve FDE \hat{p}''-(\lambda+\alpha) \hat{p} = \hat{f}
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
                ! if (ibc == 0) then ! Dirichlet BCs
                select case (g(2)%mode_fdm)
                case (FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA)
                    call INT_C2N6_LHS_E(ny, g(2)%jac, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2N6_RHS(ny, i2, g(2)%jac, p_wrk1d(1, 9), p_wrk1d(1, 11))

                case (FDM_COM6_DIRECT)
                    p_wrk1d(:, 1:7) = 0.0_wp
                    call INT_C2N6N_LHS_E(ny, g(2)%lu2(1, 8), g(2)%lu2(1, 4), lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2N6N_RHS(ny, i2, g(2)%lu2(1, 8), p_wrk1d(1, 9), p_wrk1d(1, 11))

                end select

                call PENTADFS(ny - 2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5))

                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 6))
                call PENTADSS(ny - 2, i1, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(2, 7))

                call PENTADSS(ny - 2, i2, p_wrk1d(2, 1), p_wrk1d(2, 2), p_wrk1d(2, 3), p_wrk1d(2, 4), p_wrk1d(2, 5), p_wrk1d(3, 11))

                c_wrk1d(:, 6) = c_wrk1d(:, 6) + bcs(1)*p_wrk1d(:, 6) + bcs(2)*p_wrk1d(:, 7)

                ! end if

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
            call OPR_FOURIER_B_Z_EXEC(c_tmp1, c_wrk3d)
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, a, c_tmp1)
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1, a, c_wrk3d)
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
        call c_f_pointer(c_loc(p_wrk1d(1,9)), aux_n, shape=[nfield, ny, 2]) ! lines of forcing and solution 

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        ! #######################################################################
        ! Fourier transform of forcing term; output of this section in array tmp1
        ! #######################################################################
        do ifield = 1, nfield
            if (g(3)%size > 1) then
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, &
                                          bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), c_tmp2) !, c_tmp1_n(:, :, ifield), c_wrk3d)
                call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1_n(:, :, ifield)) ! tmp2 might be overwritten
            else
                call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, &
                                          bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), c_tmp1_n(:, :, ifield)) !, c_tmp2, c_wrk3d)
            end if
        end do

        ! ###################################################################
        ! Solve FDE \hat{p}''-(\lambda+\alpha) \hat{p} = \hat{f}
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
                select case (g(2)%mode_fdm)
                case (FDM_COM6_JACOBIAN, FDM_COM6_JACPENTA)
                    call INT_C2N6_LHS_E(ny, g(2)%jac, lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2N6_RHS(ny, i2*nfield, g(2)%jac, p_wrk1d(1, 9), p_wrk1d(1, ip_sol))

                case (FDM_COM6_DIRECT)
                    p_wrk1d(:, 1:7) = 0.0_wp
                    call INT_C2N6N_LHS_E(ny, g(2)%lu2(1, 8), g(2)%lu2(1, 4), lambda, &
                            p_wrk1d(1, 1), p_wrk1d(1, 2), p_wrk1d(1, 3), p_wrk1d(1, 4), p_wrk1d(1, 5), p_wrk1d(1, 6), p_wrk1d(1, 7))
                    call INT_C2N6N_RHS(ny, i2*nfield, g(2)%lu2(1, 8), p_wrk1d(1, 9), p_wrk1d(1, ip_sol))

                end select

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
                call OPR_FOURIER_B_Z_EXEC(c_tmp1_n(:, :, ifield), c_wrk3d)
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_wrk3d, a(ifield)%field, c_tmp1_n(:, :, ifield))
            else
                call OPR_FOURIER_B_X_EXEC(nx, ny, nz, c_tmp1_n(:, :, ifield), a(ifield)%field, c_wrk3d)
            end if
        end do

        nullify (c_tmp1_n, c_tmp2, aux_n)

        return
    end subroutine OPR_HELMHOLTZ_FXZ_D_N

end module OPR_ELLIPTIC
