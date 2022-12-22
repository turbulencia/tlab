#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!#
!# Solve Lap p = a using Fourier in xOz planes, to rewritte the problem
!# as
!#     \hat{p}''-\lambda \hat{p} = \hat{a}
!#
!# where \lambda = kx^2+kz^2
!#
!# The reference value of p at the lower boundary is set to zero
!#
!########################################################################
!# ARGUMENTS
!#
!# ibc     In    BCs at j1/jmax: 0, for Dirichlet & Dirichlet
!#                               1, for Neumann   & Dirichlet
!#                               2, for Dirichlet & Neumann
!#                               3, for Neumann   & Neumann
!#
!# The global variable isize_txc_field defines the size of array txc equal
!# to (imax+2)*(jmax+2)*kmax, or larger if PARALLEL mode
!#
!########################################################################
subroutine OPR_POISSON_FXZ(flag, nx, ny, nz, g, ibc, &
                           a, dpdy, tmp1, tmp2, bcs_hb, bcs_ht, aux, wrk1d, wrk3d)

    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: ivfilter, istagger, vfilter_param
    use OPR_FOURIER
    use FDE_BVP
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

    implicit none

    logical, intent(IN) :: flag
    integer(wi), intent(IN) :: nx, ny, nz, ibc
    type(grid_dt), intent(IN) :: g(3)
    real(wp), dimension(nx, ny, nz), intent(INOUT) :: a    ! Forcing term, ans solution field p
    real(wp), dimension(nx, ny, nz), intent(INOUT) :: dpdy ! Derivative, flag .TRUE.
    real(wp), dimension(nx, nz), intent(IN) :: bcs_hb, bcs_ht   ! Boundary-condition fields
    complex(wp), dimension(isize_txc_dimz/2, nz), intent(INOUT) :: tmp1, tmp2, wrk3d
    complex(wp), dimension(ny, 2), intent(INOUT) :: aux
    complex(wp), dimension(ny, 7), intent(INOUT) :: wrk1d

    target aux, wrk1d

! -----------------------------------------------------------------------
    integer(wi) i, j, k, iglobal, kglobal, ip, isize_line
    real(wp) lambda, norm
    complex(wp), target :: bcs(3)

    real(wp), pointer :: r_wrk1d(:, :) => null(), r_aux(:, :) => null(), r_bcs(:) => null()

    ! #######################################################################
    call c_f_pointer(c_loc(wrk1d), r_wrk1d, shape=[ny*2, 7])
    call c_f_pointer(c_loc(aux), r_aux, shape=[ny*2, 2])
    call c_f_pointer(c_loc(bcs), r_bcs, shape=[3*2])

    norm = 1.0_wp/real(g(1)%size*g(3)%size,wp)

    isize_line = nx/2 + 1

! #######################################################################
! Fourier transform of forcing term; output of this section in array tmp1
! #######################################################################
    if (g(3)%size > 1) then
        call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, tmp2, tmp1, wrk3d)
        call OPR_FOURIER_F_Z_EXEC(tmp2, tmp1) ! tmp2 might be overwritten
    else
        call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a, bcs_hb, bcs_ht, tmp1, tmp2, wrk3d)
    end if

! ###################################################################
! Solve FDE \hat{p}''-\lambda \hat{p} = \hat{a}
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
            if (g(3)%size > 1) then; lambda = g(1)%mwn(iglobal, 1) + g(3)%mwn(kglobal, 1)
            else; lambda = g(1)%mwn(iglobal, 1); end if

! forcing term
            do j = 1, ny
                ip = (j - 1)*isize_line + i; aux(j, 1) = tmp1(ip, k)
            end do

! BCs
            j = ny + 1; ip = (j - 1)*isize_line + i; bcs(1) = tmp1(ip, k) ! Dirichlet or Neumann
            j = ny + 2; ip = (j - 1)*isize_line + i; bcs(2) = tmp1(ip, k) ! Dirichlet or Neumann

! -----------------------------------------------------------------------
! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
! -----------------------------------------------------------------------
            select case (ibc)

            case (3) ! Neumann   & Neumann   BCs
                if (istagger == 0) then
                    if (kglobal == 1 .and. (iglobal == 1 .or. iglobal == g(1)%size/2 + 1) .or. &
                        kglobal == g(3)%size/2 + 1 .and. (iglobal == 1 .or. iglobal == g(1)%size/2 + 1)) then
                        call FDE_BVP_SINGULAR_NN(g(2)%mode_fdm, ny, 2, &
                                                 g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 3))
                    else
                        call FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 2))
                    end if
                else ! In case of staggering only one singular mode + different modified wavenumbers
                    if (kglobal == 1 .and. iglobal == 1) then
                        call FDE_BVP_SINGULAR_NN(g(2)%mode_fdm, ny, 2, &
                                                 g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 3))
                    else
                        call FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 2))
                    end if
                end if

            case (0) ! Dirichlet & Dirichlet BCs
                if (istagger == 0) then
                    if (kglobal == 1 .and. (iglobal == 1 .or. iglobal == g(1)%size/2 + 1) .or. &
                        kglobal == g(3)%size/2 + 1 .and. (iglobal == 1 .or. iglobal == g(1)%size/2 + 1)) then
                        call FDE_BVP_SINGULAR_DD(g(2)%mode_fdm, ny, 2, &
                                                 g(2)%nodes, g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 3))
                    else
                        call FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 2))
                    end if
                else ! In case of staggering only one singular mode + different modified wavenumbers
                    if (kglobal == 1 .and. iglobal == 1) then
                        call FDE_BVP_SINGULAR_DD(g(2)%mode_fdm, ny, 2, &
                                                 g(2)%nodes, g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 3))
                    else
                        call FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 2))
                    end if
                end if

            end select

            ! Vertical filtering of p and dpdy in case of staggering
            if (ivfilter == 1) then
                call FILTER_VERTICAL_PRESSURE(aux(1, 2), wrk1d(1, 1), ny, vfilter_param, wrk1d(1, 2))
            end if

            ! Normalize
            do j = 1, ny
                ip = (j - 1)*isize_line + i
                tmp1(ip, k) = aux(j, 2)*norm ! solution
                tmp2(ip, k) = wrk1d(j, 1)*norm ! Oy derivative
            end do

        end do
    end do

! ###################################################################
! Fourier field p (based on array tmp1)
! ###################################################################
    if (g(3)%size > 1) then
        call OPR_FOURIER_B_Z_EXEC(tmp1, wrk3d)
        call OPR_FOURIER_B_X_EXEC(nx, ny, nz, wrk3d, a, tmp1)
    else
        call OPR_FOURIER_B_X_EXEC(nx, ny, nz, tmp1, a, wrk3d)
    end if

! ###################################################################
! Fourier derivatives (based on array tmp2)
! ###################################################################
    if (flag) then
        if (g(3)%size > 1) then
            call OPR_FOURIER_B_Z_EXEC(tmp2, wrk3d)
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, wrk3d, dpdy, tmp2)
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, tmp2, dpdy, wrk3d)
        end if
    end if

    return
end subroutine OPR_POISSON_FXZ
