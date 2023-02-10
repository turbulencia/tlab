#include "dns_error.h"
#ifdef USE_MPI
#include "dns_error.h"
#endif

module OPR_ELLIPTIC
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: ivfilter, istagger, vfilter_param
    use TLAB_POINTERS_3D, only: p_wrk1d
    use TLAB_POINTERS_C, only: c_wrk1d, c_wrk3d
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

    public :: OPR_POISSON_FXZ
    ! NEED TO ADD HELMHOLTZ HERE
contains
!########################################################################
!#
!# Solve Lap p = a using Fourier in xOz planes, to rewritte the problem as
!#     \hat{p}''-\lambda \hat{p} = \hat{a}
!#
!# where \lambda = kx^2+kz^2
!#
!# The reference value of p at the lower boundary is set to zero
!#
!# The global variable isize_txc_field defines the size of array txc equal
!# to (imax+2)*(jmax+2)*kmax, or larger if PARALLEL mode
!#
!########################################################################
    subroutine OPR_POISSON_FXZ(nx, ny, nz, g, ibc, p, tmp1, tmp2, bcs_hb, bcs_ht, dpdy)
        integer(wi), intent(in) :: nx, ny, nz
        integer,     intent(in) :: ibc   ! BCs at j1/jmax:  0, for Dirichlet & Dirichlet
        !                                                   1, for Neumann   & Dirichlet
        !                                                   2, for Dirichlet & Neumann
        !                                                   3, for Neumann   & Neumann
        type(grid_dt), intent(in)    :: g(3)
        real(wp),      intent(inout) :: p(nx, ny, nz)                       ! Forcing term, and solution field p
        real(wp),      intent(inout) :: tmp1(isize_txc_dimz, nz), tmp2(isize_txc_dimz, nz)
        real(wp),      intent(in)    :: bcs_hb(nx, nz), bcs_ht(nx, nz)      ! Boundary-condition fields
        real(wp),      intent(out),  optional :: dpdy(nx, ny, nz)           ! Vertical derivative of solution

        target tmp1, tmp2

! -----------------------------------------------------------------------
        integer(wi) i, j, k, iglobal, kglobal, ip, isize_line
        integer(wi) i_sing(2), k_sing(2)    ! singular global modes
        real(wp) lambda, norm

        ! #######################################################################
        call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(tmp2), c_tmp2, shape=[isize_txc_dimz/2, nz])
        call c_f_pointer(c_loc(bcs), r_bcs, shape=[3*2])

        norm = 1.0_wp/real(g(1)%size*g(3)%size, wp)

        isize_line = nx/2 + 1

        if (istagger == 0) then
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
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp2, c_tmp1, c_wrk3d)
            call OPR_FOURIER_F_Z_EXEC(c_tmp2, c_tmp1) ! tmp2 might be overwritten
        else
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, p, bcs_hb, bcs_ht, c_tmp1, c_tmp2, c_wrk3d)
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
                case (3) ! Neumann   & Neumann   BCs
                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                        call FDE_BVP_SINGULAR_NN(g(2)%mode_fdm, ny, 2, &
                                                 g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    else
                        call FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    end if

                case (0) ! Dirichlet & Dirichlet BCs
                    if (any(i_sing == iglobal) .and. any(k_sing == kglobal)) then
                        call FDE_BVP_SINGULAR_DD(g(2)%mode_fdm, ny, 2, &
                                            g(2)%nodes, g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    else
                        call FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny, 2, lambda, &
                                                g(2)%jac, p_wrk1d(:, 3), p_wrk1d(:, 1), r_bcs, p_wrk1d(:, 5), p_wrk1d(:, 7))
                    end if

                end select

                ! Vertical filtering of p and dpdy in case of staggering
                if (ivfilter == 1) then
                    call FILTER_VERTICAL_PRESSURE(ny, vfilter_param, c_wrk1d(:, 2), c_wrk1d(:, 3))
                end if

                ! Rearrange in memory and normalize
                do j = 1, ny
                    ip = (j - 1)*isize_line + i
                    c_tmp1(ip, k) = c_wrk1d(j, 2)*norm ! solution
                    c_tmp2(ip, k) = c_wrk1d(j, 3)*norm ! Oy derivative
                end do

            end do
        end do

! Fourier field p (based on array tmp1)
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

end module OPR_ELLIPTIC
