#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

!########################################################################
!#
!# Solve Lap p + \alpha p = a using Fourier in xOz planes, to rewritte
!# the problem as
!#
!#      \hat{p}''-(\lambda-\alpha) \hat{p} = \hat{a}
!#
!# where \lambda = kx^2+kz^2
!#
!########################################################################
!# ARGUMENTS
!#
!# a       In    Forcing term
!#         Out   Solution field p
!# ibc     In    BCs at j1/jmax: 0, for Dirichlet & Dirichlet
!#                               1, for Neumann   & Dirichlet
!#                               2, for Dirichlet & Neumann
!#                               3, for Neumann   & Neumann
!# bcs_??  In    BCs at j1/jmax
!#
!# The global variable isize_txc_field defines the size of array txc equal
!# to (imax+2)*(jmax+2)*kmax, or larger if PARALLEL mode
!#
!########################################################################
subroutine OPR_HELMHOLTZ_FXZ(nx, ny, nz, g, ibc, alpha, &
                             a, tmp1, tmp2, bcs_hb, bcs_ht, aux, wrk1d, wrk3d)

    use TLAB_CONSTANTS, only: efile, wp, wi
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_PROCS
    use OPR_FOURIER
    use FDE_BVP
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

    implicit none

    integer(wi) nx, ny, nz, ibc
    type(grid_dt), intent(IN) :: g(3)
    real(wp) alpha
    real(wp), dimension(nx, ny, nz) :: a
    real(wp), dimension(nx, nz) :: bcs_hb, bcs_ht
    complex(wp), dimension(isize_txc_dimz/2, nz) :: tmp1, tmp2, wrk3d
    complex(wp), dimension(ny, 2) :: aux
    complex(wp), dimension(ny, 7) :: wrk1d

    target aux, wrk1d
! -----------------------------------------------------------------------
    integer(wi) i, j, k, iglobal, kglobal, ip, isize_line
    real(wp) lambda, norm
    complex(wp), target :: bcs(3)

    real(wp), pointer :: r_wrk1d(:,:) => null(), r_aux(:,:) => null(), r_bcs(:) => null()

    ! #######################################################################
    call c_f_pointer(c_loc(wrk1d), r_wrk1d, shape=[ny*2,7])
    call c_f_pointer(c_loc(aux), r_aux, shape=[ny*2,2])
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
! Solve FDE (\hat{p}')'-\lambda \hat{p} = \hat{a}
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

            lambda = lambda - alpha

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
            if (ibc == 3) then ! Neumman BCs
               call FDE_BVP_REGULAR_NN(g(2)%mode_fdm, ny, 2, lambda, g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 2))
            else if (ibc == 0) then ! Dirichlet BCs
               call FDE_BVP_REGULAR_DD(g(2)%mode_fdm, ny, 2, lambda, g(2)%jac, r_aux(:, 2), r_aux(:, 1), r_bcs, r_wrk1d(:, 1), r_wrk1d(:, 2))
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Undeveloped BCs.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

! normalize
            do j = 1, ny
                ip = (j - 1)*isize_line + i
                tmp1(ip, k) = aux(j, 2)*norm ! solution
!        tmp2(ip,k) = wrk1d(j,1)*norm ! Oy derivative; not used in this routine
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

    return
end subroutine OPR_HELMHOLTZ_FXZ

!########################################################################
!########################################################################
! Same, but using the direct mode of FDM
subroutine OPR_HELMHOLTZ_FXZ_2(nx, ny, nz, g, ibc, alpha, &
                               a, tmp1, tmp2, bcs_hb, bcs_ht, aux, wrk1d, wrk3d)

    use TLAB_CONSTANTS, only: efile, wp, wi
    use TLAB_TYPES, only: grid_dt
    use TLAB_VARS, only: isize_field, isize_txc_dimz
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use OPR_FOURIER

    implicit none

    integer(wi) nx, ny, nz, ibc
    real(wp) alpha
    type(grid_dt), intent(IN) :: g(3)
    real(wp), dimension(isize_field) :: a
    real(wp), dimension(nx, nz) :: bcs_hb, bcs_ht
    complex(wp), dimension(isize_txc_dimz/2, nz) :: tmp1, tmp2, wrk3d
    complex(wp), dimension(ny, 2) :: aux
    real(wp), dimension(ny, 7) :: wrk1d

! -----------------------------------------------------------------------
    integer(wi) i, j, k, iglobal, kglobal, ip, isize_line
    real(wp) lambda, norm
    complex(wp) bcs(3)
    
    integer, parameter :: i1 = 1, i2 =2

! #######################################################################
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
            if (g(3)%size > 1) then; lambda = g(1)%mwn(iglobal, 2) + g(3)%mwn(kglobal, 2)
            else; lambda = g(1)%mwn(iglobal, 2); end if

            lambda = lambda - alpha

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
            if (ibc == 0) then ! Dirichlet BCs
                if (g(2)%mode_fdm == FDM_COM6_JACOBIAN .or. g(2)%mode_fdm == FDM_COM6_JACPENTA) then
                    call INT_C2N6_LHS_E(ny, g(2)%jac, lambda, &
                                        wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
                    call INT_C2N6_RHS(ny, i2, g(2)%jac, aux(1, 1), aux(1, 2))
                else if (g(2)%mode_fdm == FDM_COM6_DIRECT) then
                    wrk1d = 0.0_wp
                    call INT_C2N6N_LHS_E(ny, g(2)%lu2(1, 8), g(2)%lu2(1, 4), lambda, &
                                         wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
                    call INT_C2N6N_RHS(ny, i2, g(2)%lu2(1, 8), aux(1, 1), aux(1, 2))
                end if

                call PENTADFS(ny - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))

                call PENTADSS(ny - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
                call PENTADSS(ny - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))

                call PENTADSS(ny - 2, i2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), aux(2, 2))

                aux(:, 2) = aux(:, 2) + bcs(1)*wrk1d(:, 6) + bcs(2)*wrk1d(:, 7)

            else
                call TLAB_WRITE_ASCII(efile, 'OPR_HELMHOLT_FXZ_2. Undeveloped BCs.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

! normalize
            do j = 1, ny
                ip = (j - 1)*isize_line + i
                tmp1(ip, k) = aux(j, 2)*norm ! solution
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

    return
end subroutine OPR_HELMHOLTZ_FXZ_2

!########################################################################
!########################################################################
! Same, but for n fields
subroutine OPR_HELMHOLTZ_FXZ_2_N(nx, ny, nz, nfield, ibc, alpha, &
                                 a, tmp1, tmp2, bcs_hb, bcs_ht, aux, wrk1d, wrk3d)

    use TLAB_CONSTANTS, only: efile, wp, wi
    use TLAB_TYPES, only: pointers_dt
    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: g
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use OPR_FOURIER

    implicit none

    integer(wi) ibc, nx, ny, nz, nfield
    real(wp) alpha
    type(pointers_dt), dimension(nfield) :: a
    real(wp), dimension(nx, nz, nfield) :: bcs_hb, bcs_ht
    complex(wp), dimension(isize_txc_dimz/2, nz, nfield) :: tmp1
    complex(wp), dimension(isize_txc_dimz/2, nz) :: tmp2, wrk3d
    complex(wp), dimension(nfield, ny, 2) :: aux
    real(wp), dimension(ny, 7) :: wrk1d
! -----------------------------------------------------------------------
    integer(wi) i, j, k, iglobal, kglobal, ip, isize_line, ifield
    real(wp) lambda, norm
    complex(wp) bcs(3, nfield)

    integer, parameter :: i1 = 1, i2 =2

! #######################################################################
    norm = 1.0_wp/real(g(1)%size*g(3)%size,wp)

    isize_line = nx/2 + 1

! #######################################################################
! Fourier transform of forcing term; output of this section in array tmp1
! #######################################################################
    do ifield = 1, nfield
        if (g(3)%size > 1) then
            call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, &
                                      bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), tmp2, tmp1(1, 1, ifield), wrk3d)
            call OPR_FOURIER_F_Z_EXEC(tmp2, tmp1(1, 1, ifield)) ! tmp2 might be overwritten
        else
 call OPR_FOURIER_F_X_EXEC(nx, ny, nz, a(ifield)%field, bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield), tmp1(1, 1, ifield), tmp2, wrk3d)
        end if
    end do

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
            if (g(3)%size > 1) then; lambda = g(1)%mwn(iglobal, 2) + g(3)%mwn(kglobal, 2)
            else; lambda = g(1)%mwn(iglobal, 2); end if

            lambda = lambda - alpha

! forcing term
            do ifield = 1, nfield
                do j = 1, ny
                    ip = (j - 1)*isize_line + i; aux(ifield, j, 1) = tmp1(ip, k, ifield)
                end do
! BCs
                j = ny + 1; ip = (j - 1)*isize_line + i; 
                bcs(ifield, 1) = tmp1(ip, k, ifield)   ! Dirichlet or Neumann
                j = ny + 2; ip = (j - 1)*isize_line + i; 
                bcs(ifield, 2) = tmp1(ip, k, ifield)   ! Dirichlet or Neumann
            end do

! -----------------------------------------------------------------------
! Solve for each (kx,kz) a system of 1 complex equation as 2 independent real equations
! -----------------------------------------------------------------------
            if (ibc == 0) then ! Dirichlet BCs
                if (g(2)%mode_fdm == FDM_COM6_JACOBIAN .or. g(2)%mode_fdm == FDM_COM6_JACPENTA) then
                    call INT_C2N6_LHS_E(ny, g(2)%jac, lambda, &
                                        wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
                    call INT_C2N6_RHS(ny, i2, g(2)%jac, aux(1, 1, 1), aux(1, 1, 2))
                else if (g(2)%mode_fdm == FDM_COM6_DIRECT) then
                    wrk1d = 0.0_wp
                    call INT_C2N6N_LHS_E(ny, g(2)%lu2(1, 8), g(2)%lu2(1, 4), lambda, &
                                         wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6), wrk1d(1, 7))
                    call INT_C2N6N_RHS(ny, i2, g(2)%lu2(1, 8), aux(1, 1, 1), aux(1, 1, 2))
                end if

                call PENTADFS(ny - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))

                call PENTADSS(ny - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
                call PENTADSS(ny - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))
                call PENTADSS(ny - 2, i2*nfield, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), aux(1, 2, 2))

                do ifield = 1, nfield
                    aux(ifield, 1:ny, 2) = aux(ifield, 1:ny, 2) &
                                           + bcs(ifield, 1)*wrk1d(1:ny, 6) + bcs(ifield, 2)*wrk1d(1:ny, 7)
                end do

            else
                call TLAB_WRITE_ASCII(efile, 'OPR_HELMHOLT_FXZ_2_N. Undeveloped BCs.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

! normalize
            do ifield = 1, nfield
                do j = 1, ny
                    ip = (j - 1)*isize_line + i
                    tmp1(ip, k, ifield) = aux(ifield, j, 2)*norm ! solution
                end do
            end do
        end do
    end do

! ###################################################################
! Fourier field p (based on array tmp1)
! ###################################################################
    do ifield = 1, nfield
        if (g(3)%size > 1) then
            call OPR_FOURIER_B_Z_EXEC(tmp1(1, 1, ifield), wrk3d)
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, wrk3d, a(ifield)%field, tmp1(1, 1, ifield))
        else
            call OPR_FOURIER_B_X_EXEC(nx, ny, nz, tmp1(1, 1, ifield), a(ifield)%field, wrk3d)
        end if
    end do

    return
end subroutine OPR_HELMHOLTZ_FXZ_2_N
