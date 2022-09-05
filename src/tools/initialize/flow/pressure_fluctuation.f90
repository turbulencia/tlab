#include "dns_const.h"

!########################################################################
!#
!# solve Poisson equation for p',
!# nabla^2 p' = d/dx_i d/dx_j (rho_0 u_i u_j), assuming
!# p/\rho^\gamma0 constant
!#
!# Homentropic conditions
!#
!########################################################################
!# ARGUMENTS
!#
!# rho    In   Mean density field
!#        Out  Density mean+fluctuation field
!# p      In   Mean pressure field
!#        Out  Pressure mean+fluctuation field
!#
!########################################################################
subroutine PRESSURE_FLUCTUATION(u, v, w, rho, p, pprime, txc1, txc2, txc3, txc4, wrk1d, wrk2d, wrk3d)
    use TLAB_TYPES, only: cp, ci
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imax, jmax, kmax, isize_wrk1d
    use THERMO_VARS, only: gama0
    use FLOW_LOCAL, only: norm_ini_p

    implicit none

    real(cp), dimension(imax, jmax, kmax) :: u, v, w, rho, p, pprime
    real(cp), dimension(imax, jmax, kmax) :: txc1, txc2, txc3, txc4, wrk3d
    real(cp), dimension(imax, kmax, *) :: wrk2d
    real(cp), dimension(isize_wrk1d, *) :: wrk1d

! -------------------------------------------------------------------
    integer(ci) bcs(2, 2)

! ###################################################################

! -------------------------------------------------------------------
! Calculate RHS d/dx_i d/dx_j (u_i u_j), stored in txc4
! -------------------------------------------------------------------
! terms with u
    txc1 = rho*u*u; txc2 = rho*u*v; txc3 = rho*u*w
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc3, txc4, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc2, txc3, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc1, txc2, wrk3d, wrk2d, wrk3d)
    txc2 = 2.0_cp*(txc4 + txc3) + txc2

! rhs in txc4
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc2, txc4, wrk3d, wrk2d, wrk3d)

! terms with v
    txc1 = rho*v*v; txc2 = rho*v*w
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc2, txc3, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc1, txc2, wrk3d, wrk2d, wrk3d)
    txc2 = txc2 + 2.0_cp*txc3

! rhs in txc4
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc2, txc1, wrk3d, wrk2d, wrk3d)
    txc4 = txc4 + txc1

! terms with w
    txc1 = rho*w*w
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc1, txc2, wrk3d, wrk2d, wrk3d)

! rhs in txc4
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc2, txc1, wrk3d, wrk2d, wrk3d)
    txc4 = txc4 + txc1

! -------------------------------------------------------------------
! Solve Poisson equation
! Array pprime contains fluctuating p' (BCs are equal to zero!)
! -------------------------------------------------------------------
    if (g(1)%periodic .and. g(3)%periodic) then ! Doubly periodic in xOz
        wrk2d(:, :, 1:2) = 0.0_cp  ! bcs
        pprime = -txc4          ! change of forcing term sign
        call OPR_POISSON_FXZ(.false., imax, jmax, kmax, g, 0, &
                             pprime, wrk3d, txc1, txc2, wrk2d(1, 1, 1), wrk2d(1, 1, 2), wrk1d, wrk1d(1, 5), wrk3d)

    else                                      ! General treatment
#ifdef USE_CGLOC
! Need to define global variable with ipos,jpos,kpos,ci,cj,ck,
        call CGPOISSON(i1, imax, jmax, kmax, g(3)%size, pprime, txc4, txc3, txc2, ipos, jpos, kpos, ci, cj, ck, wrk2d)
#endif
    end if

! -------------------------------------------------------------------
! An amplification factor norm_ini_p is allowed as in previous versions
! -------------------------------------------------------------------
    rho = (norm_ini_p*pprime/p/gama0 + 1.0_cp)*rho
    p = norm_ini_p*pprime + p

    return
end subroutine PRESSURE_FLUCTUATION
