#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
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
SUBROUTINE PRESSURE_FLUCTUATION(u,v,w,rho,p,pprime, txc1,txc2,txc3,txc4, wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, isize_wrk1d
  USE THERMO_VARS, ONLY : gama0
  USE FLOW_LOCAL,    ONLY : norm_ini_p

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax) :: u,v,w,rho, p,pprime
  TREAL, DIMENSION(imax,jmax,kmax) :: txc1,txc2,txc3,txc4, wrk3d
  TREAL, DIMENSION(imax,kmax,*)    :: wrk2d
  TREAL, DIMENSION(isize_wrk1d,*)  :: wrk1d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  
! -------------------------------------------------------------------
! Calculate RHS d/dx_i d/dx_j (u_i u_j), stored in txc4
! -------------------------------------------------------------------
! terms with u
  txc1 = rho*u*u; txc2 = rho*u*v; txc3 = rho*u*w
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc3, txc4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc2, txc3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), txc1, txc2, wrk3d, wrk2d,wrk3d)
  txc2 = C_2_R*( txc4 + txc3 ) + txc2

! rhs in txc4
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), txc2, txc4, wrk3d, wrk2d,wrk3d)

! terms with v
  txc1 = rho*v*v; txc2 = rho*v*w
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc2, txc3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc1, txc2, wrk3d, wrk2d,wrk3d)
  txc2 = txc2 + C_2_R*txc3

! rhs in txc4
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), txc2, txc1, wrk3d, wrk2d,wrk3d)
  txc4 = txc4 + txc1

! terms with w
  txc1 = rho*w*w
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc1, txc2, wrk3d, wrk2d,wrk3d)

! rhs in txc4
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), txc2, txc1, wrk3d, wrk2d,wrk3d)
  txc4 = txc4 + txc1

! -------------------------------------------------------------------
! Solve Poisson equation
! Array pprime contains fluctuating p' (BCs are equal to zero!)
! -------------------------------------------------------------------
  IF ( g(1)%periodic .AND. g(3)%periodic ) THEN ! Doubly periodic in xOz
     wrk2d(:,:,1:2) = C_0_R  ! bcs
     pprime = -txc4          ! change of forcing term sign
     CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, i0, &
          pprime,wrk3d, txc1,txc2, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  ELSE                                      ! General treatment
#ifdef USE_CGLOC
! Need to define global variable with ipos,jpos,kpos,ci,cj,ck,
     CALL CGPOISSON(i1, imax,jmax,kmax, g(3)%size, pprime, txc4,txc3,txc2, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
  ENDIF


! -------------------------------------------------------------------
! An amplification factor norm_ini_p is allowed as in previous versions
! -------------------------------------------------------------------
  rho = (norm_ini_p*pprime/p/gama0 + C_1_R)*rho
  p   = norm_ini_p*pprime + p

  RETURN
END SUBROUTINE PRESSURE_FLUCTUATION
