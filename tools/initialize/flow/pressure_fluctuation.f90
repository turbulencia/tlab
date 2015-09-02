!########################################################################
!# Tool/Library INIT/FLOW
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!# 2007/03/17 - J.P. Mellado
!#              Calculation of RHS re-written
!# 2007/04/30 - J.P. Mellado
!#              Homentropic formulation
!# 2007/08/01 - J.P. Mellado
!#              Back to pressure formulation
!#
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
#include "types.h"
#include "dns_const.h"

SUBROUTINE PRESSURE_FLUCTUATION(y,dx,dy,dz, u,v,w,rho,p,pprime, &
     txc1,txc2,txc3,txc4, ipos,jpos,kpos,ci,cj,ck, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL , ONLY : imode_fdm, imax,jmax,kmax,kmax_total, isize_field, i1bc,j1bc,k1bc, isize_wrk1d
  USE THERMO_GLOBAL , ONLY : gama0
  USE FLOW_LOCAL, ONLY : norm_ini_p
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(*)              :: y, dx, dy, dz
  TREAL, DIMENSION(imax,jmax,kmax) :: u, v, w, rho, p, pprime
  TREAL, DIMENSION(imax,jmax,kmax) :: txc1, txc2, txc3, txc4, wrk3d
  TREAL, DIMENSION(*)              :: ipos, jpos, kpos, ci, cj, ck
  TREAL, DIMENSION(imax,kmax,*)    :: wrk2d
  TREAL, DIMENSION(isize_wrk1d,*)  :: wrk1d

! -------------------------------------------------------------------

! ###################################################################

! -------------------------------------------------------------------
! Calculate RHS d/dx_i d/dx_j (u_i u_j), stored in txc4
! -------------------------------------------------------------------
! terms with u
  txc1 = rho*u*u; txc2 = rho*u*v; txc3 = rho*u*w
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, txc3, txc4, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, txc2, txc3, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, txc1, txc2, i0,i0, wrk1d,wrk2d,wrk3d)
  txc2 = C_2_R*( txc4 + txc3 ) + txc2

! rhs in txc4
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, txc2, txc4, i0,i0, wrk1d,wrk2d,wrk3d)

! terms with v
  txc1 = rho*v*v; txc2 = rho*v*w
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, txc2, txc3, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, txc1, txc2, i0,i0, wrk1d,wrk2d,wrk3d)
  txc2 = txc2 + C_2_R*txc3

! rhs in txc4
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, txc2, txc1, i0,i0, wrk1d,wrk2d,wrk3d)
  txc4 = txc4 + txc1

! terms with w
  txc1 = rho*w*w
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, txc1, txc2, i0,i0, wrk1d,wrk2d,wrk3d)

! rhs in txc4
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, txc2, txc1, i0,i0, wrk1d,wrk2d,wrk3d)
  txc4 = txc4 + txc1

! -------------------------------------------------------------------
! Solve Poisson equation
! Array pprime contains fluctuating p' (BCs are equal to zero!)
! -------------------------------------------------------------------
  IF ( i1bc .EQ. 0 .AND. k1bc .EQ. 0 ) THEN ! Doubly periodic in xOz
     wrk2d(:,:,1:2) = C_0_R  ! bcs
     pprime = -txc4          ! change of forcing term sign
     CALL OPR_POISSON_FXZ(imode_fdm,i1,i0, imax,jmax,kmax, &
          y,dx,dy,dz, pprime,wrk3d, txc1,txc2, &
          wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

  ELSE                                      ! General treatment
#ifdef USE_CGLOC
     CALL CGPOISSON(i1, imax,jmax,kmax,kmax_total, i1bc,j1bc,k1bc, &
          dx,dy,dz, pprime, txc4,txc3,txc2, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
  ENDIF


! -------------------------------------------------------------------
! An amplification factor norm_ini_p is allowed as in previous versions
! -------------------------------------------------------------------
  rho = (norm_ini_p*pprime/p/gama0 + C_1_R)*rho
  p   = norm_ini_p*pprime + p

  RETURN
END SUBROUTINE PRESSURE_FLUCTUATION
