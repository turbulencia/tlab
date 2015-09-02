!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2007/10/28 - J.P. Mellado
!#              Reorganized
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the fields rho(x,y), u(x,y) and v(x,y) s.t. the axial momentum flux
!# is conserved and the continuity equation is satisfied.
!#
!# It assumes passive incompressible mixing of the temperature.
!#
!# The inputs are rho_vi, u_vi and the output of the iteration is rho_vo.
!# Array u_vo and tem_vo are auxiliar.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE FLOW_JET_SPATIAL_DENSITY(imax, jmax, iprof_tem, thick_tem, delta_tem, mean_tem, &
     ycoor_tem, diam_tem, jet_tem, iprof_u, thick_u, delta_u, mean_u, ycoor_u, diam_u, &
     jet_u, scaley, x, y, z1, p, rho_vi, u_vi, tem_vi, rho_vo, u_vo, tem_vo, wrk1d)
  
  USE DNS_CONSTANTS, ONLY : wfile

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imax, jmax, iprof_tem, iprof_u
  TREAL delta_tem, mean_tem, thick_tem, ycoor_tem, diam_tem
  TREAL delta_u, mean_u, thick_u, ycoor_u, diam_u
  TREAL scaley, jet_u(*), jet_tem(*)
  TREAL x(imax)
  TREAL y(jmax)

  TREAL z1(*), p(*)
  TREAL rho_vi(jmax), u_vi(jmax), tem_vi(jmax)
  TREAL rho_vo(imax,jmax), u_vo(jmax)
  TREAL tem_vo(jmax)
  TREAL wrk1d(jmax,*)

! -------------------------------------------------------------------
  TREAL ycenter
  TREAL tol, err, param(2)
  TREAL FLOW_JET_TEMPORAL

  TINTEGER i, j, n, nmax, ier

! ###################################################################
  param = C_0_R

! Convergence parameters
  nmax = 30
  tol = C_1EM6_R

  tem_vi(1:jmax) = C_0_R
  ycenter = y(1) + scaley*ycoor_tem
  DO j = 1,jmax
     tem_vi(j) = FLOW_JET_TEMPORAL&
          (iprof_tem, thick_tem, delta_tem, mean_tem, diam_tem, ycenter, param, y(j))
  ENDDO

#define rho_aux(j) wrk1d(j,1)
#define aux(j)     wrk1d(j,2)
! ###################################################################
! Begining loop over x-planes
! ###################################################################
! The initial condition for rho_vo is rho_vi
  DO j = 1,jmax
     rho_aux(j) = rho_vi(j)
  ENDDO
  DO i = 1,imax
     ier = 1

! Begining iteration for a fixed x-plane
     DO n = 1,nmax
! Velocity profile. Array tem_vo used as auxiliar array
        ycenter = y(1) + scaley*ycoor_u
        CALL FLOW_JET_SPATIAL_VELOCITY&
             (i1, jmax, iprof_u, thick_u, delta_u, mean_u, diam_u, ycenter,&
             jet_u(1), jet_u(2), jet_u(3), x(i), y, &
             rho_vi, u_vi, rho_aux(1), u_vo, tem_vo, aux(1), wrk1d(1,3))
! Normalized temperature and density profiles 
        ycenter = y(1) + scaley*ycoor_tem
        CALL FLOW_JET_SPATIAL_SCALAR&
             (i1, jmax, iprof_tem, thick_tem, delta_tem, mean_tem, diam_tem, diam_u, ycenter,&
             jet_tem(1), jet_tem(2), jet_tem(3), x(i), y, &
             rho_vi, u_vi, tem_vi, rho_aux(1), u_vo, tem_vo, aux(1))
        CALL THERMO_THERMAL_DENSITY(i1, jmax, i1, z1, p, tem_vo, wrk1d(1,2))
! Convergence criteria (infinity norm)
        DO j = 1,jmax
           wrk1d(j,3) = ABS(wrk1d(j,2)-rho_aux(j))
           rho_aux(j) = wrk1d(j,2)
        ENDDO
        err = MAXVAL(wrk1d(1:jmax,3),jmax)/MAXVAL(wrk1d(1:jmax,2))
        IF ( err .LT. tol ) THEN
           ier = 0
           EXIT
        ENDIF
     ENDDO

! Final check
     IF ( ier .EQ. 1 ) THEN
        CALL IO_WRITE_ASCII(wfile, 'FLOW_JET_SPATIAL: nmax reached.')
     ENDIF
     DO j = 1,jmax
        rho_vo(i,j) = rho_aux(j)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FLOW_JET_SPATIAL_DENSITY
