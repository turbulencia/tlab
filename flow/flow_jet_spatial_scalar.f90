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
!# Calculate the fields z1(x,y) s.t. the axial scalar flux
!# is conserved. 
!#
!# In this subroutine, the density and velocity field are a given input.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE FLOW_JET_SPATIAL_SCALAR(imax, jmax, iprof_z, thick_z, delta_z, mean_z, &
     diam_z, diam_u, ycenter, jet_z_a, jet_z_b, jet_z_flux, &
     x, y, rho_vi, u_vi, z_vi, rho, u, z1, wrk1d)
      
#include "types.h"

  USE DNS_CONSTANTS, ONLY : wfile
  
  IMPLICIT NONE

#include "integers.h"

  TINTEGER imax, jmax, iprof_z
  TREAL delta_z, mean_z, thick_z, diam_u, diam_z, ycenter
  TREAL jet_z_a, jet_z_b, jet_z_flux
  TREAL x(imax)
  TREAL y(jmax)

  TREAL rho_vi(jmax), u_vi(jmax), z_vi(jmax)
  TREAL rho(imax,jmax), u(imax,jmax), z1(imax,jmax)
  TREAL wrk1d(jmax,*)

! -------------------------------------------------------------------
  TREAL c1, c2
  TREAL delta, eta, ExcMom_vi, Q1, Z2, ZC, flux_aux
  TREAL dummy, diam_loc, param(2)
  TREAL SIMPSON_NU, FLOW_JET_TEMPORAL
  TREAL xi_tr, dxi_tr
  TINTEGER i, j

! ###################################################################
  param = C_0_R

! Ramaprian85 for the scalar ?
#ifdef SINGLE_PREC
  c1 =-0.6749e+0 
  c2 = 0.027e+0
#else 
  c1 =-0.6749d+0
  c2 = 0.027d+0
#endif

  Z2 = mean_z - C_05_R*delta_z

! -------------------------------------------------------------------
! Transition as a tanh profile around xi_tr between (0,2xi_tr) 
! in the slope of the half width. 
! It is set s.t. the vorticity thickness associated with this 
! tanh is half of the distance between xi_tr and the estimated end of 
! transition, 2 xi_tr.
! -------------------------------------------------------------------
  xi_tr = C_05_R*diam_z/diam_u/jet_z_a - jet_z_b
  IF ( xi_tr .LT. C_0_R ) THEN
     CALL IO_WRITE_ASCII(wfile, 'FLOW_SPATIAL_JET_VELOCITY. xi_tr negative.')
  ENDIF
  dxi_tr=xi_tr/C_8_R

! -------------------------------------------------------------------
! Normalized scalar, f(eta)
! -------------------------------------------------------------------
  DO i = 1,imax
     delta = (dxi_tr*LOG(EXP((x(i)/diam_u-xi_tr)/dxi_tr)+C_1_R)*jet_z_a+C_05_R)*diam_u
     diam_loc = C_2_R*delta

! inflow profile
     wrk1d(1:jmax,1) = C_0_R
     DO j = 1,jmax
        wrk1d(j,1) = FLOW_JET_TEMPORAL&
             (iprof_z, thick_z, delta_z, mean_z, diam_loc, ycenter, param, y(j))
     ENDDO
     ZC = wrk1d(jmax/2,1)-Z2

! Z-Z2=f(y) reference profile, stored in array z1.
     DO j = 1,jmax
        eta = (y(j)-ycenter)/delta
        z1(i,j) = EXP( c1*eta**2*(C_1_R+c2*eta**4) )
        dummy = C_05_R*(C_1_R+TANH(C_05_R*(x(i)/diam_u-xi_tr)/dxi_tr))
        z1(i,j) = dummy*z1(i,j)*ZC + (C_1_R-dummy)*(wrk1d(j,1)-Z2)
     ENDDO
  ENDDO

! -------------------------------------------------------------------
! Magnitude ZC for conserving the axial flux (w/ correction jet_u_flux)
! -------------------------------------------------------------------
! Reference momentum excess at the inflow
  DO j = 1,jmax
     wrk1d(j,2) = rho_vi(j)*u_vi(j)*(z_vi(j)-Z2)
  ENDDO
  ExcMom_vi = SIMPSON_NU(jmax, wrk1d(1,2), y)

  DO i = 1,imax
! Correction factor varying between 1 at the inflow and jet_z_flux
! at the outflow
     dummy = C_05_R*(C_1_R+TANH(C_05_R*(x(i)/diam_u-xi_tr)/dxi_tr))
     flux_aux = dummy*jet_z_flux + (C_1_R-dummy)*C_1_R

! Calculating ZC. 
! Solve second order equation Q1*ZC-J=0, positive root.
     DO j = 1,jmax
        wrk1d(j,1) = rho(i,j)*u(i,j)*z1(i,j)
     ENDDO
     Q1 = SIMPSON_NU(jmax, wrk1d(1,1), y)
     ZC = flux_aux*ExcMom_vi/Q1
     DO j = 1,jmax
        z1(i,j) = Z2 + ZC*z1(i,j)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FLOW_JET_SPATIAL_SCALAR
