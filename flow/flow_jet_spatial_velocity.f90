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
!# Calculate the fields u(x,y) and v(x,y) s.t. the axial momentum flux
!# is conserved and the continuity equation is satisfied.
!#
!# In this subroutine, the density field is a given input.
!#
!# The calculation of v assumes ycoor_u equal to 0.5. This section
!# should also be rewritten in terms of OPR_PARTIAL_ and QUAD routines.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE FLOW_JET_SPATIAL_VELOCITY&
     (imax, jmax, iprof_u, thick_u, delta_u, mean_u, diam_u, ycenter,&
     jet_u_a, jet_u_b, jet_u_flux, x, y, rho_vi, u_vi, rho, u, v, wrk1d, wrk2d)

#include "types.h"

  USE DNS_CONSTANTS, ONLY : efile, wfile

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imax, jmax, iprof_u
  TREAL delta_u, mean_u, thick_u, diam_u, ycenter
  TREAL jet_u_a, jet_u_b, jet_u_flux
  TREAL x(imax)
  TREAL y(jmax)

  TREAL rho_vi(jmax), u_vi(jmax)
  TREAL rho(imax,jmax), u(imax,jmax), v(imax,jmax)
  TREAL wrk1d(jmax,*), wrk2d(imax,jmax)

! -------------------------------------------------------------------
  TREAL c1, c2
  TREAL delta, eta, ExcMom_vi, Q1, Q2, U2, UC
  TREAL dummy, param(2), flux_aux, diam_loc
  TREAL SIMPSON_NU, FLOW_JET_TEMPORAL
  TREAL xi_tr, dxi_tr

  TINTEGER i, j, jsym

! ###################################################################
  param = C_0_R

! Bradbury profile
#ifdef SINGLE_PREC
  c1 =-0.6749e+0
  c2 = 0.027e+0
#else 
  c1 =-0.6749d+0
  c2 = 0.027d+0
#endif

  U2 = mean_u - C_05_R*delta_u

! ###################################################################
! Axial velocity, U_c*f(eta)
! ###################################################################
! -------------------------------------------------------------------
! Transition as a tanh profile around xi_tr between (0,2xi_tr) 
! in the slope of the half width. 
! It is set s.t. the vorticity thickness associated with this 
! tanh is half of the distance between xi_tr and the estimated end of 
! transition, 2 xi_tr.
! -------------------------------------------------------------------
  xi_tr = C_05_R/jet_u_a - jet_u_b
  IF ( xi_tr .LT. C_0_R ) THEN
     CALL IO_WRITE_ASCII(wfile, 'FLOW_SPATIAL_JET_VELOCITY. xi_tr negative.')
  ENDIF
  dxi_tr=xi_tr/C_8_R

! -------------------------------------------------------------------
! Normalized velocity, f(eta)
! -------------------------------------------------------------------
  DO i = 1,imax
     delta = (dxi_tr*LOG(EXP((x(i)/diam_u-xi_tr)/dxi_tr)+C_1_R)*jet_u_a+C_05_R)*diam_u
     diam_loc = C_2_R*delta

! inflow profile
     wrk1d(1:jmax,1) = C_0_R
     DO j = 1,jmax
        wrk1d(j,1) = FLOW_JET_TEMPORAL&
             (iprof_u, thick_u, delta_u, mean_u, diam_loc, ycenter, param, y(j))
     ENDDO
     UC = wrk1d(jmax/2,1)-U2

! U-U2=f(y) reference profile, stored in array u.
     DO j = 1,jmax
        eta = (y(j)-ycenter)/delta
        u(i,j) = EXP( c1*eta**2*(C_1_R+c2*eta**4) )
        dummy = C_05_R*(C_1_R+TANH(C_05_R*(x(i)/diam_u-xi_tr)/dxi_tr))
        u(i,j) = dummy*u(i,j)*UC + (C_1_R-dummy)*(wrk1d(j,1)-U2)
     ENDDO

  ENDDO

! -------------------------------------------------------------------
! Magnitude UC for conserving the axial flux (w/ correction jet_u_flux)
! -------------------------------------------------------------------
! Reference momentum excess at the inflow
  DO j = 1,jmax
     wrk1d(j,1) = rho_vi(j)*u_vi(j)*(u_vi(j)-U2)
  ENDDO
  ExcMom_vi = SIMPSON_NU(jmax, wrk1d, y)

  DO i = 1,imax
! Correction factor varying between 1 at the inflow and jet_u_flux
! at the outflow
     dummy = C_05_R*(C_1_R+TANH(C_05_R*(x(i)/diam_u-xi_tr)/dxi_tr))
     flux_aux = dummy*jet_u_flux + (C_1_R-dummy)*C_1_R

! Calculating UC. 
! Solve second order equation Q1*UC^2+Q2*UC-J=0, positive root.
     DO j = 1,jmax
        wrk1d(j,1) = rho(i,j)*u(i,j)*u(i,j)
        wrk1d(j,2) = rho(i,j)*u(i,j)
     ENDDO
     Q1 =    SIMPSON_NU(jmax, wrk1d(1,1), y)
     Q2 = U2*SIMPSON_NU(jmax, wrk1d(1,2), y)
     UC = (-Q2+SQRT(Q2*Q2+C_4_R*Q1*ExcMom_vi*flux_aux))/C_2_R/Q1

! Scaled velocity
     DO j = 1,jmax
        u(i,j) = U2 + UC*u(i,j)
     ENDDO
  ENDDO

! ###################################################################
! Lateral velocity, d(rho*v)/dy=-d(rho*u)/dx
! ###################################################################
  IF ( imax .GT. 1 ) THEN
! Backwards 1st-order derivative, -d(rho*u)/dx (w used as aux array)
     DO i = 2,imax
        DO j = 1,jmax
           Q1 = rho(i,  j)*u(i,  j)
           Q2 = rho(i-1,j)*u(i-1,j)
           wrk2d(i,j) =-(Q1-Q2)/(x(i)-x(i-1))
        ENDDO
     ENDDO
! Backwards 1st-order derivative for i=1
     i = 1
     DO j = 1,jmax
        wrk2d(i,j) = wrk2d(i+1,j)
     ENDDO

! Midpoint integration, rho*v
     DO i = 1,imax
        Q1 = wrk2d(i,jmax/2) + wrk2d(i,jmax/2+1)
        Q2 = y(jmax/2+1)-y(jmax/2)
        v(i,jmax/2+1) = C_05_R*( C_05_R*Q1*Q2 )

        DO j = jmax/2+2, jmax
           Q1 = wrk2d(i,j) + wrk2d(i,j-1)
           Q2 = y(j)-y(j-1)
           v(i,j) = v(i,j-1) + C_05_R*Q1*Q2
        ENDDO
     ENDDO

! Division by rho and antisymmetric extension
     DO i = 1,imax
        DO j = jmax/2+1,jmax
           v(i,j) = v(i,j)/rho(i,j)
           jsym = jmax - j + 1
           v(i,jsym) =-v(i,j)
        ENDDO
     ENDDO

! If only one plane, set lateral velocity to zero
  ELSE
     DO j = 1,jmax
        v(1,j) = C_0_R
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE FLOW_JET_SPATIAL_VELOCITY
