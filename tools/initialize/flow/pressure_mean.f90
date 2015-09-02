#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library FLOW/INIT
!#
!########################################################################
!# HISTORY
!#
!# 2003/01/01 - J.P. Mellado
!#              Created
!# 2007/03/19 - J.P. Mellado
!#              Modified to include all cases, not only RTI.
!#
!########################################################################
!# DESCRIPTION
!#
!# Equilibrium solution for the case of a body force acting on y coordinate.
!# Solving the equation
!#    \partial [p(x,y,z) - p_m(x,y,z)]/\partial y = rho(x,y,z)*g
!#
!# If volumetric forces are present, two types are possible:
!# 1) T/W is given, which leads to an exponential variation of p.
!# 2) rho is given, which leads to a linear variation of p.
!#
!# AIRWATER mode, needs T and q_i, to avoid recalculation
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE PRESSURE_MEAN(y,dy, p,T,s, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*)                :: y,dy, wrk2d
  TREAL, DIMENSION(imax,jmax,kmax)   :: p,T, wrk3d
  TREAL, DIMENSION(imax,jmax,kmax,*) :: s
  TREAL, DIMENSION(jmax,8)           :: wrk1d

! -------------------------------------------------------------------
  TINTEGER j, iprof_loc
  TREAL pmin,pmax, ycenter
  TREAL FLOW_SHEAR_TEMPORAL

! ###################################################################
! Constant pressure
! ###################################################################
  IF ( ibodyforce .EQ. EQNS_NONE ) THEN
     p = p_init

! ###################################################################
! Hydrostatic equilibrium
! ###################################################################
  ELSE

#define p_loc(i)       wrk1d(i,1)
#define r_loc(i)       wrk1d(i,2)
#define t_loc(i)       wrk1d(i,3)
#define e_loc(i)       wrk1d(i,4)
#define h_loc(i)       wrk1d(i,5)
#define z1_loc(i)      wrk1d(i,6)
#define z2_loc(i)      wrk1d(i,7)
#define wrk1d_loc(i)   wrk1d(i,8)

! -------------------------------------------------------------------
! Temperature profile is given
! -------------------------------------------------------------------
     IF ( iprof_rho .EQ. PROFILE_NONE ) THEN

! AIRWATER case: temperature/mixture profile is given
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. iprof_tem .GT. 0 ) THEN
           DO j = 1,jmax
              ycenter = y(1) + scaley*ycoor_tem
              t_loc(j) = FLOW_SHEAR_TEMPORAL&
                   (iprof_tem, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y(j))

              ycenter = y(1) + scaley*ycoor_i(1)
              z1_loc(j) =  FLOW_SHEAR_TEMPORAL&
                   (iprof_i(1), thick_i(1), delta_i(1), mean_i(1), ycenter, prof_i, y(j))

           ENDDO
           CALL FI_HYDROSTATIC_AIRWATER_T&
                (y, dy, z1_loc(1), t_loc(1), p_loc(1), r_loc(1), wrk1d_loc(1), wrk2d, wrk3d)
           DO j = 1,jmax
              s(:,j,:,1) = z1_loc(j)
              s(:,j,:,2) = z2_loc(j)
              T(:,j,:)   =  t_loc(j)
           ENDDO

! AIRWATER case: enthalpy/mixture profile is given
        ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. iprof_tem .LT. 0 ) THEN
           DO j = 1,jmax
              ycenter = y(1) + scaley*ycoor_tem
              iprof_loc =-iprof_tem
              h_loc(j) = FLOW_SHEAR_TEMPORAL&
                   (iprof_loc, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y(j))

              ycenter = y(1) + scaley*ycoor_i(1)
              z1_loc(j) =  FLOW_SHEAR_TEMPORAL&
                   (iprof_i(1), thick_i(1), delta_i(1), mean_i(1), ycenter, prof_i, y(j))

           ENDDO
           CALL FI_HYDROSTATIC_AIRWATER_H(jmax, y, z1_loc(1), h_loc(1), t_loc(1), p_loc(1), wrk1d_loc(1))
           DO j = 1,jmax
              s(:,j,:,1) = z1_loc(j)
              s(:,j,:,2) = z2_loc(j)
              T(:,j,:)   =  t_loc(j)
           ENDDO

! General case: temperature/mixture profile is given
        ELSE
           ycenter = y(1) + ycoor_tem*scaley
           CALL FI_HYDROSTATIC(i1, jmax, i1, ycenter, y, p_loc(1))
           DO j = 1,jmax
              p_loc(j) = p_init*EXP(p_loc(j))
           ENDDO

        ENDIF

! -------------------------------------------------------------------
! Density profile is given
! -------------------------------------------------------------------
     ELSE
        CALL IO_WRITE_ASCII(efile, 'PRESSURE_MEAN. Density case undeveloped')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

! -------------------------------------------------------------------
! 3D array. Simple case of g parallel and opposite to OY 
! -------------------------------------------------------------------
     DO j = 1,jmax
        p(:,j,:) = p_loc(j)
     ENDDO

  ENDIF

! ###################################################################
! Control
! ###################################################################
  CALL MINMAX(imax,jmax,kmax, p, pmin,pmax)

  IF ( pmin .LT. C_0_R .OR. pmax .LT. C_0_R ) THEN
     CALL IO_WRITE_ASCII(efile, 'PRESSURE_MEAN. Negative pressure.')
     CALL DNS_STOP(DNS_ERROR_NEGPRESS)
  ENDIF

  RETURN
END SUBROUTINE PRESSURE_MEAN

         
      
