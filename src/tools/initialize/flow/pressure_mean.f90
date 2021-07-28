#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE PRESSURE_MEAN(p,T,s, wrk1d,wrk2d,wrk3d)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : imax,jmax,kmax
  USE TLAB_VARS,    ONLY : rbg, pbg, tbg, sbg
  USE TLAB_VARS,    ONLY : buoyancy
  USE TLAB_PROCS
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(OUT)   :: p
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: T
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(INOUT) :: s
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: wrk1d, wrk2d, wrk3d

! -------------------------------------------------------------------
  TINTEGER j, iprof_loc
  TREAL pmin,pmax, ycenter
  TREAL PROFILES

  TREAL, DIMENSION(:), POINTER :: y,dy

! ###################################################################
! Define pointers
  y => g(2)%nodes; dy => g(2)%jac(:,1)

! ###################################################################
! Constant pressure
! ###################################################################
  IF ( buoyancy%type .EQ. EQNS_NONE ) THEN
     p = pbg%mean

! ###################################################################
! Hydrostatic equilibrium
! ###################################################################
  ELSE

#define p_loc(i)       wrk1d(i,1)
#define r_loc(i)       wrk1d(i,2)
#define t_loc(i)       wrk1d(i,3)
#define ep_loc(i)      wrk1d(i,4)
#define h_loc(i)       wrk1d(i,5)
#define z1_loc(i)      wrk1d(i,6)
#define z2_loc(i)      wrk1d(i,7)
#define z3_loc(i)      wrk1d(i,8)
#define wrk1d_loc(i)   wrk1d(i,9)

! -------------------------------------------------------------------
! Temperature profile is given
! -------------------------------------------------------------------
     IF ( rbg%type .EQ. PROFILE_NONE ) THEN

! AIRWATER case: temperature/mixture profile is given
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. tbg%type .GT. 0 ) THEN
           DO j = 1,jmax
              ycenter = y(1) + g(2)%scale *tbg%ymean
              t_loc(j) = PROFILES(tbg%type, tbg%thick, tbg%delta, tbg%mean, ycenter, tbg%parameters, y(j))

              ycenter = y(1) + g(2)%scale *sbg(1)%ymean
              z1_loc(j) = PROFILES(sbg(1)%type, sbg(1)%thick, sbg(1)%delta, sbg(1)%mean, ycenter, sbg(1)%parameters, g(2)%nodes(j))

           ENDDO
           ! CALL FI_HYDROSTATIC_AIRWATER_T&
           !      (y, dy, z1_loc(1), t_loc(1), p_loc(1), r_loc(1), wrk1d_loc(1), wrk2d, wrk3d)
           CALL TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Hydrostatic equilibrium 1 undeveloped')
           CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
           DO j = 1,jmax
              s(:,j,:,1) = z1_loc(j)
              s(:,j,:,2) = z2_loc(j)
              T(:,j,:)   =  t_loc(j)
           ENDDO

! AIRWATER case: enthalpy/mixture profile is given
        ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. tbg%type .LT. 0 ) THEN
           DO j = 1,jmax
              ycenter = y(1) + g(2)%scale *tbg%ymean
              iprof_loc =-tbg%type
              z1_loc(j) = PROFILES(iprof_loc, tbg%thick, tbg%delta, tbg%mean, ycenter, tbg%parameters, y(j))

              ycenter = y(1) + g(2)%scale *sbg(1)%ymean
              z2_loc(j) = PROFILES(sbg(1)%type, sbg(1)%thick, sbg(1)%delta, sbg(1)%mean, ycenter, sbg(1)%parameters, g(2)%nodes(j))

           ENDDO
!           CALL FI_HYDROSTATIC_H_OLD(jmax, y, z1_loc(1), ep_loc(1), t_loc(1), p_loc(1), wrk1d_loc(1))
           CALL FI_HYDROSTATIC_H(g(2), z1_loc(1), ep_loc(1), t_loc(1), p_loc(1), wrk1d_loc(1))
           DO j = 1,jmax
              s(:,j,:,1) = z2_loc(j)
              s(:,j,:,2) = z3_loc(j)
              T(:,j,:)   =  t_loc(j)
           ENDDO

! General case: temperature/mixture profile is given
        ELSE
           ycenter = y(1) + tbg%ymean*g(2)%scale
!           CALL FI_HYDROSTATIC(i1, jmax, i1, ycenter, y, p_loc(1))
           CALL TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Hydrostatic equilibrium 2 undeveloped')
           CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)

           DO j = 1,jmax
              p_loc(j) = pbg%mean*EXP(p_loc(j))
           ENDDO

        ENDIF

! -------------------------------------------------------------------
! Density profile is given
! -------------------------------------------------------------------
     ELSE
        CALL TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Density case undeveloped')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
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
     CALL TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Negative pressure.')
     CALL TLAB_STOP(DNS_ERROR_NEGPRESS)
  ENDIF

  RETURN
END SUBROUTINE PRESSURE_MEAN
