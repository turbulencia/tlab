#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Setting up a perturbation of the thermodynamic fields by a
!# displacement of the reference center plane.
!#
!# Array s enters with the scalar total field, including fluctuations.
!#
!########################################################################
SUBROUTINE DENSITY_FLUCTUATION(code, s, p, rho, T, h, disp, wrk3d)

  USE DNS_GLOBAL,    ONLY : rbg, tbg
  USE THERMO_GLOBAL, ONLY : imixture
  USE FLOW_LOCAL

  IMPLICIT NONE

  TINTEGER code

  TREAL, DIMENSION(imax,jmax,kmax)   :: T, h, rho, p, wrk3d
  TREAL, DIMENSION(imax,jmax,kmax,*) :: s
  TREAL, DIMENSION(imax,kmax)        :: disp

  ! -------------------------------------------------------------------
  TINTEGER idummy, iprof_loc
  TREAL dummy, ycenter, mean_loc, delta_loc
  TREAL AVG1V2D, PROFILES
  TREAL xcenter, amplify

  TREAL, DIMENSION(:), POINTER :: x,y,z

  ! ###################################################################
  ! Define pointers
  x => g(1)%nodes
  y => g(2)%nodes
  z => g(3)%nodes

#ifdef USE_MPI
  idsp = ims_offset_i; kdsp = ims_offset_k
#else
  idsp = 0; kdsp = 0
#endif

  ! ###################################################################
  ! Center plane displacement
  ! ###################################################################
  disp = C_0_R

  ! -------------------------------------------------------------------
  ! Broadband case
  ! -------------------------------------------------------------------
  IF ( code .EQ. 4 ) THEN
    idummy = g(2)%size; g(2)%size = 1
    CALL DNS_READ_FIELDS('scal.rand', i1, imax,i1,kmax, i1,i0, isize_field, disp, wrk3d)
    g(2)%size = idummy
    dummy = AVG1V2D(imax,i1,kmax, i1, i1, disp)     ! remove mean
    disp = disp - dummy

  ! -------------------------------------------------------------------
  ! Discrete case
  ! -------------------------------------------------------------------
  ELSE IF ( code .EQ. 5 ) THEN
    wx_1 = C_2_R * C_PI_R / g(1)%scale ! Fundamental wavelengths
    wz_1 = C_2_R * C_PI_R / g(3)%scale

    DO im = 1,fp%size
      wx = M_REAL( fp%modex(im) ) *wx_1
      wz = M_REAL( fp%modez(im) ) *wz_1

      DO k = 1,kmax
        disp(:,k) = disp(:,k) + fp%amplitude(im) *COS( wx *x(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *z(kdsp+k) +fp%phasez(im) )
      ENDDO

    ENDDO

  ENDIF

  ! -------------------------------------------------------------------
  ! Modulation
  ! -------------------------------------------------------------------
  IF ( fp%type .EQ. 1 .AND. fp%parameters(1) .GT. C_0_R ) THEN
    DO k = 1,kmax
      DO i = 1,imax
        xcenter   = x(i) - g(1)%scale *C_05_R - x(1)
        amplify   = EXP(-C_05_R*(xcenter/fp%parameters(1))**2)
        disp(i,k) = disp(i,k)*amplify
      ENDDO
    ENDDO
  ENDIF

  ! ###################################################################
  ! Perturbation in the thermodynamic fields
  ! ###################################################################
  IF ( rbg%type .EQ. PROFILE_NONE ) THEN

    IF ( tbg%type .GT. 0 ) THEN ! temperature/mixture profile is given
      DO k = 1,kmax
        DO i = 1,imax
          delta_loc = tbg%delta + (tbg%parameters(2)-tbg%parameters(1))*disp(i,k) *g(2)%scale
          mean_loc  = tbg%mean  + C_05_R*(tbg%parameters(2)+tbg%parameters(1))*disp(i,k) *g(2)%scale
          ycenter   = y(1) + g(2)%scale *tbg%ymean + disp(i,k)
          DO j = 1,jmax
            T(i,j,k) =  PROFILES(tbg%type, tbg%thick, delta_loc, mean_loc, ycenter, tbg%parameters, y(j))
          ENDDO
        ENDDO
      ENDDO

      IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL THERMO_AIRWATER_PT(imax, jmax, kmax, s, p, T)
      ENDIF

    ELSE IF ( tbg%type .LT. 0 ) THEN ! enthalpy/mixture profile is given
      DO k = 1,kmax
        DO i = 1,imax
          delta_loc = tbg%delta + (tbg%parameters(2)-tbg%parameters(1))*disp(i,k) *g(2)%scale
          mean_loc  = tbg%mean  + C_05_R*(tbg%parameters(2)+tbg%parameters(1))*disp(i,k) *g(2)%scale
          ycenter   = y(1) + g(2)%scale *tbg%ymean + disp(i,k)
          iprof_loc =-tbg%type
          DO j = 1,jmax
            h(i,j,k) =  PROFILES(iprof_loc, tbg%thick, delta_loc, mean_loc, ycenter, tbg%parameters, y(j))
          ENDDO
        ENDDO
      ENDDO

      IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL THERMO_AIRWATER_PH_RE(imax,jmax,kmax, s, p, h, T)
      ENDIF

    ENDIF

    ! compute perturbation in density
    CALL THERMO_THERMAL_DENSITY(imax,jmax,kmax, s, p, T, rho)

  ELSE ! Defined in terms of the density, to be developed

  ENDIF

  RETURN
END SUBROUTINE DENSITY_FLUCTUATION
