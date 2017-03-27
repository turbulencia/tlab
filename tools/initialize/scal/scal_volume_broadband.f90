#include "types.h"
#include "dns_const.h"

SUBROUTINE SCAL_VOLUME_BROADBAND(is, s, tmp, wrk3d)

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, imode_flow, inb_scal
  USE DNS_GLOBAL, ONLY : sbg
  USE SCAL_LOCAL, ONLY : thick_ini, ycoor_ini, norm_ini_s
  IMPLICIT NONE

#include "integers.h"

  TINTEGER is
  TREAL, DIMENSION(imax,jmax,kmax) :: s, tmp, wrk3d

! -------------------------------------------------------------------
  TREAL AVG1V2D
  TREAL ycenter, amplify, dummy
  TINTEGER j

! ###################################################################
! Read initial random field
! ###################################################################
  CALL DNS_READ_FIELDS('scal.rand', i1, imax,jmax,kmax, inb_scal,is, isize_field, tmp, wrk3d)

! Remove mean
  DO j = 1,jmax
     dummy = AVG1V2D(imax,jmax,kmax, j, i1, tmp)
     tmp(:,j,:) = tmp(:,j,:) - dummy
  ENDDO

! -------------------------------------------------------------------
! Crop to Shear Layer
! -------------------------------------------------------------------
  IF      ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
     DO j = 1, jmax
        ycenter = g(2)%nodes(j) - g(2)%scale *ycoor_ini(is) - g(2)%nodes(1)
        IF ( thick_ini(is) .eq. C_0_R ) THEN; amplify = C_1_R
        ELSE;                                 amplify = EXP(-(C_05_R*ycenter/thick_ini(is))**2); ENDIF

        tmp(:,j,:) = tmp(:,j,:)*amplify
     ENDDO

! -------------------------------------------------------------------
! Crop to Jet
! -------------------------------------------------------------------
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET   ) THEN
     DO j = 1, jmax

        ycenter =   g(2)%nodes(j) - g(2)%scale *ycoor_ini(is) - sbg(is)%diam *C_05_R - g(2)%nodes(1)
        IF ( thick_ini(is) .eq. C_0_R ) THEN; amplify = C_1_R
        ELSE;                                 amplify = EXP(-(C_05_R*ycenter/thick_ini(is))**2); ENDIF

        ycenter =-( g(2)%nodes(j) - g(2)%scale *ycoor_ini(is) + sbg(is)%diam *C_05_R - g(2)%nodes(1) )
        IF ( thick_ini(is) .eq. C_0_R ) THEN; amplify = C_1_R
        ELSE;                                 amplify = amplify + EXP(-(C_05_R*ycenter/thick_ini(is))**2); ENDIF

        tmp(:,j,:) = tmp(:,j,:)*amplify
     ENDDO

  ENDIF

! -------------------------------------------------------------------
! Scale
! -------------------------------------------------------------------
  amplify = C_0_R
  DO j = 1,jmax
     dummy = AVG1V2D(imax,jmax,kmax, j, i2, tmp)
     amplify = MAX(dummy,amplify)
  ENDDO
  amplify = norm_ini_s(is) /SQRT(amplify)

  s = s + tmp*amplify

  RETURN
END SUBROUTINE SCAL_VOLUME_BROADBAND
