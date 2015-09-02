#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2003/01/01 - J.P. Mellado
!#              Modified
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE SCAL_VOLUME_BROADBAND(is, y, dx,dz, s, tmp, wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, imode_flow
  USE DNS_GLOBAL, ONLY : diam_i
  USE DNS_GLOBAL, ONLY : area, scaley
  USE SCAL_LOCAL, ONLY : thick_ini, ycoor_ini, norm_ini_s
  IMPLICIT NONE

#include "integers.h"

  TINTEGER is
  TREAL, DIMENSION(*)              :: y, dx, dz
  TREAL, DIMENSION(imax,jmax,kmax) :: s, tmp, wrk3d

! -------------------------------------------------------------------
  TREAL AVG_IK
  TREAL dummy, savg
  TREAL ycenter, amplify
  TINTEGER j

! ###################################################################
  CALL DNS_READ_FIELDS('scal.rand', i1, imax,jmax,kmax, i1,i1, isize_field, tmp, wrk3d)

! Remove mean
  DO j = 1,jmax
     dummy = AVG_IK(imax, jmax, kmax, j, tmp, dx, dz, area)
     tmp(:,j,:) = tmp(:,j,:) - dummy
  ENDDO

! -------------------------------------------------------------------
! Crop to Shear Layer
! -------------------------------------------------------------------
  IF      ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
     DO j = 1, jmax
        ycenter = y(j) - scaley*ycoor_ini(is) - y(1)
        IF ( thick_ini(is) .eq. C_0_R ) THEN; amplify = C_1_R
        ELSE;                                 amplify = EXP(-(C_05_R*ycenter/thick_ini(is))**2); ENDIF

        tmp(:,j,:) = tmp(:,j,:)*amplify
     ENDDO

! -------------------------------------------------------------------
! Crop to Jet
! -------------------------------------------------------------------
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET   ) THEN
     DO j = 1, jmax

        ycenter =   y(j) - scaley*ycoor_ini(is) - diam_i(is)*C_05_R - y(1)
        IF ( thick_ini(is) .eq. C_0_R ) THEN; amplify = C_1_R
        ELSE;                                 amplify = EXP(-(C_05_R*ycenter/thick_ini(is))**2); ENDIF

        ycenter =-( y(j) - scaley*ycoor_ini(is) + diam_i(is)*C_05_R - y(1) )
        IF ( thick_ini(is) .eq. C_0_R ) THEN; amplify = C_1_R
        ELSE;                                 amplify = amplify + EXP(-(C_05_R*ycenter/thick_ini(is))**2); ENDIF

        tmp(:,j,:) = tmp(:,j,:)*amplify
     ENDDO

  ENDIF

! -------------------------------------------------------------------
! Scale
! -------------------------------------------------------------------
  dummy = C_0_R
  DO j = 1,jmax
     wrk3d(:,j,:) = tmp(:,j,:)**2
     savg = AVG_IK(imax,jmax,kmax, j, wrk3d, dx,dz, area)    
     dummy = MAX(savg,dummy)
  ENDDO
  dummy = norm_ini_s(is)/SQRT(dummy)

  s = s + tmp*dummy

  RETURN
END SUBROUTINE SCAL_VOLUME_BROADBAND
