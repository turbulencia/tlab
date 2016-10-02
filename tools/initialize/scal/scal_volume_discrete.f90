#include "types.h"
#include "dns_const.h"

SUBROUTINE SCAL_VOLUME_DISCRETE(is, s)

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, imode_flow
  USE SCAL_LOCAL

  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT) :: s

! -------------------------------------------------------------------
  TINTEGER i, j, k
  TREAL xi, fy, wx, wz, wxloc, wzloc
  TREAL ycenter, zprime
  TINTEGER inx2d, inx3d, inz3d

! ###################################################################
  wx = C_2_R * C_PI_R /g(1)%scale
  wz = C_2_R * C_PI_R /g(3)%scale

  ycenter = g(2)%nodes(1) +g(2)%scale *ycoor_ini(is)

! ###################################################################
! Shear
! ###################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN

     DO j = 1, jmax
        xi = g(2)%nodes(j) - ycenter
        fy = EXP(-(xi/(C_2_R*thick_ini(is)))**2)

! 2D perturbation
        DO inx2d = 1,nx2d
           wxloc = M_REAL(inx2d) *wx
           DO i = 1,imax
              zprime = A2d(inx2d) *fy &
                      *COS( wxloc*g(1)%nodes(i) +Phix2d(inx2d) )
              s(i,j,:) = s(i,j,:) + zprime
           ENDDO
        ENDDO

! 3D perturbation
        IF (kmax .GT. 1) THEN
           DO inx3d = 1,nx3d
              DO inz3d = 1,nz3d
                 wxloc = M_REAL(inx3d) *wx
                 wzloc = M_REAL(inz3d) *wz
                 DO k = 1,kmax
                    DO i = 1,imax
                       zprime = A3d(inx3d) *fy &
                               *COS( wxloc*g(1)%nodes(i) +Phix3d(inx3d) ) &
                               *COS( wzloc*g(3)%nodes(k) +Phiz3d(inz3d) )
                       s(i,j,k) = s(i,j,k) + zprime
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

     ENDDO

! ###################################################################
! Jet
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN

  ENDIF

  RETURN
END SUBROUTINE SCAL_VOLUME_DISCRETE
