#include "types.h"
#include "dns_const.h"

SUBROUTINE VELOCITY_DISCRETE(iflag, u,v,w)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : imode_flow, diam_u
  USE DNS_GLOBAL, ONLY : g
  USE FLOW_LOCAL

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER iflag
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT) :: u,v,w  

! -------------------------------------------------------------------
  TINTEGER i, j, k, jsim
  TREAL ycenter, fy, fyp, wx, wz, wxloc, wzloc
  TREAL u2d, v2d, u3d, v3d, w3d
  TINTEGER inx2d, inx3d, inz3d, idsp

  TREAL, DIMENSION(:), POINTER :: x,y,z

! ###################################################################
! Define pointers
  x => g(1)%nodes
  y => g(2)%nodes
  z => g(3)%nodes

  wx = C_2_R * C_PI_R / g(1)%scale
  wz = C_2_R * C_PI_R / g(3)%scale

! ###################################################################
! Forcing for shear
! ###################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN 
     DO j = 1,jmax
        ycenter = y(j) - g(2)%scale*ycoor_ini - y(1)
        fy  = EXP(-(ycenter/(C_2_R*thick_ini))**2)
        fyp =-ycenter*fy/(C_2_R * thick_ini**2)
        IF ( iflag .EQ. 2 ) THEN ! no-slip condition at wall
           fy = fy * (ycenter/(C_2_R*thick_ini))**2
           fyp= fyp*((ycenter/(C_2_R*thick_ini))**2 - C_1_R) 
        ENDIF

! 2D perturbation
        DO inx2d = 1,nx2d
           wxloc = M_REAL(inx2d)*wx

           DO k = 1,kmax
              DO i = 1,imax
                 v2d = A2d(inx2d) * SIN(wxloc*x(i)+Phix2d(inx2d)) * fy
                 u2d = A2d(inx2d) * COS(wxloc*x(i)+Phix2d(inx2d)) *(fyp/wxloc)
                 v(i,j,k) = v(i,j,k) + v2d
                 u(i,j,k) = u(i,j,k) + u2d
              ENDDO
           ENDDO

        ENDDO

! 3D perturbation
        IF ( g(3)%size .GT. 1 ) THEN
#ifdef USE_MPI
           idsp = ims_offset_k 
#else
           idsp = 0
#endif
           DO inx3d = 1,nx3d
              DO inz3d = 1,nz3d

                 wxloc = M_REAL(inx3d)*wx
                 wzloc = M_REAL(inz3d)*wz

                 DO k = 1,kmax
                    DO i = 1,imax
                       v3d = A3d(inx3d) * SIN(wxloc*x(i)+Phix3d(inx3d)) * &
                            SIN(wzloc*z(idsp+k)+Phiz3d(inz3d)) *  fy 
                       u3d = A3d(inx3d) * COS(wxloc*x(i)+Phix3d(inx3d)) * &
                            SIN(wzloc*z(idsp+k)+Phiz3d(inz3d)) * (fyp*C_05_R/wxloc)
                       w3d = A3d(inx3d) * SIN(wxloc*x(i)+Phix3d(inx3d)) * &
                            COS(wzloc*z(idsp+k)+Phiz3d(inz3d)) * (fyp*C_05_R/wzloc)
                       u(i,j,k) = u(i,j,k) + u3d
                       v(i,j,k) = v(i,j,k) + v3d
                       w(i,j,k) = w(i,j,k) + w3d
                    ENDDO
                 ENDDO

              ENDDO
           ENDDO

        ENDIF

     ENDDO

! ###################################################################
! Forcing for jet
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
     DO j = 1,jmax/2

        jsim = jmax - j + 1

        ycenter = y(j) - g(2)%scale*ycoor_ini + diam_u/C_2_R - y(1)
        fy  = EXP(-(ycenter/(C_2_R*thick_ini))**2)*thick_ini
        fyp =-ycenter*fy/(C_2_R*thick_ini**2)

! 2D perturbation
        DO inx2d = 1,nx2d
           wxloc = M_REAL(inx2d)*wx

           DO k = 1,kmax
              DO i = 1,imax
                 u2d = A2d(inx2d) *         SIN(wxloc*x(i)+Phix2d(inx2d)) * fyp 
                 v2d =-A2d(inx2d) * wxloc * COS(wxloc*x(i)+Phix2d(inx2d)) * fy
                 u(i,j,k) = u(i,j,k) + u2d
                 v(i,j,k) = v(i,j,k) + v2d
! varicose
                 IF (ifrcdsc_mode .EQ. 1) THEN
                    u(i,jsim,k) = u(i,jsim,k) + u2d
                    v(i,jsim,k) = v(i,jsim,k) - v2d
! sinuous
                 ELSE
                    u(i,jsim,k) = u(i,jsim,k) - u2d
                    v(i,jsim,k) = v(i,jsim,k) + v2d
                 ENDIF

              ENDDO
           ENDDO

        ENDDO

! 3D perturbation
        IF (kmax .GT. 1) THEN
           DO inx3d = 1, nx3d
              DO inz3d = 1, nz3d

                 wxloc = M_REAL(inx3d)*wx
                 wzloc = M_REAL(inz3d)*wz

                 DO k=1, kmax
                    DO i = 1,imax
                       u3d =-A3d(inx3d) * COS(wxloc*x(i)+Phix3d(inx3d)) * &
                            SIN(wzloc*z(k)+Phiz3d(inz3d)) * fyp
                       v3d = A3d(inx3d) * SIN(wxloc*x(i)+Phix3d(inx3d)) * &
                            SIN(wzloc*z(k)+Phiz3d(inz3d)) * fy * &
                            (wxloc+wzloc)
                       w3d =-A3d(inx3d) * SIN(wxloc*x(i)+Phix3d(inx3d)) * &
                            COS(wzloc*z(k)+Phiz3d(inz3d)) * fyp
                       u(i,j,k) = u(i,j,k) + u3d
                       v(i,j,k) = v(i,j,k) + v3d
                       w(i,j,k) = w(i,j,k) + w3d
! varicose
                       IF (ifrcdsc_mode .EQ. 1) THEN
                          u(i,jsim,k) = u(i,jsim,k) + u3d
                          v(i,jsim,k) = v(i,jsim,k) - v3d
                          w(i,jsim,k) = w(i,jsim,k) + w3d
! sinuous
                       ELSE
                          u(i,jsim,k) = u(i,jsim,k) - u3d
                          v(i,jsim,k) = v(i,jsim,k) + v3d
                          w(i,jsim,k) = w(i,jsim,k) - w3d
                       ENDIF

                    ENDDO
                 ENDDO

              ENDDO
           ENDDO

        ENDIF

     ENDDO

! ###################################################################
! Zero forcing
! ###################################################################
  ELSE

  ENDIF

  RETURN
END SUBROUTINE VELOCITY_DISCRETE
