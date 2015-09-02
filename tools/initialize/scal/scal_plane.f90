!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/08/16 - J.P. Mellado
!#              Created. Extracted from old SCRTI.
!#
!########################################################################
!# DESCRIPTION
!#
!# An array disp is created and added either to the mean centerplane or to
!# the mean thickness
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

SUBROUTINE SCAL_PLANE(iflag, is, x,y,z,dx,dz, s, disp)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, jmax_total,kmax_total, isize_field, inb_scal
  USE DNS_GLOBAL, ONLY : iprof_i, thick_i, delta_i, mean_i, ycoor_i, prof_i
  USE DNS_GLOBAL, ONLY : area, scalex,scaley,scalez
  USE SCAL_LOCAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER iflag, is
  TREAL, DIMENSION(*)              :: x,y,z, dx,dz
  TREAL, DIMENSION(imax,jmax,kmax) :: s
  TREAL, DIMENSION(imax,kmax)      :: disp

! -------------------------------------------------------------------
  TINTEGER i,j,k, inx2d,inx3d,inz3d
  TINTEGER idummy, idsp, kdsp
  TREAL wx,wz, wxloc,wzloc, dummy, ycenter, thick_loc,delta_loc,mean_loc
  TREAL AVG_IK, FLOW_SHEAR_TEMPORAL
  TREAL xcenter,zcenter,rcenter, amplify

! ###################################################################
  disp = C_0_R

#ifdef USE_MPI
  idsp = ims_offset_i; kdsp = ims_offset_k 
#else
  idsp = 0; kdsp = 0
#endif

! ###################################################################
! Displacement array
! ###################################################################
! -------------------------------------------------------------------
! Broadband case
! -------------------------------------------------------------------
  IF      ( iflag .EQ. 4 .OR. iflag .EQ. 6 .OR. iflag .EQ. 8 ) THEN ! use s as aux array
     idummy = jmax_total; jmax_total = 1
     CALL DNS_READ_FIELDS('scal.rand', i1, imax,i1,kmax, inb_scal,is, isize_field, disp, s)
     jmax_total = idummy
! remove mean
     dummy = AVG_IK(imax,i1,kmax, i1, disp, dx,dz, area)
     disp = disp - dummy

! -------------------------------------------------------------------
! Discrete case
! -------------------------------------------------------------------
  ELSE IF ( iflag .EQ. 5 .OR. iflag .EQ. 7 .OR. iflag .EQ. 9 ) THEN
     IF      ( imode_discrete .EQ. 1 .OR. imode_discrete .EQ. 2 ) THEN ! sinusoidal
        wx = C_2_R*C_PI_R/scalex; wz = C_2_R*C_PI_R/scalez
     
! 1D perturbation along X
        DO inx2d = 1,nx2d
           wxloc = M_REAL(inx2d)*wx
           DO k = 1,kmax; DO i = 1,imax
              disp(i,k) = disp(i,k) + A2d(inx2d)*COS(wxloc*x(idsp+i)+Phix2d(inx2d))
           ENDDO; ENDDO
        ENDDO

! 2D perturbation along X and Z
        IF (kmax_total .GT. 1) THEN
        DO inx3d = 1,nx3d; DO inz3d = 1,nz3d
           wxloc = M_REAL(inx3d)*wx; wzloc = M_REAL(inz3d)*wz
           DO k = 1,kmax; DO i = 1,imax
              disp(i,k) = disp(i,k) + &
                   A3d(inx3d)*COS(wxloc*x(idsp+i)+Phix3d(inx3d))*COS(wzloc*z(kdsp+k)+Phiz3d(inz3d))
           ENDDO; ENDDO
        ENDDO; ENDDO
        ENDIF

     ELSE IF ( imode_discrete .EQ. 3 ) THEN ! Gaussian; completed below with modulation
        disp = A2d(1)

     ENDIF

  ENDIF

! Modulation
  IF ( delta_discrete .GT. C_0_R ) THEN
     DO k = 1,kmax; DO i = 1,imax
        xcenter   = x(i+idsp) - scalex*Phix2d(1) - x(1)
        IF ( kmax_total .GT. 1 ) THEN; zcenter = z(k+kdsp) - scalez*Phiz3d(1) - z(1)
        ELSE;                          zcenter = C_0_R; ENDIF
        rcenter   = SQRT(xcenter**2+zcenter**2)
        amplify   = EXP(-(C_05_R*rcenter/delta_discrete)**2)
        disp(i,k) = disp(i,k)*amplify
     ENDDO; ENDDO
  ENDIF

! ###################################################################
! Perturbation in the scalar field
! ###################################################################
! -------------------------------------------------------------------
! Perturbation in the centerplane
! -------------------------------------------------------------------
  IF      ( iflag .EQ. 4 .OR. iflag .EQ. 5 ) THEN
     DO k = 1,kmax; DO i = 1,imax
        ycenter = y(1) + scaley*ycoor_i(is) + disp(i,k)
        DO j = 1,jmax
           s(i,j,k) =  FLOW_SHEAR_TEMPORAL&
                (iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i(1,is),y(j))
        ENDDO
     ENDDO; ENDDO

! -------------------------------------------------------------------
! Perturbation in the thickness
! -------------------------------------------------------------------
  ELSE IF ( iflag .EQ. 6 .OR. iflag .EQ. 7 ) THEN
     DO k = 1,kmax; DO i = 1,imax
        ycenter   = y(1) + scaley*ycoor_i(is)
        thick_loc = thick_i(is) + disp(i,k)
        DO j = 1,jmax
           s(i,j,k) =  FLOW_SHEAR_TEMPORAL&
                (iprof_i(is), thick_loc, delta_i(is), mean_i(is), ycenter, prof_i(1,is),y(j))
        ENDDO
     ENDDO; ENDDO

! -------------------------------------------------------------------
! Perturbation in the magnitude (constant derivative)
! -------------------------------------------------------------------
  ELSE IF ( iflag .EQ. 8 .OR. iflag .EQ. 9 ) THEN
     DO k = 1,kmax; DO i = 1,imax
        ycenter   = y(1) + scaley*ycoor_i(is)
        delta_loc = delta_i(is) + disp(i,k)
        mean_loc  =(delta_loc)*C_05_R
        thick_loc = delta_loc/delta_i(is)*thick_i(is)
        DO j = 1,jmax
           s(i,j,k) =  FLOW_SHEAR_TEMPORAL&
                (iprof_i(is), thick_loc, delta_loc, mean_loc, ycenter, prof_i(1,is),y(j))
        ENDDO
     ENDDO; ENDDO

  ENDIF

  RETURN
END SUBROUTINE SCAL_PLANE
