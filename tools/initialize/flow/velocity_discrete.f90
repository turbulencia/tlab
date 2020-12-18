#include "types.h"
#include "dns_const.h"

SUBROUTINE VELOCITY_DISCRETE(u,v,w, wrk1d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : imode_flow
  USE DNS_GLOBAL, ONLY : g, qbg
  USE FLOW_LOCAL

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_k
#endif

  IMPLICIT NONE

  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: u,v,w
  TREAL, DIMENSION(jmax,2),         INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------
  TINTEGER j, k, jsim
  TREAL wx, wz, wxloc, wzloc
  TINTEGER inx2d, inx3d, inz3d, idsp

  TREAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL, ycenter, yr
  EXTERNAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL

  TREAL, DIMENSION(:), POINTER :: x,y,z

! ###################################################################
! Define pointers
  x => g(1)%nodes
  y => g(2)%nodes
  z => g(3)%nodes

#ifdef USE_MPI
  idsp = ims_offset_k
#else
  idsp = 0
#endif

! ###################################################################
! Shape function
! ###################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
    ycenter = y(1) +g(2)%scale *ycoor_ini
    IF ( flag_wall .EQ. 3 ) THEN
      DO j = 1,jmax
        yr = y(j)-ycenter
        wrk1d(j,1) = FLOW_SHEAR_TEMPORAL( PROFILE_PARABOLIC, thick_ini, C_1_R, C_0_R, ycenter, C_0_R, y(j) )
        wrk1d(j,2) = yr /( thick_ini **2 ) *wrk1d(j,1) ! Derivative of f
        wrk1d(j,1) = wrk1d(j,1) **C_2_R
      ENDDO

    ELSE
      DO j = 1,jmax
        yr = y(j)-ycenter
        wrk1d(j,1) = FLOW_SHEAR_TEMPORAL( PROFILE_GAUSSIAN, thick_ini, C_1_R, C_0_R, ycenter, C_0_R, y(j) )
        wrk1d(j,2) = yr /( thick_ini **2 ) *wrk1d(j,1) ! Derivative of f
        IF ( flag_wall .EQ. 1 ) THEN
           wrk1d(j,1) = wrk1d(j,1) *  (yr/thick_ini)**2
           wrk1d(j,2) = wrk1d(j,2) *( (yr/thick_ini)**2 - C_2_R )
        ENDIF
      ENDDO

    ENDIF

! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
    ycenter = y(1) +g(2)%scale *ycoor_ini -C_05_R *qbg(1)%diam
    DO j = 1,jmax/2
       yr = y(j) - ycenter
       wrk1d(j,1) = FLOW_SHEAR_TEMPORAL( PROFILE_GAUSSIAN, thick_ini, C_1_R, C_0_R, ycenter, C_0_R, y(j) )
       wrk1d(j,2) =-yr /( thick_ini **2 ) *wrk1d(j,1)

       jsim = jmax - j + 1
       IF (imode_discrete .EQ. 1) THEN ! varicose
         wrk1d(jsim,1) =-wrk1d(j,1)
         wrk1d(jsim,2) = wrk1d(j,2)
       ELSE                            ! sinuous
         wrk1d(jsim,1) = wrk1d(j,1)
         wrk1d(jsim,2) =-wrk1d(j,2)
       ENDIF

     ENDDO
  ENDIF

! ###################################################################
! Forcing
! ###################################################################
  wx = C_2_R * C_PI_R / g(1)%scale
  wz = C_2_R * C_PI_R / g(3)%scale

  DO j = 1,jmax

    DO inx2d = 1,nx2d       ! 2D perturbation
      wxloc = M_REAL(inx2d)*wx

      DO k = 1,kmax
        v(:,j,k) = v(:,j,k) &
                 + A2d(inx2d) *COS( wxloc*x(:) +Phix2d(inx2d) ) *wrk1d(j,1)
        u(:,j,k) = u(:,j,k) &
                 + A2d(inx2d) *SIN( wxloc*x(:) +Phix2d(inx2d) ) *wrk1d(j,2) /wxloc
      ENDDO

    ENDDO

    DO inx3d = 1,nx3d       ! 3D perturbation
      wxloc = M_REAL(inx3d)*wx
      DO inz3d = 1,nz3d
        wzloc = M_REAL(inz3d)*wz

        DO k = idsp+1,idsp+kmax
          v(:,j,k) = v(:,j,k) &
                   + A3d(inx3d) *COS( wxloc*x(:) +Phix3d(inx3d) ) *COS( wzloc*z(k) +Phiz3d(inz3d) ) *wrk1d(j,1)
          u(:,j,k) = u(:,j,k) &
                   + A3d(inx3d) *SIN( wxloc*x(:) +Phix3d(inx3d) ) *COS( wzloc*z(k) +Phiz3d(inz3d) ) *wrk1d(j,2) *C_05_R /wxloc
          w(:,j,k) = w(:,j,k) &
                   + A3d(inx3d) *COS( wxloc*x(:) +Phix3d(inx3d) ) *SIN( wzloc*z(k) +Phiz3d(inz3d) ) *wrk1d(j,2) *C_05_R /wzloc
        ENDDO

      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE VELOCITY_DISCRETE
