#include "types.h"
#include "dns_const.h"

SUBROUTINE VELOCITY_DISCRETE(u,v,w, wrk1d, wrk2d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : g, qbg
  USE FLOW_LOCAL

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: u,v,w
  TREAL, DIMENSION(jmax,2),         INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(imax,kmax,3),    INTENT(INOUT) :: wrk2d

  ! -------------------------------------------------------------------
  TINTEGER j, k, im, idsp, kdsp
  TREAL wx, wz, wx_1, wz_1, factorx, factorz, dummy

  TREAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL, ycenter, yr
  EXTERNAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL

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
  ! Shape function
  ! ###################################################################
  SELECT CASE ( Kini%type )
  CASE ( PROFILE_PARABOLIC_SURFACE )
    ycenter = y(1) +g(2)%scale *Kini%ymean
    DO j = 1,jmax
      yr = y(j)-ycenter
      wrk1d(j,1) = FLOW_SHEAR_TEMPORAL( PROFILE_PARABOLIC, Kini%thick, C_1_R, C_0_R, ycenter, C_0_R, y(j) )
      wrk1d(j,2) = yr /( Kini%thick **2 ) *wrk1d(j,1) ! Derivative of f
      wrk1d(j,1) = wrk1d(j,1) **C_2_R
    ENDDO

  CASE ( PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_SURFACE )
    ycenter = y(1) +g(2)%scale *Kini%ymean
    DO j = 1,jmax
        yr = y(j)-ycenter
        wrk1d(j,1) = FLOW_SHEAR_TEMPORAL( PROFILE_GAUSSIAN, Kini%thick, C_1_R, C_0_R, ycenter, C_0_R, y(j) )
        wrk1d(j,2) = yr /( Kini%thick **2 ) *wrk1d(j,1) ! Derivative of f
        IF ( Kini%type .EQ. PROFILE_GAUSSIAN_SURFACE ) THEN
          wrk1d(j,1) = wrk1d(j,1) *  (yr/Kini%thick)**2
          wrk1d(j,2) = wrk1d(j,2) *( (yr/Kini%thick)**2 - C_2_R )
        ENDIF
      ENDDO

  CASE (PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_ANTISYM)
    ycenter = y(1) +g(2)%scale *Kini%ymean -C_05_R *qbg(1)%diam
    DO j = 1,jmax
      yr = y(j) - ycenter
      wrk1d(j,1) = FLOW_SHEAR_TEMPORAL( PROFILE_GAUSSIAN, Kini%thick, C_1_R, C_0_R, ycenter, C_0_R, y(j) )
      wrk1d(j,2) =-yr /( Kini%thick **2 ) *wrk1d(j,1)
    ENDDO

    ycenter = y(1) +g(2)%scale *Kini%ymean +C_05_R *qbg(1)%diam
    IF     ( fp%type .EQ. PROFILE_GAUSSIAN_ANTISYM ) THEN; factorx =-C_1_R ! varicose
    ELSEIF ( fp%type .EQ. PROFILE_GAUSSIAN_SYM     ) THEN; factorx = C_1_R ! Sinuous
    ENDIF
    DO j = 1,jmax
      yr = y(j) - ycenter
      dummy = factorx *FLOW_SHEAR_TEMPORAL( PROFILE_GAUSSIAN, Kini%thick, C_1_R, C_0_R, ycenter, C_0_R, y(j) )
      wrk1d(j,1) = wrk1d(j,1) +dummy
      wrk1d(j,2) = wrk1d(j,2) +yr /( Kini%thick **2 ) *dummy
    ENDDO

  END SELECT

  ! ###################################################################
  ! Fourier series
  ! ###################################################################
  wx_1 = C_2_R * C_PI_R / g(1)%scale ! Fundamental wavelengths
  wz_1 = C_2_R * C_PI_R / g(3)%scale

  wrk2d = C_0_R
  DO im = 1,fp%size
    wx = M_REAL( fp%modex(im) ) *wx_1
    wz = M_REAL( fp%modez(im) ) *wz_1

    ! Factor to impose solenoidal constraint
    IF     ( fp%modex(im) .EQ. 0 .AND. fp%modez(im) .EQ. 0 ) THEN; EXIT
    ELSEIF (                           fp%modez(im) .EQ. 0 ) THEN; factorx= C_1_R /wx; factorz= C_0_R
    ELSEIF ( fp%modex(im) .EQ. 0                           ) THEN; factorx= C_0_R;     factorz= C_1_R /wz
    ELSE;                                                          factorx= C_05_R/wx; factorz= C_05_R/wz
    ENDIF

    DO k = 1,kmax
      wrk2d(:,k,2) = wrk2d(:,k,2) + fp%amplitude(im) *COS( wx *x(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *z(kdsp+k) +fp%phasez(im) )
      wrk2d(:,k,1) = wrk2d(:,k,1) + fp%amplitude(im) *SIN( wx *x(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *z(kdsp+k) +fp%phasez(im) ) *factorx
      wrk2d(:,k,3) = wrk2d(:,k,3) + fp%amplitude(im) *COS( wx *x(idsp+1:idsp+imax) +fp%phasex(im) ) *SIN( wz *z(kdsp+k) +fp%phasez(im) ) *factorz
    ENDDO

  ENDDO

  ! ###################################################################
  ! Forcing
  ! ###################################################################
  DO k = 1,kmax
    DO j = 1,jmax
      u(:,j,k) = u(:,j,k) + wrk2d(:,k,1) *wrk1d(j,2)
      v(:,j,k) = v(:,j,k) + wrk2d(:,k,2) *wrk1d(j,1)
      w(:,j,k) = w(:,j,k) + wrk2d(:,k,3) *wrk1d(j,2)
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE VELOCITY_DISCRETE
