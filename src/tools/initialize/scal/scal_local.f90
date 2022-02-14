#include "types.h"
#include "dns_const.h"

MODULE SCAL_LOCAL

  USE TLAB_TYPES,  ONLY : background_dt, discrete_dt
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field, inb_scal, MAX_NSP
  USE TLAB_VARS, ONLY : g, sbg
  USE IO_FIELDS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i,ims_offset_k
#endif

  IMPLICIT NONE
  SAVE

  ! -------------------------------------------------------------------
  TINTEGER :: flag_s, flag_mixture

  TYPE(background_dt) :: Sini(MAX_NSP)                            ! Geometry of perturbation of initial boundary condition
  TREAL               :: norm_ini_s(MAX_NSP), norm_ini_radiation  ! Scaling of perturbation
  TYPE(discrete_dt)   :: fp                                       ! Discrete perturbation

  ! -------------------------------------------------------------------
  TINTEGER i, j, k

  TINTEGER im, idsp, kdsp
  TREAL wx, wz, wx_1, wz_1

  TREAL, DIMENSION(:), POINTER :: xn,zn

#include "integers.h"

CONTAINS

! ###################################################################
SUBROUTINE SCAL_SHAPE( is, wrk1d )
  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(jmax,1), INTENT(INOUT) :: wrk1d

  ! -------------------------------------------------------------------
  TREAL PROFILES, ycenter, yr
  EXTERNAL PROFILES

  TREAL, DIMENSION(:), POINTER :: yn

  ! ###################################################################
  yn => g(2)%nodes

  ycenter = yn(1) +g(2)%scale *Sini(is)%ymean
  DO j = 1,jmax
    wrk1d(j,1) = PROFILES( Sini(is)%type, Sini(is)%thick, C_1_R, C_0_R, ycenter, Sini(is)%parameters, yn(j) )
  ENDDO

  SELECT CASE ( Sini(is)%type )
  CASE ( PROFILE_GAUSSIAN_SURFACE ) ! set perturbation and its normal derivative to zero at the boundaries
    DO j = 1, jmax
      yr = C_05_R*(yn(j)-yn(1)   )/Sini(is)%thick
      wrk1d(j,1)  = wrk1d(j,1)  *TANH(yr) **2
      yr =-C_05_R*(yn(j)-yn(jmax))/Sini(is)%thick
      wrk1d(j,1)  = wrk1d(j,1)  *TANH(yr) **2
    ENDDO

  END SELECT

  RETURN
END SUBROUTINE SCAL_SHAPE

! ###################################################################
SUBROUTINE SCAL_FLUCTUATION_VOLUME( is, s, tmp, wrk1d, wrk2d, wrk3d )
  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: s
  TREAL, DIMENSION(jmax,1),         INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(imax,kmax,1),    INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(INOUT) :: tmp, wrk3d

  ! -------------------------------------------------------------------
  TREAL AVG1V2D, dummy, amplify
  EXTERNAL AVG1V2D

! ###################################################################
#ifdef USE_MPI
  idsp = ims_offset_i; kdsp = ims_offset_k
#else
  idsp = 0; kdsp = 0
#endif

  xn => g(1)%nodes
  zn => g(3)%nodes

  CALL SCAL_SHAPE( is, wrk1d )

  SELECT CASE( flag_s )
  CASE( 1 )   ! Broadband case
    CALL IO_READ_FIELDS('scal.rand', IO_SCAL, imax,jmax,kmax, inb_scal,is, tmp, wrk3d)

    amplify = C_0_R
    DO j = 1,jmax
      dummy = AVG1V2D(imax,jmax,kmax, j, i1, tmp)       ! Calculate mean
      wrk3d(:,j,:) = ( tmp(:,j,:) - dummy )*wrk1d(j,1)  ! Remove mean and apply shape function
    ENDDO

  CASE( 2 )   ! Discrete case
    wx_1 = C_2_R * C_PI_R / g(1)%scale ! Fundamental wavelengths
    wz_1 = C_2_R * C_PI_R / g(3)%scale

    wrk2d = C_0_R
    DO im = 1,fp%size
      wx = M_REAL( fp%modex(im) ) *wx_1
      wz = M_REAL( fp%modez(im) ) *wz_1
      DO k = 1,kmax
        wrk2d(:,k,1) = wrk2d(:,k,1) + fp%amplitude(im) *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) )
      ENDDO
    ENDDO

    DO k = 1,kmax; DO j = 1,jmax
      wrk3d(:,j,k) = wrk2d(:,k,1) *wrk1d(j,1)
    ENDDO; ENDDO

  END SELECT

  IF ( norm_ini_s(is) .GT. C_0_R ) CALL SCAL_NORMALIZE( is, wrk3d )

  s = s +wrk3d

  RETURN
END SUBROUTINE SCAL_FLUCTUATION_VOLUME

! ###################################################################
SUBROUTINE SCAL_FLUCTUATION_PLANE(is, s, disp)
  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(imax,jmax,kmax) :: s
  TREAL, DIMENSION(imax,kmax)      :: disp

  ! -------------------------------------------------------------------
  TINTEGER idummy, io_sizes(5)
  TREAL dummy, ycenter, thick_loc,delta_loc,mean_loc
  TREAL AVG1V2D, PROFILES
  TREAL xcenter,zcenter,rcenter, amplify

  CHARACTER*32 varname

  ! ###################################################################
#ifdef USE_MPI
  idsp = ims_offset_i; kdsp = ims_offset_k
#else
  idsp = 0; kdsp = 0
#endif

  xn => g(1)%nodes
  zn => g(3)%nodes

  ! ###################################################################
  disp = C_0_R
  SELECT CASE( flag_s )
  CASE( 4,6,8 )   ! Broadband case
    WRITE(varname,*) is; varname = TRIM(ADJUSTL(varname))
    idummy=imax*kmax; io_sizes = (/idummy,1,idummy,1,1/)
    CALL IO_READ_SUBARRAY8(i1, 'scal.rand', varname, disp, io_sizes, s) ! using array s as aux array
    dummy = AVG1V2D(imax,i1,kmax, i1, i1, disp)     ! remove mean
    disp = disp - dummy

  CASE( 5,7,9 )   ! Discrete case
      wx_1 = C_2_R * C_PI_R / g(1)%scale ! Fundamental wavelengths
      wz_1 = C_2_R * C_PI_R / g(3)%scale

      DO im = 1,fp%size
        wx = M_REAL( fp%modex(im) ) *wx_1
        wz = M_REAL( fp%modez(im) ) *wz_1

        IF ( fp%type .EQ. 2 ) THEN ! Smoothed step funtion Tanh(a*Cos(\xi/b))
          IF ( fp%parameters(2) .LE. C_0_R ) THEN; dummy = C_BIG_R;
          ELSE; dummy = C_05_R /( wx *fp%parameters(2) ); ENDIF
          DO k = 1,kmax
            disp(:,k) = disp(:,k) + fp%amplitude(im) *TANH( dummy *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) ) )
          ENDDO

        ELSE
          DO k = 1,kmax
            disp(:,k) = disp(:,k) + fp%amplitude(im) *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) )
          ENDDO

        ENDIF

      ENDDO

  END SELECT

  ! Modulation
  IF ( fp%type .EQ. 1 .AND. fp%parameters(1) .GT. C_0_R ) THEN
    DO k = 1,kmax; DO i = 1,imax
      xcenter   = g(1)%nodes(i+idsp) - g(1)%scale *fp%phasex(1) - g(1)%nodes(1)
      IF ( g(3)%size .GT. 1 ) THEN; zcenter = g(3)%nodes(k+kdsp) - g(3)%scale *fp%phasez(1) - g(3)%nodes(1)
      ELSE;                         zcenter = C_0_R; ENDIF
      rcenter   = SQRT(xcenter**2+zcenter**2)
      amplify   = EXP(-C_05_R*(rcenter/fp%parameters(1))**2)
      disp(i,k) = disp(i,k)*amplify
    ENDDO; ENDDO
  ENDIF

  ! ###################################################################
  SELECT CASE( flag_s )
  CASE( 4,5 )           ! Perturbation in the centerplane
    DO k = 1,kmax; DO i = 1,imax
      ycenter = g(2)%nodes(1) + g(2)%scale *sbg(is)%ymean + disp(i,k)
      DO j = 1,jmax
        s(i,j,k) = PROFILES(sbg(is)%type, sbg(is)%thick, sbg(is)%delta, sbg(is)%mean, ycenter, sbg(is)%parameters, g(2)%nodes(j))
      ENDDO
    ENDDO; ENDDO

  CASE( 6,7 )           ! Perturbation in the thickness
    DO k = 1,kmax; DO i = 1,imax
      ycenter   = g(2)%nodes(1) + g(2)%scale *sbg(is)%ymean
      thick_loc = sbg(is)%thick + disp(i,k)
      DO j = 1,jmax
        s(i,j,k) = PROFILES(sbg(is)%type, thick_loc, sbg(is)%delta, sbg(is)%mean, ycenter, sbg(is)%parameters, g(2)%nodes(j))
      ENDDO
    ENDDO; ENDDO

  CASE( 8,9 )           ! Perturbation in the magnitude (constant derivative)
    DO k = 1,kmax; DO i = 1,imax
      ycenter   = g(2)%nodes(1) + g(2)%scale *sbg(is)%ymean
      delta_loc = sbg(is)%delta + disp(i,k)
      mean_loc  =(delta_loc)*C_05_R
      IF ( sbg(is)%delta .GT. 0 ) THEN; thick_loc = delta_loc /sbg(is)%delta *sbg(is)%thick;
      ELSE;                             thick_loc = sbg(is)%thick; ENDIF
      DO j = 1,jmax
        s(i,j,k) = PROFILES(sbg(is)%type, thick_loc, delta_loc, mean_loc, ycenter, sbg(is)%parameters, g(2)%nodes(j))
      ENDDO
    ENDDO; ENDDO

  END SELECT

  RETURN
END SUBROUTINE SCAL_FLUCTUATION_PLANE

! ###################################################################
SUBROUTINE SCAL_NORMALIZE(is, s)
  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(INOUT) :: s

  ! -------------------------------------------------------------------
  TREAL AVG1V2D, dummy, amplify
  EXTERNAL AVG1V2D

  ! ###################################################################
  amplify= C_0_R                                      ! Maximum across the layer
  DO j = 1,jmax
    dummy = AVG1V2D(imax,jmax,kmax, j, i2, s)
    amplify = MAX(dummy,amplify)
  ENDDO

  amplify = norm_ini_s(is) /SQRT( amplify )           ! Scaling factor to normalize to maximum rms

  s = s *amplify

  RETURN
END SUBROUTINE SCAL_NORMALIZE

END MODULE SCAL_LOCAL
