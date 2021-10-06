#include "types.h"
#include "dns_const.h"

MODULE FLOW_LOCAL

  USE TLAB_TYPES,  ONLY : background_dt, discrete_dt
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field
  USE TLAB_VARS, ONLY : g, qbg

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i,ims_offset_k
#endif

  IMPLICIT NONE
  SAVE
  ! -------------------------------------------------------------------
  TINTEGER :: flag_u, flag_t, flag_dilatation, flag_mixture

  TYPE(background_dt) :: Kini                   ! Geometry of perturbation of initial boundary condition
  TREAL               :: norm_ini_u, norm_ini_p ! Scaling of perturbation
  TYPE(discrete_dt)   :: fp                     ! Discrete perturbation

  TINTEGER :: flag_wall ! Boundary conditions: 0  Free-Slip/Free-Slip
  ! 1  No-Slip/Free-Slip
  ! 2  Free-Slip/No-Slip
  ! 3  No-Slip/No-Slip
  ! -------------------------------------------------------------------
  TINTEGER i, j, k

  TINTEGER im, idsp, kdsp
  TREAL wx, wz, wx_1, wz_1

  TREAL, DIMENSION(:), POINTER :: xn,zn

#include "integers.h"

CONTAINS

  ! ###################################################################
  SUBROUTINE FLOW_SHAPE( wrk1d )
    IMPLICIT NONE

    TREAL, DIMENSION(jmax,5), INTENT(INOUT) :: wrk1d

    ! -------------------------------------------------------------------
    TINTEGER bcs(2,2)
    TREAL PROFILES, ycenter, yr
    EXTERNAL PROFILES

    TREAL, DIMENSION(:), POINTER :: yn

    ! ###################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    yn => g(2)%nodes

    ycenter = yn(1) +g(2)%scale *Kini%ymean
    DO j = 1,jmax                               ! Wall-normal velocity
       wrk1d(j,1) = PROFILES( Kini%TYPE, Kini%thick, C_1_R, C_0_R, ycenter, Kini%parameters, yn(j) )
    ENDDO
    CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), wrk1d(1,1), wrk1d(1,2), wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
    wrk1d(:,2) =-wrk1d(:,2)                     ! Negative of the derivative of f, wall-parallel velocity

    SELECT CASE ( Kini%TYPE )
    CASE ( PROFILE_PARABOLIC_SURFACE )
       ! Zero wall-parallel velocity for no-slip condition, multiply by parabolic again, f=f*f
       wrk1d(:,2) = C_2_R* wrk1d(:,2) *wrk1d(:,1)          ! Wall-parallel velocity
       wrk1d(:,1) = wrk1d(:,1) **C_2_R                     ! Wall-normal velocity

    CASE ( PROFILE_GAUSSIAN_SURFACE )
       ! Zero wall-normal derivative of wall-parallel velocity for free-slip and potentialvelocity mode, f=f*tanh
       IF ( flag_wall .EQ. 1 .OR. flag_wall .EQ. 3 ) THEN  ! jmin
          DO j = 1,jmax
             yr = C_05_R *( yn(j)-yn(1)    )/ Kini%thick
             wrk1d(j,2) = wrk1d(j,2) *TANH(yr) **3 - &       ! Wall-parallel velocity
                  wrk1d(j,1) *TANH(yr) **2 /COSH(yr) **2 *C_1_5_R /Kini%thick
             wrk1d(j,1) = wrk1d(j,1) *TANH(yr) **3           ! Wall-normal velocity
          ENDDO
       ENDIF

       IF ( flag_wall .EQ. 2 .OR. flag_wall .EQ. 3 ) THEN  ! jmax
          DO j = 1,jmax
             yr = C_05_R *( yn(jmax)-yn(j) )/ Kini%thick
             wrk1d(j,2) = wrk1d(j,2) *TANH(yr) **3 + &       ! Wall-parallel velocity
                  wrk1d(j,1) *TANH(yr) **2 /COSH(yr) **2 *C_1_5_R /Kini%thick
             wrk1d(j,1) = wrk1d(j,1) *TANH(yr) **3           ! Wall-normal velocity
          ENDDO
       ENDIF

    END SELECT

    RETURN
  END SUBROUTINE FLOW_SHAPE

  ! ###################################################################
  SUBROUTINE VELOCITY_DISCRETE(u,v,w, wrk1d, wrk2d)
    IMPLICIT NONE

    TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: u,v,w
    TREAL, DIMENSION(jmax,2),         INTENT(INOUT) :: wrk1d
    TREAL, DIMENSION(imax,kmax,3),    INTENT(INOUT) :: wrk2d

    ! -------------------------------------------------------------------
    TREAL factorx, factorz

    ! ###################################################################
#ifdef USE_MPI
    idsp = ims_offset_i; kdsp = ims_offset_k
#else
    idsp = 0; kdsp = 0
#endif

    xn => g(1)%nodes
    zn => g(3)%nodes

    CALL FLOW_SHAPE( wrk1d )

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
          wrk2d(:,k,2) = wrk2d(:,k,2) + fp%amplitude(im) *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) )
          wrk2d(:,k,1) = wrk2d(:,k,1) + fp%amplitude(im) *SIN( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) ) *factorx
          wrk2d(:,k,3) = wrk2d(:,k,3) + fp%amplitude(im) *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *SIN( wz *zn(kdsp+k) +fp%phasez(im) ) *factorz
       ENDDO

    ENDDO

    DO k = 1,kmax
       DO j = 1,jmax
          u(:,j,k) = wrk2d(:,k,1) *wrk1d(j,2)
          v(:,j,k) = wrk2d(:,k,2) *wrk1d(j,1)
          w(:,j,k) = wrk2d(:,k,3) *wrk1d(j,2)
       ENDDO
    ENDDO

    IF ( norm_ini_u .GE. C_0_R ) CALL FLOW_NORMALIZE(u,v,w)

    RETURN
  END SUBROUTINE VELOCITY_DISCRETE

  ! ###################################################################
  SUBROUTINE VELOCITY_BROADBAND(u,v,w, ax,ay,az,tmp4,tmp5, wrk1d,wrk2d,wrk3d)
    USE TLAB_VARS, ONLY : visc
    IMPLICIT NONE

    TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: u,v,w
    TREAL, DIMENSION(imax,jmax,kmax), INTENT(INOUT) :: ax,ay,az,tmp4,tmp5, wrk3d
    TREAL, DIMENSION(imax,kmax,*),    INTENT(INOUT) :: wrk2d
    TREAL, DIMENSION(jmax,*),         INTENT(INOUT) :: wrk1d

    ! -------------------------------------------------------------------
    TINTEGER ibc, bcs(2,2), bcs2(2,2)
    TREAL AVG1V2D, dummy
    EXTERNAL AVG1V2D

    ! ###################################################################
    bcs = 0

    CALL FLOW_SHAPE( wrk1d )

    dummy = visc
    CALL DNS_READ_FIELDS('flow.rand', i2, imax,jmax,kmax, i3,i1, isize_field, u, wrk3d)
    CALL DNS_READ_FIELDS('flow.rand', i2, imax,jmax,kmax, i3,i2, isize_field, v, wrk3d)
    CALL DNS_READ_FIELDS('flow.rand', i2, imax,jmax,kmax, i3,i3, isize_field, w, wrk3d)
    visc = dummy

    DO j = 1,jmax   ! Remove mean
       dummy = AVG1V2D(imax,jmax,kmax, j, i1, u)
       u(:,j,:) = u(:,j,:) -dummy
       dummy = AVG1V2D(imax,jmax,kmax, j, i1, v)
       v(:,j,:) = v(:,j,:) -dummy
       dummy = AVG1V2D(imax,jmax,kmax, j, i1, w)
       w(:,j,:) = w(:,j,:) -dummy
    ENDDO

    ! ###################################################################
    SELECT CASE( flag_u )
    CASE( 2 ) ! Velocity given
       DO j = 1,jmax
          u(:,j,:) = u(:,j,:) *wrk1d(j,2)
          v(:,j,:) = v(:,j,:) *wrk1d(j,1)
          w(:,j,:) = w(:,j,:) *wrk1d(j,2)
       ENDDO

    CASE( 3 )  ! Vorticity given, solve lap(u) = - rot(vort), vort = rot(u)
       CALL FI_CURL(imax,jmax,kmax, u,v,w, ax,ay,az, tmp4, wrk2d,wrk3d)
       DO j = 1,jmax
          ax(:,j,:) = -ax(:,j,:) *wrk1d(j,2)
          ay(:,j,:) = -ay(:,j,:) *wrk1d(j,1)
          az(:,j,:) = -az(:,j,:) *wrk1d(j,2)
       ENDDO
       CALL FI_CURL(imax,jmax,kmax, ax,ay,az, u,v,w, tmp4, wrk2d,wrk3d)

       ! Solve lap(u) = - (rot(vort))_x
       IF ( flag_wall .EQ. 0 ) THEN; ibc = 3         ! FreeSlip
       ELSE;                         ibc = 0
       ENDIF  ! NoSlip
       IF ( g(1)%periodic .AND. g(3)%periodic ) THEN
          wrk2d(:,:,1:2) = C_0_R                      ! bcs
          CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, ibc,&
               u,wrk3d, tmp4,tmp5, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)
       ELSE                                          ! General treatment, need global variable ipos,jpos,kpos,ci,cj,ck
#ifdef USE_CGLOC
          CALL CGPOISSON(i1, imax,jmax,kmax,g(3)%size, u, ax,ay,az, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
       ENDIF

       ! Solve lap(v) = - (rot(vort))_y
       ibc = 0                                       ! No penetration
       IF ( g(1)%periodic .AND. g(3)%periodic ) THEN
          wrk2d(:,:,1:2) = C_0_R                      ! bcs
          CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, ibc,&
               v,wrk3d, tmp4,tmp5, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)
       ELSE                                          ! General treatment
#ifdef USE_CGLOC
          CALL CGPOISSON(i1, imax,jmax,kmax,g(3)%size, v, ax,ay,az, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
       ENDIF

       ! Solve lap(w) = - (rot(vort))_z
       IF ( g(3)%size .GT. 1 ) THEN
          IF ( flag_wall .EQ. 0 ) THEN; ibc = 3         ! FreeSlip
          ELSE;                         ibc = 0
          ENDIF  ! NoSlip
          IF ( g(1)%periodic .AND. g(3)%periodic ) THEN
             wrk2d(:,:,1:2) = C_0_R                      ! bcs
             CALL OPR_POISSON_FXZ(.FALSE., imax,jmax,kmax, g, ibc,&
                  w,wrk3d, tmp4,tmp5, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)
          ELSE                                          ! General treatment
#ifdef USE_CGLOC
             CALL CGPOISSON(i1, imax,jmax,kmax,g(3)%size, w, ax,ay,az, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
          ENDIF
       ENDIF

       ! ###################################################################
    CASE( 4 ) ! Vector potential given
       DO j = 1,jmax
          ax(:,j,:) = u(:,j,:) *wrk1d(j,1) ! Horizontal components of vector potential give vertical velocity
          ay(:,j,:) = v(:,j,:) *wrk1d(j,2)
          az(:,j,:) = w(:,j,:) *wrk1d(j,1)
       ENDDO

       bcs2 = 0
       IF ( flag_wall.EQ.1 .OR. flag_wall.EQ.3 ) bcs2(1,1) = 1 ! bcs at ymin = 1
       IF ( flag_wall.EQ.2 .OR. flag_wall.EQ.3 ) bcs2(2,1) = 1 ! bcs at ymax = 1
       ! Cannot use fi_curl. I need to impose BCs to zero to get zero velocity there
       CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs2,g(2), az,u,   wrk3d, wrk2d,wrk3d)
       CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), ay,tmp4, wrk3d, wrk2d,wrk3d)
       u = u-tmp4
       CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), ax,v,   wrk3d, wrk2d,wrk3d)
       CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), az,tmp4, wrk3d, wrk2d,wrk3d)
       v = v-tmp4
       IF ( g(3)%size .GT. 1 ) THEN
          CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), ay,w,   wrk3d, wrk2d,wrk3d)
          CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs2,g(2), ax,tmp4, wrk3d, wrk2d,wrk3d)
          w = w-tmp4
       ENDIF

    END SELECT

    ! ###################################################################
    ! Remove dilatation (vort was not really a vorticity field because it was not solenoidal)
    IF ( flag_dilatation .EQ. 1 ) THEN
       CALL FI_SOLENOIDAL(flag_wall, imax,jmax,kmax, u,v,w, ax,ay,az,tmp4,tmp5, wrk1d,wrk2d,wrk3d)
    ENDIF

    IF ( g(3)%size .EQ. 1 ) w = C_0_R       ! Impose zero spanwise velocity in 2D case

    IF ( norm_ini_u .GE. C_0_R ) CALL FLOW_NORMALIZE(u,v,w)

    RETURN
  END SUBROUTINE VELOCITY_BROADBAND

  ! ###################################################################
  SUBROUTINE FLOW_NORMALIZE(u,v,w)
    IMPLICIT NONE

    TREAL, DIMENSION(imax,jmax,kmax) :: u, v, w

    ! -------------------------------------------------------------------
    TREAL AVG1V2D, dummy, amplify
    EXTERNAL AVG1V2D

    ! ###################################################################
    amplify= C_0_R                                      ! Maximum across the layer
    DO j = 1,jmax
       dummy = AVG1V2D(imax,jmax,kmax, j, i2, u) + AVG1V2D(imax,jmax,kmax, j, i2, v) +AVG1V2D(imax,jmax,kmax, j, i2, w)
       amplify = MAX(dummy,amplify)
    ENDDO
    amplify = C_05_R *amplify

    amplify = SQRT( norm_ini_u /amplify )           ! Scaling factor to normalize to maximum TKE

    u = u *amplify
    v = v *amplify
    w = w *amplify

    RETURN
  END SUBROUTINE FLOW_NORMALIZE

END MODULE FLOW_LOCAL
