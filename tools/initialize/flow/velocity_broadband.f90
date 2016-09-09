#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/01/01 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Given a random velocity field, impose a shape function in the 
!# corresponding vorticity or velocity field and remove solenoidal part.
!#
!########################################################################
!# ARGUMENTS 
!#
!# iflag     In   2  cropping in velocity
!#                3  cropping in vorticity
!#                4  cropping in velocity potential
!#
!########################################################################
SUBROUTINE VELOCITY_BROADBAND(iflag, u,v,w, tmp1,tmp2,tmp3,tmp4,tmp5, &
     ipos,jpos,kpos,ci,cj,ck, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imode_fdm, imax,jmax,kmax,kmax_total, i1bc,j1bc,k1bc, isize_wrk1d, isize_field
  USE DNS_GLOBAL, ONLY : imode_flow, visc, area
  USE DNS_GLOBAL, ONLY : thick_u, diam_u
  USE FLOW_LOCAL, ONLY : flag_wall, flag_dilatation, thick_ini,  ycoor_ini

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iflag

  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: u,v,w
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5, wrk3d
  TREAL, DIMENSION(*)              :: ipos,jpos,kpos, ci,cj,ck
  TREAL, DIMENSION(imax,kmax,*)    :: wrk2d
  TREAL, DIMENSION(isize_wrk1d,*)  :: wrk1d

  TARGET :: tmp1, tmp2, tmp3

! -------------------------------------------------------------------
  TINTEGER j, idsp, ibc, ibcmin, ibcmax
  TREAL ycenter, ymin, ymax, amplify, dummy

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: wx, wy, wz
  TREAL, DIMENSION(:),     POINTER :: x,y,z, dx,dy,dz

! ###################################################################
! Define pointers
  wx => tmp1
  wy => tmp2
  wz => tmp3

  x => g(1)%nodes; dx => g(1)%aux(:,1)
  y => g(2)%nodes; dy => g(2)%aux(:,1)
  z => g(3)%nodes; dz => g(3)%aux(:,1)

! ###################################################################
! Read initial random field
! ###################################################################
! read 3 first fields (u, v, w are not contiguous in memory !)
  dummy = visc
  CALL DNS_READ_FIELDS('flow.rand', i2, imax,jmax,kmax, i3,i1, isize_field, u, wrk3d)
  CALL DNS_READ_FIELDS('flow.rand', i2, imax,jmax,kmax, i3,i2, isize_field, v, wrk3d)
  CALL DNS_READ_FIELDS('flow.rand', i2, imax,jmax,kmax, i3,i3, isize_field, w, wrk3d)
  visc = dummy

! remove average
  CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, u)
  CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, v)
  CALL REYFLUCT2D(imax,jmax,kmax, dx,dz, area, w)

! change sign in the lateral velocity in the upper layer
  IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
     DO j = 1,jmax
        ycenter = y(j) - g(2)%scale*C_05_R - y(1)
        amplify = TANH(-C_05_R*ycenter/thick_u)
        v(:,j,:) = amplify*v(:,j,:)
     ENDDO
  ENDIF

! ###################################################################
! Crop field to specific geometry.
!
! If vorticity, this operation implies that div(w) is not zero and 
! that is the reason for the step at the end removing the dilatation. 
! Strictly, the field w after this step is not a vorticity because it 
! is not solenoidal.
!
! Note also that shear and jet cases are easy enough (only f(y)) and the 
! zero dilatation condition could be easily imposed
!
! ###################################################################
! -------------------------------------------------------------------
! Set field for cropping
! -------------------------------------------------------------------
  IF      ( iflag .EQ. 2 ) THEN ! velocity given
     wx = u; wy = v; wz = w
     IF ( kmax_total .EQ. 1 ) wz = 0
     
  ELSE IF ( iflag .EQ. 3 ) THEN ! vorticity given
     CALL FI_CURL(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          dx,dy,dz, u,v,w, wx,wy,wz, tmp4, wrk1d,wrk2d,wrk3d)
     IF ( kmax_total .EQ. 1 ) THEN; wx = C_0_R; wy = C_0_R; ENDIF ! exactly zero

  ELSE IF ( iflag .EQ. 4 ) THEN ! velocity potential given
     wx = u; wy = v; wz = w

  ENDIF

! -------------------------------------------------------------------
! Shear layer
! -------------------------------------------------------------------
  IF      ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
     DO j = 1,jmax
        ycenter = y(j) - g(2)%scale*ycoor_ini - y(1)
        IF ( thick_ini .eq. C_0_R ) THEN; amplify = C_0_R
        ELSE;                             amplify = EXP(-(C_05_R*ycenter/thick_ini)**2); ENDIF

        ymin = (y(j)-y(1)   )/thick_ini
        ymax = (y(j)-y(jmax))/thick_ini
        IF ( flag_wall.EQ.1 .OR. flag_wall.EQ.3 ) amplify = amplify *TANH( C_05_R*ymin) !**2
        IF ( flag_wall.EQ.2 .OR. flag_wall.EQ.3 ) amplify = amplify *TANH(-C_05_R*ymax) !**2
        v(:,j,:) = wy(:,j,:)*amplify
        IF ( flag_wall.EQ.1 .OR. flag_wall.EQ.3 ) amplify = amplify *TANH( C_05_R*ymin) !**2
        IF ( flag_wall.EQ.2 .OR. flag_wall.EQ.3 ) amplify = amplify *TANH(-C_05_R*ymax) !**2
        u(:,j,:) = wx(:,j,:)*amplify
        w(:,j,:) = wz(:,j,:)*amplify

     ENDDO

! -------------------------------------------------------------------
! Jet
! -------------------------------------------------------------------
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET   ) THEN
     DO j = 1, jmax           
        ycenter =   y(j) - g(2)%scale*ycoor_ini - diam_u*C_05_R - y(1)
        IF ( thick_ini .eq. C_0_R ) THEN; amplify = C_0_R
        ELSE;                             amplify = EXP(-(C_05_R*ycenter/thick_ini)**2); ENDIF

        ycenter =-( y(j) - g(2)%scale*ycoor_ini + diam_u*C_05_R - y(1) )
        IF ( thick_ini .eq. C_0_R ) THEN; amplify = C_0_R
        ELSE;                             amplify = amplify + EXP(-(C_05_R*ycenter/thick_ini)**2); ENDIF

        u(:,j,:) = wx(:,j,:)*amplify
        v(:,j,:) = wy(:,j,:)*amplify
        w(:,j,:) = wz(:,j,:)*amplify

     ENDDO

  ENDIF

! ###################################################################
! If vorticity case, solve lap(u) = - rot vort.
! ###################################################################
  IF ( iflag .EQ. 3 ) THEN
! -------------------------------------------------------------------
! Calculate rot vort
! -------------------------------------------------------------------
     CALL FI_CURL(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          dx,dy,dz, u,v,w, wx,wy,wz, tmp4, wrk1d,wrk2d,wrk3d)
     IF ( kmax_total .EQ. 1 ) THEN; wz = C_0_R; ENDIF ! exactly zero

! -------------------------------------------------------------------
! Solve lap(u) = - (rot vort)_x
! -------------------------------------------------------------------
     IF ( flag_wall .EQ. 0 ) THEN; ibc = 3        ! FreeSlip
     ELSE;                         ibc = 0; ENDIF ! NoSlip

! v, w are aux arrays
     IF ( i1bc .EQ. 0 .AND. k1bc .EQ. 0 ) THEN ! Doubly periodic in xOz 
        wrk2d(:,:,1:2) = C_0_R  ! bcs
        u = -wx                 ! change of forcing term sign
        CALL OPR_POISSON_FXZ(imode_fdm,i1,ibc, imax,jmax,kmax, g,&
             u,wrk3d, tmp4,tmp5, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)
        
     ELSE                                      ! General treatment
#ifdef USE_CGLOC
        CALL CGPOISSON(i1, imax,jmax,kmax,kmax_total, i1bc,j1bc,k1bc, &
             dx,dy,dz, u, wx,v,w, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
     ENDIF
        
! -------------------------------------------------------------------
! Solve lap(v) = - (rot vort)_y
! -------------------------------------------------------------------
     ibc = 0 ! No penetration

! w, wx are aux arrays
     IF ( i1bc .EQ. 0 .AND. k1bc .EQ. 0 ) THEN ! Doubly periodic in xOz
        wrk2d(:,:,1:2) = C_0_R  ! bcs
        v = -wy                 ! change sign of forcing term
        CALL OPR_POISSON_FXZ(imode_fdm,i1,ibc, imax,jmax,kmax, g,&
             v,wrk3d, tmp4,tmp5, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)
        
     ELSE                                      ! General treatment
#ifdef USE_CGLOC
        CALL CGPOISSON(i1, imax,jmax,kmax,kmax_total, i1bc,j1bc,k1bc, &
             dx,dy,dz, v, wy,w,wx, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
     ENDIF

! -------------------------------------------------------------------
! Solve lap(w) = - (rot vort)_z
! -------------------------------------------------------------------
     IF ( kmax_total .GT. 1 ) THEN
        IF ( flag_wall .EQ. 0 ) THEN; ibc = 3        ! FreeSlip
        ELSE;                         ibc = 0; ENDIF ! NoSlip

! wx, wy are aux arrays
        IF ( i1bc .EQ. 0 .AND. k1bc .EQ. 0 ) THEN ! Doubly periodic in xOz
           wrk2d(:,:,1:2) = C_0_R  ! bcs
           w = -wz                 ! change sign of forcing term
           CALL OPR_POISSON_FXZ(imode_fdm,i1,ibc, imax,jmax,kmax, g,&
                w,wrk3d, tmp4,tmp5, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),wrk3d)
        
        ELSE                                      ! General treatment
#ifdef USE_CGLOC
           CALL CGPOISSON(i1, imax,jmax,kmax,kmax_total, i1bc,j1bc,k1bc, &
                dx,dy,dz, w, wz,wx,wy, ipos,jpos,kpos,ci,cj,ck, wrk2d)
#endif
        ENDIF
        
     ELSE
        w = C_0_R
        
     ENDIF

  ENDIF

! ###################################################################
! Vector potential
! ###################################################################
  IF ( iflag .EQ. 4 ) THEN
     ibcmin = 0; ibcmax = 0
     IF ( flag_wall.EQ.1 .OR. flag_wall.EQ.3 ) ibcmin = 1
     IF ( flag_wall.EQ.2 .OR. flag_wall.EQ.3 ) ibcmax = 1
! I need to impose BCs to zero to get zero velocity there
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, w,wx,   ibcmin,ibcmax, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, v,tmp4, i0,i0,         wrk1d,wrk2d,wrk3d)
     wx = wx-tmp4
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, u,wy,   i0,i0,         wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, w,tmp4, i0,i0,         wrk1d,wrk2d,wrk3d)
     wy = wy-tmp4
     IF ( kmax_total .GT. 1 ) THEN
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v,wz,   i0,i0,         wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, u,tmp4, ibcmin,ibcmax, wrk1d,wrk2d,wrk3d)
     wz = wz-tmp4
     ELSE
     wz = C_0_R
     ENDIF

     u = wx; v = wy; w = wz

  ENDIF

! ###################################################################
! Remove Dilatational part
! 
! As explained before, due to the fact that w was really not a 
! vorticity field because it was not solenoidal.
! ###################################################################
  IF ( flag_dilatation .EQ. 1 ) THEN
     CALL SOLENOIDAL(flag_wall, u,v,w, tmp1,tmp2,tmp3,tmp4,tmp5, ipos,jpos,kpos,ci,cj,ck, wrk1d,wrk2d,wrk3d)
  ENDIF

  RETURN
END SUBROUTINE VELOCITY_BROADBAND

