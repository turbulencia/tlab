#include "types.h"
#include "dns_const.h"

SUBROUTINE VELOCITY_MEAN(u,v,w, wrk1d)

  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imode_sim, imax,jmax,kmax
  USE TLAB_VARS, ONLY : qbg
  USE TLAB_VARS, ONLY : coriolis

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: u, v, w
  TREAL, DIMENSION(jmax,*),         INTENT(INOUT) :: wrk1d

  ! -------------------------------------------------------------------
  TINTEGER iq, j
  TREAL PROFILES, ycenter, calpha, salpha
  EXTERNAL PROFILES

  !########################################################################
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN

    ! Construct reference profiles into array wrk1d
    DO iq = 1,3
      ycenter = g(2)%nodes(1) + g(2)%scale *qbg(iq)%ymean
      DO j = 1,jmax
        wrk1d(j,iq) = PROFILES( qbg(iq)%type, qbg(iq)%thick, qbg(iq)%delta, qbg(iq)%mean, ycenter, qbg(iq)%parameters, g(2)%nodes(j) )
      ENDDO
    ENDDO

    ! Construct velocity field
    IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
      calpha = COS(coriolis%parameters(1)); salpha = SIN(coriolis%parameters(1))
      wrk1d(:,3) = wrk1d(:,3) *SIGN(C_1_R,coriolis%vector(2)) ! right angular velocity vector (Garratt, 1992, p.42)

      DO j = 1,jmax
        u(:,j,:) = u(:,j,:) + wrk1d(j,1)*calpha + wrk1d(j,3)*salpha
        v(:,j,:) = v(:,j,:) + wrk1d(j,2)
        w(:,j,:) = w(:,j,:) - wrk1d(j,1)*salpha + wrk1d(j,3)*calpha
      ENDDO

    ELSE
      DO j = 1,jmax
        u(:,j,:) = u(:,j,:) + wrk1d(j,1)
        v(:,j,:) = v(:,j,:) + wrk1d(j,2)
        w(:,j,:) = w(:,j,:) + wrk1d(j,3)
      ENDDO

    ENDIF

    ! -------------------------------------------------------------------
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
! I need to pass rho; need to check
! #define rho_vi(j) wrk1d(j,1)
! #define u_vi(j)   wrk1d(j,2)
! #define aux(j)    wrk1d(j,3)
!     ycenter = g(2)%nodes(1) + g(2)%scale *qbg(iq)%ymean
!     DO j = 1,jmax
!       u_vi(j) = PROFILES( qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, ycenter, qbg(1)%parameters, g(2)%nodes(j) )
!     ENDDO
!     rho_vi(:) = rho(1,:,1)
!
!     ycenter = g(2)%nodes(1) + g(2)%scale *qbg(1)%ymean
!     CALL FLOW_SPATIAL_VELOCITY(imax,jmax, &
!     qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, qbg(1)%diam, ycenter, &
!     qbg(1)%parameters(2), qbg(1)%parameters(3), qbg(1)%parameters(4), &
!     g(1)%nodes, g(2)%nodes, rho_vi(1), u_vi(1), rho, u, v, aux(1), wrk3d)
!     IF ( g(3)%size .GT. 1 ) THEN
!       DO k = 2,kmax
!         u(:,:,k) = u(:,:,1)
!         v(:,:,k) = v(:,:,1)
!       ENDDO
!       w = w + qbg(3)%mean
!     ENDIF
! #undef rho_vi
! #undef u_vi
! #undef aux

    ENDIF

    ! -------------------------------------------------------------------
    IF ( g(3)%size .EQ. 1 ) w = C_0_R

    RETURN
  END SUBROUTINE VELOCITY_MEAN
