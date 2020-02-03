#include "types.h"
#include "dns_const.h"

SUBROUTINE VELOCITY_MEAN(rho, u,v,w, wrk1d,wrk3d)

  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imode_flow, imode_sim, imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : qbg
  USE DNS_GLOBAL, ONLY : coriolis

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax), INTENT(IN)    :: rho
  TREAL, DIMENSION(imax,jmax,kmax), INTENT(OUT)   :: u, v, w
  TREAL, DIMENSION(jmax,*),         INTENT(INOUT) :: wrk1d, wrk3d

! -------------------------------------------------------------------
  TINTEGER iq, j, k
  TREAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL, ycenter, calpha, salpha
  EXTERNAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL

! ###################################################################
! Isotropic case
! ###################################################################
  IF      ( imode_flow .EQ. DNS_FLOW_ISOTROPIC   ) THEN
     u = u + qbg(1)%mean
     v = v + qbg(2)%mean
     IF ( g(3)%size .GT. 1 ) THEN; w = w + qbg(3)%mean
     ELSE;                         w = C_0_R     ; ENDIF

! ###################################################################
! Shear layer case
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR       ) THEN
! -------------------------------------------------------------------
! Temporal
! -------------------------------------------------------------------
     IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN

! Construct reference profiles into array wrk1d
        DO iq = 1,3
           ycenter = g(2)%nodes(1) + g(2)%scale *qbg(iq)%ymean
           DO j = 1,jmax
              wrk1d(j,iq) = FLOW_SHEAR_TEMPORAL( &
                            qbg(iq)%type, qbg(iq)%thick, qbg(iq)%delta, qbg(iq)%mean, ycenter, qbg(iq)%parameters, g(2)%nodes(j) &
                            )
           ENDDO
        ENDDO

! Construct velocity field
        IF ( coriolis%type .NE. EQNS_NONE ) THEN
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

! Consistency check
        IF ( g(3)%size .EQ. 1 ) w = C_0_R

! -------------------------------------------------------------------
! Spatial
! -------------------------------------------------------------------
     ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     ENDIF

! ###################################################################
! Jet case
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
! -------------------------------------------------------------------
! Temporal
! -------------------------------------------------------------------
     IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN

! Construct reference profile into array wrk1d
! pilot to be added: ijet_pilot, rjet_pilot_thickness, rjet_pilot_velocity
        ycenter = g(2)%nodes(1) + g(2)%scale *qbg(1)%ymean
        DO j = 1,jmax
           wrk1d(j,1) =  FLOW_JET_TEMPORAL&
                (qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, qbg(1)%diam, ycenter, qbg(1)%parameters, g(2)%nodes(j))
        ENDDO

! Construct velocity field
        DO j = 1,jmax
           u(:,j,:) = u(:,j,:) + wrk1d(j,1)
        ENDDO

        v = v + qbg(2)%mean

        IF ( g(3)%size .GT. 1 ) THEN; w = w + qbg(3)%mean
        ELSE;                         w = C_0_R     ; ENDIF

! -------------------------------------------------------------------
! Spatial
! -------------------------------------------------------------------
     ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
#define rho_vi(j) wrk1d(j,1)
#define u_vi(j)   wrk1d(j,2)
#define aux(j)    wrk1d(j,3)
        DO j = 1,jmax
           rho_vi(j) = rho(1,j,1)
           u_vi(j)   = u(1,j,1)
        ENDDO
        ycenter = g(2)%nodes(1) + g(2)%scale *qbg(1)%ymean
        CALL FLOW_JET_SPATIAL_VELOCITY(imax,jmax, &
             qbg(1)%type, qbg(1)%thick, qbg(1)%delta, qbg(1)%mean, qbg(1)%diam, ycenter, &
             qbg(1)%parameters(1), qbg(1)%parameters(2), qbg(1)%parameters(3), &
             g(1)%nodes, g(2)%nodes, rho_vi(1), u_vi(1), rho, u, v, aux(1), wrk3d)
        IF ( kmax .GT. 1 ) THEN
           DO k = 2,kmax
              u(:,:,k) = u(:,:,1)
              v(:,:,k) = v(:,:,1)
           ENDDO
        ENDIF
#undef rho_vi
#undef u_vi
#undef aux

        IF ( g(3)%size .GT. 1 ) THEN; w = w + qbg(3)%mean
        ELSE;                         w = C_0_R     ; ENDIF

     ENDIF

  ENDIF

  RETURN
END SUBROUTINE VELOCITY_MEAN
