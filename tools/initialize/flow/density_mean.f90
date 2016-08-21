#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Setting the density field.
!# The default case is when temperature field is set similar to a passive 
!# scalar, and the density is then derived taking pressure constant. 
!# The inputs are however density maximum/minimum:
!#
!# T = T_0 + s(T_1-T_0), s.t. 0<=s<=1, =>
!#           => rho = W / [ W_0/rho_0 + s(W_1/rho_1-W_0/rho_0) ]
!#
!# This allow to compare the evolution of the passive scalar with that 
!# of the temperature.
!#
!# If volumetric forces are present, then rho is computed from
!# dp/dy = rho g, to ensure hydrostatic equilibrium including numerical
!# errors.
!#
!# Spatial case with mulstispecies is not jet done
!#
!########################################################################
SUBROUTINE DENSITY_MEAN(rho, p,T,s, txc, wrk1d,wrk2d,wrk3d)
  
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE
  
#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(IN)    :: p,T
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(OUT)   :: rho
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(OUT)   :: s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: txc, wrk3d
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: wrk1d,wrk2d

! -------------------------------------------------------------------
  TREAL ycenter, dummy
  TINTEGER i, j, k, ij, is
  TREAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL
  EXTERNAL FLOW_SHEAR_TEMPORAL, FLOW_JET_TEMPORAL

  TREAL, DIMENSION(:), POINTER :: x,y,dy
  
! ###################################################################
! Define pointers
  x => g(1)%nodes
  y => g(2)%nodes; dy => g(2)%aux(:,1)
   
! ###################################################################
! Isotropic case
! ###################################################################
  IF      ( imode_flow .EQ. DNS_FLOW_ISOTROPIC ) THEN
     rho =  rho + mean_rho

! ###################################################################
! Shear layer case
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR     ) THEN

! -------------------------------------------------------------------
! Temporal shear layer case without volumetric force: 
! Calculate density from equation of state
! -------------------------------------------------------------------
     IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN
        IF ( ibodyforce .EQ. EQNS_NONE ) THEN

#define TEM_MEAN_LOC(i,j,k) wrk3d(i,j,k)
#define RHO_MEAN_LOC(i,j,k) txc(i,j,k)

! temperature/mixture profile are given
           IF ( iprof_rho .EQ. PROFILE_NONE ) THEN
              ycenter = y(1) + scaley*ycoor_tem
              DO j = 1,jmax
                 dummy =  FLOW_SHEAR_TEMPORAL&
                      (iprof_tem, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y(j))
                 DO k = 1,kmax
                    DO i = 1,imax
                       TEM_MEAN_LOC(i,j,k) = dummy
                    ENDDO
                 ENDDO
              ENDDO

              DO is = 1,inb_scal
                 ycenter = y(1) + scaley*ycoor_i(is)
                 DO j = 1,jmax
                    dummy =  FLOW_SHEAR_TEMPORAL&
                         (iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i,y(j))
                    DO k = 1,kmax
                       DO i = 1,imax
                          s(i,j,k,is) = dummy
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO

! define liquid content in AirWater case: (p,T) given
              IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
                 CALL THERMO_AIRWATER_PT(imax, jmax, kmax, s, p, TEM_MEAN_LOC(1,1,1))
              ENDIF

              CALL THERMO_THERMAL_DENSITY&
                   (imax, jmax, kmax, s, p, TEM_MEAN_LOC(1,1,1), RHO_MEAN_LOC(1,1,1))
              DO k = 1,kmax
                 DO ij = 1,imax*jmax
                    rho(ij,1,k) = rho(ij,1,k) + RHO_MEAN_LOC(ij,1,k)
                 ENDDO
              ENDDO

! density profile itself is given
           ELSE
              ycenter = y(1) + scaley*ycoor_rho
              DO j = 1,jmax
                 dummy =  FLOW_SHEAR_TEMPORAL&
                      (iprof_rho, thick_rho, delta_rho, mean_rho, ycenter, prof_rho, y(j))
                 DO k = 1,kmax
                    DO i = 1,imax
                       rho(i,j,k) = rho(i,j,k) + dummy
                    ENDDO
                 ENDDO
              ENDDO

           ENDIF

#undef TEM_MEAN_LOC
#undef RHO_MEAN_LOC

! -------------------------------------------------------------------
! Temporal shear layer case with volumetric force:
! Calculate density from hydrostatic equilibrium
! assuming a volumetric force along OY
! -------------------------------------------------------------------
        ELSE
! AIRWATER case. Routine PARTIAL_Y introduces small errors in equilibrium
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
              CALL THERMO_THERMAL_DENSITY(imax, jmax, kmax, s, p, T, rho)

! General case
           ELSE
              CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
                   dy, p, txc, i0, i0, wrk1d, wrk2d, wrk3d)
              DO k = 1,kmax
                 DO ij = 1,imax*jmax
                    rho(ij,1,k) = rho(ij,1,k) + txc(ij,1,k)/body_vector(2)
                 ENDDO
              ENDDO
           ENDIF

        ENDIF

! -------------------------------------------------------------------
! Spatial shear layer case
! -------------------------------------------------------------------
     ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
        CALL IO_WRITE_ASCII(efile, 'DENSITY_MEAN. Spatial shear layer undeveloped')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)

     ENDIF

! ###################################################################
! Jet case
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN

! -------------------------------------------------------------------
! Temporal jet
! -------------------------------------------------------------------

! temperature/mixture profile are given
     IF ( iprof_rho .EQ. PROFILE_NONE ) THEN
        ycenter = y(1) + scaley*ycoor_tem
        DO j = 1,jmax
           dummy = FLOW_JET_TEMPORAL&
                (iprof_tem, thick_tem, delta_tem, mean_tem, diam_tem, ycenter, prof_tem, y(j))
           DO k = 1,kmax
              DO i = 1,imax
                 wrk3d(i,j,k) = dummy
              ENDDO
           ENDDO
        ENDDO

        DO is = 1,inb_scal
           ycenter = y(1) + scaley*ycoor_i(is)
           DO j = 1,jmax
              dummy =  FLOW_JET_TEMPORAL&
                   (iprof_i(is), thick_i(is), delta_i(is), mean_i(is), diam_i(is), ycenter, prof_i, y(j))
              DO k = 1,kmax
                 DO i = 1,imax
                    s(i,j,k,is) = dummy
                 ENDDO
              ENDDO
           ENDDO
        ENDDO

        CALL THERMO_THERMAL_DENSITY(imax, jmax, kmax, s, p, wrk3d, rho)

! density profile itself is given
     ELSE
     ENDIF

! -------------------------------------------------------------------
! Spatial jet
! Only if there is a density variation. Constant density is already
! initialized in previous routine segment.
! -------------------------------------------------------------------
     IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. delta_rho .NE. C_0_R ) THEN

! temperature/mixture profile are given
        IF ( iprof_rho .EQ. PROFILE_NONE ) THEN
#define rho_vi(j) wrk1d(j,1)
#define u_vi(j)   wrk1d(j,2)
#define aux1(j)   wrk1d(j,3)
#define aux2(j)   wrk1d(j,4)
#define aux3(j)   wrk1d(j,5)
#define aux4(j)   wrk1d(j,6)
! Inflow profile of density
           DO j = 1,jmax
              rho_vi(j) = rho(1,j,1)
           ENDDO

! Inflow profile of axial velocity
           u_vi(1:jmax) = C_0_R
           ycenter = y(1) + scaley*ycoor_u
           DO j = 1,jmax
              u_vi(j) = FLOW_JET_TEMPORAL&
                   (iprof_u, thick_u, delta_u, mean_u, diam_u, ycenter, prof_u, y(j))
           ENDDO

! 2D distribution of density
           CALL FLOW_JET_SPATIAL_DENSITY(imax, jmax, iprof_tem, thick_tem, delta_tem, mean_tem, &
                ycoor_tem, diam_tem, jet_tem, iprof_u, thick_u, delta_u, mean_u, ycoor_u, diam_u,&
                jet_u, scaley, x, y, s,p,rho_vi(1),u_vi(1),aux1(1),rho,aux2(1),aux3(1),aux4(1)) 
                
           IF ( kmax .GT. 1 ) THEN
              DO k = 2,kmax
                 DO ij = 1,imax*jmax
                    rho(ij,1,k) = rho(ij,1,1)
                 ENDDO
              ENDDO
           ENDIF

! density profile itself is given
        ELSE
           ycenter = y(1) + scaley*ycoor_rho
           DO j = 1,jmax
              dummy =  FLOW_JET_TEMPORAL&
                   (iprof_rho, thick_rho, delta_rho, mean_rho, diam_rho, ycenter, prof_rho, y(j))
              DO k = 1,kmax
                 DO i = 1,imax
                    rho(i,j,k) = rho(i,j,k) + dummy
                 ENDDO
              ENDDO
           ENDDO

        ENDIF
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE DENSITY_MEAN
