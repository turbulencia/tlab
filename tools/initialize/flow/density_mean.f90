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
  USE DNS_GLOBAL,    ONLY : g, j1bc
  USE DNS_GLOBAL,    ONLY : imode_sim, imode_flow, imode_fdm, inb_scal, imax,jmax,kmax
  USE DNS_GLOBAL,    ONLY : iprof_tem, mean_tem, delta_tem, thick_tem, ycoor_tem, prof_tem, diam_tem, jet_tem
  USE DNS_GLOBAL,    ONLY : iprof_rho, mean_rho, delta_rho, thick_rho, ycoor_rho, prof_rho, diam_rho, jet_rho
  USE DNS_GLOBAL,    ONLY : iprof_u, mean_u, delta_u, thick_u, ycoor_u, prof_u, diam_u, jet_u
  USE DNS_GLOBAL,    ONLY : iprof_i, mean_i, delta_i, thick_i, ycoor_i, prof_i, diam_i, jet_i
  USE DNS_GLOBAL,    ONLY : buoyancy
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
        IF ( buoyancy%type .EQ. EQNS_NONE ) THEN

#define TEM_MEAN_LOC(i,j,k) wrk3d(i,j,k)
#define RHO_MEAN_LOC(i,j,k) txc(i,j,k)

! temperature/mixture profile are given
           IF ( iprof_rho .EQ. PROFILE_NONE ) THEN
              ycenter = y(1) + g(2)%scale*ycoor_tem
              DO j = 1,jmax
                 dummy =  FLOW_SHEAR_TEMPORAL&
                      (iprof_tem, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y(j))
                 TEM_MEAN_LOC(:,j,:) = dummy
              ENDDO

              DO is = 1,inb_scal
                 ycenter = y(1) + g(2)%scale*ycoor_i(is)
                 DO j = 1,jmax
                    dummy =  FLOW_SHEAR_TEMPORAL&
                         (iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i,y(j))
                    s(:,j,:,is) = dummy
                 ENDDO
              ENDDO

! define liquid content in AirWater case: (p,T) given
              IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
                 CALL THERMO_AIRWATER_PT(imax, jmax, kmax, s, p, TEM_MEAN_LOC(1,1,1))
              ENDIF

              CALL THERMO_THERMAL_DENSITY&
                   (imax, jmax, kmax, s, p, TEM_MEAN_LOC(1,1,1), RHO_MEAN_LOC(1,1,1))
              rho(:,:,:) = rho(:,:,:) + RHO_MEAN_LOC(:,:,:)

! density profile itself is given
           ELSE
              ycenter = y(1) + g(2)%scale*ycoor_rho
              DO j = 1,jmax
                 dummy =  FLOW_SHEAR_TEMPORAL&
                      (iprof_rho, thick_rho, delta_rho, mean_rho, ycenter, prof_rho, y(j))
                 rho(:,j,:) = rho(:,j,:) + dummy
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
              dummy = C_1_R /buoyancy%vector(2)
              rho(:,:,:) = rho(:,:,:) + txc(:,:,:) *dummy
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
        ycenter = y(1) + g(2)%scale*ycoor_tem
        DO j = 1,jmax
           dummy = FLOW_JET_TEMPORAL&
                (iprof_tem, thick_tem, delta_tem, mean_tem, diam_tem, ycenter, prof_tem, y(j))
           wrk3d(:,j,:) = dummy
        ENDDO

        DO is = 1,inb_scal
           ycenter = y(1) + g(2)%scale*ycoor_i(is)
           DO j = 1,jmax
              dummy =  FLOW_JET_TEMPORAL&
                   (iprof_i(is), thick_i(is), delta_i(is), mean_i(is), diam_i(is), ycenter, prof_i, y(j))
              s(:,j,:,is) = dummy
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
           rho_vi(:) = rho(1,:,1)

! Inflow profile of axial velocity
           ycenter = y(1) + g(2)%scale*ycoor_u
           DO j = 1,jmax
              u_vi(j) = FLOW_JET_TEMPORAL&
                   (iprof_u, thick_u, delta_u, mean_u, diam_u, ycenter, prof_u, y(j))
           ENDDO

! 2D distribution of density
           CALL FLOW_JET_SPATIAL_DENSITY(imax, jmax, iprof_tem, thick_tem, delta_tem, mean_tem, &
                ycoor_tem, diam_tem, jet_tem, iprof_u, thick_u, delta_u, mean_u, ycoor_u, diam_u,&
                jet_u, g(2)%scale, x, y, s,p,rho_vi(1),u_vi(1),aux1(1),rho,aux2(1),aux3(1),aux4(1)) 
                
           DO k = 2,kmax
              rho(:,:,k) = rho(:,:,1)
           ENDDO

! density profile itself is given
        ELSE
           ycenter = y(1) + g(2)%scale*ycoor_rho
           DO j = 1,jmax
              dummy =  FLOW_JET_TEMPORAL&
                   (iprof_rho, thick_rho, delta_rho, mean_rho, diam_rho, ycenter, prof_rho, y(j))
              rho(:,j,:) = rho(:,j,:) + dummy
           ENDDO

        ENDIF
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE DENSITY_MEAN
