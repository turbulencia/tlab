#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE AVG_SCAL_XZ(is, q,s, s_local, dsdx,dsdy,dsdz, tmp1,tmp2,tmp3, mean2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL, ONLY : g, area
  USE DNS_GLOBAL, ONLY : imode_eqns, imode_flow, idiffusion, itransport
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, inb_scal, inb_scal_array
  USE DNS_GLOBAL, ONLY : buoyancy, radiation, transport
  USE DNS_GLOBAL, ONLY : rbg, sbg, qbg
  USE DNS_GLOBAL, ONLY : rbackground, ribackground
  USE DNS_GLOBAL, ONLY : visc, schmidt, froude
  USE THERMO_GLOBAL, ONLY : imixture, thermo_param
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                           INTENT(IN)    :: is
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)    :: q, s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(IN)    :: s_local
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: dsdx,dsdy,dsdz, tmp1,tmp2,tmp3, wrk3d
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: mean2d, wrk1d
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk2d

  TARGET q, s

! -----------------------------------------------------------------------
  TINTEGER, PARAMETER :: MAX_VARS_GROUPS = 10
  TINTEGER i,j,k, bcs(2,2)
  TREAL diff, SIMPSON_NU, UPPER_THRESHOLD, LOWER_THRESHOLD
  TREAL delta_m, delta_w, delta_s, delta_s_area, delta_s_position, delta_s_value 
  TREAL smin_loc, smax_loc
  TREAL delta_hb01, delta_ht01, delta_h01
  TREAL delta_sb01, delta_st01, delta_s01
  TREAL mixing1, mixing2
  TINTEGER jloc_max(1)

  TINTEGER ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nmax
  
  TREAL dummy, coefT, coefR, coefQ

  TREAL VAUXPOS(11)
  TINTEGER ivauxpos
  CHARACTER*32 name, groupname(MAX_VARS_GROUPS)
  CHARACTER*250 line1, varname(MAX_VARS_GROUPS)
  CHARACTER*850 line2

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u,v,w, rho, vis

! ###################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero
  
 ! Define pointers
  u => q(:,:,:,1)
  v => q(:,:,:,2)
  w => q(:,:,:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL   )THEN
     rho => q(:,:,:,5)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)
  ENDIF

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

! Variable definition and memory management
! -----------------------------------------------------------------------
! Independent variables
  ig(1) = 1; ng = 1
  IF      ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
#define VAUXPRE1 mean2d(j,ig(1))
#define VAUXPRE2 mean2d(j,ig(1)+1)
#define VAUXPRE3 mean2d(j,ig(1)+2)
#define VAUXPRE4 mean2d(j,ig(1)+3)
#define VAUXPRE5 mean2d(j,ig(1)+4)
     sg(1) = 5

     varname(1) = 'Y SM SW SS SR'

  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
#define VAUXPRE1 mean2d(j,ig(1))
     sg(1) = 1

     varname(1) = 'Y'
     
  ENDIF

! -----------------------------------------------------------------------
! Dependent variables
  ig(2) = ig(1)+ sg(1); ng = ng + 1
#define rS(j)     mean2d(j,ig(2)  )
#define fS(j)     mean2d(j,ig(2)+1)
#define rS_y(j)   mean2d(j,ig(2)+2)
#define fS_y(j)   mean2d(j,ig(2)+3)
#define rQ(j)     mean2d(j,ig(2)+4)
#define fQ(j)     mean2d(j,ig(2)+5)
  sg(2) = 6
     
  groupname(2) = 'Mean'
  varname(2)   = 'rS fS rS_y fS_y rQ fQ'
  IF ( radiation%active(is) ) THEN
     varname(2) = TRIM(ADJUSTL(varname(2)))//' rQrad rQradC'
     sg(2) = sg(2) + 2
  ENDIF
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR .OR. imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     varname(2) = TRIM(ADJUSTL(varname(2)))//' rQeva'
     sg(2) = sg(2) + 1     
  ENDIF  
  IF ( transport%active(is) ) THEN
     varname(2) = TRIM(ADJUSTL(varname(2)))//' rQtra rQtraC'
     sg(2) = sg(2) + 2
  ENDIF
  
! -----------------------------------------------------------------------
  ig(3) = ig(2)+ sg(2); ng = ng + 1
#define Rsu(j)    mean2d(j,ig(3)  )
#define Rsv(j)    mean2d(j,ig(3)+1)
#define Rsw(j)    mean2d(j,ig(3)+2)
#define fS2(j)    mean2d(j,ig(3)+3)
#define fS3(j)    mean2d(j,ig(3)+4)
#define fS4(j)    mean2d(j,ig(3)+5)
#define fS5(j)    mean2d(j,ig(3)+6)
#define fS6(j)    mean2d(j,ig(3)+7)
#define rS2(j)    mean2d(j,ig(3)+8)
#define rS3(j)    mean2d(j,ig(3)+9)
#define rS4(j)    mean2d(j,ig(3)+10)
#define rS5(j)    mean2d(j,ig(3)+11)
#define rS6(j)    mean2d(j,ig(3)+12)
#define fSRHO(j)  mean2d(j,ig(3)+13)
#define rSRHO(j)  mean2d(j,ig(3)+14)
  sg(3) = 15

  groupname(3) = 'Fluctuations'
  varname(3)   = 'Rsu Rsv Rsw fS2 fS3 fS4 fS5 fS6 rS2 rS3 rS4 rS5 rS6 fRS rRS'
  
! -----------------------------------------------------------------------
  ig(4) = ig(3)+ sg(3); ng = ng + 1
#define Rss_t(j)  mean2d(j,ig(4)  )
#define Pss(j)    mean2d(j,ig(4)+1)
#define Ess(j)    mean2d(j,ig(4)+2)
#define Tssy(j)   mean2d(j,ig(4)+3)
#define Tssy_y(j) mean2d(j,ig(4)+4)
#define Dss(j)    mean2d(j,ig(4)+5)
#define Qss(j)    mean2d(j,ig(4)+6)
  sg(4) = 7

  groupname(4) = 'fS2Budget'
  varname(4)   = 'Rss_t Pss Ess Tssy Tssy_y Dss Qss'
  
! -----------------------------------------------------------------------
  ig(5) = ig(4)+ sg(4); ng = ng + 1
#define var_x(j)  mean2d(j,ig(5)  )
#define var_y(j)  mean2d(j,ig(5)+1)
#define var_z(j)  mean2d(j,ig(5)+2)
#define skew_x(j) mean2d(j,ig(5)+3)
#define skew_y(j) mean2d(j,ig(5)+4)
#define skew_z(j) mean2d(j,ig(5)+5)
#define flat_x(j) mean2d(j,ig(5)+6)
#define flat_y(j) mean2d(j,ig(5)+7)
#define flat_z(j) mean2d(j,ig(5)+8)
  sg(5) = 9

  groupname(5) = 'DerivativeFluctuations'
  varname(5)   = 'VarSx VarSy VarSz SkewSx SkewSy SkewSz FlatSx FlatSy FlatSz'

#define L_AVGMAX 47
  
! -----------------------------------------------------------------------
! Auxiliary variables depending on y and t; this last group is not written
  ig(6) = ig(5)+ sg(5); ng = ng + 1
#define rR(j)     mean2d(j,ig(6))
#define fU(j)     mean2d(j,ig(6)+1)
#define fV(j)     mean2d(j,ig(6)+2)
#define fW(j)     mean2d(j,ig(6)+3)
#define fU_y(j)   mean2d(j,ig(6)+4)
#define fV_y(j)   mean2d(j,ig(6)+5)
#define fW_y(j)   mean2d(j,ig(6)+6)
#define F1(j)     mean2d(j,ig(6)+7)
#define F2(j)     mean2d(j,ig(6)+8)
#define F3(j)     mean2d(j,ig(6)+9)
#define F2_y(j)   mean2d(j,ig(6)+10)
#define rR2(j)    mean2d(j,ig(6)+11)
 sg(6) = 12

! -----------------------------------------------------------------------
  nmax = ig(ng) +sg(ng) -1
  IF ( MAX_AVG_TEMPORAL .LT. nmax ) THEN
     CALL IO_WRITE_ASCII(efile,'AVERAGES_SCAL_XZ. Not enough space in local arrays.')
     CALL DNS_STOP(LES_ERROR_AVGTMP)
  ENDIF
  mean2d(:,1:nmax) = C_0_R

  ng   = ng -1
  nmax = ig(ng) +sg(ng) -1 ! the last group is not written out

  IF ( L_AVGMAX .LT. nmax ) THEN
     CALL IO_WRITE_ASCII(efile,'AVERAGES_SCAL_XZ. Not enough space in format definition.')
     CALL DNS_STOP(LES_ERROR_AVGTMP)
  ENDIF

! #######################################################################
  WRITE(line1,*) itime; line1 = 'Calculating scal statistics at It'//TRIM(ADJUSTL(line1))//'...'
  CALL IO_WRITE_ASCII(lfile,line1)

! #######################################################################
! Preliminary data
! #######################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     CALL AVG_IK_V(imax,jmax,kmax, jmax, u, g(1)%jac,g(3)%jac, fU(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, v, g(1)%jac,g(3)%jac, fV(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, w, g(1)%jac,g(3)%jac, fW(1), wrk1d(1,1), area)
     rR(:) = C_1_R
     
  ELSE
     dsdx = rho *u
     dsdy = rho *v
     dsdz = rho *w
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdx, g(1)%jac,g(3)%jac, fU(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdy, g(1)%jac,g(3)%jac, fV(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdz, g(1)%jac,g(3)%jac, fW(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, rho,  g(1)%jac,g(3)%jac, rR(1), wrk1d(1,1), area)
     fU(:) = fU(:) /rR(:)
     fV(:) = fV(:) /rR(:)
     fW(:) = fW(:) /rR(:)
     fS(:) = fS(:) /rR(:)
     
  ENDIF

! #######################################################################
! Moments
! #######################################################################
  CALL AVG_IK_V(imax,jmax,kmax, jmax, s_local, g(1)%jac,g(3)%jac, rS(1), wrk1d(1,1), area)
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     fS(:) = rS(:)
  ELSE
     wrk3d = rho *u
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, fS(1), wrk1d(1,1), area)
     fS(:) = fS(:) /rR(:)
  ENDIF

! -----------------------------------------------------------------------
  DO j = 1,jmax
     wrk3d(:,j,:) = s_local(:,j,:) - rS(j)
  ENDDO
  tmp1 = wrk3d *wrk3d
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS2(1), wrk1d(1,1), area)
  tmp1 = wrk3d *tmp1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS3(1), wrk1d(1,1), area)
  tmp1 = wrk3d *tmp1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS4(1), wrk1d(1,1), area)
  tmp1 = wrk3d *tmp1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS5(1), wrk1d(1,1), area)
  tmp1 = wrk3d *tmp1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS6(1), wrk1d(1,1), area)

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     fS2(:) = rS2(:)
     fS3(:) = rS3(:)
     fS4(:) = rS4(:)
     fS5(:) = rS5(:)
     fS6(:) = rS6(:)
     
  ELSE
     DO j = 1,jmax
        wrk3d(:,j,:) = s_local(:,j,:) - fS(j)
     ENDDO
     tmp1 = wrk3d *wrk3d *rho
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS2(1), wrk1d(1,1), area)
     tmp1 = wrk3d *tmp1
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS3(1), wrk1d(1,1), area)
     tmp1 = wrk3d *tmp1
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS4(1), wrk1d(1,1), area)
     tmp1 = wrk3d *tmp1
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS5(1), wrk1d(1,1), area)
     tmp1 = wrk3d *tmp1
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS6(1), wrk1d(1,1), area)
     fS2(:) = fS2(:) /rR(:)
     fS3(:) = fS3(:) /rR(:)
     fS4(:) = fS4(:) /rR(:)
     fS5(:) = fS5(:) /rR(:)
     fS6(:) = fS6(:) /rR(:)

  ENDIF

! #######################################################################
! Correlations
! #######################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     DO j = 1,jmax
        wrk3d(:,j,:) = s_local(:,j,:) - rS(j)
     ENDDO
     DO j = 1,jmax
        tmp1(:,j,:) = wrk3d(:,j,:) *(u(:,j,:)-fU(j))
        tmp2(:,j,:) = wrk3d(:,j,:) *(v(:,j,:)-fV(j))
        tmp3(:,j,:) = wrk3d(:,j,:) *(w(:,j,:)-fW(j))
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Rsu(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, Rsv(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, Rsw(1), wrk1d(1,1), area)

  ELSE
     DO j = 1,jmax
        wrk3d(:,j,:) =(s_local(:,j,:) - fS(j)) *rho(:,j,:)
     ENDDO
     DO j = 1,jmax
        tmp1(:,j,:) = wrk3d(:,j,:) *(u(:,j,:)-fU(j))
        tmp2(:,j,:) = wrk3d(:,j,:) *(v(:,j,:)-fV(j))
        tmp3(:,j,:) = wrk3d(:,j,:) *(w(:,j,:)-fW(j))
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Rsu(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, Rsv(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, Rsw(1), wrk1d(1,1), area)
     Rsu(:) = Rsu(:) /rR(:)
     Rsv(:) = Rsv(:) /rR(:)
     Rsw(:) = Rsw(:) /rR(:)

  ENDIF

! -----------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     rSRHO(:) = C_0_R
     fSRHO(:) = C_0_R
     
  ELSE
     DO j = 1,jmax
        tmp1(:,j,:) = (rho(:,j,:) - rR(j)) *(s_local(:,j,:) - rS(j))
        tmp2(:,j,:) = (rho(:,j,:) - rR(j)) *(s_local(:,j,:) - fS(j))
        tmp3(:,j,:) = (rho(:,j,:) - rR(j)) *(rho(:,j,:) - rR(j))
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rSRHO(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, fSRHO(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, rR2(1),   wrk1d(1,1), area)

     DO j = 1,jmax
        IF ( rR2(j) .GT. C_0_R ) THEN
           IF ( rS2(j) .GT. C_0_R ) THEN; rSRHO(j) = rSRHO(j)/sqrt(rR2(j)*rS2(j))
           ELSE;                          rSRHO(j) = C_2_R; ENDIF
           
           IF ( fS2(j) .GT. C_0_R ) THEN; fSRHO(j) = fSRHO(j)/sqrt(rR2(j)*fS2(j))
           ELSE;                          fSRHO(j) = C_2_R; ENDIF
           
        ELSE
           rSRHO(j) = C_BIG_R
           fSRHO(j) = C_BIG_R
        
        ENDIF
     ENDDO
  ENDIF

! #######################################################################
! Source terms
! #######################################################################
  dsdx = C_0_R; dsdy = C_0_R; dsdz = C_0_R; tmp1 = C_0_R; tmp2 = C_0_R; tmp3 = C_0_R

  IF ( radiation%active(is) ) THEN ! Radiation in tmp1
     IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, s(1,1,1,radiation%scalar(is)), tmp2)
        CALL OPR_RADIATION     (radiation, imax,jmax,kmax, g(2), tmp2,                          tmp1, wrk1d,wrk3d)
        CALL OPR_RADIATION_FLUX(radiation, imax,jmax,kmax, g(2), tmp2,                          dsdx, wrk1d,wrk3d)
        CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, ribackground, tmp1)
        tmp2 = C_0_R
        
     ELSE
        CALL OPR_RADIATION     (radiation, imax,jmax,kmax, g(2), s(1,1,1,radiation%scalar(is)), tmp1, wrk1d,wrk3d)
        CALL OPR_RADIATION_FLUX(radiation, imax,jmax,kmax, g(2), s(1,1,1,radiation%scalar(is)), dsdx, wrk1d,wrk3d)
     ENDIF
  ENDIF

  IF ( transport%active(is) ) THEN ! Transport in tmp3
     CALL FI_TRANSPORT     (transport, i1, imax,jmax,kmax, is, s,tmp3, dsdy, wrk2d,wrk3d)
     CALL FI_TRANSPORT_FLUX(transport,     imax,jmax,kmax, is, s,dsdz)
  ENDIF

  IF ( is .GT. inb_scal ) THEN     ! Diagnostic variables; I overwrite tmp1 and dsdx and recalculate them.
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        coefQ = C_1_R                        ! Coefficient in the evaporation term
        coefR = C_0_R                        ! Coefficient in the radiation term
        coefT = C_0_R                        ! Coefficient in the transport term
        IF ( is .EQ. inb_scal_array+1 ) THEN ! Default values are for liquid; defining them for buoyancy
           coefQ = buoyancy%parameters(inb_scal_array) /froude
           coefR = buoyancy%parameters(inb_scal) /froude
           DO i = 1,inb_scal; coefT = coefT + transport%parameters(i) *buoyancy%parameters(i) /froude; ENDDO
        ENDIF
     
        CALL THERMO_AIRWATER_LINEAR_SOURCE(imax,jmax,kmax, s, dsdx,dsdy,dsdz) ! calculate xi in dsdx
        CALL FI_GRADIENT(imax,jmax,kmax, dsdx,tmp2, tmp1, wrk2d,wrk3d)
        
        dummy=-diff *coefQ
        tmp2 = dsdz *tmp2 *dummy         ! evaporation source
        
        IF ( transport%active(is) .OR. radiation%active(is) ) THEN ! preparing correction terms into dsdz
           CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), dsdx,tmp1, wrk3d, wrk2d,wrk3d)
           dsdz = dsdz *tmp1
        ENDIF
        
        IF ( radiation%active(is) ) THEN ! radiation source; needs dsdy
           CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,1,1,radiation%scalar(is)), tmp1, wrk1d,wrk3d)
           dummy= thermo_param(2) *coefQ
           tmp1 = tmp1 *( coefR + dsdy *dummy  )
! Correction term needs dsdz
           CALL OPR_RADIATION_FLUX(radiation, imax,jmax,kmax, g(2), s(1,1,1,radiation%scalar(is)), dsdx, wrk1d,wrk3d)
           dsdx = dsdx *dsdz *dummy
        ELSE
           tmp1 = C_0_R; dsdx = C_0_R
        ENDIF

        IF ( transport%active(is) ) THEN ! transport source; needs dsdy
           dummy= coefQ
           tmp3 = tmp3 *( coefT + dsdy *dummy )
! Correction term needs dsdz
           dummy= dummy
           dsdy =-s(:,:,:,transport%scalar(1))**(C_1_R+transport%auxiliar(1))
           dsdz = dsdy *dsdz *dummy
        ELSE
           dsdz = C_0_R
        ENDIF
        
     ELSE
        IF ( buoyancy%type .NE. EQNS_EXPLICIT ) THEN
           CALL FI_GRADIENT(imax,jmax,kmax, s,dsdx, dsdy, wrk2d,wrk3d)
           CALL FI_BUOYANCY_SOURCE(buoyancy, isize_field, s, dsdx, wrk3d) ! dsdx contains gradient
           tmp1 = wrk3d* diff /froude
        ENDIF
        
     ENDIF

  ENDIF

! -----------------------------------------------------------------------
! Calculating averages
  k = ig(2)+5
  IF ( radiation%active(is) ) THEN
     k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d(1,1), area)
     k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdx, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d(1,1), area) ! correction term or flux
  ENDIF
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR .OR. imixture .EQ. MIXT_TYPE_AIRWATER ) THEN ! evaporation
     k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d(1,1), area)
  ENDIF
  IF ( transport%active(is) ) THEN
     k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d(1,1), area)
     k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdz, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d(1,1), area) ! correction term or flux
  ENDIF

  wrk3d = tmp1 + tmp2 + tmp3 ! total
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, rQ(1), wrk1d(1,1), area)
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     fQ(:) = rQ(:)

     DO j = 1,jmax
        wrk3d(:,j,:) =(s_local(:,j,:) - rS(j)) *wrk3d(:,j,:)
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Qss(1), wrk1d(1,1), area) ! transport equation
     Qss(:) = Qss(:) *C_2_R

  ELSE     
     DO j = 1,jmax
        wrk3d(:,j,:) = wrk3d(:,j,:) *rho(:,j,:)
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, fQ(1), wrk1d(1,1), area)
     fQ(:) = fQ(:) /rR(:)

     DO j = 1,jmax
        wrk3d(:,j,:) =(s_local(:,j,:) - fS(j)) *wrk3d(:,j,:)
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Qss(1), wrk1d(1,1), area)
     Qss(:) = Qss(:) *C_2_R /rR(:)

  ENDIF

! #######################################################################
! Derivative terms
! #######################################################################
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s_local, dsdx,    wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s_local, dsdy,    wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s_local, dsdz,    wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1,     bcs, g(2), rS(1),   rS_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1,     bcs, g(2), fS(1),   fS_y(1), wrk3d, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Derivatives Fluctuations
  tmp1 = dsdx *dsdx
  DO j = 1,jmax
     tmp2(:,j,:) =(dsdy(:,j,:) - rS_y(j)) *(dsdy(:,j,:) - rS_y(j))
  ENDDO
  tmp3 = dsdz *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, var_x(1), wrk1d(1,1), area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, var_y(1), wrk1d(1,1), area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, var_z(1), wrk1d(1,1), area)

  tmp1 = tmp1 *dsdx
  DO j = 1,jmax
     tmp2(:,j,:) =  tmp2(:,j,:) *(dsdy(:,j,:) - rS_y(j))
  ENDDO
  tmp3 = tmp3 *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, skew_x(1), wrk1d(1,1), area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, skew_y(1), wrk1d(1,1), area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, skew_z(1), wrk1d(1,1), area)
  
  tmp1 = tmp1 *dsdx
  DO j = 1,jmax
     tmp2(:,j,:) =  tmp2(:,j,:) *(dsdy(:,j,:) - rS_y(j))
  ENDDO
  tmp3 = tmp3 *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, flat_x(1), wrk1d(1,1), area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, flat_y(1), wrk1d(1,1), area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, flat_z(1), wrk1d(1,1), area)

  DO j = 1,jmax
     IF ( var_x(j) .GT. C_SMALL_R ) THEN
        skew_x(j) = skew_x(j) / var_x(j)**C_1_5_R
        flat_x(j) = flat_x(j) / var_x(j)**C_2_R
     ELSE
        skew_x(j) = C_BIG_R
        flat_x(j) = C_BIG_R
     ENDIF
     IF ( var_y(j) .GT. C_SMALL_R ) THEN
        skew_y(j) = skew_y(j) / var_y(j)**C_1_5_R
        flat_y(j) = flat_y(j) / var_y(j)**C_2_R
     ELSE
        skew_y(j) = C_BIG_R
        flat_y(j) = C_BIG_R
     ENDIF
     IF ( var_z(j) .GT. C_SMALL_R ) THEN
        skew_z(j) = skew_z(j) / var_z(j)**C_1_5_R
        flat_z(j) = flat_z(j) / var_z(j)**C_2_R
     ELSE
        skew_z(j) = C_BIG_R
        flat_z(j) = C_BIG_R
     ENDIF
  ENDDO

! -----------------------------------------------------------------------
! Variance transport equation; dissipation
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     tmp1 = vis *dsdx
     tmp2 = vis *dsdy
     tmp3 = vis *dsdz
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, F1(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, F2(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, F3(1), wrk1d(1,1), area)
  ELSE
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdx, g(1)%jac,g(3)%jac, F1(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdy, g(1)%jac,g(3)%jac, F2(1), wrk1d(1,1), area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdz, g(1)%jac,g(3)%jac, F3(1), wrk1d(1,1), area)
  ENDIF

  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
     DO j = 1,jmax
        wrk3d(:,j,:) = (tmp1(:,j,:) - F1(j)) * dsdx(:,j,:)          &
                     + (tmp2(:,j,:) - F2(j)) *(dsdy(:,j,:)-fS_y(j)) &
                     + (tmp3(:,j,:) - F3(j)) * dsdz(:,j,:)
     ENDDO
  ELSE
     DO j = 1,jmax
        wrk3d(:,j,:) = (dsdx(:,j,:) - F1(j)) * dsdx(:,j,:)          &
                     + (dsdy(:,j,:) - F2(j)) *(dsdy(:,j,:)-fS_y(j)) &
                     + (dsdz(:,j,:) - F3(j)) * dsdz(:,j,:)
     ENDDO
  ENDIF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Ess(1), wrk1d(1,1), area)
  Ess(:) = Ess(:) *C_2_R *diff /rR(:)

! -----------------------------------------------------------------------
! Variance transport equation; turbulent transport and diffusion
  DO j = 1,jmax
     wrk3d(:,j,:) = s_local(:,j,:)-fS(j)
  ENDDO
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     DO j = 1,jmax
        tmp1(:,j,:) = wrk3d(:,j,:)* ( (v(:,j,:)-fV(j))*wrk3d(:,j,:) - C_2_R *diff *(dsdy(:,j,:)-F2(j)) )
     ENDDO     
  ELSE
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        DO j = 1,jmax
           tmp1(:,j,:) = wrk3d(:,j,:)* ( rho(:,j,:)*(v(:,j,:)-fV(j))*wrk3d(:,j,:) - C_2_R *diff *(vis(:,j,:)*dsdy(:,j,:)-F2(j)) )
        ENDDO
        
     ELSE
        DO j = 1,jmax
           tmp1(:,j,:) = wrk3d(:,j,:)* ( rho(:,j,:)*(v(:,j,:)-fV(j))*wrk3d(:,j,:) - C_2_R *diff *(           dsdy(:,j,:)-F2(j)) )
        ENDDO
     ENDIF
     
  ENDIF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1,  g(1)%jac,g(3)%jac, Tssy(1), wrk1d(1,1), area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Dss(1),  wrk1d(1,1), area)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tssy(1), Tssy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), F2(1),   F2_y(1),   wrk3d, wrk2d,wrk3d)
  Dss(:) = Dss(:) *F2_y(j) *C_2_R *diff

! -----------------------------------------------------------------------
  Pss(:)   =-C_2_R *Rsv(:) *fS_y(:)
  Rss_t(:) = Pss(:) - Ess(:) + (Dss(:) + Qss(:) -Tssy_y(:)) /rR(:)

! #######################################################################
! Global quantites
! #######################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN

! -----------------------------------------------------------------------
! Based on delta_u
! -----------------------------------------------------------------------
     IF ( ABS(qbg(1)%delta) .GT. C_SMALL_R ) THEN
        CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fU(1), fU_y(1), wrk3d, wrk2d,wrk3d)
        delta_w = ABS(qbg(1)%delta)/MAXVAL(ABS(fU_y(1:jmax)))

        DO j=1, jmax
           wrk1d(j,1) = rR(j)*(C_025_R - (fU(j)/qbg(1)%delta)**2)
        ENDDO
        delta_m = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)/rbg%mean

        DO j=1, jmax
           wrk1d(j,1) = ( diff*F2(j) -  rR(j)*Rsv(j) )*fS_y(j)
        ENDDO
        delta_s_area = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)*C_2_R/(rbg%mean*qbg(1)%delta)

     ELSE
        delta_m      = C_1_R
        delta_w      = C_1_R
        delta_s_area = C_1_R

     ENDIF

! -----------------------------------------------------------------------
! Based on rbg%delta 
! -----------------------------------------------------------------------
     IF ( ABS(rbg%delta) .GT. C_SMALL_R ) THEN
        dummy = rbg%mean + (C_05_R-C_1EM2_R)*rbg%delta
        delta_hb01 = LOWER_THRESHOLD(jmax, dummy, rR(1), g(2)%nodes)
        dummy = rbg%mean - (C_05_R-C_1EM2_R)*rbg%delta
        delta_ht01 = UPPER_THRESHOLD(jmax, dummy, rR(1), g(2)%nodes)
        delta_h01 = delta_ht01 - delta_hb01

     ELSE
        delta_h01 = C_1_R

     ENDIF

! -----------------------------------------------------------------------
! Based on scalar
! -----------------------------------------------------------------------
     IF ( ABS(sbg(is)%delta) .GT. C_SMALL_R ) THEN
        jloc_max = MAXLOC(ABS(fS_y(1:jmax))); j = jloc_max(1)

        wrk1d(1:jmax,1) = fS(jmax) - fS(1:jmax)
        delta_s_area    = SIMPSON_NU(jmax-j+1, wrk1d(j,1), g(2)%nodes(j)) / ABS(sbg(is)%delta)
 
        delta_s          = ABS(sbg(is)%delta)/ABS(fS_y(j))
        delta_s_position = g(2)%nodes(j)
        delta_s_value    = fS(j)

! Thickness
        dummy      = sbg(is)%mean + (C_05_R-C_1EM3_R) *sbg(is)%delta
        delta_sb01 = LOWER_THRESHOLD(jmax, dummy, rS(1), g(2)%nodes)
        dummy      = sbg(is)%mean - (C_05_R-C_1EM3_R) *sbg(is)%delta
        delta_st01 = UPPER_THRESHOLD(jmax, dummy, rS(1), g(2)%nodes)
        
        delta_sb01 =              (g(2)%nodes(1) + sbg(is)%ymean *g(2)%scale) - delta_sb01
        delta_st01 = delta_st01 - (g(2)%nodes(1) + sbg(is)%ymean *g(2)%scale)
        delta_s01  = delta_st01 + delta_sb01

!  Mixing, Youngs' definition
        smin_loc = sbg(is)%mean - C_05_R*ABS(sbg(is)%delta)
        smax_loc = sbg(is)%mean + C_05_R*ABS(sbg(is)%delta)
        wrk3d = (s_local - smin_loc) *(smax_loc -s_local)
        CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d, wrk1d(1,2), area)
        mixing1 = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)
        DO j = 1,jmax
           wrk1d(j,1)=(rS(j)-smin_loc)*(smax_loc-rS(j))
        ENDDO
        mixing1 = mixing1/SIMPSON_NU(jmax, wrk1d, g(2)%nodes)

! Mixing, Cook's definition
        smin_loc = sbg(is)%mean - C_05_R*ABS(sbg(is)%delta)
        smax_loc = sbg(is)%mean + C_05_R*ABS(sbg(is)%delta)
        DO k = 1,kmax
           DO i = 1,imax*jmax
              wrk3d(i,1,k) = MIN(s_local(i,1,k)-smin_loc,smax_loc-s_local(i,1,k))
           ENDDO
        ENDDO
        CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d, wrk1d(1,2), area)
        mixing2 = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)
        DO j = 1,jmax
           wrk1d(j,1) = MIN(rS(j)-smin_loc,smax_loc-rS(j))
        ENDDO
        mixing2 = mixing2/SIMPSON_NU(jmax, wrk1d, g(2)%nodes)

     ELSE
        delta_sb01 = C_1_R
        delta_st01 = C_1_R
        delta_s01  = C_1_R
        delta_s    = C_1_R
        mixing1    = C_1_R
        mixing2    = C_1_R

     ENDIF

! Independent variables
     DO j = 1,jmax            
        VAUXPRE1 =  g(2)%nodes(j)
        VAUXPRE2 = (g(2)%nodes(j)-g(2)%scale *qbg(1)%ymean  - g(2)%nodes(1))/delta_m
        VAUXPRE3 = (g(2)%nodes(j)-g(2)%scale *qbg(1)%ymean  - g(2)%nodes(1))/delta_w
        VAUXPRE4 = (g(2)%nodes(j)-g(2)%scale *sbg(is)%ymean - g(2)%nodes(1))/delta_s01
        VAUXPRE5 = (g(2)%nodes(j)-g(2)%scale *rbg%ymean     - g(2)%nodes(1))/delta_h01
     ENDDO
     
! -----------------------------------------------------------------------
! Jet
! -----------------------------------------------------------------------
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN ! Not developed yet
     DO j = 1,jmax            
        VAUXPRE1 = g(2)%nodes(j)
     ENDDO
     
  ENDIF

! #######################################################################
! Output
! #######################################################################
#define LOC_UNIT_ID 23
#define LOC_STATUS 'unknown'

  WRITE(line1,*) is; line1='avg'//TRIM(ADJUSTL(line1))//'s'
  WRITE(name,*) itime; name=TRIM(ADJUSTL(line1))//TRIM(ADJUSTL(name))

! -----------------------------------------------------------------------
! TkStat file; header
! -----------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

! #include "dns_open_file.h"
#ifdef USE_RECLEN
     OPEN(UNIT=LOC_UNIT_ID, RECL=1050, FILE=name, STATUS='unknown') ! this is probably outdated
#else
     OPEN(UNIT=LOC_UNIT_ID, FILE=name, STATUS='unknown')
#endif

! Header
     WRITE(LOC_UNIT_ID, '(A8,E14.7E3)') 'RTIME = ', rtime
     ! WRITE(LOC_UNIT_ID, '(A7,I8)') 'IMAX = ', i1
     ! WRITE(LOC_UNIT_ID, '(A7,I8)') 'JMAX = ', jmax
     ! WRITE(LOC_UNIT_ID) 'RTIME = ', rtime
     ! WRITE(LOC_UNIT_ID) 'IMAX = ', i1
     ! WRITE(LOC_UNIT_ID) 'JMAX = ', jmax

! Independent variables     
     line2 = 'I J '//TRIM(ADJUSTL(varname(1)))
     
! Dependent variables depending on y and t
     DO k = 2,ng
        WRITE(LOC_UNIT_ID,1010) 'GROUP = '//TRIM(ADJUSTL(groupname(k)))//' '//TRIM(ADJUSTL(varname(k)))
        ! WRITE(LOC_UNIT_ID) 'GROUP = '//TRIM(ADJUSTL(groupname(k)))//' '//TRIM(ADJUSTL(varname(k)))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(varname(k)))
     ENDDO

! Dependent variables dependent on t only
     IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        line1 = 'Delta_m Delta_w Delta_s Delta_s_Area Delta_s_Position Delta_s_Value'
        WRITE(LOC_UNIT_ID,1010) 'GROUP = ShearThicknesses '//TRIM(ADJUSTL(line1))
        ! WRITE(LOC_UNIT_ID) 'GROUP = ShearThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

        line1 = 'Delta_sb01 Delta_st01 Delta_s01 mixing_Youngs mixing_Cook'
        WRITE(LOC_UNIT_ID,1010) 'GROUP = MixingThicknesses '//TRIM(ADJUSTL(line1))
        ! WRITE(LOC_UNIT_ID) 'GROUP = MixingThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     ENDIF

     WRITE(LOC_UNIT_ID,1010) TRIM(ADJUSTL(line2))
     ! WRITE(LOC_UNIT_ID) TRIM(ADJUSTL(line2))

! Body     
     DO j = 1,jmax            

        ivauxpos = 0        
        IF ( j .EQ. jmax/2 ) THEN
           IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
              ivauxpos = 11
              VAUXPOS(1) = delta_m
              VAUXPOS(2) = delta_w
              VAUXPOS(3) = delta_s
              VAUXPOS(4) = delta_s_area
              VAUXPOS(5) = delta_s_position
              VAUXPOS(6) = delta_s_value 
              VAUXPOS(7) = delta_sb01
              VAUXPOS(8) = delta_st01
              VAUXPOS(9) = delta_s01
              VAUXPOS(10)= mixing1
              VAUXPOS(11)= mixing2
           ENDIF
        ENDIF

        WRITE(LOC_UNIT_ID,1020) 1, j, (mean2d(j,k),k=1,nmax), (VAUXPOS(k),k=1,ivauxpos)
        ! WRITE(LOC_UNIT_ID) mean2d(1:jmax,1:nmax)

     ENDDO

     CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1010 FORMAT(A)
1020 FORMAT(I5,(1X,I5),L_AVGMAX(1X,G_FORMAT_R),11(1X,G_FORMAT_R))

END SUBROUTINE AVG_SCAL_XZ

