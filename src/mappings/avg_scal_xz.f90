#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# Assumes statistical homogeneity in xOz, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# In the incompressible case, the array p has been
!# pointed to dudz and the pressure field is stored there; do not
!# use array tmp3 until pressure block
!#
!########################################################################

SUBROUTINE AVG_SCAL_XZ(is, q,s, s_local, dsdx,dsdy,dsdz, tmp1,tmp2,tmp3, mean2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture, thermo_param
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                           INTENT(IN   ) :: is
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN   ) :: q, s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(IN   ) :: s_local
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: dsdx,dsdy,dsdz, tmp1,tmp2,tmp3, wrk3d
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: mean2d, wrk1d
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk2d

  TARGET q, s, tmp3

  ! -----------------------------------------------------------------------
  TINTEGER, PARAMETER :: MAX_VARS_GROUPS = 10
  TINTEGER i,j,k, bcs(2,2)
  TREAL diff, dummy, coefT, coefR, coefQ, c23

  TINTEGER ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nmax

  CHARACTER*32 name, groupname(MAX_VARS_GROUPS)
  CHARACTER*250 line1, varname(MAX_VARS_GROUPS)
  CHARACTER*850 line2

  TINTEGER nvar
  CHARACTER(LEN=32) varname2(MAX_AVG_TEMPORAL)

  ! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u,v,w,p, rho, vis

  ! ###################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  ! Define pointers
  u => q(:,:,:,1)
  v => q(:,:,:,2)
  w => q(:,:,:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
    rho => q(:,:,:,5)
    p   => q(:,:,:,6)
    IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)
  ELSE
    p   => tmp3
  END IF

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc /schmidt(is)
  END IF

  c23 = C_2_R /C_3_R

  ! -----------------------------------------------------------------------
  ! Dependent variables
  ng = 1; ig(ng) = 1
#define rS(j)     mean2d(j,ig(1)  )
#define fS(j)     mean2d(j,ig(1)+1)
#define rS_y(j)   mean2d(j,ig(1)+2)
#define fS_y(j)   mean2d(j,ig(1)+3)
#define rQ(j)     mean2d(j,ig(1)+4)
#define fQ(j)     mean2d(j,ig(1)+5)
  sg(ng) = 6

  groupname(ng) = 'Mean'
  varname(ng)   = 'rS fS rS_y fS_y rQ fQ'
  IF ( radiation%active(is) ) THEN
    varname(ng) = TRIM(ADJUSTL(varname(ng)))//' rQrad rQradC'
    sg(ng) = sg(ng) + 2
  END IF
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR .OR. imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
    varname(ng) = TRIM(ADJUSTL(varname(ng)))//' rQeva'
    sg(ng) = sg(ng) + 1
  END IF
  IF ( transport%active(is) ) THEN
    varname(ng) = TRIM(ADJUSTL(varname(ng)))//' rQtra rQtraC'
    sg(ng) = sg(ng) + 2
  END IF

  ! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rsu(j)    mean2d(j,ig(2)  )
#define Rsv(j)    mean2d(j,ig(2)+1)
#define Rsw(j)    mean2d(j,ig(2)+2)
#define fS2(j)    mean2d(j,ig(2)+3)
#define fS3(j)    mean2d(j,ig(2)+4)
#define fS4(j)    mean2d(j,ig(2)+5)
#define rS2(j)    mean2d(j,ig(2)+6)
#define rS3(j)    mean2d(j,ig(2)+7)
#define rS4(j)    mean2d(j,ig(2)+8)
  sg(ng) = 9

  groupname(ng) = 'Fluctuations'
  varname(ng)   = 'Rsu Rsv Rsw fS2 fS3 fS4 rS2 rS3 rS4'

  ! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rss_t(j)  mean2d(j,ig(3)  )
#define Css(j)    mean2d(j,ig(3)+1)
#define Pss(j)    mean2d(j,ig(3)+2)
#define Ess(j)    mean2d(j,ig(3)+3)
#define Tssy1(j)  mean2d(j,ig(3)+4)
#define Tssy2(j)  mean2d(j,ig(3)+5)
#define Tssy_y(j) mean2d(j,ig(3)+6)
#define Dss(j)    mean2d(j,ig(3)+7)
#define Qss(j)    mean2d(j,ig(3)+8)
  sg(ng) = 9

  groupname(ng) = 'RssBudget'
  varname(ng)   = 'Rss_t Css Pss Ess Tssy1 Tssy2 Tssy_y Dss Qss'

  ! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rsv_t(j)  mean2d(j,ig(4)  )
#define Csv(j)    mean2d(j,ig(4)+1)
#define Psv(j)    mean2d(j,ig(4)+2)
#define Esv(j)    mean2d(j,ig(4)+3)
#define PIsv(j)   mean2d(j,ig(4)+4)
#define Tsvy1(j)  mean2d(j,ig(4)+5)
#define Tsvy2(j)  mean2d(j,ig(4)+6)
#define Tsvy3(j)  mean2d(j,ig(4)+7)
#define Tsvy_y(j) mean2d(j,ig(4)+8)
#define Dsv(j)    mean2d(j,ig(4)+9)
#define Gsv(j)    mean2d(j,ig(4)+10)
#define Bsv(j)    mean2d(j,ig(4)+11)
#define Fsv(j)    mean2d(j,ig(4)+12)
#define Qsv(j)    mean2d(j,ig(4)+13)
  sg(ng) = 14

  groupname(ng) = 'RsvBudget'
  varname(ng)   = 'Rsv_t Csv Psv Esv PIsv Tsvy1 Tsvy2 Tsvy3 Tsvy_y Dsv Gsv Bsv Fsv Qsv'

  ! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define S_x2(j)  mean2d(j,ig(5)  )
#define S_y2(j)  mean2d(j,ig(5)+1)
#define S_z2(j)  mean2d(j,ig(5)+2)
#define S_x3(j)  mean2d(j,ig(5)+3)
#define S_y3(j)  mean2d(j,ig(5)+4)
#define S_z3(j)  mean2d(j,ig(5)+5)
#define S_x4(j)  mean2d(j,ig(5)+6)
#define S_y4(j)  mean2d(j,ig(5)+7)
#define S_z4(j)  mean2d(j,ig(5)+8)
  sg(ng) = 9

  groupname(ng) = 'DerivativeFluctuations'
  varname(ng)   = 'S_x2 S_y2 S_z2 S_x3 S_y3 S_z3 S_x4 S_y4 S_z4'

#define L_AVGMAX 38

  ! -----------------------------------------------------------------------
  ! Auxiliary variables depending on y and t; this last group is not written
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define rR(j)       mean2d(j,ig(6))
#define fU(j)       mean2d(j,ig(6)+1)
#define fV(j)       mean2d(j,ig(6)+2)
#define fW(j)       mean2d(j,ig(6)+3)
#define fU_y(j)     mean2d(j,ig(6)+4)
#define fV_y(j)     mean2d(j,ig(6)+5)
#define fW_y(j)     mean2d(j,ig(6)+6)
#define rV(j)       mean2d(j,ig(6)+7)
#define rV_y(j)     mean2d(j,ig(6)+8)
#define Rvv(j)      mean2d(j,ig(6)+9)
#define Rss_y(j)    mean2d(j,ig(6)+10)
#define Rsv_y(j)    mean2d(j,ig(6)+11)
#define Fy(j)       mean2d(j,ig(6)+12)
#define Fy_y(j)     mean2d(j,ig(6)+13)
#define Tau_yy(j)   mean2d(j,ig(6)+14)
#define Tau_yy_y(j) mean2d(j,ig(6)+15)
#define rP(j)       mean2d(j,ig(6)+16)
#define aux(j)      mean2d(j,ig(6)+17)
  sg(ng) = 18

  ! -----------------------------------------------------------------------
  nmax = ig(ng) +sg(ng) -1
  IF ( MAX_AVG_TEMPORAL .LT. nmax ) THEN
    CALL IO_WRITE_ASCII(efile,'AVERAGES_SCAL_XZ. Not enough space in local arrays.')
    CALL DNS_STOP(LES_ERROR_AVGTMP)
  END IF
  mean2d(:,1:nmax) = C_0_R

  ng   = ng -1
  nmax = ig(ng) +sg(ng) -1 ! the last group is not written out

  ! #######################################################################
  WRITE(line1,*) itime; line1 = 'Calculating scal statistics at It'//TRIM(ADJUSTL(line1))//'...'
  CALL IO_WRITE_ASCII(lfile,line1)

  ! #######################################################################
  ! Preliminary data
  ! #######################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
    CALL AVG_IK_V(imax,jmax,kmax, jmax, u, g(1)%jac,g(3)%jac, fU(1), wrk1d, area)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, v, g(1)%jac,g(3)%jac, fV(1), wrk1d, area)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, w, g(1)%jac,g(3)%jac, fW(1), wrk1d, area)
    rR(:) = C_1_R

  ELSE
    dsdx = rho *u
    dsdy = rho *v
    dsdz = rho *w
    CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdx, g(1)%jac,g(3)%jac, fU(1), wrk1d, area)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdy, g(1)%jac,g(3)%jac, fV(1), wrk1d, area)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdz, g(1)%jac,g(3)%jac, fW(1), wrk1d, area)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, rho,  g(1)%jac,g(3)%jac, rR(1), wrk1d, area)
    fU(:) = fU(:) /rR(:)
    fV(:) = fV(:) /rR(:)
    fW(:) = fW(:) /rR(:)

  END IF

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fU(1), fU_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fU(1), fV_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fU(1), fW_y(1), wrk3d, wrk2d,wrk3d)

  CALL AVG_IK_V(imax,jmax,kmax, jmax, v, g(1)%jac,g(3)%jac, aux(1), wrk1d, area)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), aux(1), rV_y(1), wrk3d, wrk2d,wrk3d)

  ! #######################################################################
  ! Moments (Reynolds and Favre averages)
  ! #######################################################################
  CALL AVG_IK_V(imax,jmax,kmax, jmax, s_local, g(1)%jac,g(3)%jac, rS(1), wrk1d, area)

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
    fS(:) = rS(:)
  ELSE
    wrk3d = rho *s_local
    CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, fS(1), wrk1d, area)
    fS(:) = fS(:) /rR(:)
  END IF

  ! -----------------------------------------------------------------------
  DO j = 1,jmax
    wrk3d(:,j,:) = s_local(:,j,:) - rS(j)
  END DO
  tmp1 = wrk3d *wrk3d
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS2(1), wrk1d, area)
  tmp1 = wrk3d *tmp1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS3(1), wrk1d, area)
  tmp1 = wrk3d *tmp1
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS4(1), wrk1d, area)

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
    fS2(:) = rS2(:)
    fS3(:) = rS3(:)
    fS4(:) = rS4(:)

  ELSE
    DO j = 1,jmax
      wrk3d(:,j,:) = s_local(:,j,:) - fS(j)
    END DO
    tmp1 = wrk3d *wrk3d *rho
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS2(1), wrk1d, area)
    tmp1 = wrk3d *tmp1
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS3(1), wrk1d, area)
    tmp1 = wrk3d *tmp1
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, fS4(1), wrk1d, area)
    fS2(:) = fS2(:) /rR(:)
    fS3(:) = fS3(:) /rR(:)
    fS4(:) = fS4(:) /rR(:)

  END IF

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fS2(1), Rss_y(1), wrk3d, wrk2d,wrk3d)

  DO j = 1,jmax
    wrk3d(:,j,:) = ( v(:,j,:) - fV(j) ) *( v(:,j,:) - fV(j) )
  END DO
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) wrk3d = wrk3d *rho
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Rvv(1), wrk1d, area)
  Rvv(:) = Rvv(:) /rR(:)

  ! #######################################################################
  ! Correlations
  ! #######################################################################
  DO j = 1,jmax
    wrk3d(:,j,:) = s_local(:,j,:) -rS(j)
  END DO
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) wrk3d = wrk3d *rho

  DO j = 1,jmax
    dsdx(:,j,:) = wrk3d(:,j,:) *(u(:,j,:)-fU(j))
    dsdy(:,j,:) = wrk3d(:,j,:) *(v(:,j,:)-fV(j))
    dsdz(:,j,:) = wrk3d(:,j,:) *(w(:,j,:)-fW(j))
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdx, g(1)%jac,g(3)%jac, Rsu(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdy, g(1)%jac,g(3)%jac, Rsv(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdz, g(1)%jac,g(3)%jac, Rsw(1), wrk1d, area)
  Rsu(:) = Rsu(:) /rR(:)
  Rsv(:) = Rsv(:) /rR(:)
  Rsw(:) = Rsw(:) /rR(:)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Rsv(1), Rsv_y(1), wrk3d, wrk2d,wrk3d)

  ! -----------------------------------------------------------------------
  ! turbulent transport terms
  DO j = 1,jmax
    tmp1(:,j,:) = dsdy(:,j,:) *(s_local(:,j,:) -fS(j))
    tmp2(:,j,:) = dsdy(:,j,:) *(v(:,j,:)       -fV(j))
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Tssy1(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, Tsvy1(1), wrk1d, area)

  ! -----------------------------------------------------------------------
  ! Pressure terms in transport equations
  CALL AVG_IK_V(imax,jmax,kmax, jmax, p, g(1)%jac,g(3)%jac, rP(1), wrk1d, area)

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s_local,dsdy, wrk3d, wrk2d,wrk3d)
  DO j = 1,jmax
    tmp1(:,j,:) = (p(:,j,:) -rP(j)) *(s_local(:,j,:) - fS(j)  )
    tmp2(:,j,:) = (p(:,j,:) -rP(j)) *(dsdy(:,j,:)    - fS_y(j))
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Tsvy3(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, PIsv(1), wrk1d, area)
  PIsv(:) = PIsv(:) /rR(:)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rP(1), aux(1), wrk3d, wrk2d,wrk3d)
  Gsv(:) = (rS(:) -fS(:)) *aux(:)

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
    END IF
  END IF

  IF ( transport%active(is) ) THEN ! Transport in tmp3
    CALL FI_TRANSPORT     (transport, i1, imax,jmax,kmax, is, s,tmp3, dsdy, wrk2d,wrk3d)
    CALL FI_TRANSPORT_FLUX(transport,     imax,jmax,kmax, is, s,dsdz)
  END IF

  IF ( is .GT. inb_scal ) THEN     ! Diagnostic variables; I overwrite tmp1 and dsdx and recalculate them.
    IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
      coefQ = C_1_R                        ! Coefficient in the evaporation term
      coefR = C_0_R                        ! Coefficient in the radiation term
      coefT = C_0_R                        ! Coefficient in the transport term
      IF ( is .EQ. inb_scal_array+1 ) THEN ! Default values are for liquid; defining them for buoyancy
        coefQ = buoyancy%parameters(inb_scal_array) /froude
        coefR = buoyancy%parameters(inb_scal) /froude
        DO i = 1,inb_scal
          coefT = coefT + transport%parameters(i) /settling *buoyancy%parameters(i) /froude
        END DO
      END IF

      CALL THERMO_AIRWATER_LINEAR_SOURCE(imax,jmax,kmax, s, dsdx,dsdy,dsdz) ! calculate xi in dsdx
      CALL FI_GRADIENT(imax,jmax,kmax, dsdx,tmp2, tmp1, wrk2d,wrk3d)

      dummy=-diff *coefQ
      tmp2 = dsdz *tmp2 *dummy         ! evaporation source

      IF ( transport%active(is) .OR. radiation%active(is) ) THEN ! preparing correction terms into dsdz
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), dsdx,tmp1, wrk3d, wrk2d,wrk3d)
        dsdz = dsdz *tmp1
      END IF

      IF ( radiation%active(is) ) THEN ! radiation source; needs dsdy
        CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,1,1,radiation%scalar(is)), tmp1, wrk1d,wrk3d)
        dummy= thermo_param(2) *coefQ
        tmp1 = tmp1 *( coefR + dsdy *dummy  )
        ! Correction term needs dsdz
        CALL OPR_RADIATION_FLUX(radiation, imax,jmax,kmax, g(2), s(1,1,1,radiation%scalar(is)), dsdx, wrk1d,wrk3d)
        dsdx = dsdx *dsdz *dummy
      ELSE
        tmp1 = C_0_R; dsdx = C_0_R
      END IF

      IF ( transport%active(is) ) THEN ! transport source; needs dsdy
        dummy= coefQ
        tmp3 = tmp3 *( coefT + dsdy *dummy )
        ! Correction term needs dsdz
        CALL FI_TRANSPORT_FLUX(transport, imax,jmax,kmax, is, s,dsdy)
        dsdz = dsdy *dsdz *dummy
      ELSE
        tmp3 = C_0_R; dsdz = C_0_R
      END IF

    ELSE
      IF ( buoyancy%TYPE .NE. EQNS_EXPLICIT ) THEN
        CALL FI_GRADIENT(imax,jmax,kmax, s,dsdx, dsdy, wrk2d,wrk3d)
        CALL FI_BUOYANCY_SOURCE(buoyancy, isize_field, s, dsdx, wrk3d) ! dsdx contains gradient
        tmp1 = wrk3d* diff /froude
      END IF

    END IF

  END IF

  ! -----------------------------------------------------------------------
  ! Calculating averages
  k = ig(1)+5
  IF ( radiation%active(is) ) THEN
    k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d, area)
    k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdx, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d, area) ! correction term or flux
  END IF
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR .OR. imixture .EQ. MIXT_TYPE_AIRWATER ) THEN           ! evaporation
    k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d, area)
  END IF
  IF ( transport%active(is) ) THEN
    k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d, area)
    k = k + 1; CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdz, g(1)%jac,g(3)%jac, mean2d(1,k), wrk1d, area) ! correction term or flux
  END IF

  wrk3d = tmp1 + tmp2 + tmp3 ! total
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, rQ(1), wrk1d, area)
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
    fQ(:) = rQ(:)

    DO j = 1,jmax
      tmp1(:,j,:) =(s_local(:,j,:) - rS(j)) *wrk3d(:,j,:)
      tmp2(:,j,:) =(v(:,j,:) - fV(j))       *wrk3d(:,j,:)
    END DO
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Qss(1), wrk1d, area)                ! transport equation
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, Qsv(1), wrk1d, area)                ! transport equation
    Qss(:) = Qss(:) *C_2_R

  ELSE
    DO j = 1,jmax
      wrk3d(:,j,:) = wrk3d(:,j,:) *rho(:,j,:)
    END DO
    CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, fQ(1), wrk1d, area)
    fQ(:) = fQ(:) /rR(:)

    DO j = 1,jmax
      tmp1(:,j,:) =(s_local(:,j,:) -fS(j)) *wrk3d(:,j,:)
      tmp2(:,j,:) =(v(:,j,:)       -fV(j)) *wrk3d(:,j,:)
    END DO
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Qss(1), wrk1d, area)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, Qsv(1), wrk1d, area)
    Qss(:) = Qss(:) *C_2_R /rR(:)
    Qsv(:) = Qsv(:) /rR(:)

  END IF

  ! #######################################################################
  ! Derivatives
  ! #######################################################################
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s_local, dsdx,    wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s_local, dsdy,    wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s_local, dsdz,    wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1,     bcs, g(2), rS(1),   rS_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1,     bcs, g(2), fS(1),   fS_y(1), wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)

  ! -----------------------------------------------------------------------
  wrk3d = tmp2 *C_2_R -tmp1 -tmp3
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Tau_yy(1), wrk1d, area)
  Tau_yy(:) = Tau_yy(:) *visc *c23
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tau_yy(1), Tau_yy_y(1), wrk3d, wrk2d,wrk3d)

  ! Transport term
  DO j = 1,jmax
    wrk3d(:,j,:) = (wrk3d(:,j,:) -Tau_yy(j)) *(s_local(:,j,:) -fS(j))
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Tsvy2(1), wrk1d, area)
  Tsvy2(:)  =-Tsvy2(:)  *visc *c23

  ! Dissipation terms; mean terms substracted below
  wrk3d = dsdy *( ( tmp2 *C_2_R -tmp1 -tmp3 ) *c23 *visc + tmp2 *diff )
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, Esv(1), wrk1d, area)

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, tmp1, wrk3d, wrk2d,wrk3d)
  tmp2 = dsdx *( ( tmp1 +tmp2 ) *visc + tmp2 *diff )
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) tmp2 = tmp2 *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, aux(1), wrk1d, area)
  Esv(:) = Esv(:) +aux(:)

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  tmp2 = dsdz *( ( tmp3 +tmp2 ) *visc + tmp2 *diff )
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) tmp2 = tmp2 *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, aux(1), wrk1d, area)
  Esv(:) = Esv(:) +aux(:)

  ! Dissipation terms; mean terms substracted below
  wrk3d = dsdx*dsdx +dsdy*dsdy +dsdz*dsdz
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Ess(1), wrk1d, area)
  Ess(:) = Ess(:) *diff *C_2_R

  ! -----------------------------------------------------------------------
  ! Moments
  DO j = 1,jmax
    wrk3d(:,j,:) = dsdy(:,j,:) - rS_y(j)
  END DO

  tmp1 = dsdx *dsdx
  tmp2 = wrk3d*wrk3d
  tmp3 = dsdz *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, S_x2(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, S_y2(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, S_z2(1), wrk1d, area)

  tmp1 = tmp1 *dsdx
  tmp2 = tmp2*wrk3d
  tmp3 = tmp3 *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, S_x3(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, S_y3(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, S_z3(1), wrk1d, area)

  tmp1 = tmp1 *dsdx
  tmp2 = tmp2*wrk3d
  tmp3 = tmp3 *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, S_x4(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, S_y4(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, S_z4(1), wrk1d, area)

  ! -----------------------------------------------------------------------
  ! Molecular fluxes
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) dsdy = dsdy *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dsdy, g(1)%jac,g(3)%jac, Fy(1), wrk1d, area)
  Fy(:) = Fy(:) *diff
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Fy(1),   Fy_y(1),   wrk3d, wrk2d,wrk3d)

  ! Contribution to turbulent transport
  DO j = 1,jmax
    tmp1(:,j,:) = dsdy(:,j,:) *(s_local(:,j,:) -fS(j))
    tmp2(:,j,:) = dsdy(:,j,:) *(v(:,j,:)       -fV(j))
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Tssy2(1), wrk1d, area)
  Tssy2(:) =-Tssy2(:) *diff *C_2_R
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, aux(1), wrk1d, area)
  Tsvy2(:) = Tsvy2(:) -aux(:) *diff

  ! Contribution to dissipation
  Ess(:) = ( Ess(:) -Fy(:)     *rS_y(:) -Fy(:) *rS_y(:) ) /rR(:)
  Esv(:) = ( Esv(:) -Tau_yy(:) *rS_y(:) -Fy(:) *rV_y(:) ) /rR(:)

  ! #######################################################################
  ! Source terms in transport equations
  ! #######################################################################
  IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
    CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, wrk3d)
  ELSE
    wrk1d(1:jmax,1) = C_0_R
    CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
  ENDIF
  dummy =  C_1_R /froude
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) wrk3d = wrk3d *rho
  DO j = 1,jmax
    wrk3d(:,j,:) = (s_local(:,j,:)-fS(j)) *wrk3d(:,j,:) *dummy
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Bsv(1), wrk1d, area)
  Bsv(:) = Bsv(:) /rR(:)

  ! #######################################################################
  ! Complete budget equations
  ! #######################################################################
  ! Transport terms
  aux(:)   = Tssy1(:) +Tssy2(:)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), aux(1), Tssy_y(1), wrk3d, wrk2d,wrk3d)
  aux(:)   = Tsvy1(:) +Tsvy2(:) +Tsvy3(:)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), aux(1), Tsvy_y(1), wrk3d, wrk2d,wrk3d)

  ! Convective terms
  Css(:)   =-fV(:) *Rss_y(:)
  Csv(:)   =-fV(:) *Rsv_y(:)

  ! Production terms
  Pss(:)   =-Rsv(:) *fS_y(:) *C_2_R
  Psv(:)   =-Rsv(:) *fV_y(:) -Rvv(:) *fS_y(:)

  ! Diffusion variable-density terms
  Dss(:)   = (rS(:)-fS(:)) *Fy_y(:) *C_2_R
  Dsv(:)   = (rS(:)-fS(:)) *Tau_yy_y(:) +(rV(:)-fV(:)) *Fy_y(:)

  ! Transient terms
  Rss_t(:) = Css(:) +Pss(:) -Ess(:)         +Qss(:) +( Dss(:)         -Tssy_y(:) ) /rR(:)
  Rsv_t(:) = Csv(:) +Psv(:) -Esv(:) +Bsv(:) +Qsv(:) +( Dsv(:) -Gsv(:) -Tsvy_y(:) ) /rR(:)

  ! ###################################################################
  ! Output
  ! #######################################################################
#ifdef USE_NETCDF
  line2 = ''
  DO k = 1,ng
    line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(varname(k)))
  END DO
  nvar = MAX_AVG_TEMPORAL
  CALL LIST_STRING( line2, nvar, varname2 )

  WRITE(line1,*) is; line1='avg'//TRIM(ADJUSTL(line1))//'s'
  WRITE(name,*) itime; name=TRIM(ADJUSTL(line1))//TRIM(ADJUSTL(name))
  CALL IO_WRITE_AVERAGES( name, itime,rtime, nvar,jmax, g(2)%nodes, varname2, mean2d )

#else
  ! -----------------------------------------------------------------------
  ! TkStat file
  ! -----------------------------------------------------------------------
  IF ( L_AVGMAX .LT. nmax ) THEN
    CALL IO_WRITE_ASCII(efile,'AVERAGES_SCAL_XZ. Not enough space in format definition.')
    CALL DNS_STOP(LES_ERROR_AVGTMP)
  END IF

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

    WRITE(line1,*) is; line1='avg'//TRIM(ADJUSTL(line1))//'s'
    WRITE(name,*) itime; name=TRIM(ADJUSTL(line1))//TRIM(ADJUSTL(name))

#define LOC_UNIT_ID 23
#define LOC_STATUS 'unknown'
#ifdef USE_RECLEN
    OPEN(UNIT=LOC_UNIT_ID, RECL=1050, FILE=name, STATUS=LOC_STATUS) ! this is probably outdated
#else
    OPEN(UNIT=LOC_UNIT_ID, FILE=name, STATUS=LOC_STATUS)
#endif

    WRITE(LOC_UNIT_ID, '(A8,E14.7E3)') 'RTIME = ', rtime
    line2 = 'I J Y'     ! Independent variables
    DO k = 1,ng         ! Dependent variables depending on y and t
      WRITE(LOC_UNIT_ID,1010) 'GROUP = '//TRIM(ADJUSTL(groupname(k)))//' '//TRIM(ADJUSTL(varname(k)))
      line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(varname(k)))
    END DO
    DO k = 1,11         ! Dependent variables depending on t, for consistency with old format
      line2 = TRIM(ADJUSTL(line2))//' '//'dummy'
    END DO
    WRITE(LOC_UNIT_ID,1010) 'GROUP = '//TRIM(ADJUSTL(line1))
    WRITE(LOC_UNIT_ID,1010) 'GROUP = '//TRIM(ADJUSTL(line1))
    WRITE(LOC_UNIT_ID,1010) TRIM(ADJUSTL(line2))

    DO j = 1,jmax
      WRITE(LOC_UNIT_ID,1020) 1, j, g(2)%nodes(j), (mean2d(j,k),k=1,nmax)
    END DO

    CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  END IF
#endif

1010 FORMAT(A)
1020 FORMAT(I5,(1X,I5),L_AVGMAX(1X,G_FORMAT_R),11(1X,G_FORMAT_R))

#endif

  RETURN
END SUBROUTINE AVG_SCAL_XZ
