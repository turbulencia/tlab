#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE AVG_SCAL_XZ(is, q,s, s_local, dsdx,dsdy,dsdz, tmp1,tmp2,tmp3, mean2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL, ONLY : g, area
  USE DNS_GLOBAL, ONLY : imode_eqns, idiffusion, itransport
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, inb_scal, inb_scal_array
  USE DNS_GLOBAL, ONLY : buoyancy, radiation, transport
  USE DNS_GLOBAL, ONLY : rbackground, ribackground
  USE DNS_GLOBAL, ONLY : visc, schmidt, froude, settling
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

  TARGET q, s

  ! -----------------------------------------------------------------------
  TINTEGER, PARAMETER :: MAX_VARS_GROUPS = 10
  TINTEGER i,j,k, bcs(2,2)
  TREAL diff, dummy, coefT, coefR, coefQ

  TINTEGER ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nmax

  CHARACTER*32 name, groupname(MAX_VARS_GROUPS)
  CHARACTER*250 line1, varname(MAX_VARS_GROUPS)
  CHARACTER*850 line2

  TINTEGER nvar
  CHARACTER(LEN=32) varname2(MAX_AVG_TEMPORAL)

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
  END IF

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc /schmidt(is)
  END IF

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
    varname(2) = TRIM(ADJUSTL(varname(2)))//' rQrad rQradC'
    sg(ng) = sg(ng) + 2
  END IF
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR .OR. imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
    varname(2) = TRIM(ADJUSTL(varname(2)))//' rQeva'
    sg(ng) = sg(ng) + 1
  END IF
  IF ( transport%active(is) ) THEN
    varname(2) = TRIM(ADJUSTL(varname(2)))//' rQtra rQtraC'
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
#define fSR(j)    mean2d(j,ig(2)+9)
#define rSR(j)    mean2d(j,ig(2)+10)
  sg(ng) = 11

  groupname(ng) = 'Fluctuations'
  varname(ng)   = 'Rsu Rsv Rsw fS2 fS3 fS4 rS2 rS3 rS4 fRS rRS'

  ! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rss_t(j)  mean2d(j,ig(3)  )
#define Pss(j)    mean2d(j,ig(3)+1)
#define Ess(j)    mean2d(j,ig(3)+2)
#define Tssy(j)   mean2d(j,ig(3)+3)
#define Tssy_y(j) mean2d(j,ig(3)+4)
#define Dss(j)    mean2d(j,ig(3)+5)
#define Qss(j)    mean2d(j,ig(3)+6)
  sg(ng) = 7

  groupname(ng) = 'fS2Budget'
  varname(ng)   = 'Rss_t Pss Ess Tssy Tssy_y Dss Qss'

  ! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define S_x2(j)  mean2d(j,ig(4)  )
#define S_y2(j)  mean2d(j,ig(4)+1)
#define S_z2(j)  mean2d(j,ig(4)+2)
#define S_x3(j)  mean2d(j,ig(4)+3)
#define S_y3(j)  mean2d(j,ig(4)+4)
#define S_z3(j)  mean2d(j,ig(4)+5)
#define S_x4(j)  mean2d(j,ig(4)+6)
#define S_y4(j)  mean2d(j,ig(4)+7)
#define S_z4(j)  mean2d(j,ig(4)+8)
  sg(ng) = 9

  groupname(ng) = 'DerivativeFluctuations'
  varname(ng)   = 'Sx2 Sy2 Sz2 Sx3 Sy3 Sz3 Sx4 Sy4 Sz4'

#define L_AVGMAX 38

  ! -----------------------------------------------------------------------
  ! Auxiliary variables depending on y and t; this last group is not written
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define rR(j)     mean2d(j,ig(5))
#define fU(j)     mean2d(j,ig(5)+1)
#define fV(j)     mean2d(j,ig(5)+2)
#define fW(j)     mean2d(j,ig(5)+3)
#define fU_y(j)   mean2d(j,ig(5)+4)
#define fV_y(j)   mean2d(j,ig(5)+5)
#define fW_y(j)   mean2d(j,ig(5)+6)
#define F1(j)     mean2d(j,ig(5)+7)
#define F2(j)     mean2d(j,ig(5)+8)
#define F3(j)     mean2d(j,ig(5)+9)
#define F2_y(j)   mean2d(j,ig(5)+10)
  sg(ng) = 11

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

  ! #######################################################################
  ! Correlations
  ! #######################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
    DO j = 1,jmax
      wrk3d(:,j,:) = s_local(:,j,:) - rS(j)
    END DO
  ELSE
    DO j = 1,jmax
      wrk3d(:,j,:) =(s_local(:,j,:) - fS(j)) *rho(:,j,:)
    END DO
  END IF

  DO j = 1,jmax
    tmp1(:,j,:) = wrk3d(:,j,:) *(u(:,j,:)-fU(j))
    tmp2(:,j,:) = wrk3d(:,j,:) *(v(:,j,:)-fV(j))
    tmp3(:,j,:) = wrk3d(:,j,:) *(w(:,j,:)-fW(j))
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, Rsu(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, Rsv(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, Rsw(1), wrk1d, area)

  Rsu(:) = Rsu(:) /rR(:)
  Rsv(:) = Rsv(:) /rR(:)
  Rsw(:) = Rsw(:) /rR(:)

  ! -----------------------------------------------------------------------
  IF ( .NOT. (imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC) )THEN
    DO j = 1,jmax
      tmp1(:,j,:) = (rho(:,j,:) - rR(j)) *(s_local(:,j,:) - rS(j))
      tmp2(:,j,:) = (rho(:,j,:) - rR(j)) *(s_local(:,j,:) - fS(j))
    END DO
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, rS(1), wrk1d, area)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, fS(1), wrk1d, area)
  END IF

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
      wrk3d(:,j,:) =(s_local(:,j,:) - rS(j)) *wrk3d(:,j,:)
    END DO
    CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Qss(1), wrk1d, area)                ! transport equation
    Qss(:) = Qss(:) *C_2_R

  ELSE
    DO j = 1,jmax
      wrk3d(:,j,:) = wrk3d(:,j,:) *rho(:,j,:)
    END DO
    CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, fQ(1), wrk1d, area)
    fQ(:) = fQ(:) /rR(:)

    DO j = 1,jmax
      wrk3d(:,j,:) =(s_local(:,j,:) - fS(j)) *wrk3d(:,j,:)
    END DO
    CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Qss(1), wrk1d, area)
    Qss(:) = Qss(:) *C_2_R /rR(:)

  END IF

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
  END DO
  tmp3 = dsdz *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, S_x2(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, S_y2(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, S_z2(1), wrk1d, area)

  tmp1 = tmp1 *dsdx
  DO j = 1,jmax
    tmp2(:,j,:) =  tmp2(:,j,:) *(dsdy(:,j,:) - rS_y(j))
  END DO
  tmp3 = tmp3 *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, S_x3(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, S_y3(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, S_z3(1), wrk1d, area)

  tmp1 = tmp1 *dsdx
  DO j = 1,jmax
    tmp2(:,j,:) =  tmp2(:,j,:) *(dsdy(:,j,:) - rS_y(j))
  END DO
  tmp3 = tmp3 *dsdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, S_x4(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, S_y4(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, S_z4(1), wrk1d, area)

  ! -----------------------------------------------------------------------
  ! Variance transport equation; dissipation
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
    tmp1 = vis *dsdx
    tmp2 = vis *dsdy
    tmp3 = vis *dsdz
  ELSE
    tmp1 = dsdx
    tmp2 = dsdy
    tmp3 = dsdz
  END IF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1, g(1)%jac,g(3)%jac, F1(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp2, g(1)%jac,g(3)%jac, F2(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp3, g(1)%jac,g(3)%jac, F3(1), wrk1d, area)

  DO j = 1,jmax
    wrk3d(:,j,:) = (tmp1(:,j,:) - F1(j)) * dsdx(:,j,:)          &
                 + (tmp2(:,j,:) - F2(j)) *(dsdy(:,j,:)-fS_y(j)) &
                 + (tmp3(:,j,:) - F3(j)) * dsdz(:,j,:)
  END DO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Ess(1), wrk1d, area)
  Ess(:) = Ess(:) *C_2_R *diff /rR(:)

  ! -----------------------------------------------------------------------
  ! Variance transport equation; turbulent transport
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
    DO j = 1,jmax
      wrk3d(:,j,:) = s_local(:,j,:)-fS(j)
      tmp1(:,j,:) = wrk3d(:,j,:)* (            (v(:,j,:)-fV(j))*wrk3d(:,j,:) - C_2_R *diff *(tmp2(:,j,:)-F2(j)) )
    END DO
  ELSE
    DO j = 1,jmax
      wrk3d(:,j,:) = s_local(:,j,:)-fS(j)
      tmp1(:,j,:) = wrk3d(:,j,:)* ( rho(:,j,:)*(v(:,j,:)-fV(j))*wrk3d(:,j,:) - C_2_R *diff *(tmp2(:,j,:)-F2(j)) )
    END DO
  END IF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, tmp1,  g(1)%jac,g(3)%jac, Tssy(1), wrk1d, area)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tssy(1), Tssy_y(1), wrk3d, wrk2d,wrk3d)

  ! -----------------------------------------------------------------------
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), F2(1),   F2_y(1),   wrk3d, wrk2d,wrk3d)
  Dss(:)   = (rS(:)-fS(:)) *F2_y(:) *C_2_R *diff

  Pss(:)   =-C_2_R *Rsv(:) *fS_y(:)

  Rss_t(:) = Pss(:) - Ess(:) + (Dss(:) + Qss(:) -Tssy_y(:)) /rR(:)

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
#endif

  ! -----------------------------------------------------------------------
  ! TkStat file; header
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

  RETURN

1010 FORMAT(A)
1020 FORMAT(I5,(1X,I5),L_AVGMAX(1X,G_FORMAT_R),11(1X,G_FORMAT_R))

END SUBROUTINE AVG_SCAL_XZ
