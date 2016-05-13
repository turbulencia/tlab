#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define rR(j)     mean2d(j,1)
#define fU(j)     mean2d(j,2)
#define fV(j)     mean2d(j,3)
#define fW(j)     mean2d(j,4)
#define fS(j)     mean2d(j,5)
#define rS(j)     mean2d(j,6)
#define Rsv(j)    mean2d(j,7)
#define Rsu(j)    mean2d(j,8)
#define Rsw(j)    mean2d(j,9)
#define fS2(j)    mean2d(j,10)
#define fS3(j)    mean2d(j,11)
#define fS4(j)    mean2d(j,12)
#define fS5(j)    mean2d(j,13)
#define fS6(j)    mean2d(j,14)
#define rS2(j)    mean2d(j,15)
#define rS3(j)    mean2d(j,16)
#define rS4(j)    mean2d(j,17)
#define rS5(j)    mean2d(j,18)
#define rS6(j)    mean2d(j,19)
#define Tssy(j)   mean2d(j,20)
#define Ess(j)    mean2d(j,21)
#define Dss(j)    mean2d(j,22)
#define F1(j)     mean2d(j,23)
#define F2(j)     mean2d(j,24)
#define F3(j)     mean2d(j,25)
#define rSRHO(j)  mean2d(j,26)
#define fSRHO(j)  mean2d(j,27)
#define fS_y(j)   mean2d(j,28)
#define fU_y(j)   mean2d(j,29)
#define fV_y(j)   mean2d(j,30)
#define fW_y(j)   mean2d(j,31)
#define F2_y(j)   mean2d(j,32)
#define Tssy_y(j) mean2d(j,33)
#define rS_y(j)   mean2d(j,34)
#define fQ(j)     mean2d(j,35)
#define rQ(j)     mean2d(j,36)
#define Qss(j)    mean2d(j,37)

! we need space at the end for additional calculations
#define L_AVGMAX 42

!########################################################################
!# Tool/Library AVERAGES
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!# 2007/08/23 - J.P. Mellado
!#              Cleaned. Adding derivative statistics.
!# 2013/12/16 - J.P. Mellado
!#              Passing whole array s and computing sources inside.
!#              Cleaning the buoyancy case.
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE AVG_SCAL_TEMPORAL_LAYER(is, y,dx,dy,dz, q,s, s_local, dsdx,dsdy,dsdz, tmp1,tmp2,tmp3, mean2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  USE DNS_GLOBAL, ONLY : imode_eqns, imode_flow, idiffusion, iradiation, itransport, ibodyforce
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, inb_scal, inb_scal_array, imode_fdm, i1bc,j1bc,k1bc, area, scaley
  USE DNS_GLOBAL, ONLY : body_param, rad_param, trans_param,irad_scalar
  USE DNS_GLOBAL, ONLY : mean_rho, delta_rho, ycoor_rho, delta_u, ycoor_u, mean_i, delta_i, ycoor_i
  USE DNS_GLOBAL, ONLY : visc, schmidt, froude,settling
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
  TREAL, DIMENSION(*),                INTENT(IN)    :: y,dx,dy,dz
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)    :: q, s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(IN)    :: s_local
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: dsdx,dsdy,dsdz, tmp1,tmp2,tmp3, wrk3d
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: mean2d, wrk1d
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk2d

  TARGET q, s

! -----------------------------------------------------------------------
  TINTEGER i,j,k
  TREAL diff, AVG_IK, SIMPSON_NU, UPPER_THRESHOLD, LOWER_THRESHOLD
  TREAL delta_m, delta_w, delta_s, delta_s_area, delta_s_position, delta_s_value 
  TREAL smin_loc, smax_loc
  TREAL delta_hb01, delta_ht01, delta_h01
  TREAL delta_sb01, delta_st01, delta_s01
  TREAL Pss, Rss_t
  TREAL rR2, zprim, fz1, fz2, fz3
  TREAL mixing1, mixing2
  TREAL var_x, var_y, var_z, skew_x, skew_y, skew_z, flat_x, flat_y, flat_z
  TINTEGER jloc_max(1)

  TREAL dummy, dummy1, dummy2, dummy3

  TREAL VAUXPRE(5), VAUXPOS(11)
  TINTEGER ivauxpre, ivauxpos
  CHARACTER*32 fname
  CHARACTER*250 line1
  CHARACTER*850 line2

  TINTEGER iavgpos
  
! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u,v,w, rho, vis

! #######################################################################
  WRITE(line1,*) itime; line1 = 'Calculating scal statistics at It'//TRIM(ADJUSTL(line1))//'...'
  CALL IO_WRITE_ASCII(lfile,line1)

  IF ( MAX_AVG_TEMPORAL .LT. L_AVGMAX ) THEN
     CALL IO_WRITE_ASCII(efile,'AVG_SCAL_TEMPORAL_LAYER. Not enough space.')
     CALL DNS_STOP(LES_ERROR_AVGTMP)
  ENDIF
  mean2d(:,1:L_AVGMAX) = C_0_R

! Define pointers
  u => q(:,:,:,1)
  v => q(:,:,:,2)
  w => q(:,:,:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL   )THEN
     rho => q(:,:,:,5)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)
  ENDIF

! -----------------------------------------------------------------------
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF

  iavgpos = 0
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
     iavgpos = iavgpos + 3
  ENDIF
  IF ( itransport .NE. EQNS_NONE ) THEN
     iavgpos = iavgpos + 2
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        trans_param(inb_scal_array    ) = C_1_R ! liquid
        trans_param(inb_scal_array + 1) = C_1_R ! buoyancy
     ENDIF
  ENDIF

! -----------------------------------------------------------------------
! TkStat file; header
! -----------------------------------------------------------------------
#ifdef USE_MPI
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)

  IF ( ims_pro .EQ. 0 ) THEN
#endif
  WRITE(line1,*) is;    line1='avg'//TRIM(ADJUSTL(line1))//'s'
  WRITE(fname,*) itime; fname=TRIM(ADJUSTL(line1))//TRIM(ADJUSTL(fname))

#ifdef USE_RECLEN
     OPEN(UNIT=23, RECL=1050, FILE=fname, STATUS='unknown') ! this is probably outdated
#else
     OPEN(UNIT=23, FILE=fname, STATUS='unknown')
#endif

     WRITE(23, '(A8,E14.7E3)') 'RTIME = ', rtime

! Independent variables
     line2='I J Y SM SW SS SR'

! Dependent variables depending on y and t
     line1 = 'rS fS rS_y fS_y rQ fQ'
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN ! Source term partition
        line1 = TRIM(ADJUSTL(line1))//' rQrad rQradC rQeva'
     ENDIF
     IF ( itransport .NE. EQNS_NONE ) THEN
        line1 = TRIM(ADJUSTL(line1))//' rQtra rQtraC'
     ENDIF
     WRITE(i23,1010) 'GROUP = Mean '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Rsu Rsv Rsw fS2 fS3 fS4 fS5 fS6 rS2 rS3 rS4 rS5 rS6 fRS rRS'
     WRITE(i23,1010) 'GROUP = Fluctuations '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Rss_t Pss Ess Tssy Tssy_y Dss Qss'
     WRITE(i23,1010) 'GROUP = fS2Budget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'VarSx VarSy VarSz SkewSx SkewSy SkewSz FlatSx FlatSy FlatSz'
     WRITE(i23,1010) 'GROUP = DerivativeFluctuations '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

! dependent variables dependent on t only
     IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        line1 = 'Delta_m Delta_w Delta_s Delta_s_Area Delta_s_Position Delta_s_Value'
        WRITE(i23,1010) 'GROUP = ShearThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

        line1 = 'Delta_sb01 Delta_st01 Delta_s01 mixing_Youngs mixing_Cook'
        WRITE(i23,1010) 'GROUP = MixingThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     ENDIF

     WRITE(i23,1010) TRIM(ADJUSTL(line2))

1010 FORMAT(A)

#ifdef USE_MPI
  ENDIF
#endif

! #######################################################################
! Calculating source terms
! #######################################################################
  dsdx = C_0_R; dsdy = C_0_R; dsdz = C_0_R; tmp1 = C_0_R; tmp2 = C_0_R; tmp3 = C_0_R

  IF ( is .EQ. irad_scalar ) THEN      ! Radiation in tmp1
     CALL OPR_RADIATION     (iradiation, imax,jmax,kmax, dy, rad_param, s(1,1,1,inb_scal_array), tmp1, wrk1d,wrk3d)
     CALL OPR_RADIATION_FLUX(iradiation, imax,jmax,kmax, dy, rad_param, s(1,1,1,inb_scal_array), dsdx, wrk1d,wrk3d)
  ENDIF

  IF ( itransport .NE. EQNS_NONE) THEN ! Transport in tmp3
     CALL FI_TRANS_FLUX(itransport, i1, imax,jmax,kmax, is, inb_scal_array, trans_param, settling, dy, s,tmp3, dsdy, wrk1d,wrk2d,wrk3d)
  ENDIF

  IF ( is .GT. inb_scal ) THEN         ! Diagnostic variables
     IF      ( is .EQ. inb_scal_array   .AND. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN ! liquid
        dummy3 = C_1_R; dummy2 = C_0_R; dummy1 = C_0_R

     ELSE IF ( is .EQ. inb_scal_array+1 .AND. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN ! buoyancy
        dummy3 = body_param(inb_scal_array) /froude
        dummy2 = C_0_R; IF ( inb_scal .EQ. 2 ) dummy2 = dummy2 + body_param(2) /froude;
        dummy1 = trans_param(1) *body_param(1) /froude; IF ( inb_scal .EQ. 2 ) dummy1 = dummy1 + trans_param(2) *body_param(2) /froude
     ENDIF

     IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        CALL THERMO_AIRWATER_LINEAR_SOURCE(imax,jmax,kmax, s, dsdy,dsdz, dsdx) ! calculate xi in dsdx
        CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, dx,dy,dz, dsdx,tmp2, tmp1, wrk1d,wrk2d,wrk3d)
        
        dummy=-diff *dummy3
        tmp2 = dsdz *tmp2 *dummy ! evaporation source in tmp2
        
        IF ( itransport .NE. EQNS_NONE ) THEN ! transport source; needs dsdy
           tmp3 = tmp3 *( dummy1 + dsdy *dummy3 )
        ENDIF
        
        IF ( iradiation .NE. EQNS_NONE .AND. inb_scal .EQ. 2 ) THEN  ! radiation source; needs dsdy from before
           CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(1,1,1,inb_scal_array), tmp1, wrk1d,wrk3d)
           dummy= thermo_param(2) *dummy3
           tmp1 = tmp1 *( dummy2 + dsdy *dummy  )
           
! Correction term
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, dsdx,dsdy, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL OPR_RADIATION_FLUX(iradiation, imax,jmax,kmax, dy, rad_param, s(1,1,1,inb_scal_array), dsdx, wrk1d,wrk3d)
           dsdx = dsdx *dsdy *dsdz *dummy
           
        ELSE
           tmp1 = C_0_R; dsdx = C_0_R

        ENDIF

     ELSE
        CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, dx,dy,dz, s,dsdx, dsdy, wrk1d,wrk2d,wrk3d)
        CALL FI_BUOYANCY_SOURCE(ibodyforce, isize_field, body_param, s, dsdx, wrk3d) ! dsdx contains gradient
        tmp1 = wrk3d* diff /froude
        
     ENDIF

  ENDIF

! #######################################################################
! Main loop
! #######################################################################
  DO j = 1,jmax

! -----------------------------------------------------------------------
! Mean
! -----------------------------------------------------------------------
! Reynolds
     IF     ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
        rR(j) = C_1_R
     ELSEIF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        rR(j) = C_1_R ! not yet developed
     ELSE
        rR(j) = AVG_IK(imax,jmax,kmax, j, rho, dx,dz, area)
     ENDIF
     rS(j) = AVG_IK(imax,jmax,kmax, j, s_local,  dx,dz, area)
! Favre
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
        fU(j) = AVG_IK(imax,jmax,kmax, j, u, dx,dz, area)
        fV(j) = AVG_IK(imax,jmax,kmax, j, v, dx,dz, area)
        fW(j) = AVG_IK(imax,jmax,kmax, j, w, dx,dz, area)
        fS(j) = rS(j)

     ELSE
        DO k=1,kmax; DO i=1,imax
           wrk3d(i,1,k) = rho(i,j,k)*u(i,j,k)
           wrk3d(i,2,k) = rho(i,j,k)*v(i,j,k)
           wrk3d(i,3,k) = rho(i,j,k)*w(i,j,k)
           wrk3d(i,4,k) = rho(i,j,k)*s_local(i,j,k)
        ENDDO; ENDDO
        fU(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)/rR(j)
        fV(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)/rR(j)
        fW(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)/rR(j)
        fS(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)/rR(j)

     ENDIF

! -----------------------------------------------------------------------
! Fluctuations
! -----------------------------------------------------------------------
! Reynolds
     DO k = 1,kmax; DO i = 1,imax
        zprim = s_local(i,j,k)-rS(j)
        wrk3d(i,1,k) = zprim**2
        wrk3d(i,2,k) = wrk3d(i,1,k)*zprim
        wrk3d(i,3,k) = wrk3d(i,2,k)*zprim
        wrk3d(i,4,k) = wrk3d(i,3,k)*zprim
        wrk3d(i,5,k) = wrk3d(i,4,k)*zprim
     ENDDO; ENDDO
     rS2(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     rS3(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     rS4(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     rS5(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     rS6(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)

! Favre
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
        fS2(j) = rS2(j)
        fS3(j) = rS3(j)
        fS4(j) = rS4(j)
        fS5(j) = rS5(j)
        fS6(j) = rS6(j)

     ELSE
        DO k = 1,kmax; DO i = 1,imax
           zprim = s_local(i,j,k)-fS(j)
           wrk3d(i,1,k) =   rho(i,j,k)*zprim**2
           wrk3d(i,2,k) = wrk3d(i,1,k)*zprim
           wrk3d(i,3,k) = wrk3d(i,2,k)*zprim
           wrk3d(i,4,k) = wrk3d(i,3,k)*zprim
           wrk3d(i,5,k) = wrk3d(i,4,k)*zprim
        ENDDO; ENDDO
        fS2(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)/rR(j)
        fS3(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)/rR(j)
        fS4(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)/rR(j)
        fS5(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)/rR(j)
        fS6(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)/rR(j)

     ENDIF

! -----------------------------------------------------------------------
! Turbulent fluxes
! -----------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
        DO k = 1,kmax; DO i = 1,imax
           wrk3d(i,3,k) = (w(i,j,k)-fW(j))*(s_local(i,j,k)-fS(j))
           wrk3d(i,4,k) = (v(i,j,k)-fV(j))*(s_local(i,j,k)-fS(j))
           wrk3d(i,5,k) = (u(i,j,k)-fU(j))*(s_local(i,j,k)-fS(j))
        ENDDO; ENDDO
        Rsw(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
        Rsv(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
        Rsu(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)

     ELSE
        DO k = 1,kmax; DO i = 1,imax
           wrk3d(i,3,k) = rho(i,j,k)*(w(i,j,k)-fW(j))*(s_local(i,j,k)-fS(j))
           wrk3d(i,4,k) = rho(i,j,k)*(v(i,j,k)-fV(j))*(s_local(i,j,k)-fS(j))
           wrk3d(i,5,k) = rho(i,j,k)*(u(i,j,k)-fU(j))*(s_local(i,j,k)-fS(j))
        ENDDO; ENDDO
        Rsw(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)/rR(j)
        Rsv(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)/rR(j)
        Rsu(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)/rR(j)

     ENDIF

! -----------------------------------------------------------------------
! Production terms
! -----------------------------------------------------------------------
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = tmp1(i,j,k) + tmp2(i,j,k) + tmp3(i,j,k)
        wrk3d(i,2,k) = wrk3d(i,1,k)*(s_local(i,j,k)-fS(j)) 
     ENDDO; ENDDO
     rQ(j)  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     Qss(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)*C_2_R ! transport equation

     fQ(j)  = rQ(j)

     k = L_AVGMAX-iavgpos
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        k = k + 1; mean2d(j,k) = AVG_IK(imax,jmax,kmax, j, tmp1, dx,dz, area) ! radiation
        k = k + 1; mean2d(j,k) = AVG_IK(imax,jmax,kmax, j, dsdx, dx,dz, area) ! correction term
        k = k + 1; mean2d(j,k) = AVG_IK(imax,jmax,kmax, j, tmp2, dx,dz, area) ! evaporation
     ENDIF
     IF ( itransport .NE. EQNS_NONE ) THEN
        k = k + 1; mean2d(j,k) = AVG_IK(imax,jmax,kmax, j, tmp3, dx,dz, area)
        k = k + 1; mean2d(j,k) = C_0_R ! Correction term; to be developed.
     ENDIF
     
! -----------------------------------------------------------------------
! density-scalar correlations
! -----------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
        rSRHO(j) = C_0_R
        fSRHO(j) = C_0_R

     ELSE
        DO k = 1,kmax; DO i = 1,imax
           wrk3d(i,1,k) =   (s_local(i,j,k)-rS(j))*(rho(i,j,k)-rR(j))
           wrk3d(i,2,k) =   (s_local(i,j,k)-fS(j))*(rho(i,j,k)-rR(j))
           wrk3d(i,3,k) = (rho(i,j,k)-rR(j))*(rho(i,j,k)-rR(j))
        ENDDO; ENDDO
        rSRHO(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
        fSRHO(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
        rR2      = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)

        IF ( rR2 .GT. C_0_R ) THEN
           IF ( rS2(j) .GT. C_0_R ) THEN; rSRHO(j) = rSRHO(j)/sqrt(rR2*rS2(j))
           ELSE;                          rSRHO(j) = C_2_R; ENDIF
           
           IF ( fS2(j) .GT. C_0_R ) THEN; fSRHO(j) = fSRHO(j)/sqrt(rR2*fS2(j))
           ELSE;                          fSRHO(j) = C_2_R; ENDIF

        ELSE
           rSRHO(j) = C_BIG_R
           fSRHO(j) = C_BIG_R

        ENDIF

     ENDIF

  ENDDO

! #######################################################################
! Derivative terms
! #######################################################################
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, s_local, dsdx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, s_local, dsdy, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, s_local, dsdz, i0,i0, wrk1d,wrk2d,wrk3d)

  DO j = 1,jmax

     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        DO k = 1,kmax; DO i = 1,imax
           wrk3d(i,1,k) = vis(i,j,k)*dsdx(i,j,k)
           wrk3d(i,2,k) = vis(i,j,k)*dsdy(i,j,k)
           wrk3d(i,3,k) = vis(i,j,k)*dsdz(i,j,k)
        ENDDO; ENDDO
        F1(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)*diff
        F2(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)*diff
        F3(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)*diff

     ELSE
        F1(j) = AVG_IK(imax,jmax,kmax, j, dsdx, dx,dz, area)*diff
        F2(j) = AVG_IK(imax,jmax,kmax, j, dsdy, dx,dz, area)*diff
        F3(j) = AVG_IK(imax,jmax,kmax, j, dsdz, dx,dz, area)*diff

     ENDIF

! transport term
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
        DO k = 1,kmax; DO i = 1,imax
           wrk3d(i,1,k) = (s_local(i,j,k)-fS(j))*&
                ( (v(i,j,k)-fV(j))*(s_local(i,j,k)-fS(j)) - C_2_R*(diff*dsdy(i,j,k)-F2(j)) )
        ENDDO; ENDDO

     ELSE
        IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
           DO k = 1,kmax; DO i = 1,imax
              wrk3d(i,1,k) = (s_local(i,j,k)-fS(j))*&
                   ( rho(i,j,k)*(v(i,j,k)-fV(j))*(s_local(i,j,k)-fS(j))&
                   - C_2_R*(diff*vis(i,j,k)*dsdy(i,j,k)-F2(j)) )
           ENDDO; ENDDO
        ELSE
           DO k = 1,kmax; DO i = 1,imax
              wrk3d(i,1,k) = (s_local(i,j,k)-fS(j))*&
                   ( rho(i,j,k)*(v(i,j,k)-fV(j))*(s_local(i,j,k)-fS(j))&
                   - C_2_R*(diff*dsdy(i,j,k)-F2(j)) )
           ENDDO; ENDDO
        ENDIF
        
     ENDIF
     Tssy(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     
  ENDDO

! Oy Mean Derivatives
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rS(1),   rS_y(1),   i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, fS(1),   fS_y(1),   i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Tssy(1), Tssy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, F2(1),   F2_y(1),   i0,i0, wrk1d,wrk2d,wrk3d)

! #######################################################################
! Global quantites
! #######################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN

! -----------------------------------------------------------------------
! Based on delta_u 
! -----------------------------------------------------------------------
     IF ( ABS(delta_u) .GT. C_SMALL_R ) THEN
        CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, fU(1), fU_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
        delta_w = ABS(delta_u)/MAXVAL(ABS(fU_y(1:jmax)))

        DO j=1, jmax
           wrk1d(j,1) = rR(j)*(C_025_R - (fU(j)/delta_u)**2)
        ENDDO
        delta_m = SIMPSON_NU(jmax, wrk1d, y)/mean_rho

        DO j=1, jmax
           wrk1d(j,1) = ( F2(j) -  rR(j)*Rsv(j) )*fS_y(j)
        ENDDO
        delta_s_area = SIMPSON_NU(jmax, wrk1d, y)*C_2_R/(mean_rho*delta_u)

     ELSE
        delta_m      = C_1_R
        delta_w      = C_1_R
        delta_s_area = C_1_R

     ENDIF

! -----------------------------------------------------------------------
! Based on delta_rho 
! -----------------------------------------------------------------------
     IF ( ABS(delta_rho) .GT. C_SMALL_R ) THEN
        dummy = mean_rho + (C_05_R-C_1EM2_R)*delta_rho
        delta_hb01 = LOWER_THRESHOLD(jmax, dummy, rR(1), y)
        dummy = mean_rho - (C_05_R-C_1EM2_R)*delta_rho
        delta_ht01 = UPPER_THRESHOLD(jmax, dummy, rR(1), y)
        delta_h01 = delta_ht01 - delta_hb01

     ELSE
        delta_h01 = C_1_R

     ENDIF

! -----------------------------------------------------------------------
! Based on scalar
! -----------------------------------------------------------------------
     IF ( ABS(delta_i(is)) .GT. C_SMALL_R ) THEN
        jloc_max = MAXLOC(ABS(fS_y(1:jmax))); j = jloc_max(1)

        wrk1d(1:jmax,1) = fS(jmax) - fS(1:jmax)
        delta_s_area     = SIMPSON_NU(jmax-j+1, wrk1d(j,1), y(j)) / ABS(delta_i(is))
 
        delta_s          = ABS(delta_i(is))/ABS(fS_y(j))
        delta_s_position = y(j)
        delta_s_value    = fS(j)

! Thickness
        dummy = mean_i(is) + (C_05_R-C_1EM3_R)*delta_i(is)
        delta_sb01 = LOWER_THRESHOLD(jmax, dummy, rS(1), y)
        dummy = mean_i(is) - (C_05_R-C_1EM3_R)*delta_i(is)
        delta_st01 = UPPER_THRESHOLD(jmax, dummy, rS(1), y)
        
        delta_sb01 = (y(1) + ycoor_i(is)*scaley) - delta_sb01
        delta_st01 = delta_st01 - (y(1) + ycoor_i(is)*scaley)
        delta_s01 = delta_st01 + delta_sb01
        
!  Mixing, Youngs' definition
        smin_loc = mean_i(is) - C_05_R*ABS(delta_i(is))
        smax_loc = mean_i(is) + C_05_R*ABS(delta_i(is))
        DO k = 1,kmax
           DO i = 1,imax*jmax
              wrk3d(i,1,k) = (s_local(i,1,k)-smin_loc)*(smax_loc-s_local(i,1,k))
           ENDDO
        ENDDO
        DO j = 1,jmax
           wrk1d(j,1) = AVG_IK(imax, jmax, kmax, j, wrk3d, dx, dz, area)
        ENDDO
        mixing1 = SIMPSON_NU(jmax, wrk1d, y)
        DO j = 1,jmax
           wrk1d(j,1)=(rS(j)-smin_loc)*(smax_loc-rS(j))
        ENDDO
        mixing1 = mixing1/SIMPSON_NU(jmax, wrk1d, y)

! Mixing, Cook's definition
        smin_loc = mean_i(is) - C_05_R*ABS(delta_i(is))
        smax_loc = mean_i(is) + C_05_R*ABS(delta_i(is))
        DO k = 1,kmax
           DO i = 1,imax*jmax
              wrk3d(i,1,k) = MIN(s_local(i,1,k)-smin_loc,smax_loc-s_local(i,1,k))
           ENDDO
        ENDDO
        DO j = 1,jmax
           wrk1d(j,1)=AVG_IK(imax, jmax, kmax, j, wrk3d, dx, dz, area)
        ENDDO
        mixing2 = SIMPSON_NU(jmax, wrk1d, y)
        DO j = 1,jmax
           wrk1d(j,1) = MIN(rS(j)-smin_loc,smax_loc-rS(j))
        ENDDO
        mixing2 = mixing2/SIMPSON_NU(jmax, wrk1d, y)

     ELSE
        delta_sb01 = C_1_R
        delta_st01 = C_1_R
        delta_s01  = C_1_R
        delta_s    = C_1_R
        mixing1    = C_1_R
        mixing2    = C_1_R

     ENDIF

! -----------------------------------------------------------------------
! Jet
! -----------------------------------------------------------------------
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
! Not developed yet

  ENDIF

! #######################################################################
! Variance transport equation
! #######################################################################
  DO j = 1,jmax
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        DO k = 1,kmax; DO i = 1,imax
           dummy = diff*vis(i,j,k) 
           fz1 = dummy*dsdx(i,j,k)-F1(j)
           fz2 = dummy*dsdy(i,j,k)-F2(j)
           fz3 = dummy*dsdz(i,j,k)-F3(j)
           wrk3d(i,1,k) = fz1*dsdx(i,j,k) + fz2*(dsdy(i,j,k)-fS_y(j)) + fz3*dsdz(i,j,k)
           wrk3d(i,2,k) = s_local(i,j,k)-fS(j)
        ENDDO; ENDDO
     ELSE
        DO k = 1,kmax; DO i = 1,imax
           fz1 = diff*dsdx(i,j,k)-F1(j)
           fz2 = diff*dsdy(i,j,k)-F2(j)
           fz3 = diff*dsdz(i,j,k)-F3(j)
           wrk3d(i,1,k) = fz1*dsdx(i,j,k) + fz2*(dsdy(i,j,k)-fS_y(j)) + fz3*dsdz(i,j,k)
           wrk3d(i,2,k) = s_local(i,j,k)-fS(j)
        ENDDO; ENDDO
     ENDIF
     Ess(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)*C_2_R/rR(j)
     Dss(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)*C_2_R*F2_y(j)

  ENDDO

! #######################################################################
! Final loop
! #######################################################################
  DO j = 1,jmax            

! -----------------------------------------------------------------------
! Derivatives Fluctuations
! -----------------------------------------------------------------------
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) =  dsdx(i,j,k)**2
        wrk3d(i,2,k) = (dsdy(i,j,k)-rS_y(j))**2
        wrk3d(i,3,k) =  dsdz(i,j,k)**2
        wrk3d(i,4,k) = wrk3d(i,1,k)* dsdx(i,j,k)
        wrk3d(i,5,k) = wrk3d(i,2,k)*(dsdy(i,j,k)-rS_y(j))
        wrk3d(i,6,k) = wrk3d(i,3,k)* dsdz(i,j,k)
        wrk3d(i,7,k) = wrk3d(i,4,k)* dsdx(i,j,k)
        wrk3d(i,8,k) = wrk3d(i,5,k)*(dsdy(i,j,k)-rS_y(j))
        wrk3d(i,9,k) = wrk3d(i,6,k)* dsdz(i,j,k)
     ENDDO; ENDDO
     var_x  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     var_y  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     var_z  = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     skew_x = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     skew_y = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     skew_z = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area)
     flat_x = AVG_IK(imax,jmax,kmax, i7, wrk3d, dx,dz, area)
     flat_y = AVG_IK(imax,jmax,kmax, i8, wrk3d, dx,dz, area)
     flat_z = AVG_IK(imax,jmax,kmax, i9, wrk3d, dx,dz, area)
     IF ( var_x .GT. C_SMALL_R ) THEN
        skew_x = skew_x / var_x**C_1_5_R
        flat_x = flat_x / var_x**C_2_R
     ELSE
        skew_x = C_BIG_R
        flat_x = C_BIG_R
     ENDIF
     IF ( var_y .GT. C_SMALL_R ) THEN
        skew_y = skew_y / var_y**C_1_5_R
        flat_y = flat_y / var_y**C_2_R
     ELSE
        skew_y = C_BIG_R
        flat_y = C_BIG_R
     ENDIF
     IF ( var_z .GT. C_SMALL_R ) THEN
        skew_z = skew_z / var_z**C_1_5_R
        flat_z = flat_z / var_z**C_2_R
     ELSE
        skew_z = C_BIG_R
        flat_z = C_BIG_R
     ENDIF

! -----------------------------------------------------------------------
! Final calculations
! -----------------------------------------------------------------------
     Pss   =-C_2_R*Rsv(j)*fS_y(j)
     Rss_t = Pss - Ess(j) + (-Tssy_y(j) + Dss(j) + Qss(j))/rR(j)

! #######################################################################
! Output
! #######################################################################
#ifdef USE_MPI
     IF ( ims_pro .EQ. 0 ) THEN
#endif

        IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
           ivauxpre = 5
           VAUXPRE(1) =  y(j)
           VAUXPRE(2) = (y(j)-scaley*ycoor_u     - y(1))/delta_m
           VAUXPRE(3) = (y(j)-scaley*ycoor_u     - y(1))/delta_w
           VAUXPRE(4) = (y(j)-scaley*ycoor_i(is) - y(1))/delta_s01
           VAUXPRE(5) = (y(j)-scaley*ycoor_rho   - y(1))/delta_h01

        ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
! Not developed yet; for TkStat compatibility with previous files
           ivauxpre = 5
           VAUXPRE(1) = y(j)
           VAUXPRE(2) = y(j)
           VAUXPRE(3) = y(j)
           VAUXPRE(4) = y(j)
           VAUXPRE(5) = y(j)

        ENDIF

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
           ELSE
              ivauxpos = 0
           ENDIF

        ELSE
           ivauxpos = 0

        ENDIF

        WRITE(23,1020) 1, j, (VAUXPRE(k),k=1,ivauxpre),&
             rS(j), fS(j), rS_y(j), fS_y(j), rQ(j), fQ(j), (mean2d(j,L_AVGMAX-iavgpos+k),k=1,iavgpos), &
             Rsu(j), Rsv(j), Rsw(j), &
             fS2(j), fS3(j), fS4(j), fS5(j), fS6(j),&
             rS2(j), rS3(j), rS4(j), rS5(j), rS6(j),&
             fSRHO(j), rSRHO(j), &
             Rss_t, Pss, Ess(j), Tssy(j), Tssy_y(j), Dss(j), Qss(j), &
             var_x, var_y, var_z, skew_x, skew_y, skew_z, flat_x, flat_y, flat_z,&
             (VAUXPOS(k),k=1,ivauxpos)

! using the maximum 37+5+11
1020    FORMAT(I5,(1X,I5),5(1X,G_FORMAT_R),L_AVGMAX(1X,G_FORMAT_R),11(1X,G_FORMAT_R))
        
#ifdef USE_MPI
     ENDIF
#endif

  ENDDO

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CLOSE(23)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE AVG_SCAL_TEMPORAL_LAYER

