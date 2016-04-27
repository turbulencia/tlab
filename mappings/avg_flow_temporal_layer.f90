#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define rU(j)     mean2d(j,1)
#define rV(j)     mean2d(j,2)
#define rW(j)     mean2d(j,3)
#define rP(j)     mean2d(j,4)
#define rR(j)     mean2d(j,5)
#define rT(j)     mean2d(j,6)

#define fU(j)     mean2d(j,7)
#define fV(j)     mean2d(j,8)
#define fW(j)     mean2d(j,9)
#define fT(j)     mean2d(j,10)

#define rUf(j)    mean2d(j,11)
#define rVf(j)    mean2d(j,12)
#define rWf(j)    mean2d(j,13)

#define rU_y(j)   mean2d(j,14)
#define rV_y(j)   mean2d(j,15)
#define rW_y(j)   mean2d(j,16)
#define fU_y(j)   mean2d(j,17)
#define fV_y(j)   mean2d(j,18)
#define fW_y(j)   mean2d(j,19)
#define rP_y(j)   mean2d(j,20)
#define rR_y(j)   mean2d(j,21)
#define rT_y(j)   mean2d(j,22)

#define Rxx(j)    mean2d(j,23)
#define Ryy(j)    mean2d(j,24)
#define Rzz(j)    mean2d(j,25)
#define Rxy(j)    mean2d(j,26)
#define Rxz(j)    mean2d(j,27)
#define Ryz(j)    mean2d(j,28)
#define Rxx_y(j)  mean2d(j,29)
#define Ryy_y(j)  mean2d(j,30)
#define Rzz_y(j)  mean2d(j,31)
#define Rxy_y(j)  mean2d(j,32)
#define Rxz_y(j)  mean2d(j,33)
#define Ryz_y(j)  mean2d(j,34)

#define Txxy(j)   mean2d(j,35)
#define Txyy(j)   mean2d(j,36)
#define Tyyy(j)   mean2d(j,37)
#define Tzzy(j)   mean2d(j,38)
#define Txzy(j)   mean2d(j,39)
#define Tyzy(j)   mean2d(j,40)
#define Txxy_y(j) mean2d(j,41)
#define Txyy_y(j) mean2d(j,42)
#define Tyyy_y(j) mean2d(j,43)
#define Tzzy_y(j) mean2d(j,44)
#define Txzy_y(j) mean2d(j,45)
#define Tyzy_y(j) mean2d(j,46)

#define rP2(j)    mean2d(j,47)
#define rR2(j)    mean2d(j,48)
#define rT2(j)    mean2d(j,49)
#define fT2(j)    mean2d(j,50)
#define rRT(j)    mean2d(j,51)
#define T_ac(j)   mean2d(j,52)
#define T_en(j)   mean2d(j,53)

#define PIxx(j)   mean2d(j,54)
#define PIyy(j)   mean2d(j,55)
#define PIzz(j)   mean2d(j,56)
#define PIxy(j)   mean2d(j,57)
#define PIxz(j)   mean2d(j,58)
#define PIyz(j)   mean2d(j,59)

#define rGamma(j) mean2d(j,60)
#define c2(j)     mean2d(j,61)
#define ugradp(j) mean2d(j,62)

#define Tau_xx(j) mean2d(j,63)
#define Tau_yy(j) mean2d(j,64)
#define Tau_zz(j) mean2d(j,65)
#define Tau_xy(j) mean2d(j,66)
#define Tau_xz(j) mean2d(j,67)
#define Tau_yz(j) mean2d(j,68)

#define Tau_xy_y(j) mean2d(j,69)
#define Tau_yy_y(j) mean2d(j,70)
#define Tau_yz_y(j) mean2d(j,71)

#define rR2_y(j)    mean2d(j,72)

#define Ty1(j)      mean2d(j,73)
#define Ty2(j)      mean2d(j,74)
#define Ty3(j)      mean2d(j,75)
#define Ty1_y(j)    mean2d(j,76)
#define Ty2_y(j)    mean2d(j,77)
#define Ty3_y(j)    mean2d(j,78)

#define vortx(j)    mean2d(j,79)
#define vorty(j)    mean2d(j,80)
#define vortz(j)    mean2d(j,81)
#define vortx2(j)   mean2d(j,82)
#define vorty2(j)   mean2d(j,83)
#define vortz2(j)   mean2d(j,84)

#define re(j)       mean2d(j,85)
#define fe(j)       mean2d(j,86)
#define rh(j)       mean2d(j,87)
#define fh(j)       mean2d(j,88)
#define rs(j)       mean2d(j,89)
#define fs(j)       mean2d(j,90)
#define rs2(j)      mean2d(j,91)
#define fs2(j)      mean2d(j,92)

#define psat(j)     mean2d(j,93)
#define bfreq_fr(j) mean2d(j,94)
#define bfreq_eq(j) mean2d(j,95)
#define lapse_fr(j) mean2d(j,96)
#define lapse_eq(j) mean2d(j,97)
#define potem_fr(j) mean2d(j,98)
#define potem_eq(j) mean2d(j,99)
#define p_mod(j)    mean2d(j,100)

#define rB(j)       mean2d(j,101)
#define Bxx(j)      mean2d(j,102)
#define Byy(j)      mean2d(j,103)
#define Bzz(j)      mean2d(j,104)
#define Bxy(j)      mean2d(j,105)
#define Bxz(j)      mean2d(j,106)
#define Byz(j)      mean2d(j,107)

#define Fxx(j)      mean2d(j,108)
#define Fyy(j)      mean2d(j,109)
#define Fzz(j)      mean2d(j,110)
#define Fxy(j)      mean2d(j,111)
#define Fxz(j)      mean2d(j,112)
#define Fyz(j)      mean2d(j,113)

#define rSb(j)      mean2d(j,114)
#define rB_y(j)     mean2d(j,115)

#define L_AVGMAX 115

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
!# 2007/05/11 - J.P. Mellado
!#              Cleaned. Adding vorticity and buoyancy frequency.
!# 2007/05/11 - J.P. Mellado
!#              Adding lateral derivative skewness and variances, 
!#              redefining (isotropic) Re_\lambda and adding 
!#              longitudinal Re_\lambda
!# 2007/06/22 - J.P. Mellado
!#              Adding energy, enthalpy and entropy.
!# 2007/07/31 - J.P. Mellado
!#              Adding lapse rates, potential temperatures and 
!#              Rxx transport equation.
!# 2007/08/30 - J.P. Mellado
!#              Adding derivative flatness.
!# 2011/03/23 - C. Ansorge (cedrick@gmx.net)
!#              modified stress budgets to include horizontal stresses,
!#              Rxz and Ryz stress budgets and Coriolis and Buoyancy Tensor
!# 
!########################################################################
!# DESCRIPTION
!#
!# Calculate statistics along horizontal planes.
!# 
!# Assumes statistical homogeneity in xOz, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# To be used in the incompressible case, the array p has been 
!# pointed to dudz and the pressure field is stored there; do not
!# use array dudz until pressure block
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE AVG_FLOW_TEMPORAL_LAYER(y,dx,dy,dz, q,s,&
     dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, mean2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL, MAX_PROF
  USE DNS_GLOBAL, ONLY : imode_eqns, imode_flow, itransport, ibodyforce
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, inb_scal, imode_fdm, i1bc,j1bc,k1bc, area, scaley
  USE DNS_GLOBAL, ONLY : froude, visc, rossby
  USE DNS_GLOBAL, ONLY : body_vector, body_param
  USE DNS_GLOBAL, ONLY : rotn_vector, icoriolis_y
  USE DNS_GLOBAL, ONLY : iprof_i, mean_i, delta_i, thick_i, ycoor_i, prof_i
  USE DNS_GLOBAL, ONLY : delta_u, ycoor_u
  USE DNS_GLOBAL, ONLY : mean_rho, delta_rho, ycoor_rho
  USE THERMO_GLOBAL, ONLY : imixture, MRATIO, GRATIO
  USE THERMO_GLOBAL, ONLY : THERMO_AI, WGHT_INV
#ifdef TRACE_ON
  USE DNS_CONSTANTS, ONLY : tfile
#endif
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(*),                INTENT(IN)    :: y,dx,dy,dz
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)    :: q, s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, wrk3d
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: mean2d, wrk1d
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk2d

  TARGET q, dudz

! -------------------------------------------------------------------
  TINTEGER i, j, k, flag !, flag_buoyancy
  TREAL AVG_IK, SIMPSON_NU, FLOW_SHEAR_TEMPORAL, UPPER_THRESHOLD, LOWER_THRESHOLD
  TREAL delta_m, delta_m_p, delta_w
  TINTEGER iprof
  TREAL ycenter, thick, delta, mean, param(MAX_PROF)
  TREAL rho_min, rho_max, delta_hb01, delta_ht01, delta_h01, mixing1, mixing2
  TREAL delta_hb25, delta_ht25, delta_h25
  TREAL u_friction, d_friction, a_friction
  TREAL eta, lambda_x,lambda_y,lambda_z, relambda_x,relambda_y,relambda_z, relambda_iso
  TREAL dil, dummy
  TREAL Cxx, Pxx, Exx, Gxx, Dxx, Rxx_t
  TREAL Cxy, Pxy, Exy, Gxy, Dxy, Rxy_t
  TREAL Cyy, Pyy, Eyy, Gyy, Dyy, Ryy_t
  TREAL Czz, Pzz, Ezz, Gzz, Dzz, Rzz_t
  TREAL Cyz, Pyz, Eyz, Gyz, Dyz, Ryz_t
  TREAL Cxz, Pxz, Exz, Gxz, Dxz, Rxz_t
  TREAL Con, Prd, Eps, Gkin,Dkin,Tke_t, Tke, Pi, Buo, Ty_y
  TREAL Phi, M_t, Pot, SourcePot
  TREAL var_dil, var_uy,  var_uz,  var_vx,  var_vz,  var_wx,  var_wy,  var_x,   var_y,   var_z
  TREAL skew_ux, skew_uy, skew_uz, skew_vx, skew_vy, skew_vz, skew_wx, skew_wy, skew_wz
  TREAL flat_ux, flat_uy, flat_uz, flat_vx, flat_vy, flat_vz, flat_wx, flat_wy, flat_wz
  TREAL rRP, re2, fe2, rh2, fh2, rho_ac, rho_en
  TREAL tau11, tau22, tau33, tau12, tau23, tau13
  TREAL up, vp, wp, upy, vpy, wpy
  TREAL c23, prefactor
  TREAL r_prime, p_prime, T_prime, u_prime, v_prime, w_prime
  TREAL rey_flux_x, rey_flux_y, rey_flux_z
  TREAL rey_dil1, rey_dil2, rey_trp, rey_prod, rey_conv
  TREAL eddy_diff, eddy_visc, eddy_prandtl, ri_f, ri_g

  TREAL VAUXPRE(4), VAUXPOS(14), L_RATIO, Q_RATIO, WMEAN_INV, C_RATIO
  TINTEGER ivauxpre, ivauxpos
  CHARACTER*32 fname
  CHARACTER*250 line1
  CHARACTER*1300 line2

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u,v,w,p, e,rho, vis

! ###################################################################
  WRITE(line1,*) itime; line1 = 'Calculating flow statistics at It'//TRIM(ADJUSTL(line1))//'...'
  CALL IO_WRITE_ASCII(lfile,line1)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 1')
#endif

  c23 = C_2_R/C_3_R

  prefactor = GRATIO*MRATIO ! (gama0-C_1_R)*mach*mach

  IF ( MAX_AVG_TEMPORAL .LT. L_AVGMAX ) THEN
     CALL IO_WRITE_ASCII(efile,'AVG_FLOW_TEMPORAL_LAYER. Not enough space.')
     CALL DNS_STOP(LES_ERROR_AVGTMP)
  ENDIF
  mean2d(:,1:L_AVGMAX) = C_0_R

! Define pointers
  u   => q(:,:,:,1)
  v   => q(:,:,:,2)
  w   => q(:,:,:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     e   => q(:,:,:,4)
     rho => q(:,:,:,5)
     p   => q(:,:,:,6)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)
  ELSE
     p   => dudz
  ENDIF

! in case we need the buoyancy statistics
  ! IF ( ibodyforce .EQ. EQNS_BOD_BILINEAR           .OR. &
  !      ibodyforce .EQ. EQNS_BOD_QUADRATIC          .OR. &
  !      imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
  !    flag_buoyancy = 1
  ! ELSE 
  !    flag_buoyancy = 0   
  ! ENDIF

! -------------------------------------------------------------------
! TkStat file; header
! -------------------------------------------------------------------
#ifdef USE_MPI
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)

  IF ( ims_pro .EQ. 0 ) THEN
#endif
     WRITE(fname,*) itime; fname='avg'//TRIM(ADJUSTL(fname))

#ifdef USE_RECLEN
     OPEN(UNIT=23, RECL=1613, FILE=fname, STATUS='unknown') ! this is probably outdated
#else
     OPEN(UNIT=23, FILE=fname, STATUS='unknown')
#endif

     WRITE(23, '(A8,E14.7E3)') 'RTIME = ', rtime

! Independent variables
     line2='I J Y SM SW SR'

! Dependent variables depending on y and t
     line1 = 'rR rU rV rW rP rT re rh rs rB fU fV fW fT fe fh fs'
     WRITE(i23,1010) 'GROUP = Mean '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Tke Rxx Ryy Rzz Rxy Rxz Ryz rP2 rR2 rT2 fT2 re2 fe2 rh2 fh2 rs2 fs2'
     WRITE(i23,1010) 'GROUP = Fluctuations '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Wx Wy Wz Wx2 Wy2 Wz2'
     WRITE(i23,1010) 'GROUP = Vorticity '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Rxx_t Bxx Cxx Pxx Exx PIxx Fxx Txxy_y Txxy Gxx Dxx'
     WRITE(i23,1010) 'GROUP = RxxBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Ryy_t Byy Cyy Pyy Eyy PIyy Fyy Tyyy_y Tyyy Gyy Dyy '
     WRITE(i23,1010) 'GROUP = RyyBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Rzz_t Bzz Czz Pzz Ezz PIzz Fzz Tzzy_y Tzzy Gzz Dzz '
     WRITE(i23,1010) 'GROUP = RzzBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Rxy_t Bxy Cxy Pxy Exy PIxy Fxy Txyy_y Txyy Gxy Dxy '
     WRITE(i23,1010) 'GROUP = RxyBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Rxz_t Bxz Cxz Pxz Exz PIxz Fxz Txzy_y Txzy Gxz Dxz '
     WRITE(i23,1010) 'GROUP = RxzBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Ryz_t Byz Cyz Pyz Eyz PIyz Fyz Tyzy_y Tyzy Gyz Dyz '
     WRITE(i23,1010) 'GROUP = RyzBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Tke_t Buo Con Prd Eps Pi Trp Trp1 Trp2 Trp3 Trp1_y Trp2_y Trp3_y G D Phi UgradP'
     WRITE(i23,1010) 'GROUP = TkeBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Eta LambdaUx LambdaVy LambdaWz ReLambdaUx ReLambdaVy ReLambdaWz ReLambdaIso '
     WRITE(i23,1010) 'GROUP = Scales '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'VarDilatation VarUx VarUy VarUz VarVx VarVy VarVz VarWx VarWy VarWz '&
          //'SkewUx SkewUy SkewUz SkewVx SkewVy SkewVz SkewWx SkewWy SkewWz '&
          //'FlatUx FlatUy FlatUz FlatVx FlatVy FlatVz FlatWx FlatWy FlatWz '
     WRITE(i23,1010) 'GROUP = DerivativeFluctuations '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'gamma C2 Rho_ac Rho_en T_ac T_en M_t rRP rRT'
     WRITE(i23,1010) 'GROUP = Acoustics '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'RhoFluxX RhoFluxY RhoFluxZ RhoDil1 RhoDil2 RhoTrp RhoProd RhoConv '
     WRITE(i23,1010) 'GROUP = RhoBudget '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     line1 = 'Pot Source rSb BuoyFreq_fr BuoyFreq_eq LapseRate_fr LapseRate_eq '&
          //'PotTemp_fr PotTemp_eq SaturationPressure rP0 rPmod Ri_f Ri_g'
     WRITE(i23,1010) 'GROUP = Stratification '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

! other dependent variables, to be sorted
     line1 = 'EddyDiff EddyVisc TurbPrandtl'
     WRITE(i23,1010) 'GROUP = TurbDiffusivities '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

! dependent variables function of t only
     IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        line1 = 'Delta_m Delta_m_p Delta_w'
        WRITE(i23,1010) 'GROUP = ShearThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

        line1 = 'Delta_hb01 Delta_ht01 Delta_h01 Delta_hb25 Delta_ht25 Delta_h25 '&
             //'mixing_Youngs mixing_Cook'
        WRITE(i23,1010) 'GROUP = MixingThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

        line1 = 'FrictionVelocity FrictionThickness FrictionAngle'
        WRITE(i23,1010) 'GROUP = FrictionTerms '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))
     ENDIF

     WRITE(i23,1010) TRIM(ADJUSTL(line2))

1010 FORMAT(A)

#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
! Main data
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 2')
#endif

  DO j = 1,jmax
     rU(j) = AVG_IK(imax,jmax,kmax, j, u, dx, dz, area)
     rV(j) = AVG_IK(imax,jmax,kmax, j, v, dx, dz, area)
     rW(j) = AVG_IK(imax,jmax,kmax, j, w, dx, dz, area)

     IF     ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
        rR(j) = C_1_R; fU(j) = rU(j); fV(j) = rV(j); fW(j) = rW(j)

     ELSEIF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC      ) THEN ! not yet developed
        rR(j) = C_1_R; fU(j) = rU(j); fV(j) = rV(j); fW(j) = rW(j)

     ELSE
        rR(j) = AVG_IK(imax,jmax,kmax, j, rho, dx, dz, area)
        DO k = 1,kmax; DO i = 1,imax
           wrk3d(i,1,k) = rho(i,j,k)*u(i,j,k)
           wrk3d(i,2,k) = rho(i,j,k)*v(i,j,k)
           wrk3d(i,3,k) = rho(i,j,k)*w(i,j,k)
        ENDDO; ENDDO
        fU(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx, dz, area)/rR(j)
        fV(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx, dz, area)/rR(j)
        fW(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx, dz, area)/rR(j)

     ENDIF

     rUf(j) = rU(j) - fU(j)
     rVf(j) = rV(j) - fV(j)
     rWf(j) = rW(j) - fW(j)

     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        DO k = 1,kmax; DO i = 1,imax
           up = u(i,j,k) - fU(j)
           vp = v(i,j,k) - fV(j)
           wp = w(i,j,k) - fW(j)
           
           wrk3d(i,1,k) = up**2
           wrk3d(i,2,k) = vp**2
           wrk3d(i,3,k) = wp**2
           wrk3d(i,4,k) = up*vp
           wrk3d(i,5,k) = up*wp
           wrk3d(i,6,k) = vp*wp
           
        ENDDO; ENDDO

     ELSE
        DO k = 1,kmax; DO i = 1,imax
           up = u(i,j,k) - fU(j)
           vp = v(i,j,k) - fV(j)
           wp = w(i,j,k) - fW(j)
           
           wrk3d(i,1,k) = rho(i,j,k)*up**2
           wrk3d(i,2,k) = rho(i,j,k)*vp**2
           wrk3d(i,3,k) = rho(i,j,k)*wp**2
           wrk3d(i,4,k) = rho(i,j,k)*up*vp
           wrk3d(i,5,k) = rho(i,j,k)*up*wp
           wrk3d(i,6,k) = rho(i,j,k)*vp*wp
              
           wrk3d(i,7,k) = (rho(i,j,k) - rR(j))**2
              
        ENDDO; ENDDO
        rR2(j) = AVG_IK(imax,jmax,kmax, i7, wrk3d, dx, dz, area)

     ENDIF
     Rxx(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx, dz, area)/rR(j)
     Ryy(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx, dz, area)/rR(j)
     Rzz(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx, dz, area)/rR(j)
     Rxy(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx, dz, area)/rR(j)
     Rxz(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx, dz, area)/rR(j)
     Ryz(j) = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx, dz, area)/rR(j)

  ENDDO

! ###################################################################
! Pressure; array p used only in this section
!
! dudx = du/dx
! dudy = du/dy
! dudz = 
! dvdx = dv/dx
! dvdy = dv/dy
! dvdz = 
! dwdx =       ; dp/dx
! dwdy =       ; dp/dy
! dwdz = dw/dz ; dp/dz
! ###################################################################
! mean and fluctutation
  DO j = 1,jmax
     rP(j) = AVG_IK(imax,jmax,kmax, j, p, dx, dz, area)

     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = (p(i,j,k)-rP(j))**2
     ENDDO; ENDDO
     rP2(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx, dz, area)
  ENDDO

! Pressure convection term
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, p, dwdx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, p, dwdy, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, p, dwdz, i0,i0, wrk1d,wrk2d,wrk3d)
  DO j = 1,jmax
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = u(i,j,k)*dwdx(i,j,k)+v(i,j,k)*dwdy(i,j,k)+w(i,j,k)*dwdz(i,j,k)
     ENDDO; ENDDO
     ugradp(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx, dz, area)
  ENDDO

! Pressure Strain Terms
! 9 derivatives are here recomputed; ok, this routine is not called that often
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u, dudx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, dvdy, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, w, dwdz, i0,i0, wrk1d,wrk2d,wrk3d)

  DO j = 1,jmax
     DO k = 1,kmax; DO i = 1,imax
        p_prime  = p(i,j,k)-rP(j)
        wrk3d(i,1,k) = p_prime*dudx(i,j,k)
        wrk3d(i,2,k) = p_prime*dvdy(i,j,k)               ! no need to substract rV_y
        wrk3d(i,3,k) = p_prime*dwdz(i,j,k)
     ENDDO; ENDDO

     PIxx(j) = C_2_R*AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     PIyy(j) = C_2_R*AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     PIzz(j) = C_2_R*AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)

  ENDDO

  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, u, dudy, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v, dvdx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, u, dwdz, i0,i0, wrk1d,wrk2d,wrk3d) !dudz not free
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, w, dwdx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, v, dvdz, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, w, dwdy, i0,i0, wrk1d,wrk2d,wrk3d)

  DO j = 1,jmax
     DO k = 1,kmax; DO i = 1,imax
        p_prime = p(i,j,k) - rP(j)
        wrk3d(i,1,k) = p_prime*(dudy(i,j,k)+dvdx(i,j,k)) ! no need to substract rU_y
        wrk3d(i,2,k) = p_prime*(dwdz(i,j,k)+dwdx(i,j,k))
        wrk3d(i,3,k) = p_prime*(dvdz(i,j,k)+dwdy(i,j,k)) ! no need to substract rW_y
     ENDDO; ENDDO

     PIxy(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     PIxz(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     PIyz(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     
  ENDDO

! Turbulent transport terms
! Only the pressure-velocity-correlation contribution; rest below
  DO j = 1,jmax
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = (p(i,j,k)-rP(j))*(u(i,j,k)-fU(j))
        wrk3d(i,2,k) = (p(i,j,k)-rP(j))*(v(i,j,k)-fV(j))
        wrk3d(i,3,k) = (p(i,j,k)-rP(j))*(w(i,j,k)-fW(j))
     ENDDO; ENDDO
     Txyy(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     Ty2(j)  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     Tyyy(j) = C_2_R*Ty2(j)
     Tyzy(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
  ENDDO

! ###################################################################
! Thermodynamic variables
!
! dudx = gamma
! dudy = dsdy
! dudz = dTdy
! dvdx = cp
! dvdy = drdy
! dvdz = psat
! dwdx = T
! dwdy = dpdy
! dwdz = 
! ###################################################################
#define T_LOC(i,j,k)     dwdx(i,j,k)
#define GAMMA_LOC(i,j,k) dudx(i,j,k)
#define S_LOC(i,j,k)     dwdz(i,j,k)

  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN

! -------------------------------------------------------------------
! Main fields
! -------------------------------------------------------------------
     CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T_LOC(1,1,1), wrk3d)
     CALL THERMO_GAMMA(imax,jmax,kmax, s, T_LOC(1,1,1), GAMMA_LOC(1,1,1))
     CALL THERMO_ENTROPY(imax,jmax,kmax, s, T_LOC(1,1,1), p, S_LOC(1,1,1))

     DO j = 1,jmax
! Means
        rT(j)     = AVG_IK(imax, jmax, kmax, j, T_LOC(1,1,1),     dx,dz, area)
        re(j)     = AVG_IK(imax, jmax, kmax, j, e,                dx,dz, area)
        rs(j)     = AVG_IK(imax, jmax, kmax, j, S_LOC(1,1,1),     dx,dz, area)
        rGamma(j) = AVG_IK(imax, jmax, kmax, j, GAMMA_LOC(1,1,1), dx,dz, area)

        DO k = 1,kmax
           DO i = 1,imax
              wrk3d(i,1,k) = rho(i,j,k)*e(i,j,k)
              wrk3d(i,2,k) = e(i,j,k) + prefactor*p(i,j,k)/rho(i,j,k)
              wrk3d(i,3,k) = rho(i,j,k)*e(i,j,k) + prefactor*p(i,j,k)
              wrk3d(i,4,k) = rho(i,j,k)*S_LOC(i,j,k)
              wrk3d(i,5,k) = rho(i,j,k)*T_LOC(i,j,k)
              wrk3d(i,6,k) = GAMMA_LOC(i,j,k)*p(i,j,k)/rho(i,j,k)
           ENDDO
        ENDDO
        fe(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)/rR(j)
        rh(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
        fh(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)/rR(j)
        fs(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)/rR(j)
        fT(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)/rR(j)
        c2(j) = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area) ! speed of sound
      
! Fluctuations  
        DO k = 1,kmax
           DO i = 1,imax
              wrk3d(i,1,k) =            (S_LOC(i,j,k)-rs(j))**2
              wrk3d(i,2,k) = rho(i,j,k)*(S_LOC(i,j,k)-fs(j))**2
              wrk3d(i,3,k) =            (T_LOC(i,j,k)-rT(j))**2
              wrk3d(i,4,k) = rho(i,j,k)*(T_LOC(i,j,k)-fT(j))**2
              wrk3d(i,5,k) = (rho(i,j,k)-rR(j))*(T_LOC(i,j,k)-fT(j))
           ENDDO
        ENDDO
        rs2(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
        fs2(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)/rR(j)
        rT2(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
        fT2(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)/rR(j)
        rRT(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     
! Acoustic and entropic density and temperature fluctuations
        DO k=1, kmax
           DO i=1, imax
              r_prime = rho(i,j,k)   - rR(j)
              p_prime = p(i,j,k)     - rP(j)
              T_prime = T_LOC(i,j,k) - fT(j)
              rho_ac = p_prime/c2(j)
              rho_en = r_prime - rho_ac
              T_ac(j) = fT(j)*(p_prime/rP(j)-rho_ac/rR(j))
              T_en(j) = T_prime - T_ac(j)
              
              wrk3d(i,1,k) = T_ac(j)*T_ac(j)
              wrk3d(i,2,k) = T_en(j)*T_en(j)
           ENDDO
        ENDDO

        T_ac(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
        T_en(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)

     ENDDO

! -------------------------------------------------------------------
! Buoyancy frequency & saturation pressure
! -------------------------------------------------------------------
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, rho, dvdy, i0, i0, wrk1d,wrk2d,wrk3d)
  
     CALL THERMO_POLYNOMIAL_PSAT(imax, jmax, kmax, T_LOC(1,1,1), dvdz)
     CALL THERMO_CP(imax, jmax, kmax, s, GAMMA_LOC(1,1,1), dvdx)
     
     DO j = 1,jmax
        DO k = 1,kmax
           DO i = 1,imax
              wrk3d(i,1,k) = dwdy(i,j,k)/p(i,j,k)/GAMMA_LOC(i,j,k)-dvdy(i,j,k)/rho(i,j,k)
              wrk3d(i,2,k) = C_1_R/dvdx(i,j,k)
              wrk3d(i,3,k) = T_LOC(i,j,k)*( (MRATIO*p(i,j,k))**(C_1_R/GAMMA_LOC(i,j,k)-C_1_R) )
           ENDDO
        ENDDO
        bfreq_fr(j) =-AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)*body_vector(2)
        lapse_fr(j) =-AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)*body_vector(2)*prefactor
        potem_fr(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)        
        psat(j)     = AVG_IK(imax,jmax,kmax,  j, dvdz,  dx,dz, area)
     ENDDO
     
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, T_LOC(1,1,1), dudz, &
             i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, s(1,1,1,2), dudy, &
             i0,i0, wrk1d,wrk2d,wrk3d)
        DO j = 1,jmax
           DO k=1, kmax
              DO i=1, imax
                 L_RATIO = THERMO_AI(6,1,1)-THERMO_AI(6,1,3)&
                      - (THERMO_AI(1,1,3)-THERMO_AI(1,1,1)+GRATIO*WGHT_INV(1))*T_LOC(i,j,k)
                 L_RATIO = L_RATIO/(GRATIO*WGHT_INV(1)*T_LOC(i,j,k))
                 Q_RATIO = C_1_R/(MRATIO*p(i,j,k)/dvdz(i,j,k)-C_1_R)
                 WMEAN_INV = (Q_RATIO+C_1_R)*(C_1_R-s(i,j,k,1))*WGHT_INV(2)
                 
                 wrk3d(i,1,k) = (C_1_R+Q_RATIO*L_RATIO)/WMEAN_INV/&
                      (GAMMA_LOC(i,j,k)/(GAMMA_LOC(i,j,k)-C_1_R)+Q_RATIO*L_RATIO*L_RATIO)

                 wrk3d(i,2,k) = (dudz(i,j,k)-body_vector(2)*MRATIO*wrk3d(i,1,k))/T_LOC(i,j,k)*&
                      (C_1_R+WGHT_INV(1)/WGHT_INV(2)*L_RATIO/(C_1_R-s(i,j,k,1)))
                 wrk3d(i,2,k) = wrk3d(i,2,k) - WGHT_INV(2)/WMEAN_INV*dudy(i,j,k)
                 
                 C_RATIO = THERMO_AI(1,1,2)+s(i,j,k,1)*(THERMO_AI(1,1,3)-THERMO_AI(1,1,2))
                 C_RATIO = (C_1_R-s(i,j,k,1))*GRATIO*WGHT_INV(2)/C_RATIO
                 
                 wrk3d(i,3,k) = T_LOC(i,j,k)/( (MRATIO*p(i,j,k))**C_RATIO )*exp(Q_RATIO*C_RATIO*L_RATIO)
                 wrk3d(i,3,k) = wrk3d(i,3,k)*(C_1_R+Q_RATIO)**C_RATIO&
                      /((MRATIO*p(i,j,k)/dvdz(i,j,k))**(Q_RATIO*C_RATIO))
              ENDDO
           ENDDO
           lapse_eq(j) =-AVG_IK(imax,jmax,kmax, i1, wrk3d, dx, dz, area)*body_vector(2)*MRATIO
           bfreq_eq(j) =-AVG_IK(imax,jmax,kmax, i2, wrk3d, dx, dz, area)*body_vector(2)
           potem_eq(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx, dz, area)
           
        ENDDO
     ELSE
        DO j = 1,jmax
           lapse_eq(j) = C_0_R
           bfreq_eq(j) = C_0_R
           potem_eq(j) = C_0_R
        ENDDO
     ENDIF
     
     CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rP(1), p_mod(1), i0,i0, wrk1d,wrk2d,wrk3d)
     DO j = 1,jmax
        p_mod(j) = -p_mod(j)+body_vector(2)*rR(j)
     ENDDO
     
#undef GAMMA_LOC
#undef T_LOC
#undef S_LOC

  ENDIF

! ###################################################################
! Potential energy
!
! dudx = buoyancy
! dudy = scalar gradient GiGi
!
! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN

     IF ( ibodyforce .NE. EQNS_NONE ) THEN
! buoyancy field
     ycenter = y(1) + scaley*ycoor_i(1)
     iprof   = iprof_i(1)
     thick   = thick_i(1)
     delta   = delta_i(1)
     mean    = mean_i (1)
     param(:)= prof_i (:,1)
     DO j = 1,jmax
        wrk1d(j,1) = FLOW_SHEAR_TEMPORAL(iprof, thick, delta, mean, ycenter, param, y(j))
        wrk1d(j,3) = C_0_R
     ENDDO
     flag = EQNS_BOD_LINEAR
     CALL FI_BUOYANCY(flag,       i1,jmax,i1,     body_param, wrk1d(1,1), wrk1d(1,2), wrk1d(1,3))
     CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s,          dudx,       wrk1d(1,2))

! scalar gradient for the source term; prefactor in dummy
     ! IF ( flag_buoyancy .EQ. 1 ) THEN
     !    CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, s, dvdx, i0,i0, wrk1d,wrk2d,wrk3d)
     !    CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, s, dvdy, i0,i0, wrk1d,wrk2d,wrk3d)
     !    CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, s, dvdz, i0,i0, wrk1d,wrk2d,wrk3d)
     !    dudy = dvdx*dvdx + dvdy*dvdy + dvdz*dvdz
     
     !    dummy =-(body_param(2)-body_param(1))/(C_1_R-body_param(3)) - body_param(2)/body_param(3)
     !    dummy = dummy/(C_4_R*body_param(4)) * visc/schmidt(1) / froude
     ! ENDIF

! buoyancy terms
     DO j = 1,jmax
        rB(j) = AVG_IK(imax,jmax,kmax, j, dudx, dx,dz, area)
        DO k = 1,kmax; DO i = 1,imax
           wrk3d(i,1,k) = (u(i,j,k)-rU(j))*(dudx(i,j,k)-rB(j))
           wrk3d(i,2,k) = (v(i,j,k)-rV(j))*(dudx(i,j,k)-rB(j))
           wrk3d(i,3,k) = (w(i,j,k)-rW(j))*(dudx(i,j,k)-rB(j))
        ENDDO; ENDDO

        rB(j) = (C_1_R/froude)*rB(j)

        Bxx(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
        Byy(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
        Bzz(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)

        Bxy(j) = Bxx(j)*body_vector(2) + Byy(j)*body_vector(1) ! body_vector includes the Froude
        Bxz(j) = Bxx(j)*body_vector(3) + Bzz(j)*body_vector(1)
        Byz(j) = Byy(j)*body_vector(3) + Bzz(j)*body_vector(2)

        Bxx(j) = C_2_R*Bxx(j)*body_vector(1)
        Byy(j) = C_2_R*Byy(j)*body_vector(2)
        Bzz(j) = C_2_R*Bzz(j)*body_vector(3)

        ! IF ( flag_buoyancy .EQ. 1 ) THEN
        !    DO k = 1,kmax; DO i = 1,imax
        !       wrk3d(i,4,k) =-dudy(i,j,k)*dummy / &
        !         (cosh(C_05_R*(s(i,j,k,1)-body_param(3))/body_param(4)))**2
        !    ENDDO; ENDDO
        !    rSb(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
        ! ENDIF

     ENDDO

     ENDIF

  ELSE ! Compressible case is not yet finished
     DO j = 1,jmax
        Bxx(j) =-rR(j)*rUf(j)*body_vector(1)
        Byy(j) =-rR(j)*rVf(j)*body_vector(2)
        Bzz(j) =-rR(j)*rWf(j)*body_vector(3)
        rSb(j) = C_0_R
     ENDDO

  ENDIF

! ###################################################################
! # Array storage of velocity gradient tensor
! #
! # dudx = d U / d x
! # dudy = d U / d y
! # dudz = d U / d z
! # dvdx = d V / d x
! # dvdy = d V / d y
! # dvdz = d V / d z
! # dwdx = d W / d x
! # dwdy = d W / d y
! # dwdz = d W / d z
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 3')
#endif

  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u, dudx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, u, dudy, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, u, dudz, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v, dvdx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, dvdy, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, v, dvdz, i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, w, dwdx, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, w, dwdy, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, w, dwdz, i0,i0, wrk1d,wrk2d,wrk3d)

! ###################################################################
! Vorticity
! ###################################################################
  DO j = 1,jmax
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = dwdy(i,j,k) - dvdz(i,j,k)
        wrk3d(i,2,k) = dudz(i,j,k) - dwdx(i,j,k)
        wrk3d(i,3,k) = dvdx(i,j,k) - dudy(i,j,k)
     ENDDO; ENDDO
     vortx(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     vorty(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     vortz(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)

     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = (wrk3d(i,1,k)-vortx(j))**2
        wrk3d(i,2,k) = (wrk3d(i,2,k)-vorty(j))**2
        wrk3d(i,3,k) = (wrk3d(i,3,k)-vortz(j))**2
     ENDDO; ENDDO
     vortx2(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     vorty2(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     vortz2(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)

  ENDDO

! ##################################################################
! Averaged viscous shear-stress tensor
! ##################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 4')
#endif

  DO j = 1,jmax
     IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        DO k = 1,kmax; DO i = 1,imax
           dil = (dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))*c23
              
           wrk3d(i,1,k) = vis(i,j,k)*(C_2_R*dudx(i,j,k)-dil)
           wrk3d(i,2,k) = vis(i,j,k)*(C_2_R*dvdy(i,j,k)-dil)
           wrk3d(i,3,k) = vis(i,j,k)*(C_2_R*dwdz(i,j,k)-dil)
           wrk3d(i,4,k) = vis(i,j,k)*(dudy(i,j,k)+dvdx(i,j,k))
           wrk3d(i,5,k) = vis(i,j,k)*(dudz(i,j,k)+dwdx(i,j,k))
           wrk3d(i,6,k) = vis(i,j,k)*(dvdz(i,j,k)+dwdy(i,j,k))
        ENDDO; ENDDO

     ELSE
        DO k = 1,kmax; DO i=1,imax
           dil = (dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))*c23
              
           wrk3d(i,1,k) = C_2_R*dudx(i,j,k)-dil
           wrk3d(i,2,k) = C_2_R*dvdy(i,j,k)-dil
           wrk3d(i,3,k) = C_2_R*dwdz(i,j,k)-dil
           wrk3d(i,4,k) = dudy(i,j,k)+dvdx(i,j,k)
           wrk3d(i,5,k) = dudz(i,j,k)+dwdx(i,j,k)
           wrk3d(i,6,k) = dvdz(i,j,k)+dwdy(i,j,k)
        ENDDO; ENDDO

     ENDIF

     Tau_xx(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area) *visc
     Tau_yy(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area) *visc
     Tau_zz(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area) *visc
     Tau_xy(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area) *visc
     Tau_xz(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area) *visc
     Tau_yz(j) = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area) *visc

  ENDDO

! ##################################################################
! Y Mean Derivatives
! ##################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 5')
#endif

  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Rxx(1), Rxx_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Ryy(1), Ryy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Rzz(1), Rzz_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Rxy(1), Rxy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Rxz(1), Rxz_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Ryz(1), Ryz_y(1), i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, fU(1), fU_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, fV(1), fV_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, fW(1), fW_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rU(1), rU_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rV(1), rV_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rW(1), rW_y(1), i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rP(1),  rP_y(1),  i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rT(1),  rT_y(1),  i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rR(1),  rR_y(1),  i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rR2(1), rR2_y(1), i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Tau_xy(1), Tau_xy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Tau_yy(1), Tau_yy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Tau_yz(1), Tau_yz_y(1), i0,i0, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rB(1), rB_y(1), i0,i0, wrk1d,wrk2d,wrk3d)

! ##################################################################
! Global quantites
! ##################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 6')
#endif

  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN

! -------------------------------------------------------------------
! Based on delta_u 
! -------------------------------------------------------------------
! Vorticity thickness and momentum thickness
     IF ( ABS(delta_u) .GT. C_SMALL_R ) THEN
        delta_w = delta_u/MAXVAL(ABS(fU_y(1:jmax)))

        DO j=1, jmax
           wrk1d(j,1) = rR(j)*(C_025_R-(fU(j)/delta_u)**2)
        ENDDO
        delta_m = SIMPSON_NU(jmax, wrk1d, y)/mean_rho

        DO j=1, jmax
           wrk1d(j,1) = ( Tau_xy(j) -  rR(j)*Rxy(j) )*fU_y(j)
        ENDDO
        delta_m_p = SIMPSON_NU(jmax, wrk1d, y)*C_2_R/(mean_rho*delta_u**3)

     ELSE
        delta_w   = C_1_R
        delta_m   = C_1_R
        delta_m_p = C_1_R

     ENDIF

! -------------------------------------------------------------------
! Based on delta_rho
! -------------------------------------------------------------------
! 1% and 25% thickness
     IF ( imode_eqns .NE. DNS_EQNS_INCOMPRESSIBLE .AND. &
          imode_eqns .NE. DNS_EQNS_ANELASTIC      .AND. &
          ABS(delta_rho) .GT. C_SMALL_R                 ) THEN
        dummy = mean_rho + (C_05_R-C_1EM2_R)*delta_rho
        delta_hb01 = LOWER_THRESHOLD(jmax, dummy, rR(1), y)
        dummy = mean_rho - (C_05_R-C_1EM2_R)*delta_rho
        delta_ht01 = UPPER_THRESHOLD(jmax, dummy, rR(1), y)

        delta_hb01 = (y(1) + scaley*ycoor_rho) - delta_hb01  
        delta_ht01 = delta_ht01 - (y(1) + scaley*ycoor_rho)
        delta_h01  = delta_ht01 + delta_hb01

        dummy = mean_rho + (C_05_R-C_025_R)*delta_rho
        delta_hb25 = LOWER_THRESHOLD(jmax, dummy, rR(1), y)
        dummy = mean_rho - (C_05_R-C_025_R)*delta_rho
        delta_ht25 = UPPER_THRESHOLD(jmax, dummy, rR(1), y)

        delta_hb25 = (y(1) + scaley*ycoor_rho) - delta_hb25  
        delta_ht25 = delta_ht25 - (y(1) + scaley*ycoor_rho)
        delta_h25  = delta_ht25 + delta_hb25

! Mixing, Youngs definition
        rho_min = mean_rho - C_05_R*ABS(delta_rho)
        rho_max = mean_rho + C_05_R*ABS(delta_rho)
        DO k = 1,kmax
           DO i = 1,imax*jmax
              wrk3d(i,1,k) = (rho(i,1,k)-rho_min)*(rho_max-rho(i,1,k))
           ENDDO
        ENDDO
        DO j = 1,jmax
           wrk1d(j,1) = AVG_IK(imax, jmax, kmax, j, wrk3d, dx, dz, area)
        ENDDO
        mixing1 = SIMPSON_NU(jmax, wrk1d, y)
        DO j = 1,jmax
           wrk1d(j,1)=(rR(j)-rho_min)*(rho_max-rR(j))
        ENDDO
        mixing1 = mixing1/SIMPSON_NU(jmax, wrk1d, y)

! Mixing, Cook's definition
        rho_min = mean_rho - C_05_R*ABS(delta_rho)
        rho_max = mean_rho + C_05_R*ABS(delta_rho)
        DO k = 1,kmax
           DO i = 1,imax*jmax
              wrk3d(i,1,k) = MIN(rho(i,1,k)-rho_min,rho_max-rho(i,1,k))
           ENDDO
        ENDDO
        DO j = 1,jmax
           wrk1d(j,1)=AVG_IK(imax, jmax, kmax, j, wrk3d, dx, dz, area)
        ENDDO
        mixing2 = SIMPSON_NU(jmax, wrk1d, y)
        DO j = 1,jmax
           wrk1d(j,1) = MIN(rR(j)-rho_min,rho_max-rR(j))
        ENDDO
        mixing2 = mixing2/SIMPSON_NU(jmax, wrk1d, y)

     ELSE
        delta_h01 = C_1_R
        mixing1   = C_1_R
        mixing2   = C_1_R

     ENDIF

! -------------------------------------------------------------------
! Friction velocity terms
! -------------------------------------------------------------------
     u_friction = SQRT(Tau_xy(1)**2 + Tau_yz(1)**2)
     u_friction = SQRT(u_friction)

     d_friction = u_friction*rossby

     IF ( u_friction .GT. C_SMALL_R ) THEN 
        a_friction = ATAN2(Tau_yz(1),Tau_xy(1))*C_18_R*C_10_R/C_PI_R
     ELSE
        a_friction = C_0_R
     ENDIF

! -------------------------------------------------------------------
! Jet
! -------------------------------------------------------------------
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
! not developed yet

  ENDIF

! ##################################################################
! Turbulent transport terms
! p' contribution has been calculated before
! ##################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 10')
#endif

  DO j = 1,jmax
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN 
        DO k = 1,kmax
           DO i = 1,imax
              dil = (dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))*c23
              tau11 = visc*(C_2_R*dudx(i,j,k)-dil)-Tau_xx(j)
              tau22 = visc*(C_2_R*dvdy(i,j,k)-dil)-Tau_yy(j)
              tau33 = visc*(C_2_R*dwdz(i,j,k)-dil)-Tau_zz(j)
              tau12 = visc*(dudy(i,j,k)+dvdx(i,j,k))-Tau_xy(j)
              tau13 = visc*(dudz(i,j,k)+dwdx(i,j,k))-Tau_xz(j)
              tau23 = visc*(dvdz(i,j,k)+dwdy(i,j,k))-Tau_yz(j)
              up = u(i,j,k)-fU(j)
              vp = v(i,j,k)-fV(j)
              wp = w(i,j,k)-fW(j)
              
              wrk3d(i,1,k) = (up**2)*vp -C_2_R*up*tau12       
              wrk3d(i,2,k) = (vp**3)    -C_2_R*vp*tau22 
              wrk3d(i,3,k) = (wp**2)*vp -C_2_R*wp*tau23       
              wrk3d(i,4,k) = (vp**2)*up -      up*tau22-vp*tau12
              wrk3d(i,5,k) = up*vp*wp   -      up*tau23-wp*tau12  ! new Txzy
              wrk3d(i,6,k) = (vp**2)*wp -      vp*tau23-wp*tau22  ! new Tyzy

              wrk3d(i,7,k) = C_05_R*(up**2+vp**2+wp**2)*vp 
              wrk3d(i,9,k) =-up*tau12-vp*tau22-wp*tau23
           ENDDO
        ENDDO

     ELSE
        IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
           DO k = 1,kmax
              DO i = 1,imax
                 dil = (dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))*c23
                 
                 tau11 = visc*vis(i,j,k)*(C_2_R*dudx(i,j,k)-dil)-Tau_xx(j)
                 tau22 = visc*vis(i,j,k)*(C_2_R*dvdy(i,j,k)-dil)-Tau_yy(j)
                 tau33 = visc*vis(i,j,k)*(C_2_R*dwdz(i,j,k)-dil)-Tau_zz(j)
                 tau12 = visc*vis(i,j,k)*(dudy(i,j,k)+dvdx(i,j,k))-Tau_xy(j)
                 tau13 = visc*vis(i,j,k)*(dudz(i,j,k)+dwdx(i,j,k))-Tau_xz(j)
                 tau23 = visc*vis(i,j,k)*(dvdz(i,j,k)+dwdy(i,j,k))-Tau_yz(j)
                 
                 up      = u(i,j,k)-fU(j)
                 vp      = v(i,j,k)-fV(j)
                 wp      = w(i,j,k)-fW(j)
                 
                 wrk3d(i,1,k) = rho(i,j,k)*(up**2)*vp  -C_2_R*up*tau12          !T_112
                 wrk3d(i,2,k) = rho(i,j,k)*(vp**3)     -C_2_R*vp*tau22          !T_222
                 wrk3d(i,3,k) = rho(i,j,k)*(wp**2)*vp  -C_2_R*wp*tau23          !T_332
                 wrk3d(i,4,k) = rho(i,j,k)*(vp**2)*up  -      up*tau22-vp*tau12 !T_122
                 wrk3d(i,5,k) = rho(i,j,k)*up*vp*wp    -      up*tau23-wp*tau12 !T_132 new
                 wrk3d(i,6,k) = rho(i,j,k)*(vp**2)*wp  -      vp*tau23-wp*tau22 !T_232 new

! Partition of turbulent transport term for TKE equation
                 wrk3d(i,7,k) = C_05_R*rho(i,j,k)*(up**2+vp**2+wp**2)*vp 
                 wrk3d(i,9,k) =-up*tau12-vp*tau22-wp*tau23
              ENDDO
           ENDDO

        ELSE
           DO k = 1,kmax
              DO i = 1,imax
                 dil = (dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))*c23
                 tau11 = visc*(C_2_R*dudx(i,j,k)-dil)-Tau_xx(j)
                 tau22 = visc*(C_2_R*dvdy(i,j,k)-dil)-Tau_yy(j)
                 tau33 = visc*(C_2_R*dwdz(i,j,k)-dil)-Tau_zz(j)
                 tau12 = visc*(dudy(i,j,k)+dvdx(i,j,k))-Tau_xy(j)
                 tau13 = visc*(dudz(i,j,k)+dwdx(i,j,k))-Tau_xz(j)
                 tau23 = visc*(dvdz(i,j,k)+dwdy(i,j,k))-Tau_yz(j)
                 up = u(i,j,k)-fU(j)
                 vp = v(i,j,k)-fV(j)
                 wp = w(i,j,k)-fW(j)
                 
                 wrk3d(i,1,k) = rho(i,j,k)*(up**2)*vp -C_2_R*up*tau12            ! Txxy
                 wrk3d(i,2,k) = rho(i,j,k)*(vp**3)    -C_2_R*vp*tau22            ! Tyyy
                 wrk3d(i,3,k) = rho(i,j,k)*(wp**2)*vp -C_2_R*wp*tau23            ! Tzzy
                 wrk3d(i,4,k) = rho(i,j,k)*(vp**2)*up -      up*tau22-vp*tau12   ! Txyy
                 wrk3d(i,5,k) = rho(i,j,k)*up*vp*wp   -      up*tau23-wp*tau12   ! Txzy new
                 wrk3d(i,6,k) = rho(i,j,k)*(vp**2)*wp -      vp*tau23-wp*tau22   ! Tyzy new

                 wrk3d(i,7,k) = C_05_R*rho(i,j,k)*(up**2+vp**2+wp**2)*vp 
                 wrk3d(i,9,k) =-up*tau12-vp*tau22-wp*tau23

              ENDDO
           ENDDO
           
        ENDIF
     ENDIF
     
     Txxy(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx, dz, area)
     Tyyy(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx, dz, area) + Tyyy(j) ! add p' terms
     Tzzy(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx, dz, area)
     Txyy(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx, dz, area) + Txyy(j) ! add p' terms
     Txzy(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx, dz, area)
     Tyzy(j) = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx, dz, area) + Tyzy(j) ! add p' terms


     Ty1(j) = AVG_IK(imax,jmax,kmax, i7, wrk3d, dx, dz, area)
!    Ty2(j) contains only pressure terms
     Ty3(j) = AVG_IK(imax,jmax,kmax, i9, wrk3d, dx, dz, area)

  ENDDO

  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Txxy(1), Txxy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Tyyy(1), Tyyy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Tzzy(1), Tzzy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Txyy(1), Txyy_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Txzy(1), Txzy_y(1), i0,i0, wrk1d,wrk2d,wrk3d) !new
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Tyzy(1), Tyzy_y(1), i0,i0, wrk1d,wrk2d,wrk3d) !new
  

  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Ty1(1), Ty1_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Ty2(1), Ty2_y(1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, Ty3(1), Ty3_y(1), i0,i0, wrk1d,wrk2d,wrk3d)

! ###################################################################
! Final loop
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 11')
#endif

  DO j = 1,jmax
! -------------------------------------------------------------------
! Derivatives Fluctuations. Taylor Microscales
! -------------------------------------------------------------------
! Longitudinal terms and Taylor microscales
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) =  dudx(i,j,k)         **2
        wrk3d(i,2,k) = (dvdy(i,j,k)-rV_y(j))**2
        wrk3d(i,3,k) =  dwdz(i,j,k)         **2
        wrk3d(i,4,k) =  dudx(i,j,k)         **3
        wrk3d(i,5,k) = (dvdy(i,j,k)-rV_y(j))**3
        wrk3d(i,6,k) =  dwdz(i,j,k)         **3
        wrk3d(i,7,k) =  dudx(i,j,k)         **4
        wrk3d(i,8,k) = (dvdy(i,j,k)-rV_y(j))**4
        wrk3d(i,9,k) =  dwdz(i,j,k)         **4
     ENDDO; ENDDO
     var_x  = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
     var_y  = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)
     var_z  = AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)
     skew_ux= AVG_IK(imax, jmax, kmax, i4, wrk3d, dx, dz, area)
     skew_vy= AVG_IK(imax, jmax, kmax, i5, wrk3d, dx, dz, area)
     skew_wz= AVG_IK(imax, jmax, kmax, i6, wrk3d, dx, dz, area)
     flat_ux= AVG_IK(imax, jmax, kmax, i7, wrk3d, dx, dz, area)
     flat_vy= AVG_IK(imax, jmax, kmax, i8, wrk3d, dx, dz, area)
     flat_wz= AVG_IK(imax, jmax, kmax, i9, wrk3d, dx, dz, area)

     IF ( var_x .GT. C_0_R ) THEN
        lambda_x = SQRT(Rxx(j)/var_x)
        skew_ux = skew_ux / var_x**C_1_5_R
        flat_ux = flat_ux / var_x**C_2_R
     ELSE
        lambda_x = C_BIG_R
        skew_ux = C_BIG_R
        flat_ux = C_BIG_R
     ENDIF
     IF ( var_y .GT. C_0_R ) THEN
        lambda_y = SQRT(Ryy(j)/var_y)
        skew_vy = skew_vy / var_y**C_1_5_R
        flat_vy = flat_vy / var_y**C_2_R
     ELSE
        lambda_y = C_BIG_R
        skew_vy = C_BIG_R
        flat_vy = C_BIG_R
     ENDIF
     IF ( var_z .GT. C_0_R ) THEN
        lambda_z = SQRT(Rzz(j)/var_z)
        skew_wz = skew_wz / var_z**C_1_5_R
        flat_wz = flat_wz / var_z**C_2_R
     ELSE
        lambda_z = C_BIG_R
        skew_wz = C_BIG_R
        flat_wz = C_BIG_R
     ENDIF

! Lateral terms U
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = (dudy(i,j,k)-rU_y(j))**2
        wrk3d(i,2,k) =  dudz(i,j,k)         **2
        wrk3d(i,3,k) = (dudy(i,j,k)-rU_y(j))**3
        wrk3d(i,4,k) =  dudz(i,j,k)         **3
        wrk3d(i,5,k) = (dudy(i,j,k)-rU_y(j))**4
        wrk3d(i,6,k) =  dudz(i,j,k)         **4
     ENDDO; ENDDO
     var_uy  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     var_uz  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     skew_uy = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     skew_uz = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     flat_uy = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     flat_uz = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area)
     IF ( var_uy .GT. C_0_R ) THEN
        skew_uy = skew_uy / var_uy**C_1_5_R
        flat_uy = flat_uy / var_uy**C_2_R
     ELSE
        skew_uy = C_BIG_R
        flat_uy = C_BIG_R
     ENDIF
     IF ( var_uz .GT. C_0_R ) THEN
        skew_uz = skew_uz / var_uz**C_1_5_R
        flat_uz = flat_uz / var_uz**C_2_R
     ELSE
        skew_uz = C_BIG_R
        flat_uz = C_BIG_R
     ENDIF

! Lateral terms V
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = dvdx(i,j,k)**2
        wrk3d(i,2,k) = dvdz(i,j,k)**2
        wrk3d(i,3,k) = dvdx(i,j,k)**3
        wrk3d(i,4,k) = dvdz(i,j,k)**3
        wrk3d(i,5,k) = dvdx(i,j,k)**4
        wrk3d(i,6,k) = dvdz(i,j,k)**4
     ENDDO; ENDDO
     var_vx  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     var_vz  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     skew_vx = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     skew_vz = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     flat_vx = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     flat_vz = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area)
     IF ( var_vx .GT. C_0_R ) THEN
        skew_vx = skew_vx / var_vx**C_1_5_R
        flat_vx = flat_vx / var_vx**C_2_R
     ELSE
        skew_vx = C_BIG_R
        flat_vx = C_BIG_R
     ENDIF
     IF ( var_vz .GT. C_0_R ) THEN
        skew_vz = skew_vz / var_vz**C_1_5_R
        flat_vz = flat_vz / var_vz**C_2_R
     ELSE
        skew_vz = C_BIG_R
        flat_vz = C_BIG_R
     ENDIF

! Lateral terms W
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) =  dwdx(i,j,k)         **2
        wrk3d(i,2,k) = (dwdy(i,j,k)-rW_y(j))**2
        wrk3d(i,3,k) =  dwdx(i,j,k)         **3
        wrk3d(i,4,k) = (dwdy(i,j,k)-rW_y(j))**3
        wrk3d(i,5,k) =  dwdx(i,j,k)         **4
        wrk3d(i,6,k) = (dwdy(i,j,k)-rW_y(j))**4
     ENDDO; ENDDO
     var_wx  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     var_wy  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     skew_wx = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     skew_wy = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     flat_wx = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     flat_wy = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area)
     IF ( var_wx .GT. C_0_R ) THEN
        skew_wx = skew_wx / var_wx**C_1_5_R
        flat_wx = flat_wx / var_wx**C_2_R
     ELSE
        skew_wx = C_BIG_R
        flat_wx = C_BIG_R
     ENDIF
     IF ( var_wy .GT. C_0_R ) THEN
        skew_wy = skew_wy / var_wy**C_1_5_R
        flat_wy = flat_wy / var_wy**C_2_R
     ELSE
        skew_wy = C_BIG_R
        flat_wy = C_BIG_R
     ENDIF

! Dilatation fluctuation
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = (dudx(i,j,k)+dvdy(i,j,k)-rV_y(j)+dwdz(i,j,k))**2
     ENDDO; ENDDO
     var_dil = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)

! -------------------------------------------------------------------
! Thermodynamic fluctuations
! -------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN

        DO k = 1,kmax
           DO i = 1,imax
              wrk3d(i,1,k) = (e(i,j,k)-re(j))**2
              wrk3d(i,2,k) = rho(i,j,k)*(e(i,j,k)-fe(j))**2
              wrk3d(i,3,k) = (e(i,j,k)+prefactor*p(i,j,k)/rho(i,j,k)-rh(j))**2
              wrk3d(i,4,k) = rho(i,j,k)*(e(i,j,k)+prefactor*p(i,j,k)/rho(i,j,k)-fh(j))**2
           ENDDO
        ENDDO
        re2 = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
        fe2 = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)/rR(j)
        rh2 = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
        fh2 = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)/rR(j)

! Correlations
        DO k=1, kmax
           DO i=1, imax
              wrk3d(i,1,k) = (rho(i,j,k)-rR(j))*(p(i,j,k)-rP(j))
           ENDDO
        ENDDO        
        rRP = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        
        IF ( rR2(j) .GT. C_0_R .AND. rP2(j) .GT. C_0_R ) THEN
           rRP = rRP/sqrt(rR2(j)*rP2(j))
        ELSE
           rRP = C_2_R
        ENDIF
        
        IF ( rR2(j) .GT. C_0_R .AND. rT2(j) .GT. C_0_R ) THEN
           rRT(j) = rRT(j)/sqrt(rR2(j)*rT2(j))
        ELSE
           rRT(j) = C_2_R
        ENDIF

! Acoustic and entropic density fluctuations
        DO k=1, kmax
           DO i=1, imax
              r_prime = rho(i,j,k) - rR(j)
              p_prime = p(i,j,k)   - rP(j)
              rho_ac = p_prime/c2(j)
              rho_en = r_prime - rho_ac
              
              wrk3d(i,1,k) = rho_ac*rho_ac
              wrk3d(i,2,k) = rho_en*rho_en
           ENDDO
        ENDDO

        rho_ac = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        rho_en = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)

! Turbulent Mach number        
        M_t = SQRT((Rxx(j)+Ryy(j)+Rzz(j))/c2(j))

     ELSE
        re2 = C_0_R
        fe2 = C_0_R
        rh2 = C_0_R
        fh2 = C_0_R
        rRP = C_0_R
        rho_ac = C_0_R
        rho_en = C_0_R
        M_t = C_0_R

     ENDIF

! -------------------------------------------------------------------
! Averaged Density Fluctuations Budget
! -------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        DO k = 1,kmax; DO i = 1,imax
           r_prime = rho(i,j,k)-rR(j)
           dil = (dudx(i,j,k)+dvdy(i,j,k)-rV_y(j)+dwdz(i,j,k))
           u_prime = u(i,j,k) - rU(j)
           v_prime = v(i,j,k) - rV(j)
           w_prime = w(i,j,k) - rW(j)
           
           wrk3d(i,1,k) = dil*r_prime
           wrk3d(i,2,k) = dil*r_prime**2
           wrk3d(i,3,k) = r_prime*u_prime
           wrk3d(i,4,k) = r_prime*v_prime
           wrk3d(i,5,k) = r_prime*w_prime
           wrk3d(i,6,k) = v_prime*r_prime**2
        ENDDO; ENDDO
        
        rey_dil1   = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        rey_dil2   = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)
        rey_flux_x = AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)
        rey_flux_y = AVG_IK(imax, jmax, kmax, i4, wrk3d, dx, dz, area)
        rey_flux_z = AVG_IK(imax, jmax, kmax, i5, wrk3d, dx, dz, area)
        rey_trp    = AVG_IK(imax, jmax, kmax, i6, wrk3d, dx, dz, area)
        
        rey_prod =-C_2_R*(rey_flux_y*rR_y(j)+rR2(j)*rV_y(j))
        rey_conv =-rV(j)*rR2_y(j)
        rey_dil1 = C_2_R*rR(j)*rey_dil1
        
        IF( rR_y(j) .NE. C_0_R ) THEN
           eddy_diff =-rey_flux_y/rR_y(j)
        ELSE
           eddy_diff = C_BIG_R
        ENDIF
     
     ELSE
        IF( rB_y(j) .NE. C_0_R ) THEN
           eddy_diff = C_05_R*Byy(j)/rB_y(j)
        ELSE
           eddy_diff = C_BIG_R
        ENDIF
     
     ENDIF

     dummy =  rU_y(j)**2 + rW_y(j)**2 
     IF ( dummy .NE. C_0_R ) THEN
        eddy_visc = SQRT( (Rxy(j)**2+Ryz(j)**2) / dummy )
        ri_g      = rB_y(j) / dummy
     ELSE
        eddy_visc = C_BIG_R
        ri_g      = C_BIG_R
     ENDIF

     IF ( eddy_diff .NE. C_0_R ) THEN
        eddy_prandtl = eddy_visc/eddy_diff
     ELSE
        eddy_prandtl = C_0_R
     ENDIF

! -------------------------------------------------------------------
! Dissipation Terms
! -------------------------------------------------------------------
     IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        DO k = 1,kmax
           DO i = 1,imax              
              dil = (dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))*c23

              tau11 = visc*vis(i,j,k)*(C_2_R*dudx(i,j,k)-dil)-Tau_xx(j)
              tau22 = visc*vis(i,j,k)*(C_2_R*dvdy(i,j,k)-dil)-Tau_yy(j)
              tau33 = visc*vis(i,j,k)*(C_2_R*dwdz(i,j,k)-dil)-Tau_zz(j)
              tau12 = visc*vis(i,j,k)*(dudy(i,j,k)+dvdx(i,j,k))-Tau_xy(j)
              tau13 = visc*vis(i,j,k)*(dudz(i,j,k)+dwdx(i,j,k))-Tau_xz(j)
              tau23 = visc*vis(i,j,k)*(dvdz(i,j,k)+dwdy(i,j,k))-Tau_yz(j)

              upy = dudy(i,j,k)-rU_y(j)
              vpy = dvdy(i,j,k)-rV_y(j)
              wpy = dwdy(i,j,k)-rW_y(j)
              
              wrk3d(i,1,k) = tau11*dudx(i,j,k) + tau12*upy + tau13*dudz(i,j,k)   ! Exx
              wrk3d(i,2,k) = tau12*dvdx(i,j,k) + tau22*vpy + tau23*dvdz(i,j,k)   ! Eyy
              wrk3d(i,3,k) = tau13*dwdx(i,j,k) + tau23*wpy + tau33*dwdz(i,j,k)   ! Ezz
              wrk3d(i,4,k) = tau11*dvdx(i,j,k) + tau12*vpy + tau13*dvdz(i,j,k) & ! Exy
                   + tau12*dudx(i,j,k) + tau22*upy + tau23*dudz(i,j,k)              
              
              wrk3d(i,6,k) = visc*vis(i,j,k)* & ! Mean Viscous Dissipation phi
                   ( dudx(i,j,k)**2 + dvdy(i,j,k)**2 + dwdz(i,j,k)**2 &
                   + C_05_R*( (dudy(i,j,k)+dvdx(i,j,k))**2 + (dudz(i,j,k)+dwdx(i,j,k))**2 &
                            + (dvdz(i,j,k)+dwdy(i,j,k))**2 ) &
                   - ((dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))**2)/C_3_R )

           ENDDO
        ENDDO
     
     ELSE
        DO k = 1,kmax
           DO i = 1,imax
              dil = (dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))*c23

              tau11 = visc*(C_2_R*dudx(i,j,k)-dil)-Tau_xx(j)
              tau22 = visc*(C_2_R*dvdy(i,j,k)-dil)-Tau_yy(j)
              tau33 = visc*(C_2_R*dwdz(i,j,k)-dil)-Tau_zz(j)
              tau12 = visc*(dudy(i,j,k)+dvdx(i,j,k))-Tau_xy(j)
              tau13 = visc*(dudz(i,j,k)+dwdx(i,j,k))-Tau_xz(j)
              tau23 = visc*(dvdz(i,j,k)+dwdy(i,j,k))-Tau_yz(j)

              upy = dudy(i,j,k)-rU_y(j)
              vpy = dvdy(i,j,k)-rV_y(j)
              wpy = dwdy(i,j,k)-rW_y(j)
              
              wrk3d(i,1,k) = tau11*dudx(i,j,k) + tau12*upy + tau13*dudz(i,j,k)
              wrk3d(i,2,k) = tau12*dvdx(i,j,k) + tau22*vpy + tau23*dvdz(i,j,k)
              wrk3d(i,3,k) = tau13*dwdx(i,j,k) + tau23*wpy + tau33*dwdz(i,j,k)
              wrk3d(i,4,k) = tau11*dvdx(i,j,k) + tau12*vpy + tau13*dvdz(i,j,k) &
                           + tau12*dudx(i,j,k) + tau22*upy + tau23*dudz(i,j,k)              
              wrk3d(i,7,k) = tau13*dudx(i,j,k) + tau23*upy + tau33*dudz(i,j,k) &
                           + tau11*dwdx(i,j,k) + tau12*wpy + tau13*dwdz(i,j,k) 
              wrk3d(i,8,k) = tau13*dvdx(i,j,k) + tau23*vpy + tau33*dvdz(i,j,k) &
                           + tau12*dwdx(i,j,k) + tau22*wpy + tau23*dwdz(i,j,k) 

              wrk3d(i,6,k) = visc*&
                   ( dudx(i,j,k)**2 + dvdy(i,j,k)**2 + dwdz(i,j,k)**2 &
                   + C_05_R*( (dudy(i,j,k)+dvdx(i,j,k))**2 + (dudz(i,j,k)+dwdx(i,j,k))**2 &
                            + (dvdz(i,j,k)+dwdy(i,j,k))**2 ) &
                   - ((dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))**2)/C_3_R )

           ENDDO
        ENDDO

     ENDIF

     Exx = C_2_R*AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)/rR(j)
     Eyy = C_2_R*AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)/rR(j)
     Ezz = C_2_R*AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)/rR(j)
     Exy =       AVG_IK(imax, jmax, kmax, i4, wrk3d, dx, dz, area)/rR(j)
     Exz =       AVG_IK(imax, jmax, kmax, i7, wrk3d, dx, dz, area)/rR(j) ! new
     Eyz =       AVG_IK(imax, jmax, kmax, i8, wrk3d, dx, dz, area)/rR(j) ! new

     Phi = C_2_R*AVG_IK(imax, jmax, kmax, i6, wrk3d, dx, dz, area)

! -------------------------------------------------------------------
! Convective Terms 
! -------------------------------------------------------------------
     Cxx =-fV(j)*Rxx_y(j)
     Cyy =-fV(j)*Ryy_y(j)
     Czz =-fV(j)*Rzz_y(j)
     Cxy =-fV(j)*Rxy_y(j)
     Cxz =-fV(j)*Rxz_y(j) ! new
     Cyz =-fV(j)*Ryz_y(j) ! new 

! -------------------------------------------------------------------
! Production Terms
! -------------------------------------------------------------------
     Pxx =-C_2_R*Rxy(j)*fU_y(j)
     Pyy =-C_2_R*Ryy(j)*fV_y(j)
     Pzz =-C_2_R*Ryz(j)*fW_y(j)                ! was set to zero before 
     Pxy =-( Rxy(j)*fV_y(j) + Ryy(j)*fU_y(j) ) 
     Pxz =-( Rxy(j)*fW_y(j) + Ryz(j)*fU_y(j) ) ! new
     Pyz =-( Ryy(j)*fW_y(j) + Ryz(j)*fV_y(j) ) ! new

! -------------------------------------------------------------------
! Pressure Variable-Density  Terms
! -------------------------------------------------------------------
     Gxx = C_0_R
     Gyy = C_2_R*rVf(j)*rP_y(j)
     Gzz = C_0_R
     Gxy =       rUf(j)*rP_y(j)
     Gxz = C_0_R
     Gyz =       rWf(j)*rP_y(j)

! -------------------------------------------------------------------
! Viscous Variable-Density  Terms
! -------------------------------------------------------------------
     Dxx = C_2_R*rUf(j)*Tau_xy_y(j)
     Dyy = C_2_R*rVf(j)*Tau_yy_y(j)
     Dzz = C_2_R*rWf(j)*Tau_yz_y(j)
     Dxy = rUf(j)*Tau_yy_y(j) + rVf(j)*Tau_xy_y(j)
     Dxz = rUf(j)*Tau_yz_y(j) + rWf(j)*Tau_xy_y(j) ! new
     Dyz = rVf(j)*Tau_yz_y(j) + rWf(j)*Tau_yy_y(j) ! new

! -------------------------------------------------------------------
! Coriolis Terms 
! -------------------------------------------------------------------
     IF ( icoriolis_y .NE. EQNS_NONE ) THEN ! contribution from angular velocity Oy
        dummy = rotn_vector(2)
        Fxx(j) = dummy *C_2_R * Rxz(j)
        Fyy(j) =        C_0_R
        Fzz(j) =-dummy *C_2_R * Rxz(j)
        Fxy(j) = dummy        * Ryz(j)
        Fxz(j) = dummy        *(Rzz(j)-Rxx(j))
        Fyz(j) =-dummy        * Rxy(j)
     ENDIF

! -------------------------------------------------------------------
! Buoyancy Terms 
! -------------------------------------------------------------------
! Calculated in Section Potential Energy
     
! -------------------------------------------------------------------
! Transient terms
! -------------------------------------------------------------------
     Rxx_t = -Fxx(j) + Bxx(j) + Cxx + Pxx - Exx + ( PIxx(j) - Txxy_y(j) - Gxx + Dxx ) /rR(j)
     Ryy_t = -Fyy(j) + Byy(j) + Cyy + Pyy - Eyy + ( PIyy(j) - Tyyy_y(j) - Gyy + Dyy ) /rR(j)
     Rzz_t = -Fzz(j) + Bzz(j) + Czz + Pzz - Ezz + ( PIzz(j) - Tzzy_y(j) - Gzz + Dzz ) /rR(j)
     Rxy_t = -Fxy(j) + Bxy(j) + Cxy + Pxy - Exy + ( PIxy(j) - Txyy_y(j) - Gxy + Dxy ) /rR(j)
     Rxz_t = -Fxz(j) + Bxz(j) + Cxz + Pxz - Exz + ( PIxz(j) - Txzy_y(j) - Gxz + Dxz ) /rR(j)
     Ryz_t = -Fyz(j) + Byz(j) + Cyz + Pyz - Eyz + ( PIyz(j) - Tyzy_y(j) - Gyz + Dyz ) /rR(j)

! -------------------------------------------------------------------
! Kinetic energy equation
! -------------------------------------------------------------------
     Tke  = C_05_R*(Rxx(j)    + Ryy(j)    + Rzz(j)   )

     Buo  = C_05_R*(Bxx(j)    + Byy(j)    + Bzz(j)   )
     Con  = C_05_R*(Cxx       + Cyy       + Czz      )
     Prd  = C_05_R*(Pxx       + Pyy       + Pzz      )
     Pi   = C_05_R*(PIxx(j)   + PIyy(j)   + PIzz(j)  )
     Eps  = C_05_R*(Exx       + Eyy       + Ezz      )
     Ty_y = C_05_R*(Txxy_y(j) + Tyyy_y(j) + Tzzy_y(j))
     Gkin = C_05_R*(Gxx       + Gyy       + Gzz      )
     Dkin = C_05_R*(Dxx       + Dyy       + Dzz      )

     Tke_t= Buo + Con + Prd - Eps + ( - Ty_y + Pi - Gkin + Dkin ) / rR(j)

! -------------------------------------------------------------------
! Potential energy equation
! -------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        Pot       = -rB(j)*(y(j) - y(1) - scaley*ycoor_i(inb_scal))
        SourcePot =-rSb(j)*(y(j) - y(1) - scaley*ycoor_i(inb_scal))

     ELSE
        Pot       =-rR(j)*(y(j) - y(1) - scaley*ycoor_rho)*body_vector(2)
        SourcePot = C_0_R
        
     ENDIF

     IF ( Prd .NE. C_0_R ) THEN
        ri_f =-Buo / Prd ! BuoyancyDestruction / ShearProduction
     ELSE
        ri_f = C_BIG_R
     ENDIF

! -------------------------------------------------------------------
! Kolmogorov microscale and Taylor Reynolds number
! -------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
        eta = visc

     ELSE
        IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN; wrk3d(:,1,:) = visc*vis(:,j,:)/rho(:,j,:)
        ELSE;                          wrk3d(:,1,:) = visc/rho(:,j,:)            ;ENDIF
        eta = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)

     ENDIF

     IF ( eta .GT. C_0_R ) THEN; relambda_x = SQRT(Rxx(j))*lambda_x/eta
     ELSE;                       relambda_x = C_BIG_R; ENDIF

     IF ( eta .GT. C_0_R ) THEN; relambda_y = SQRT(Ryy(j))*lambda_y/eta
     ELSE;                       relambda_y = C_BIG_R; ENDIF

     IF ( eta .GT. C_0_R ) THEN; relambda_z = SQRT(Rzz(j))*lambda_z/eta
     ELSE;                       relambda_z = C_BIG_R; ENDIF

     IF ( eta .GT. C_0_R .AND. Eps .GT. C_0_R ) THEN
        relambda_iso = ((Rxx(j)+Ryy(j)+Rzz(j))/C_3_R)* SQRT(C_15_R/(eta*Eps))
        eta = (eta**3/Eps)**C_025_R
     ELSE
        relambda_iso = C_BIG_R
        eta = C_BIG_R
     ENDIF

! ###################################################################
! Output
! ###################################################################
#ifdef USE_MPI
     IF ( ims_pro .EQ. 0 ) THEN
#endif

        IF      ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
           ivauxpre = 4
           VAUXPRE(1) =  y(j)
           VAUXPRE(2) = (y(j)-scaley*ycoor_u  -y(1))/delta_m
           VAUXPRE(3) = (y(j)-scaley*ycoor_u  -y(1))/delta_w
           VAUXPRE(4) = (y(j)-scaley*ycoor_rho-y(1))/delta_h01

        ELSE IF ( imode_flow .EQ. DNS_FLOW_JET   ) THEN
! Not developed yet; for TkStat compatibility with previous files
           ivauxpre = 4
           VAUXPRE(1) = y(j)
           VAUXPRE(2) = y(j)
           VAUXPRE(3) = y(j)
           VAUXPRE(4) = y(j)

        ENDIF

        IF ( j .EQ. jmax/2 ) THEN
           IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
              ivauxpos = 14
              VAUXPOS(1) = delta_m
              VAUXPOS(2) = delta_m_p
              VAUXPOS(3) = delta_w
              VAUXPOS(4) = delta_hb01
              VAUXPOS(5) = delta_ht01
              VAUXPOS(6) = delta_h01
              VAUXPOS(7) = delta_hb25
              VAUXPOS(8) = delta_ht25
              VAUXPOS(9) = delta_h25
              VAUXPOS(10)= mixing1
              VAUXPOS(11)= mixing2
              VAUXPOS(12)= u_friction
              VAUXPOS(13)= d_friction
              VAUXPOS(14)= a_friction
           ELSE
              ivauxpos = 0
           ENDIF

        ELSE
           ivauxpos = 0
           
        ENDIF

        WRITE(23,1020) 1, j, (VAUXPRE(k),k=1,ivauxpre),&
             rR(j), rU(j), rV(j), rW(j), rP(j), rT(j), re(j), rh(j), rs(j), rB(j),&! 10
             fU(j), fV(j), fW(j), fT(j), fe(j), fh(j), fs(j),&                     !  7
             Tke, Rxx(j), Ryy(j), Rzz(j), Rxy(j), Rxz(j), Ryz(j), &                !  7
             rP2(j), rR2(j), rT2(j), fT2(j),  re2, fe2, rh2, fh2, rs2(j), fs2(j),& ! 10
             vortx(j), vorty(j), vortz(j), vortx2(j), vorty2(j), vortz2(j),&       !  6
             Rxx_t, Bxx(j), Cxx, Pxx, Exx, PIxx(j), Fxx(j), Txxy_y(j), Txxy(j), Gxx, Dxx,&     ! 11
             Ryy_t, Byy(j), Cyy, Pyy, Eyy, PIyy(j), Fyy(j), Tyyy_y(j), Tyyy(j), Gyy, Dyy,&     ! 11
             Rzz_t, Bzz(j), Czz, Pzz, Ezz, PIzz(j), Fzz(j), Tzzy_y(j), Tzzy(j), Gzz, Dzz,&     ! 11
             Rxy_t, Bxy(j), Cxy, Pxy, Exy, PIxy(j), Fxy(j), Txyy_y(j), Txyy(j), Gxy, Dxy,&     ! 11
             Rxz_t, Bxz(j), Cxz, Pxz, Exz, PIxz(j), Fxz(j), Txzy_y(j), Txzy(j), Gxz, Dxz,&     ! 11
             Ryz_t, Byz(j), Cyz, Pyz, Eyz, PIyz(j), Fyz(j), Tyzy_y(j), Tyzy(j), Gyz, Dyz,&     ! 11
             Tke_t, Buo,    Con, Prd, Eps, Pi,              Ty_y,      Ty1(j),Ty2(j),Ty3(j),&  ! 10
             Ty1_y(j), Ty2_y(j), Ty3_y(j), Gkin, Dkin, Phi, ugradp(j), &                       !  7
             eta, lambda_x,lambda_y,lambda_z, relambda_x,relambda_y,relambda_z, relambda_iso,& !  8
             var_dil, var_x,var_uy,var_uz, var_vx,var_y,var_vz, var_wx,var_wy,var_z,&          ! 10
             skew_ux,skew_uy,skew_uz, skew_vx,skew_vy,skew_vz, skew_wx,skew_wy,skew_wz,&       !  9
             flat_ux,flat_uy,flat_uz, flat_vx,flat_vy,flat_vz, flat_wx,flat_wy,flat_wz,&       !  9
             rGamma(j), c2(j), rho_ac, rho_en, T_ac(j), T_en(j), M_t, rRP, rRT(j),&            !  9
             rey_flux_x,rey_flux_y,rey_flux_z, rey_dil1,rey_dil2, rey_trp, rey_prod, rey_conv,&!  8
             Pot, SourcePot, rSb(j),&                                        ! 3
             bfreq_fr(j), bfreq_eq(j), lapse_fr(j), lapse_eq(j),&            ! 4
             potem_fr(j), potem_eq(j), psat(j), &                            ! 3
             rP(j)-C_05_R*(rP(jmax/2)+rP(jmax/2+1)), p_mod(j), ri_f, ri_g,&  ! 4
             eddy_diff, eddy_visc, eddy_prandtl, &                           ! 3
             (VAUXPOS(k),k=1,ivauxpos)
! using the maximum 193+4+14
1020       FORMAT(I5,(1X,I5),4(1X,G_FORMAT_R),193(1X,G_FORMAT_R),14(1X,G_FORMAT_R))

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
END SUBROUTINE AVG_FLOW_TEMPORAL_LAYER

