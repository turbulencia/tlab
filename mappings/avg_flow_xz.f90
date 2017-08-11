#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

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
SUBROUTINE AVG_FLOW_XZ(q,s, dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, mean2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imode_eqns, imode_flow, itransport, inb_scal
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, area
  USE DNS_GLOBAL, ONLY : froude, visc, rossby
  USE DNS_GLOBAL, ONLY : buoyancy, coriolis
  USE DNS_GLOBAL, ONLY : pbg, rbg, sbg, qbg
  USE DNS_GLOBAL, ONLY : bbackground, epbackground, pbackground, rbackground, tbackground
  USE THERMO_GLOBAL, ONLY : imixture, MRATIO, GRATIO
  USE THERMO_GLOBAL, ONLY : THERMO_AI, WGHT_INV, gama0
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

  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)    :: q,s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, wrk3d
  TREAL, DIMENSION(jmax,*),           INTENT(INOUT) :: mean2d, wrk1d
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk2d

  TARGET q, dudz

! -------------------------------------------------------------------
  TINTEGER, PARAMETER :: MAX_VARS_GROUPS = 20
  TINTEGER i,j,k, bcs(2,2)
  TREAL SIMPSON_NU, UPPER_THRESHOLD, LOWER_THRESHOLD
  TREAL delta_m, delta_m_p, delta_w
  TREAL rho_min, rho_max, delta_hb01, delta_ht01, delta_h01, mixing1, mixing2
  TREAL delta_hb25, delta_ht25, delta_h25
  TREAL u_friction, d_friction, a_friction
  TREAL dummy
  TREAL c23, prefactor

  TINTEGER ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nmax

  TREAL VAUXPOS(14)
  TINTEGER ivauxpos
  CHARACTER*32 name, groupname(MAX_VARS_GROUPS)
  CHARACTER*250 line1, varname(MAX_VARS_GROUPS)
  CHARACTER*1300 line2

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u,v,w,p, e,rho, vis

! ###################################################################
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

! Define pointers
  u => q(:,:,:,1)
  v => q(:,:,:,2)
  w => q(:,:,:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     e   => q(:,:,:,4)
     rho => q(:,:,:,5)
     p   => q(:,:,:,6)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)
  ELSE
     p   => dudz
  ENDIF

  c23 = C_2_R/C_3_R

  prefactor = GRATIO*MRATIO ! (gama0-C_1_R)*mach*mach

! Variable definition and memory management
! -----------------------------------------------------------------------
! Independent variables
  ig(1) = 1; ng = 1
  IF      ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
#define VAUXPRE1 mean2d(j,ig(1))
#define VAUXPRE2 mean2d(j,ig(1)+1)
#define VAUXPRE3 mean2d(j,ig(1)+2)
#define VAUXPRE4 mean2d(j,ig(1)+3)
     sg(1) = 4

     varname(1) = 'Y SM SW SR'

  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
#define VAUXPRE1 mean2d(j,ig(1))
     sg(1) = 1

     varname(1) = 'Y'
     
  ENDIF

! -----------------------------------------------------------------------
! Dependent variables
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define rR(j)     mean2d(j,ig(2)  )
#define rU(j)     mean2d(j,ig(2)+1)
#define rV(j)     mean2d(j,ig(2)+2)
#define rW(j)     mean2d(j,ig(2)+3)
#define rP(j)     mean2d(j,ig(2)+4)
#define rT(j)     mean2d(j,ig(2)+5)
#define re(j)     mean2d(j,ig(2)+6)
#define rh(j)     mean2d(j,ig(2)+7)
#define rs(j)     mean2d(j,ig(2)+8)
#define rB(j)     mean2d(j,ig(2)+9)
#define fU(j)     mean2d(j,ig(2)+10)
#define fV(j)     mean2d(j,ig(2)+11)
#define fW(j)     mean2d(j,ig(2)+12)
#define fT(j)     mean2d(j,ig(2)+13)
#define fe(j)     mean2d(j,ig(2)+14)
#define fh(j)     mean2d(j,ig(2)+15)
#define fs(j)     mean2d(j,ig(2)+16)
  sg(ng) = 17
     
  groupname(ng) = 'Mean'
  varname(ng)   = 'rR rU rV rW rP rT re rh rs rB fU fV fW fT fe fh fs'
  
! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Tke(j)    mean2d(j,ig(3)  )
#define Rxx(j)    mean2d(j,ig(3)+1)
#define Ryy(j)    mean2d(j,ig(3)+2)
#define Rzz(j)    mean2d(j,ig(3)+3)
#define Rxy(j)    mean2d(j,ig(3)+4)
#define Rxz(j)    mean2d(j,ig(3)+5)
#define Ryz(j)    mean2d(j,ig(3)+6)
#define rP2(j)    mean2d(j,ig(3)+7)
#define rR2(j)    mean2d(j,ig(3)+8)
#define rT2(j)    mean2d(j,ig(3)+9)
#define fT2(j)    mean2d(j,ig(3)+10)
#define re2(j)    mean2d(j,ig(3)+11)
#define fe2(j)    mean2d(j,ig(3)+12)
#define rh2(j)    mean2d(j,ig(3)+13)
#define fh2(j)    mean2d(j,ig(3)+14)
#define rs2(j)    mean2d(j,ig(3)+15)
#define fs2(j)    mean2d(j,ig(3)+16)
  sg(ng) = 17

  groupname(ng) = 'Fluctuations'
  varname(ng)   = 'Tke Rxx Ryy Rzz Rxy Rxz Ryz rP2 rR2 rT2 fT2 re2 fe2 rh2 fh2 rs2 fs2'
  
! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define vortx(j)  mean2d(j,ig(4)  )
#define vorty(j)  mean2d(j,ig(4)+1)
#define vortz(j)  mean2d(j,ig(4)+2)
#define vortx2(j) mean2d(j,ig(4)+3)
#define vorty2(j) mean2d(j,ig(4)+4)
#define vortz2(j) mean2d(j,ig(4)+5)
  sg(ng) = 6

  groupname(ng) = 'Vorticity'
  varname(ng)   = 'Wx Wy Wz Wx2 Wy2 Wz2'
  
! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rxx_t(j)  mean2d(j,ig(5)  )
#define Bxx(j)    mean2d(j,ig(5)+1)
#define Cxx(j)    mean2d(j,ig(5)+2)
#define Pxx(j)    mean2d(j,ig(5)+3)
#define Exx(j)    mean2d(j,ig(5)+4)
#define PIxx(j)   mean2d(j,ig(5)+5)
#define Fxx(j)    mean2d(j,ig(5)+6)
#define Txxy_y(j) mean2d(j,ig(5)+7)
#define Txxy(j)   mean2d(j,ig(5)+8)
#define Gxx(j)    mean2d(j,ig(5)+9)
#define Dxx(j)    mean2d(j,ig(5)+10)
  sg(ng) = 11

  groupname(ng) = 'RxxBudget'
  varname(ng)   = 'Rxx_t Bxx Cxx Pxx Exx PIxx Fxx Txxy_y Txxy Gxx Dxx'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Ryy_t(j)  mean2d(j,ig(6)  )
#define Byy(j)    mean2d(j,ig(6)+1)
#define Cyy(j)    mean2d(j,ig(6)+2)
#define Pyy(j)    mean2d(j,ig(6)+3)
#define Eyy(j)    mean2d(j,ig(6)+4)
#define PIyy(j)   mean2d(j,ig(6)+5)
#define Fyy(j)    mean2d(j,ig(6)+6)
#define Tyyy_y(j) mean2d(j,ig(6)+7)
#define Tyyy(j)   mean2d(j,ig(6)+8)
#define Gyy(j)    mean2d(j,ig(6)+9)
#define Dyy(j)    mean2d(j,ig(6)+10)
  sg(ng) = 11

  groupname(ng) = 'RyyBudget'
  varname(ng)   = 'Ryy_t Byy Cyy Pyy Eyy PIyy Fyy Tyyy_y Tyyy Gyy Dyy'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rzz_t(j)  mean2d(j,ig(7)  )
#define Bzz(j)    mean2d(j,ig(7)+1)
#define Czz(j)    mean2d(j,ig(7)+2)
#define Pzz(j)    mean2d(j,ig(7)+3)
#define Ezz(j)    mean2d(j,ig(7)+4)
#define PIzz(j)   mean2d(j,ig(7)+5)
#define Fzz(j)    mean2d(j,ig(7)+6)
#define Tzzy_y(j) mean2d(j,ig(7)+7)
#define Tzzy(j)   mean2d(j,ig(7)+8)
#define Gzz(j)    mean2d(j,ig(7)+9)
#define Dzz(j)    mean2d(j,ig(7)+10)
  sg(ng) = 11

  groupname(ng) = 'RzzBudget'
  varname(ng)   = 'Rzz_t Bzz Czz Pzz Ezz PIzz Fzz Tzzy_y Tzzy Gzz Dzz'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rxy_t(j)  mean2d(j,ig(8)  )
#define Bxy(j)    mean2d(j,ig(8)+1)
#define Cxy(j)    mean2d(j,ig(8)+2)
#define Pxy(j)    mean2d(j,ig(8)+3)
#define Exy(j)    mean2d(j,ig(8)+4)
#define PIxy(j)   mean2d(j,ig(8)+5)
#define Fxy(j)    mean2d(j,ig(8)+6)
#define Txyy_y(j) mean2d(j,ig(8)+7)
#define Txyy(j)   mean2d(j,ig(8)+8)
#define Gxy(j)    mean2d(j,ig(8)+9)
#define Dxy(j)    mean2d(j,ig(8)+10)
  sg(ng) = 11

  groupname(ng) = 'RxyBudget'
  varname(ng)   = 'Rxy_t Bxy Cxy Pxy Exy PIxy Fxy Txyy_y Txyy Gxy Dxy'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Rxz_t(j)  mean2d(j,ig(9)  )
#define Bxz(j)    mean2d(j,ig(9)+1)
#define Cxz(j)    mean2d(j,ig(9)+2)
#define Pxz(j)    mean2d(j,ig(9)+3)
#define Exz(j)    mean2d(j,ig(9)+4)
#define PIxz(j)   mean2d(j,ig(9)+5)
#define Fxz(j)    mean2d(j,ig(9)+6)
#define Txzy_y(j) mean2d(j,ig(9)+7)
#define Txzy(j)   mean2d(j,ig(9)+8)
#define Gxz(j)    mean2d(j,ig(9)+9)
#define Dxz(j)    mean2d(j,ig(9)+10)
  sg(ng) = 11

  groupname(ng) = 'RxzBudget'
  varname(ng)   = 'Rxz_t Bxz Cxz Pxz Exz PIxz Fxz Txzy_y Txzy Gxz Dxz'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Ryz_t(j)  mean2d(j,ig(10)  )
#define Byz(j)    mean2d(j,ig(10)+1)
#define Cyz(j)    mean2d(j,ig(10)+2)
#define Pyz(j)    mean2d(j,ig(10)+3)
#define Eyz(j)    mean2d(j,ig(10)+4)
#define PIyz(j)   mean2d(j,ig(10)+5)
#define Fyz(j)    mean2d(j,ig(10)+6)
#define Tyzy_y(j) mean2d(j,ig(10)+7)
#define Tyzy(j)   mean2d(j,ig(10)+8)
#define Gyz(j)    mean2d(j,ig(10)+9)
#define Dyz(j)    mean2d(j,ig(10)+10)
  sg(ng) = 11

  groupname(ng) = 'RyzBudget'
  varname(ng)   = 'Ryz_t Byz Cyz Pyz Eyz PIyz Fyz Tyzy_y Tyzy Gyz Dyz'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Tke_t(j)  mean2d(j,ig(11)  )
#define Buo(j)    mean2d(j,ig(11)+1)
#define Con(j)    mean2d(j,ig(11)+2)
#define Prd(j)    mean2d(j,ig(11)+3)
#define Eps(j)    mean2d(j,ig(11)+4)
#define Pi(j)     mean2d(j,ig(11)+5)
#define Ty_y(j)   mean2d(j,ig(11)+6)
#define Ty1(j)    mean2d(j,ig(11)+7)
#define Ty2(j)    mean2d(j,ig(11)+8)
#define Ty3(j)    mean2d(j,ig(11)+9)
#define Ty1_y(j)  mean2d(j,ig(11)+10)
#define Ty2_y(j)  mean2d(j,ig(11)+11)
#define Ty3_y(j)  mean2d(j,ig(11)+12)
#define Gkin(j)   mean2d(j,ig(11)+13)
#define Dkin(j)   mean2d(j,ig(11)+14)
#define Phi(j)    mean2d(j,ig(11)+15)
#define ugradp(j) mean2d(j,ig(11)+16)
  sg(ng) = 17

  groupname(ng) = 'TkeBudget'
  varname(ng)   = 'Tke_t Buo Con Prd Eps Pi Trp Trp1 Trp2 Trp3 Trp1_y Trp2_y Trp3_y G D Phi UgradP'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define eta(j)    mean2d(j,ig(12)  )
#define lxx(j)    mean2d(j,ig(12)+1)
#define lyy(j)   mean2d(j,ig(12)+2)
#define lzz(j)    mean2d(j,ig(12)+3)
#define re_x(j)   mean2d(j,ig(12)+4)
#define re_y(j)   mean2d(j,ig(12)+5)
#define re_z(j)   mean2d(j,ig(12)+6)
#define re_iso(j) mean2d(j,ig(12)+7)
  sg(ng) = 8

  groupname(ng) = 'Scales'
  varname(ng)   = 'Eta LambdaUx LambdaVy LambdaWz ReLambdaUx ReLambdaVy ReLambdaWz ReLambdaIso'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define var_dil(j) mean2d(j,ig(13)  )
#define var_ux(j)  mean2d(j,ig(13)+1)
#define var_uy(j)  mean2d(j,ig(13)+2)
#define var_uz(j)  mean2d(j,ig(13)+3)
#define var_vx(j)  mean2d(j,ig(13)+4)
#define var_vy(j)  mean2d(j,ig(13)+5)
#define var_vz(j)  mean2d(j,ig(13)+6)
#define var_wx(j)  mean2d(j,ig(13)+7)
#define var_wy(j)  mean2d(j,ig(13)+8)
#define var_wz(j)  mean2d(j,ig(13)+9)
#define skew_ux(j) mean2d(j,ig(13)+10)
#define skew_uy(j) mean2d(j,ig(13)+11)
#define skew_uz(j) mean2d(j,ig(13)+12)
#define skew_vx(j) mean2d(j,ig(13)+13)
#define skew_vy(j) mean2d(j,ig(13)+14)
#define skew_vz(j) mean2d(j,ig(13)+15)
#define skew_wx(j) mean2d(j,ig(13)+16)
#define skew_wy(j) mean2d(j,ig(13)+17)
#define skew_wz(j) mean2d(j,ig(13)+18)
#define flat_ux(j) mean2d(j,ig(13)+19)
#define flat_uy(j) mean2d(j,ig(13)+20)
#define flat_uz(j) mean2d(j,ig(13)+21)
#define flat_vx(j) mean2d(j,ig(13)+22)
#define flat_vy(j) mean2d(j,ig(13)+23)
#define flat_vz(j) mean2d(j,ig(13)+24)
#define flat_wx(j) mean2d(j,ig(13)+25)
#define flat_wy(j) mean2d(j,ig(13)+26)
#define flat_wz(j) mean2d(j,ig(13)+27)
  sg(ng) = 28

  groupname(ng) = 'DerivativeFluctuations'
  varname(ng)   = 'VarDilatation '&
                //'VarUx VarUy VarUz VarVx VarVy VarVz VarWx VarWy VarWz '&
                //'SkewUx SkewUy SkewUz SkewVx SkewVy SkewVz SkewWx SkewWy SkewWz '&
                //'FlatUx FlatUy FlatUz FlatVx FlatVy FlatVz FlatWx FlatWy FlatWz'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define rGamma(j)  mean2d(j,ig(14)  )
#define c2(j)      mean2d(j,ig(14)+1)
#define rho_ac(j)  mean2d(j,ig(14)+2)
#define rho_en(j)  mean2d(j,ig(14)+3)
#define T_ac(j)    mean2d(j,ig(14)+4)
#define T_en(j)    mean2d(j,ig(14)+5)
#define M_t(j)     mean2d(j,ig(14)+6)
#define rRP(j)     mean2d(j,ig(14)+7)
#define rRT(j)     mean2d(j,ig(14)+8)
  sg(ng) = 9

  groupname(ng) = 'Acoustics'
  varname(ng)   = 'gamma C2 Rho_ac Rho_en T_ac T_en M_t rRP rRT'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define rey_flux_x(j) mean2d(j,ig(15)  )
#define rey_flux_y(j) mean2d(j,ig(15)+1)
#define rey_flux_z(j) mean2d(j,ig(15)+2)
#define rey_dil1(j)   mean2d(j,ig(15)+3)
#define rey_dil2(j)   mean2d(j,ig(15)+4)
#define rey_trp(j)    mean2d(j,ig(15)+5)
#define rey_prod(j)   mean2d(j,ig(15)+6)
#define rey_conv(j)   mean2d(j,ig(15)+7)
  sg(ng) = 8

  groupname(ng) = 'RhoBudget'
  varname(ng)   = 'RhoFluxX RhoFluxY RhoFluxZ RhoDil1 RhoDil2 RhoTrp RhoProd RhoConv'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define Pot(j)       mean2d(j,ig(16)  )
#define SourcePot(j) mean2d(j,ig(16)+1)
#define rSb(j)       mean2d(j,ig(16)+2)
#define bfreq_fr(j)  mean2d(j,ig(16)+3)
#define bfreq_eq(j)  mean2d(j,ig(16)+4)
#define lapse_fr(j)  mean2d(j,ig(16)+5)
#define lapse_eq(j)  mean2d(j,ig(16)+6)
#define potem_fr(j)  mean2d(j,ig(16)+7)
#define potem_eq(j)  mean2d(j,ig(16)+8)
#define psat(j)      mean2d(j,ig(16)+9)
#define pref(j)      mean2d(j,ig(16)+10)
#define relhum(j)    mean2d(j,ig(16)+11)
#define ri_f(j)      mean2d(j,ig(16)+12)
#define ri_g(j)      mean2d(j,ig(16)+13)
  sg(ng) = 14

  groupname(ng) = 'Stratification'
  varname(ng)   = 'Pot Source rSb BuoyFreq_fr BuoyFreq_eq LapseRate_fr LapseRate_eq '&
                //'PotTemp_fr PotTemp_eq SaturationPressure rP0 RelativeHumidity Ri_f Ri_g'

! -----------------------------------------------------------------------
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define eddy_diff(j)     mean2d(j,ig(17)  )
#define eddy_visc(j)     mean2d(j,ig(17)+1)
#define eddy_prandtl(j)  mean2d(j,ig(17)+2)
  sg(ng) = 3

  groupname(ng) = 'TurbDiffusivities'
  varname(ng)   = 'EddyDiff EddyVisc TurbPrandtl'
  
! -----------------------------------------------------------------------
! Auxiliary variables depending on y and t; this last group is not written
  ng = ng + 1; ig(ng) = ig(ng-1)+ sg(ng-1)
#define rUf(j)    mean2d(j,ig(18))
#define rVf(j)    mean2d(j,ig(18)+1)
#define rWf(j)    mean2d(j,ig(18)+2)

#define rU_y(j)   mean2d(j,ig(18)+3)
#define rV_y(j)   mean2d(j,ig(18)+4)
#define rW_y(j)   mean2d(j,ig(18)+5)
#define fU_y(j)   mean2d(j,ig(18)+6)
#define fV_y(j)   mean2d(j,ig(18)+7)
#define fW_y(j)   mean2d(j,ig(18)+8)
#define rP_y(j)   mean2d(j,ig(18)+9)
#define rR_y(j)   mean2d(j,ig(18)+10)
#define rT_y(j)   mean2d(j,ig(18)+11)
#define rB_y(j)   mean2d(j,ig(18)+12)

#define Rxx_y(j)  mean2d(j,ig(18)+13)
#define Ryy_y(j)  mean2d(j,ig(18)+14)
#define Rzz_y(j)  mean2d(j,ig(18)+15)
#define Rxy_y(j)  mean2d(j,ig(18)+16)
#define Rxz_y(j)  mean2d(j,ig(18)+17)
#define Ryz_y(j)  mean2d(j,ig(18)+18)
#define rR2_y(j)  mean2d(j,ig(18)+19)

#define Tau_yy(j)   mean2d(j,ig(18)+20)
#define Tau_xy(j)   mean2d(j,ig(18)+21)
#define Tau_yz(j)   mean2d(j,ig(18)+22)
#define Tau_xy_y(j) mean2d(j,ig(18)+23)
#define Tau_yy_y(j) mean2d(j,ig(18)+24)
#define Tau_yz_y(j) mean2d(j,ig(18)+25)
  sg(ng) = 26

#define L_AVGMAX 226

! -----------------------------------------------------------------------
  nmax = ig(ng) +sg(ng) -1
  IF ( MAX_AVG_TEMPORAL .LT. nmax ) THEN
     CALL IO_WRITE_ASCII(efile,'AVERAGES_FLOW_XZ. Not enough space in local arrays.')
     CALL DNS_STOP(LES_ERROR_AVGTMP)
  ENDIF
  mean2d(:,1:nmax) = C_0_R

  ng   = ng -1
  nmax = ig(ng) +sg(ng) -1 ! the last group is not written out

  IF ( L_AVGMAX .LT. nmax ) THEN
     CALL IO_WRITE_ASCII(efile,'AVERAGES_FLOW_XZ. Not enough space in format definition.')
     CALL DNS_STOP(LES_ERROR_AVGTMP)
  ENDIF

! #######################################################################
  WRITE(line1,*) itime; line1 = 'Calculating flow statistics at It'//TRIM(ADJUSTL(line1))//'...'
  CALL IO_WRITE_ASCII(lfile,line1)

! ###################################################################
! Averages (do not overwrite dudz; it cotains p for incompressible case)
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 2')
#endif

  CALL AVG_IK_V(imax,jmax,kmax, jmax, u, g(1)%jac,g(3)%jac, rU(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, v, g(1)%jac,g(3)%jac, rV(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, w, g(1)%jac,g(3)%jac, rW(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, p, g(1)%jac,g(3)%jac, rP(1), wrk1d, area)

  IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
     rR(:) = rbackground(:)

     fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC      ) THEN
     CALL THERMO_ANELASTIC_DENSITY(imax,jmax,kmax, s, epbackground,pbackground, dwdx)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdx, g(1)%jac,g(3)%jac, rR(1),  wrk1d, area)     
     
     fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

  ELSE
     dwdx = rho *u
     dwdy = rho *v
     dwdz = rho *w
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdx, g(1)%jac,g(3)%jac, fU(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdy, g(1)%jac,g(3)%jac, fV(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdz, g(1)%jac,g(3)%jac, fW(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, rho,  g(1)%jac,g(3)%jac, rR(1), wrk1d, area)
     fU(:) = fU(:) /rR(:)
     fV(:) = fV(:) /rR(:)
     fW(:) = fW(:) /rR(:)

  ENDIF

  rUf(:) = rU(:) - fU(:)
  rVf(:) = rV(:) - fV(:)
  rWf(:) = rW(:) - fW(:)  

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fU(1), fU_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fV(1), fV_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), fW(1), fW_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rU(1), rU_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rV(1), rV_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rW(1), rW_y(1), wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rP(1),  rP_y(1),  wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rR(1),  rR_y(1),  wrk3d, wrk2d,wrk3d)
 
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     pref(:) = pbackground(:)
  ELSE
     pref(:) = rP(:)
  ENDIF
  
! #######################################################################
! Main covariances (do not overwrite dudz; it contains p for incompressible case)
! #######################################################################
  DO j = 1,jmax
     dwdx(:,j,:) = u(:,j,:) - fU(j)
     dwdy(:,j,:) = v(:,j,:) - fV(j)
     dwdz(:,j,:) = w(:,j,:) - fW(j)
  ENDDO
  
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     dvdx = dwdx *dwdx
     dvdy = dwdy *dwdy
     dvdz = dwdz *dwdz
  ELSE
     dvdx = dwdx *dwdx *rho
     dvdy = dwdy *dwdy *rho
     dvdz = dwdz *dwdz *rho
  ENDIF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, Rxx(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, Ryy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, Rzz(1), wrk1d, area)
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     Rxx(:) = Rxx(:) /rR(:)
     Ryy(:) = Ryy(:) /rR(:)
     Rzz(:) = Rzz(:) /rR(:)
  ENDIF
  
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     dvdx = dwdx *dwdy
     dvdy = dwdx *dwdz
     dvdz = dwdy *dwdz
  ELSE
     dvdx = dwdx *dwdy *rho
     dvdy = dwdx *dwdz *rho
     dvdz = dwdy *dwdz *rho
  ENDIF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, Rxy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, Rxz(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, Ryz(1), wrk1d, area)
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     Rxy(:) = Rxy(:) /rR(:)
     Rxz(:) = Rxz(:) /rR(:)
     Ryz(:) = Ryz(:) /rR(:)
  ENDIF
  
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Rxx(1), Rxx_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Ryy(1), Ryy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Rzz(1), Rzz_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Rxy(1), Rxy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Rxz(1), Rxz_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Ryz(1), Ryz_y(1), wrk3d, wrk2d,wrk3d)

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
     rR2(:) = C_0_R

  ELSE
     IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL THERMO_ANELASTIC_DENSITY(imax,jmax,kmax, s, epbackground,pbackground, wrk3d)
        DO j = 1,jmax
           wrk3d(:,j,:)= wrk3d(:,j,:)-rR(j)
        ENDDO
     ELSE
        DO j = 1,jmax
           wrk3d(:,j,:)= rho(:,j,:)-rR(j)
        ENDDO
     ENDIF
     dvdx = wrk3d*wrk3d
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx,g(1)%jac,g(3)%jac, rR2(1), wrk1d, area)

     CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rR2(1), rR2_y(1), wrk3d, wrk2d,wrk3d)

! Density Fluctuations Budget
     DO j = 1,jmax
        dvdx(:,j,:) = u(:,j,:) - rU(j)
        dvdy(:,j,:) = v(:,j,:) - rV(j)
        dvdz(:,j,:) = w(:,j,:) - rW(j)
     ENDDO
     dvdx = dvdx *wrk3d
     dvdy = dvdy *wrk3d
     dvdz = dvdz *wrk3d
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, rey_flux_x(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, rey_flux_y(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, rey_flux_z(1), wrk1d, area)
     dvdy = dvdy *wrk3d
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, rey_trp(1),    wrk1d, area)

  ENDIF

! Triple-velocity correlations
  dvdx = dwdx *dwdx *dwdy
  dvdy = dwdy *dwdy *dwdy
  dvdz = dwdz *dwdz *dwdy
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     dvdx = dvdx *rho
     dvdy = dvdy *rho
     dvdz = dvdz *rho
  ENDIF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, Txxy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, Tyyy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, Tzzy(1), wrk1d, area)
  Ty1(:) = ( Txxy(:) + Tyyy(:) + Tzzy(:) )*C_05_R

  dvdx = dwdx *dwdy *dwdy
  dvdy = dwdx *dwdy *dwdz
  dvdz = dwdy *dwdy *dwdz
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     dvdx = dvdx *rho
     dvdy = dvdy *rho
     dvdz = dvdz *rho
  ENDIF
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, Txyy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, Txzy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, Tyzy(1), wrk1d, area)

! Pressure
  DO j = 1,jmax
     dvdz(:,j,:) = p(:,j,:) - rP(j)
  ENDDO
  wrk3d = dvdz *dvdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, rP2(1), wrk1d, area)

! Pressure-velocity correlation in TKE transport terms
  dwdx = dwdx *dvdz
  dwdy = dwdy *dvdz
  dwdz = dwdz *dvdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdx, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Txyy(:) = Txyy(:) + wrk1d(:,2)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdy, g(1)%jac,g(3)%jac, Ty2(1),     wrk1d, area)
  Tyyy(:) = Tyyy(:) + Ty2(:) *C_2_R
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdz, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Tyzy(:) = Tyzy(:) + wrk1d(:,2)

! ###################################################################
! Pressure; array dudz containing p is used only up to this section
!
! dudx = du/dx
! dudy = du/dy
! dudz = p
! dvdx = dv/dx
! dvdy = dv/dy
! dvdz = p_prime
! dwdx =       ; dp/dx
! dwdy =       ; dp/dy
! dwdz = dw/dz ; dp/dz
! ###################################################################
! Pressure convection term
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), p, dwdx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), p, dwdy, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), p, dwdz, wrk3d, wrk2d,wrk3d)
  wrk3d = u *dwdx + v*dwdy + w*dwdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, ugradp(1), wrk1d, area)

! Pressure Strain Terms
! 9 derivatives are here recomputed; ok, this routine is not called that often
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, dudx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, dvdy, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, dwdz, wrk3d, wrk2d,wrk3d)
  dudx = dvdz *dudx ! dvdz contains the pressure fluctuation
  dvdy = dvdz *dvdy ! no need to substract rV_y
  dwdz = dvdz *dwdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dudx, g(1)%jac,g(3)%jac, PIxx(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, PIyy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdz, g(1)%jac,g(3)%jac, PIzz(1), wrk1d, area)
  PIxx(:) = PIxx(:) *C_2_R
  PIyy(:) = PIyy(:) *C_2_R
  PIzz(:) = PIzz(:) *C_2_R

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, dudy, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, dvdx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, dwdz, wrk3d, wrk2d,wrk3d) !dudz not free
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, dwdx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, dudx, wrk3d, wrk2d,wrk3d) !dvdz not free
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, dvdy, wrk3d, wrk2d,wrk3d)
  dudy = dvdz *(dudy +dvdx) ! no need to substract rU_y
  dwdz = dvdz *(dwdz +dwdx)
  dudx = dvdz *(dudx +dvdy) ! no need to substract rW_y
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dudy, g(1)%jac,g(3)%jac, PIxy(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dwdz, g(1)%jac,g(3)%jac, PIxz(1), wrk1d, area)
  CALL AVG_IK_V(imax,jmax,kmax, jmax, dudx, g(1)%jac,g(3)%jac, PIyz(1), wrk1d, area)

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
#define GAMMA_LOC(i,j,k) dudx(i,j,k)
#define T_LOC(i,j,k)     dwdx(i,j,k)
#define S_LOC(i,j,k)     dwdz(i,j,k)

  IF      ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
     rT(:)   = tbackground(:)

  ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC      ) THEN
     CALL THERMO_ANELASTIC_TEMPERATURE(imax,jmax,kmax, s, epbackground, T_LOC(1,1,1))
     CALL AVG_IK_V(imax,jmax,kmax, jmax, T_LOC(1,1,1), g(1)%jac,g(3)%jac, rT(1),  wrk1d, area)

     DO j = 1,jmax
        dvdx(:,j,:) = (T_LOC(:,j,:)-rT(j))**2
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx,         g(1)%jac,g(3)%jac, rT2(1), wrk1d, area)
     
     CALL THERMO_POLYNOMIAL_PSAT(imax,jmax,kmax, T_LOC(1,1,1), dvdz)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz,         g(1)%jac,g(3)%jac, psat(1), wrk1d, area)

     CALL THERMO_ANELASTIC_RELATIVEHUMIDITY(imax,jmax,kmax, s, epbackground,pbackground, T_LOC(1,1,1), wrk3d)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d,        g(1)%jac,g(3)%jac, relhum(1), wrk1d, area)
     
     CALL THERMO_ANELASTIC_THETA  (imax,jmax,kmax, s, epbackground,pbackground, wrk3d)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d,        g(1)%jac,g(3)%jac, potem_fr(1), wrk1d, area)
     CALL THERMO_ANELASTIC_THETA_V(imax,jmax,kmax, s, epbackground,pbackground, wrk3d)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d,        g(1)%jac,g(3)%jac, potem_eq(1), wrk1d, area)
     
     dummy = C_1_R /( pbg%parameters(1) *gama0 )
     bfreq_fr(:) = -rR_y(:) /rbackground(:) -dummy *rR(:) /pbackground(:) 
     bfreq_fr(:) = bfreq_fr(:) *buoyancy%vector(2)

  ELSE

! -------------------------------------------------------------------
! Main fields
! -------------------------------------------------------------------
     CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T_LOC(1,1,1), wrk3d)
     CALL THERMO_GAMMA(imax,jmax,kmax, s, T_LOC(1,1,1), GAMMA_LOC(1,1,1))
     CALL THERMO_ENTROPY(imax,jmax,kmax, s, T_LOC(1,1,1), p, S_LOC(1,1,1))

     CALL AVG_IK_V(imax,jmax,kmax, jmax, T_LOC(1,1,1),    g(1)%jac,g(3)%jac, rT(1),     wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, e,               g(1)%jac,g(3)%jac, re(1),     wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, S_LOC(1,1,1),    g(1)%jac,g(3)%jac, rs(1),     wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, GAMMA_LOC(1,1,1),g(1)%jac,g(3)%jac, rGamma(1), wrk1d, area)

! Means
     dudy = rho *e
     dudz = e + prefactor *p /rho
     dvdx = rho *dudz
     dvdy = rho *dwdz    ! rho *S_LOC
     dvdz = rho *dwdx    ! rho *T_LOC
     wrk3d= dudx *p /rho ! GAMMA_LOC *p /rho = speed of sound
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudy, g(1)%jac,g(3)%jac, fe(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudz, g(1)%jac,g(3)%jac, rh(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, fh(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, fs(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, fT(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d,g(1)%jac,g(3)%jac, c2(1), wrk1d, area)
     fe(:) = fe(:) /rR(:)
     fh(:) = fh(:) /rR(:)
     fs(:) = fs(:) /rR(:)
     fT(:) = fT(:) /rR(:)

     CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rT(1),  rT_y(1),  wrk3d, wrk2d,wrk3d)

! Turbulent Mach number        
     M_t(:) = SQRT((Rxx(:)+Ryy(:)+Rzz(:)) /c2(:))

! Covariances
     DO j = 1,jmax
        dudy(:,j,:) =                     (S_LOC(:,j,:)-rs(j))**2
        dudz(:,j,:) = rho(:,j,:)         *(S_LOC(:,j,:)-fs(j))**2
        dvdx(:,j,:) =                     (T_LOC(:,j,:)-rT(j))**2
        dvdy(:,j,:) = rho(:,j,:)         *(T_LOC(:,j,:)-fT(j))**2
        dvdz(:,j,:) =(rho(:,j,:) -rR(j)) *(T_LOC(:,j,:)-fT(j))
        wrk3d(:,j,:)=(rho(:,j,:) -rR(j)) *(p(:,j,:)-rP(j))
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudy, g(1)%jac,g(3)%jac, rs2(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudz, g(1)%jac,g(3)%jac, fs2(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, rT2(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, fT2(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, rRT(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d,g(1)%jac,g(3)%jac, rRP(1), wrk1d, area)
     fs2(:) = fs2(:) /rR(:)
     fT2(:) = fT2(:) /rR(:)

     DO j = 1,jmax
        IF ( rR2(j) .GT. C_0_R .AND. rP2(j) .GT. C_0_R ) THEN; rRP(j) = rRP(j)/SQRT(rR2(j)*rP2(j))
        ELSE;                                                  rRP(j) = C_2_R; ENDIF
        
        IF ( rR2(j) .GT. C_0_R .AND. rT2(j) .GT. C_0_R ) THEN; rRT(j) = rRT(j)/SQRT(rR2(j)*rT2(j))
        ELSE;                                                  rRT(j) = C_2_R; ENDIF

     ENDDO

     DO j = 1,jmax
        dudy(:,j,:) =             (e(:,j,:)-re(j))**2
        dudz(:,j,:) = rho(:,j,:) *(e(:,j,:)-fe(j))**2
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudy, g(1)%jac,g(3)%jac, re2(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudz, g(1)%jac,g(3)%jac, fe2(1), wrk1d, area)
     fe2(:) = fe2(:) /rR(:)

     wrk3d = e + prefactor *p /rho
     DO j = 1,jmax
        dudy(:,j,:) =             (wrk3d(:,j,:)-rh(j))**2
        dudz(:,j,:) = rho(:,j,:) *(wrk3d(:,j,:)-fh(j))**2
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudy, g(1)%jac,g(3)%jac, rh2(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudz, g(1)%jac,g(3)%jac, fh2(1), wrk1d, area)
     fh2(:) = fh2(:) /rR(:)

! Acoustic and entropic density and temperature fluctuations
     DO j = 1,jmax
        dudy(:,j,:) = p(:,j,:) - rP(j)                 ! pprime
        dudz(:,j,:) = dudy(:,j,:) /c2(j)               ! rho_ac
        dvdx(:,j,:) = rho(:,j,:) - rR(j) - dudz(:,j,:) ! rho_en = rprime - rho_ac
        dvdy(:,j,:) =(dudy(:,j,:) /rP(j) - dudz(:,j,:) /rR(j) ) *fT(j) ! T_ac
        dvdz(:,j,:) = T_LOC(:,j,:) - fT(j) - dvdy(:,j,:)               ! T_en = Tprime - T_ac
     ENDDO
     dudz = dudz *dudz
     dvdx = dvdx *dvdx
     dvdy = dvdy *dvdy
     dvdz = dvdz *dvdz
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudz, g(1)%jac,g(3)%jac, rho_ac(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, rho_en(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, T_ac(1),   wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, T_en(1),   wrk1d, area)

! -------------------------------------------------------------------
! Buoyancy frequency & saturation pressure
! -------------------------------------------------------------------
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), rho, dvdy, i0, i0, wrk1d,wrk2d,wrk3d)
  
     CALL THERMO_POLYNOMIAL_PSAT(imax, jmax, kmax, T_LOC(1,1,1), dvdz)
     CALL THERMO_CP(imax, jmax, kmax, s, GAMMA_LOC(1,1,1), dvdx)
     
     DO j = 1,jmax
        dudy(:,j,:) = dwdy(:,j,:) /p(:,j,:) /GAMMA_LOC(:,j,:) - dvdy(:,j,:) /rho(:,j,:)
        dvdx(:,j,:) = C_1_R /dvdx(:,j,:)
        dvdy(:,j,:) = T_LOC(:,j,:) *( (MRATIO*p(:,j,:))**(C_1_R/GAMMA_LOC(:,j,:)-C_1_R) )
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dudy, g(1)%jac,g(3)%jac, bfreq_fr(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, lapse_fr(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, potem_fr(1), wrk1d, area)
     CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, psat(1),     wrk1d, area)
     bfreq_fr(:) =-bfreq_fr(:) *buoyancy%vector(2)
     lapse_fr(:) =-lapse_fr(:) *buoyancy%vector(2) *prefactor
     
#undef S_LOC

#define L_RATIO   dvdx
#define Q_RATIO   dvdy
#define WMEAN_INV dwdy
#define C_RATIO   dwdz

     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), T_LOC(1,1,1), dudz, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(1,1,1,2),   dudy, wrk3d, wrk2d,wrk3d)

        dummy     = THERMO_AI(1,1,3) -THERMO_AI(1,1,1) +GRATIO *WGHT_INV(1)
        L_RATIO   = THERMO_AI(6,1,1) -THERMO_AI(6,1,3) - dummy *dwdx ! dwdx is T_LOC
        L_RATIO   = L_RATIO /( GRATIO *WGHT_INV(1) *dwdx )
        Q_RATIO   = C_1_R /( MRATIO *p /dvdz -C_1_R )                ! dvdz is psat
        WMEAN_INV =(Q_RATIO+C_1_R) *( C_1_R -s(:,:,:,1) )*WGHT_INV(2)

        wrk3d = ( C_1_R +Q_RATIO *L_RATIO )/ WMEAN_INV /&
               ( dudx /( dudx -C_1_R ) +Q_RATIO *L_RATIO *L_RATIO )  ! dudx is GAMMA_LOC
        CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, lapse_eq(1), wrk1d, area)
        lapse_eq(:) =-lapse_eq(:) *buoyancy%vector(2) *MRATIO

        dummy = WGHT_INV(1) /WGHT_INV(2)
        wrk3d = ( dudz -buoyancy%vector(2) *MRATIO *wrk3d )/dwdx &
               *( C_1_R +dummy *L_RATIO /( C_1_R -s(:,:,:,1) ) )
        wrk3d = wrk3d - WGHT_INV(2) /WMEAN_INV *dudy
        CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, bfreq_eq(1), wrk1d, area)
        bfreq_eq(:) =-bfreq_eq(:) *buoyancy%vector(2)
        
        C_RATIO = THERMO_AI(1,1,2)+s(:,:,:,1)*(THERMO_AI(1,1,3)-THERMO_AI(1,1,2))
        C_RATIO = (C_1_R-s(:,:,:,1))*GRATIO*WGHT_INV(2)/C_RATIO
        wrk3d = dwdx /( (MRATIO*p)**C_RATIO ) *EXP( Q_RATIO *C_RATIO *L_RATIO )
        wrk3d = wrk3d *( C_1_R +Q_RATIO )**C_RATIO /( (MRATIO *p /dvdz)**(Q_RATIO*C_RATIO) )
        CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, potem_eq(1), wrk1d, area)
        
     ENDIF
     
#undef L_RATIO
#undef Q_RATIO
#undef WMEAN_INV
#undef C_RATIO

  ENDIF

#undef GAMMA_LOC
#undef T_LOC

! ###################################################################
! Potential energy
!
! dudx = buoyancy
!
! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN

     IF ( buoyancy%type .NE. EQNS_NONE ) THEN
        IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
           CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, dudx)
        ELSE
           CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, dudx, bbackground)
        ENDIF
        
        CALL AVG_IK_V(imax,jmax,kmax, jmax, dudx, g(1)%jac,g(3)%jac, rB(1), wrk1d, area)
        DO j = 1,jmax
           dvdx(:,j,:) = (u(:,j,:)-rU(j))*(dudx(:,j,:)-rB(j))
           dvdy(:,j,:) = (v(:,j,:)-rV(j))*(dudx(:,j,:)-rB(j))
           dvdz(:,j,:) = (w(:,j,:)-rW(j))*(dudx(:,j,:)-rB(j))
        ENDDO
        CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdx, g(1)%jac,g(3)%jac, Bxx(1), wrk1d, area)
        CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdy, g(1)%jac,g(3)%jac, Byy(1), wrk1d, area)
        CALL AVG_IK_V(imax,jmax,kmax, jmax, dvdz, g(1)%jac,g(3)%jac, Bzz(1), wrk1d, area)
        Bxy(:) = Bxx(:) *buoyancy%vector(2) + Byy(:) *buoyancy%vector(1) ! buoyancy%vector includes the Froude
        Bxz(:) = Bxx(:) *buoyancy%vector(3) + Bzz(:) *buoyancy%vector(1)
        Byz(:) = Byy(:) *buoyancy%vector(3) + Bzz(:) *buoyancy%vector(2)
        
        Bxx(:) = C_2_R *Bxx(:) *buoyancy%vector(1)
        Byy(:) = C_2_R *Byy(:) *buoyancy%vector(2)
        Bzz(:) = C_2_R *Bzz(:) *buoyancy%vector(3)
        
        dummy = C_1_R /froude
        rB(:) = rB(:) *dummy

!        pmod(:) =-rP_y(:) + SIGN(rB(:),buoyancy%vector(2))
        
        CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), rB(1), rB_y(1), wrk3d, wrk2d,wrk3d)

        rSb(:) = C_0_R ! not yet developed
                
     ENDIF
     
  ELSE ! Compressible case is not yet finished
     Bxx(:) =-rR(:)*rUf(:)*buoyancy%vector(1)
     Byy(:) =-rR(:)*rVf(:)*buoyancy%vector(2)
     Bzz(:) =-rR(:)*rWf(:)*buoyancy%vector(3)

!     pmod(:) =-rP_y(:) +buoyancy%vector(2) *rR(:)

     rSb(:) = C_0_R

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

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, dudx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u, dudy, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), u, dudz, wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, dvdx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, dvdy, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), v, dvdz, wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, dwdx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w, dwdy, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, dwdz, wrk3d, wrk2d,wrk3d)

! ###################################################################
! Vorticity
! ###################################################################
  wrk3d = dwdy - dvdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, vortx(1), wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = ( wrk3d(:,j,:) - vortx(j) ) **2
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, vortx2(1), wrk1d, area)

  wrk3d = dudz - dwdx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, vorty(1), wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = ( wrk3d(:,j,:) - vorty(j) ) **2
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, vorty2(1), wrk1d, area)

  wrk3d = dvdx - dudy
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, vortz(1), wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = ( wrk3d(:,j,:) - vortz(j) ) **2
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, vortz2(1), wrk1d, area)

! ###################################################################
! Derivatives Fluctuations. Taylor Microscales
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 11')
#endif

! -------------------------------------------------------------------
! Longitudinal terms and Taylor microscales
  wrk3d = dudx  *dudx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_ux(1),  wrk1d, area)
  wrk3d = wrk3d *dudx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_ux(1), wrk1d, area)
  wrk3d = wrk3d *dudx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_ux(1), wrk1d, area)

  DO j = 1,jmax
     wrk3d(:,j,:) = (dvdy(:,j,:) -rV_y(j)) *(dvdy(:,j,:)-rV_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_vy(1),  wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = wrk3d(:,j,:) *(dvdy(:,j,:)-rV_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_vy(1), wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = wrk3d(:,j,:) *(dvdy(:,j,:)-rV_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_vy(1), wrk1d, area)

  wrk3d = dwdz  *dwdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_wz(1),  wrk1d, area)
  wrk3d = wrk3d *dwdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_wz(1), wrk1d, area)
  wrk3d = wrk3d *dwdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_wz(1), wrk1d, area)

  DO j = 1,jmax
     IF ( var_ux(j) .GT. C_0_R ) THEN
        lxx(j) = SQRT(Rxx(j)/var_ux(j))
        skew_ux(j) = skew_ux(j) / var_ux(j)**C_1_5_R
        flat_ux(j) = flat_ux(j) / var_ux(j)**C_2_R
     ELSE
        lxx(j) = C_BIG_R
        skew_ux(j) = C_BIG_R
        flat_ux(j) = C_BIG_R
     ENDIF
     IF ( var_vy(j) .GT. C_0_R ) THEN
        lyy(j) = SQRT(Ryy(j)/var_vy(j))
        skew_vy(j) = skew_vy(j) / var_vy(j)**C_1_5_R
        flat_vy(j) = flat_vy(j) / var_vy(j)**C_2_R
     ELSE
        lyy(j) = C_BIG_R
        skew_vy(j) = C_BIG_R
        flat_vy(j) = C_BIG_R
     ENDIF
     IF ( var_wz(j) .GT. C_0_R ) THEN
        lzz(j) = SQRT(Rzz(j)/var_wz(j))
        skew_wz(j) = skew_wz(j) / var_wz(j)**C_1_5_R
        flat_wz(j) = flat_wz(j) / var_wz(j)**C_2_R
     ELSE
        lzz(j) = C_BIG_R
        skew_wz(j) = C_BIG_R
        flat_wz(j) = C_BIG_R
     ENDIF

  ENDDO

! -------------------------------------------------------------------
! Lateral terms U
  DO j = 1,jmax
     wrk3d(:,j,:) = (dudy(:,j,:) -rU_y(j)) *(dudy(:,j,:)-rU_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_uy(1),  wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = wrk3d(:,j,:) *(dudy(:,j,:)-rU_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_uy(1), wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = wrk3d(:,j,:) *(dudy(:,j,:)-rU_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_uy(1), wrk1d, area)

  wrk3d = dudz  *dudz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_uz(1),  wrk1d, area)
  wrk3d = wrk3d *dudz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_uz(1), wrk1d, area)
  wrk3d = wrk3d *dudz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_uz(1), wrk1d, area)

  DO j = 1,jmax
     IF ( var_uy(j) .GT. C_0_R ) THEN
        skew_uy(j) = skew_uy(j) / var_uy(j)**C_1_5_R
        flat_uy(j) = flat_uy(j) / var_uy(j)**C_2_R
     ELSE
        skew_uy(j) = C_BIG_R
        flat_uy(j) = C_BIG_R
     ENDIF
     IF ( var_uz(j) .GT. C_0_R ) THEN
        skew_uz(j) = skew_uz(j) / var_uz(j)**C_1_5_R
        flat_uz(j) = flat_uz(j) / var_uz(j)**C_2_R
     ELSE
        skew_uz(j) = C_BIG_R
        flat_uz(j) = C_BIG_R
     ENDIF
  ENDDO

! -------------------------------------------------------------------
! Lateral terms V
  wrk3d = dvdx  *dvdx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_vx(1),  wrk1d, area)
  wrk3d = wrk3d *dvdx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_vx(1), wrk1d, area)
  wrk3d = wrk3d *dvdx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_vx(1), wrk1d, area)

  wrk3d = dvdz  *dvdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_vz(1),  wrk1d, area)
  wrk3d = wrk3d *dvdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_vz(1), wrk1d, area)
  wrk3d = wrk3d *dvdz
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_vz(1), wrk1d, area)

  DO j = 1,jmax
     IF ( var_vx(j) .GT. C_0_R ) THEN
        skew_vx(j) = skew_vx(j) / var_vx(j)**C_1_5_R
        flat_vx(j) = flat_vx(j) / var_vx(j)**C_2_R
     ELSE
        skew_vx(j) = C_BIG_R
        flat_vx(j) = C_BIG_R
     ENDIF
     IF ( var_vz(j) .GT. C_0_R ) THEN
        skew_vz(j) = skew_vz(j) / var_vz(j)**C_1_5_R
        flat_vz(j) = flat_vz(j) / var_vz(j)**C_2_R
     ELSE
        skew_vz(j) = C_BIG_R
        flat_vz(j) = C_BIG_R
     ENDIF

  ENDDO

! -------------------------------------------------------------------
! Lateral terms W
  wrk3d = dwdx  *dwdx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_wx(1),  wrk1d, area)
  wrk3d = wrk3d *dwdx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_wx(1), wrk1d, area)
  wrk3d = wrk3d *dwdx
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_wx(1), wrk1d, area)

  DO j = 1,jmax
     wrk3d(:,j,:) = (dwdy(:,j,:) -rW_y(j)) *(dwdy(:,j,:)-rW_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_wy(1),  wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = wrk3d(:,j,:) *(dwdy(:,j,:)-rW_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, skew_wy(1), wrk1d, area)
  DO j = 1,jmax
     wrk3d(:,j,:) = wrk3d(:,j,:) *(dwdy(:,j,:)-rW_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, flat_wy(1), wrk1d, area)

  DO j = 1,jmax
     IF ( var_wx(j) .GT. C_0_R ) THEN
        skew_wx(j) = skew_wx(j) / var_wx(j)**C_1_5_R
        flat_wx(j) = flat_wx(j) / var_wx(j)**C_2_R
     ELSE
        skew_wx(j) = C_BIG_R
        flat_wx(j) = C_BIG_R
     ENDIF
     IF ( var_wy(j) .GT. C_0_R ) THEN
        skew_wy(j) = skew_wy(j) / var_wy(j)**C_1_5_R
        flat_wy(j) = flat_wy(j) / var_wy(j)**C_2_R
     ELSE
        skew_wy(j) = C_BIG_R
        flat_wy(j) = C_BIG_R
     ENDIF

  ENDDO

! -------------------------------------------------------------------
! Dilatation fluctuation
  wrk3d = dudx + dvdy + dwdz
  DO j = 1,jmax
     wrk3d(:,j,:) = (wrk3d(:,j,:) -rV_y(j)) *(wrk3d(:,j,:) -rV_y(j))
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, var_dil(1), wrk1d, area)

! ###################################################################
! Density Fluctuations Budget
! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     wrk3d = dudx + dvdy + dwdz
     DO j = 1,jmax
        wrk3d(:,j,:) = ( wrk3d(:,j,:) -rV_y(j) ) *( rho(:,j,:)-rR(j) )
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, rey_dil1(1), wrk1d, area)

     DO j = 1,jmax
        wrk3d(:,j,:) = wrk3d(:,j,:) *( rho(:,j,:)-rR(j) )
     ENDDO
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, rey_dil2(1), wrk1d, area)

  ENDIF

! ##################################################################
! Mean viscous dissipation rate
! ##################################################################
  wrk3d = dudx **2 + dvdy **2 + dwdz **2 + C_05_R*( (dudy+dvdx)**2 + (dudz+dwdx)**2 + (dvdz+dwdy)**2 )&
        - (dudx+dvdy+dwdz)**2 /C_3_R

  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Phi(1), wrk1d, area)
  Phi(:) = Phi(:) *visc *C_2_R

! ###################################################################
! Dissipation Terms; final operation after viscous terms below
! ###################################################################
  wrk3d = ( dudx +dvdy +dwdz ) *c23
  wrk3d = ( dudx *C_2_R -wrk3d ) *dudx + ( dudy +dvdx ) *dudy + ( dudz +dwdx ) *dudz
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Exx(1), wrk1d, area)

  wrk3d = ( dudx +dvdy +dwdz ) *c23
  wrk3d = ( dvdy *C_2_R -wrk3d ) *dvdy + ( dudy +dvdx ) *dvdx + ( dvdz +dwdy ) *dvdz
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Eyy(1), wrk1d, area)

  wrk3d = ( dudx +dvdy +dwdz ) *c23
  wrk3d = ( dwdz *C_2_R -wrk3d ) *dwdz + ( dwdy +dvdz ) *dwdy + ( dwdx +dudz ) *dwdx
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Ezz(1), wrk1d, area)

  wrk3d = ( dudx +dvdy +dwdz ) *c23
  wrk3d = ( dudx *C_2_R -wrk3d ) *dvdx + ( dudy +dvdx ) *dvdy + ( dudz +dwdx ) *dvdz &
        + ( dvdy *C_2_R -wrk3d ) *dudy + ( dudy +dvdx ) *dudx + ( dvdz +dwdy ) *dudz
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Exy(1), wrk1d, area)

  wrk3d = ( dudx +dvdy +dwdz ) *c23
  wrk3d = ( dudx *C_2_R -wrk3d ) *dwdx + ( dudy +dvdx ) *dwdy + ( dudz +dwdx ) *dwdz &
        + ( dwdz *C_2_R -wrk3d ) *dudz + ( dudz +dwdx ) *dudx + ( dvdz +dwdy ) *dudy
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Exz(1), wrk1d, area)

  wrk3d = ( dudx +dvdy +dwdz ) *c23
  wrk3d = ( dvdy *C_2_R -wrk3d ) *dwdy + ( dudy +dvdx ) *dwdx + ( dvdz +dwdy ) *dwdz &
        + ( dwdz *C_2_R -wrk3d ) *dvdz + ( dudz +dwdx ) *dvdx + ( dvdz +dwdy ) *dvdy
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Eyz(1), wrk1d, area)
  
! ##################################################################
! Viscous shear-stress tensor
! ##################################################################
  wrk3d = dvdy *C_2_R -dudx -dwdz
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Tau_yy(1), wrk1d, area)
  DO j = 1,jmax ! fluctuation tau22'
     dvdy(:,j,:) = ( wrk3d(:,j,:) - Tau_yy(j) ) *c23
  ENDDO
  Tau_yy(:) = Tau_yy(:) *visc *c23

  wrk3d = dudy +dvdx
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Tau_xy(1), wrk1d, area)
  DO j = 1,jmax ! fluctuation tau12'
     dudy(:,j,:) = wrk3d(:,j,:) - Tau_xy(j)
  ENDDO
  Tau_xy(:) = Tau_xy(:) *visc

  wrk3d = dvdz +dwdy
  IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) wrk3d = wrk3d *vis
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Tau_yz(1), wrk1d, area)
  DO j = 1,jmax ! fluctuation tau23'
     dwdy(:,j,:) = wrk3d(:,j,:) - Tau_yz(j)
  ENDDO
  Tau_yz(:) = Tau_yz(:) *visc

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tau_xy(1), Tau_xy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tau_yy(1), Tau_yy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tau_yz(1), Tau_yz_y(1), wrk3d, wrk2d,wrk3d)

! -------------------------------------------------------------------
! Contribution to turbulent transport terms
  DO j = 1,jmax
     wrk3d(:,j,:) = dudy(:,j,:) *( u(:,j,:) -fU(j) ) ! -2*u'*tau12'
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Txxy(:) = Txxy(:) - wrk1d(:,2) *visc *C_2_R
  Ty3(:)  = Ty3(:)  - wrk1d(:,2) *visc 

  DO j = 1,jmax
     wrk3d(:,j,:) = dvdy(:,j,:) *( v(:,j,:) -fV(j) ) ! -2*v'*tau22'
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Tyyy(:) = Tyyy(:) - wrk1d(:,2) *visc *C_2_R
  Ty3(:)  = Ty3(:)  - wrk1d(:,2) *visc       

  DO j = 1,jmax
     wrk3d(:,j,:) = dwdy(:,j,:) *( w(:,j,:) -fW(j) ) ! -2*w'*tau23'
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Tzzy(:) = Tzzy(:) - wrk1d(:,2) *visc *C_2_R
  Ty3(:)  = Ty3(:)  - wrk1d(:,2) *visc 

  DO j = 1,jmax
     wrk3d(:,j,:) = dvdy(:,j,:) *( u(:,j,:) -fU(j) ) + dudy(:,j,:) *( v(:,j,:) -fV(j) )! -u'*tau22' -v'*tau12'
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Txyy(:) = Txyy(:) - wrk1d(:,2) *visc

  DO j = 1,jmax
     wrk3d(:,j,:) = dwdy(:,j,:) *( u(:,j,:) -fU(j) ) + dudy(:,j,:) *( w(:,j,:) -fW(j) )! -u'*tau23' -w'*tau12'
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Txzy(:) = Txzy(:) - wrk1d(:,2) *visc

  DO j = 1,jmax
     wrk3d(:,j,:) = dwdy(:,j,:) *( v(:,j,:) -fV(j) ) + dvdy(:,j,:) *( w(:,j,:) -fW(j) )! -v'*tau23' -w'*tau22'
  ENDDO
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d(1,2), wrk1d, area)
  Tyzy(:) = Tyzy(:) - wrk1d(:,2) *visc

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Txxy(1), Txxy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tyyy(1), Tyyy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tzzy(1), Tzzy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Txyy(1), Txyy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Txzy(1), Txzy_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Tyzy(1), Tyzy_y(1), wrk3d, wrk2d,wrk3d)

  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Ty1(1), Ty1_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Ty2(1), Ty2_y(1), wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), Ty3(1), Ty3_y(1), wrk3d, wrk2d,wrk3d)

! -------------------------------------------------------------------
! Contribution to dissipation
  Exx(:) = ( Exx(:) *visc - Tau_xy(:) *rU_y(:) )*C_2_R
  Eyy(:) = ( Eyy(:) *visc - Tau_yy(:) *rV_y(:) )*C_2_R
  Ezz(:) = ( Ezz(:) *visc - Tau_yz(:) *rW_y(:) )*C_2_R
  Exy(:) =   Exy(:) *visc - Tau_xy(:) *rV_y(:) - Tau_yy(:) *rU_y(:)
  Exz(:) =   Exz(:) *visc - Tau_xy(:) *rW_y(:) - Tau_yz(:) *rU_y(:)
  Eyz(:) =   Eyz(:) *visc - Tau_yy(:) *rW_y(:) - Tau_yz(:) *rV_y(:)
  
! ###################################################################
! Complete budget equations
! ###################################################################
! Density fluctuations budget equation
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
     rey_prod(:) =-C_2_R*(rey_flux_y(:)*rR_y(:)+rR2(:)*rV_y(:))
     rey_conv(:) =-rV(:)*rR2_y(:)
     rey_dil1(:) = C_2_R*rR(:)*rey_dil1(:)
     
     DO j = 1,jmax
        IF( rR_y(j) .NE. C_0_R ) THEN
           eddy_diff(j) =-rey_flux_y(j)/rR_y(j)
        ELSE
           eddy_diff(j) = C_BIG_R
        ENDIF
     ENDDO

  ELSE
     DO j = 1,jmax
        IF( rB_y(j) .NE. C_0_R ) THEN
           eddy_diff(j) = C_05_R*Byy(j)/rB_y(j)
        ELSE
           eddy_diff(j) = C_BIG_R
        ENDIF
     ENDDO
     
  ENDIF
  
  DO j = 1,jmax
     dummy =  rU_y(j)**2 + rW_y(j)**2 
     IF ( dummy .NE. C_0_R ) THEN
        eddy_visc(j) = SQRT( (Rxy(j)**2+Ryz(j)**2) / dummy )
        ri_g(j)      = rB_y(j) / dummy
     ELSE
        eddy_visc(j) = C_BIG_R
        ri_g(j)      = C_BIG_R
     ENDIF
     
     IF ( eddy_diff(j) .NE. C_0_R ) THEN
        eddy_prandtl(j) = eddy_visc(j)/eddy_diff(j)
     ELSE
        eddy_prandtl(j) = C_0_R
     ENDIF
     
  ENDDO

! Rij Convective Terms 
  Cxx(:) =-fV(:)*Rxx_y(:)
  Cyy(:) =-fV(:)*Ryy_y(:)
  Czz(:) =-fV(:)*Rzz_y(:)
  Cxy(:) =-fV(:)*Rxy_y(:)
  Cxz(:) =-fV(:)*Rxz_y(:)
  Cyz(:) =-fV(:)*Ryz_y(:)
     
! Rij Production Terms
  Pxx(:) =-C_2_R*Rxy(:)*fU_y(:)
  Pyy(:) =-C_2_R*Ryy(:)*fV_y(:)
  Pzz(:) =-C_2_R*Ryz(:)*fW_y(:)
  Pxy(:) =-( Rxy(:)*fV_y(:) + Ryy(:)*fU_y(:) ) 
  Pxz(:) =-( Rxy(:)*fW_y(:) + Ryz(:)*fU_y(:) )
  Pyz(:) =-( Ryy(:)*fW_y(:) + Ryz(:)*fV_y(:) )

! Rij Pressure Variable-Density Terms
  Gxx(:) = C_0_R
  Gyy(:) = C_2_R*rVf(:)*rP_y(:)
  Gzz(:) = C_0_R
  Gxy(:) =       rUf(:)*rP_y(:)
  Gxz(:) = C_0_R
  Gyz(:) =       rWf(:)*rP_y(:)

! Rij Viscous Variable-Density Terms
  Dxx(:) = C_2_R*rUf(:)*Tau_xy_y(:)
  Dyy(:) = C_2_R*rVf(:)*Tau_yy_y(:)
  Dzz(:) = C_2_R*rWf(:)*Tau_yz_y(:)
  Dxy(:) = rUf(:)*Tau_yy_y(:) + rVf(:)*Tau_xy_y(:)
  Dxz(:) = rUf(:)*Tau_yz_y(:) + rWf(:)*Tau_xy_y(:)
  Dyz(:) = rVf(:)*Tau_yz_y(:) + rWf(:)*Tau_yy_y(:)

! Rij Coriolis Terms 
  IF ( coriolis%active(2) ) THEN ! contribution from angular velocity Oy
     dummy = coriolis%vector(2)
     Fxx(:) = dummy *C_2_R * Rxz(:)
     Fyy(:) =        C_0_R
     Fzz(:) =-dummy *C_2_R * Rxz(:)
     Fxy(:) = dummy        * Ryz(:)
     Fxz(:) = dummy        *(Rzz(:)-Rxx(:))
     Fyz(:) =-dummy        * Rxy(:)
  ENDIF
     
! Rij Buoyancy Terms; Calculated in Section Potential Energy
    
! Rij Transient terms
  Rxx_t(:) = -Fxx(:) + Bxx(:) + Cxx(:) + Pxx(:) - Exx(:) + ( PIxx(:) - Txxy_y(:) - Gxx(:) + Dxx(:) ) /rR(:)
  Ryy_t(:) = -Fyy(:) + Byy(:) + Cyy(:) + Pyy(:) - Eyy(:) + ( PIyy(:) - Tyyy_y(:) - Gyy(:) + Dyy(:) ) /rR(:)
  Rzz_t(:) = -Fzz(:) + Bzz(:) + Czz(:) + Pzz(:) - Ezz(:) + ( PIzz(:) - Tzzy_y(:) - Gzz(:) + Dzz(:) ) /rR(:)
  Rxy_t(:) = -Fxy(:) + Bxy(:) + Cxy(:) + Pxy(:) - Exy(:) + ( PIxy(:) - Txyy_y(:) - Gxy(:) + Dxy(:) ) /rR(:)
  Rxz_t(:) = -Fxz(:) + Bxz(:) + Cxz(:) + Pxz(:) - Exz(:) + ( PIxz(:) - Txzy_y(:) - Gxz(:) + Dxz(:) ) /rR(:)
  Ryz_t(:) = -Fyz(:) + Byz(:) + Cyz(:) + Pyz(:) - Eyz(:) + ( PIyz(:) - Tyzy_y(:) - Gyz(:) + Dyz(:) ) /rR(:)

! Kinetic energy equation
  Tke(:)  = C_05_R*(Rxx(:)    + Ryy(:)    + Rzz(:)   )
  
  Buo(:)  = C_05_R*(Bxx(:)    + Byy(:)    + Bzz(:)   )
  Con(:)  = C_05_R*(Cxx(:)    + Cyy(:)    + Czz(:)   )
  Prd(:)  = C_05_R*(Pxx(:)    + Pyy(:)    + Pzz(:)   )
  Pi(:)   = C_05_R*(PIxx(:)   + PIyy(:)   + PIzz(:)  )
  Eps(:)  = C_05_R*(Exx(:)    + Eyy(:)    + Ezz(:)   )
  Ty_y(:) = C_05_R*(Txxy_y(:) + Tyyy_y(:) + Tzzy_y(:))
  Gkin(:) = C_05_R*(Gxx(:)    + Gyy(:)    + Gzz(:)   )
  Dkin(:) = C_05_R*(Dxx(:)    + Dyy(:)    + Dzz(:)   )
  
  Tke_t(:)= Buo(:) + Con(:) + Prd(:) - Eps(:) + ( - Ty_y(:) + Pi(:) - Gkin(:) + Dkin(:) ) /rR(:)

! Potential energy equation
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     Pot(:)       = -rB(:)*(g(2)%nodes(:) - g(2)%nodes(1) - g(2)%scale *sbg(inb_scal)%ymean)
     SourcePot(:) =-rSb(:)*(g(2)%nodes(:) - g(2)%nodes(1) - g(2)%scale *sbg(inb_scal)%ymean)
     
  ELSE
     Pot(:)       =-rR(:)*(g(2)%nodes(:) - g(2)%nodes(1) - g(2)%scale*rbg%ymean)*buoyancy%vector(2)
     SourcePot(:) = C_0_R
     
  ENDIF

  DO j = 1,jmax
     IF ( Prd(j) .NE. C_0_R ) THEN
        ri_f(j) =-Buo(j) / Prd(j) ! BuoyancyDestruction / ShearProduction
     ELSE
        ri_f(j) = C_BIG_R
     ENDIF
  ENDDO

! Kolmogorov microscale and Taylor Reynolds number
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
     eta(:) = visc
     
  ELSE
     IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN; wrk3d = visc *vis /rho
     ELSE;                                            wrk3d = visc      /rho ; ENDIF
     CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, eta(1), wrk1d, area)        
  ENDIF
     
  DO j = 1,jmax
     IF ( eta(j) .GT. C_0_R ) THEN; re_x(j) = SQRT(Rxx(j))*lxx(j)/eta(j)
     ELSE;                          re_x(j) = C_BIG_R; ENDIF

     IF ( eta(j) .GT. C_0_R ) THEN; re_y(j) = SQRT(Ryy(j))*lyy(j)/eta(j)
     ELSE;                          re_y(j) = C_BIG_R; ENDIF

     IF ( eta(j) .GT. C_0_R ) THEN; re_z(j) = SQRT(Rzz(j))*lzz(j)/eta(j)
     ELSE;                          re_z(j) = C_BIG_R; ENDIF

     IF ( eta(j) .GT. C_0_R .AND. Eps(j) .GT. C_0_R ) THEN
        re_iso(j) = ((Rxx(j)+Ryy(j)+Rzz(j))/C_3_R)* SQRT(C_15_R/(eta(j)*Eps(j)))
        eta(j)    = (eta(j)**3/Eps(j))**C_025_R
     ELSE
        re_iso(j) = C_BIG_R
        eta(j)    = C_BIG_R
     ENDIF
  ENDDO

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
     IF ( ABS(qbg(1)%delta) .GT. C_SMALL_R ) THEN
        delta_w = qbg(1)%delta/MAXVAL(ABS(fU_y(1:jmax)))

        DO j=1, jmax
           wrk1d(j,1) = rR(j)*(C_025_R-(fU(j)/qbg(1)%delta)**2)
        ENDDO
        delta_m = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)/rbg%mean
           
        DO j=1, jmax
           wrk1d(j,1) = ( Tau_xy(j) -  rR(j)*Rxy(j) )*fU_y(j)
        ENDDO
        delta_m_p = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)*C_2_R/(rbg%mean*qbg(1)%delta**3)

     ELSE
        delta_w   = C_1_R
        delta_m   = C_1_R
        delta_m_p = C_1_R

     ENDIF

! -------------------------------------------------------------------
! Based on delta_rho
! -------------------------------------------------------------------
! 1% and 25% thickness
     IF ( imode_eqns .NE. DNS_EQNS_INCOMPRESSIBLE .AND. imode_eqns .NE. DNS_EQNS_ANELASTIC &
          .AND. ABS(rbg%delta) .GT. C_SMALL_R ) THEN
        dummy = rbg%mean + (C_05_R-C_1EM2_R)*rbg%delta
        delta_hb01 = LOWER_THRESHOLD(jmax, dummy, rR(1), g(2)%nodes)
        dummy = rbg%mean - (C_05_R-C_1EM2_R)*rbg%delta
        delta_ht01 = UPPER_THRESHOLD(jmax, dummy, rR(1), g(2)%nodes)

        delta_hb01 = (g(2)%nodes(1) + g(2)%scale*rbg%ymean) - delta_hb01  
        delta_ht01 = delta_ht01 - (g(2)%nodes(1) + g(2)%scale*rbg%ymean)
        delta_h01  = delta_ht01 + delta_hb01

        dummy = rbg%mean + (C_05_R-C_025_R)*rbg%delta
        delta_hb25 = LOWER_THRESHOLD(jmax, dummy, rR(1), g(2)%nodes)
        dummy = rbg%mean - (C_05_R-C_025_R)*rbg%delta
        delta_ht25 = UPPER_THRESHOLD(jmax, dummy, rR(1), g(2)%nodes)

        delta_hb25 = (g(2)%nodes(1) + g(2)%scale*rbg%ymean) - delta_hb25  
        delta_ht25 = delta_ht25 - (g(2)%nodes(1) + g(2)%scale*rbg%ymean)
        delta_h25  = delta_ht25 + delta_hb25

! Mixing, Youngs definition
        rho_min = rbg%mean - C_05_R*ABS(rbg%delta)
        rho_max = rbg%mean + C_05_R*ABS(rbg%delta)
        wrk3d = (rho - rho_min) *(rho_max -rho)
        CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d, wrk1d(1,2), area)
        mixing1 = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)
        DO j = 1,jmax
           wrk1d(j,1)=(rR(j)-rho_min)*(rho_max-rR(j))
        ENDDO
        mixing1 = mixing1/SIMPSON_NU(jmax, wrk1d, g(2)%nodes)

! Mixing, Cook's definition
        rho_min = rbg%mean - C_05_R*ABS(rbg%delta)
        rho_max = rbg%mean + C_05_R*ABS(rbg%delta)
        DO k = 1,kmax
           DO i = 1,imax*jmax
              wrk3d(i,1,k) = MIN(rho(i,1,k)-rho_min,rho_max-rho(i,1,k))
           ENDDO
        ENDDO
        CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, wrk1d, wrk1d(1,2), area)
        mixing2 = SIMPSON_NU(jmax, wrk1d, g(2)%nodes)
        DO j = 1,jmax
           wrk1d(j,1) = MIN(rR(j)-rho_min,rho_max-rR(j))
        ENDDO
        mixing2 = mixing2/SIMPSON_NU(jmax, wrk1d, g(2)%nodes)

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
! Independent variables
     DO j = 1,jmax
        VAUXPRE1 =  g(2)%nodes(j)
        VAUXPRE2 = (g(2)%nodes(j)-g(2)%scale *qbg(1)%ymean -g(2)%nodes(1))/delta_m
        VAUXPRE3 = (g(2)%nodes(j)-g(2)%scale *qbg(1)%ymean -g(2)%nodes(1))/delta_w
        VAUXPRE4 = (g(2)%nodes(j)-g(2)%scale *rbg%ymean    -g(2)%nodes(1))/delta_h01
     ENDDO

! -------------------------------------------------------------------
! Jet
! -------------------------------------------------------------------
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
! not developed yet
     DO j = 1,jmax
        VAUXPRE1 = g(2)%nodes(j)
     ENDDO

  ENDIF

! ###################################################################
! Output
! ###################################################################
#define LOC_UNIT_ID 23
#define LOC_STATUS 'unknown'

  WRITE(name,*) itime; name='avg'//TRIM(ADJUSTL(name))

! -----------------------------------------------------------------------
! TkStat file; header
! -----------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

#ifdef USE_RECLEN
     OPEN(UNIT=LOC_UNIT_ID, RECL=1050, FILE=name, STATUS='unknown') ! this is probably outdated
#else
     OPEN(UNIT=LOC_UNIT_ID, FILE=name, STATUS='unknown')
#endif

! Header
     WRITE(LOC_UNIT_ID, '(A8,E14.7E3)') 'RTIME = ', rtime

! Independent variables     
     line2 = 'I J '//TRIM(ADJUSTL(varname(1)))
     
! Dependent variables depending on y and t
     DO k = 2,ng
        WRITE(LOC_UNIT_ID,1010) 'GROUP = '//TRIM(ADJUSTL(groupname(k)))//' '//TRIM(ADJUSTL(varname(k)))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(varname(k)))
     ENDDO

! Dependent variables dependent on t only
     IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        line1 = 'Delta_m Delta_m_p Delta_w'
        WRITE(LOC_UNIT_ID,1010) 'GROUP = ShearThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

        line1 = 'Delta_hb01 Delta_ht01 Delta_h01 Delta_hb25 Delta_ht25 Delta_h25 '&
             //'mixing_Youngs mixing_Cook'
        WRITE(LOC_UNIT_ID,1010) 'GROUP = MixingThicknesses '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

        line1 = 'FrictionVelocity FrictionThickness FrictionAngle'
        WRITE(LOC_UNIT_ID,1010) 'GROUP = FrictionTerms '//TRIM(ADJUSTL(line1))
        line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))
     ENDIF

     WRITE(LOC_UNIT_ID,1010) TRIM(ADJUSTL(line2))

! Body     
     DO j = 1,jmax            

        ivauxpos = 0        
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
           ENDIF
        ENDIF

        WRITE(LOC_UNIT_ID,1020) 1, j, (mean2d(j,k),k=1,nmax), (VAUXPOS(k),k=1,ivauxpos)

     ENDDO
     
     CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

1010 FORMAT(A)
1020 FORMAT(I5,(1X,I5),L_AVGMAX(1X,G_FORMAT_R),14(1X,G_FORMAT_R))

END SUBROUTINE AVG_FLOW_XZ

