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
SUBROUTINE AVG_FLOW_XZ(y,dx,dy,dz, q,s,&
     dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, mean2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : MAX_AVG_TEMPORAL, MAX_PROF
  USE DNS_GLOBAL, ONLY : imode_eqns, imode_flow, itransport, ibodyforce
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, inb_scal, inb_scal_array, imode_fdm, i1bc,j1bc,k1bc, area, scaley
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
  TINTEGER, PARAMETER :: MAX_VARS_GROUPS = 20
  TINTEGER i,j,k, is
  TREAL AVG_IK, SIMPSON_NU, FLOW_SHEAR_TEMPORAL, UPPER_THRESHOLD, LOWER_THRESHOLD
  TREAL delta_m, delta_m_p, delta_w
  TREAL ycenter
  TREAL rho_min, rho_max, delta_hb01, delta_ht01, delta_h01, mixing1, mixing2
  TREAL delta_hb25, delta_ht25, delta_h25
  TREAL u_friction, d_friction, a_friction
  TREAL dil, dummy
  TREAL tau11, tau22, tau33, tau12, tau23, tau13
  TREAL up, vp, wp, upy, vpy, wpy
  TREAL c23, prefactor
  TREAL r_prime, p_prime, T_prime, u_prime, v_prime, w_prime

  TINTEGER ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nmax

  TREAL VAUXPOS(14), L_RATIO, Q_RATIO, WMEAN_INV, C_RATIO
  TINTEGER ivauxpos
  CHARACTER*32 name, groupname(MAX_VARS_GROUPS)
  CHARACTER*250 line1, varname(MAX_VARS_GROUPS)
  CHARACTER*1300 line2

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u,v,w,p, e,rho, vis

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 1')
#endif

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
#define pmod(j)      mean2d(j,ig(16)+11)
#define ri_f(j)      mean2d(j,ig(16)+12)
#define ri_g(j)      mean2d(j,ig(16)+13)
  sg(ng) = 14

  groupname(ng) = 'Stratification'
  varname(ng)   = 'Pot Source rSb BuoyFreq_fr BuoyFreq_eq LapseRate_fr LapseRate_eq '&
                //'PotTemp_fr PotTemp_eq SaturationPressure rP0 rPmod Ri_f Ri_g'

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

#define Tau_xx(j) mean2d(j,ig(18)+20)
#define Tau_yy(j) mean2d(j,ig(18)+21)
#define Tau_zz(j) mean2d(j,ig(18)+22)
#define Tau_xy(j) mean2d(j,ig(18)+23)
#define Tau_xz(j) mean2d(j,ig(18)+24)
#define Tau_yz(j) mean2d(j,ig(18)+25)

#define Tau_xy_y(j) mean2d(j,ig(18)+26)
#define Tau_yy_y(j) mean2d(j,ig(18)+27)
#define Tau_yz_y(j) mean2d(j,ig(18)+28)
  sg(ng) = 29

#define L_AVGMAX 226
! -----------------------------------------------------------------------
  nmax = ig(ng) +sg(ng) -1
!  print*,nmax
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
              rho_ac(j) = p_prime/c2(j)
              rho_en(j) = r_prime - rho_ac(j)
              T_ac(j) = fT(j)*(p_prime/rP(j)-rho_ac(j)/rR(j))
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
     
     CALL PARTIAL_Y(imode_fdm, i1,jmax,i1, j1bc, dy, rP(1), pmod(1), i0,i0, wrk1d,wrk2d,wrk3d)
     DO j = 1,jmax
        pmod(j) = -pmod(j)+body_vector(2)*rR(j)
     ENDDO
     
#undef GAMMA_LOC
#undef T_LOC
#undef S_LOC

  ENDIF

  DO j = 1,jmax
     pref(j) = rP(j)-C_05_R*(rP(jmax/2)+rP(jmax/2+1))
  ENDDO

! ###################################################################
! Potential energy
!
! dudx = buoyancy
!
! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN

     IF ( ibodyforce .NE. EQNS_NONE ) THEN
! buoyancy field as used in the integration of the equations (as in dns_profiles)
        DO is = 1,inb_scal
           ycenter = y(1) + scaley*ycoor_i(is)
           DO j = 1,jmax
              wrk1d(j,is) = FLOW_SHEAR_TEMPORAL(iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i(1,is), y(j))
           ENDDO
        ENDDO
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
           CALL THERMO_AIRWATER_LINEAR(i1,jmax,i1, wrk1d(1,1), wrk1d(1,inb_scal_array))
        ENDIF
        wrk1d(:,inb_scal_array+1) = C_0_R
        CALL FI_BUOYANCY(ibodyforce, i1,  jmax,i1,   body_param, wrk1d(1,1), wrk1d(1,inb_scal_array+2), wrk1d(1,inb_scal_array+1))
        CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s,          dudx,                      wrk1d(1,inb_scal_array+2))

! buoyancy terms
        DO j = 1,jmax
           rB(j) = AVG_IK(imax,jmax,kmax, j, dudx, dx,dz, area)
           DO k = 1,kmax; DO i = 1,imax
              wrk3d(i,1,k) = (u(i,j,k)-rU(j))*(dudx(i,j,k)-rB(j))
              wrk3d(i,2,k) = (v(i,j,k)-rV(j))*(dudx(i,j,k)-rB(j))
              wrk3d(i,3,k) = (w(i,j,k)-rW(j))*(dudx(i,j,k)-rB(j))
           ENDDO; ENDDO

           dummy = C_1_R /froude
           rB(j) = rB(j) *dummy
           
           Bxx(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
           Byy(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
           Bzz(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
           
           Bxy(j) = Bxx(j)*body_vector(2) + Byy(j)*body_vector(1) ! body_vector includes the Froude
           Bxz(j) = Bxx(j)*body_vector(3) + Bzz(j)*body_vector(1)
           Byz(j) = Byy(j)*body_vector(3) + Bzz(j)*body_vector(2)
           
           Bxx(j) = C_2_R*Bxx(j)*body_vector(1)
           Byy(j) = C_2_R*Byy(j)*body_vector(2)
           Bzz(j) = C_2_R*Bzz(j)*body_vector(3)
           
           rSb(j) = C_0_R ! not yet developed

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
     var_ux(j)  = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
     var_vy(j)  = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)
     var_wz(j)  = AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)
     skew_ux(j)= AVG_IK(imax, jmax, kmax, i4, wrk3d, dx, dz, area)
     skew_vy(j)= AVG_IK(imax, jmax, kmax, i5, wrk3d, dx, dz, area)
     skew_wz(j)= AVG_IK(imax, jmax, kmax, i6, wrk3d, dx, dz, area)
     flat_ux(j)= AVG_IK(imax, jmax, kmax, i7, wrk3d, dx, dz, area)
     flat_vy(j)= AVG_IK(imax, jmax, kmax, i8, wrk3d, dx, dz, area)
     flat_wz(j)= AVG_IK(imax, jmax, kmax, i9, wrk3d, dx, dz, area)

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

! Lateral terms U
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = (dudy(i,j,k)-rU_y(j))**2
        wrk3d(i,2,k) =  dudz(i,j,k)         **2
        wrk3d(i,3,k) = (dudy(i,j,k)-rU_y(j))**3
        wrk3d(i,4,k) =  dudz(i,j,k)         **3
        wrk3d(i,5,k) = (dudy(i,j,k)-rU_y(j))**4
        wrk3d(i,6,k) =  dudz(i,j,k)         **4
     ENDDO; ENDDO
     var_uy(j)  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     var_uz(j)  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     skew_uy(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     skew_uz(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     flat_uy(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     flat_uz(j) = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area)
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

! Lateral terms V
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = dvdx(i,j,k)**2
        wrk3d(i,2,k) = dvdz(i,j,k)**2
        wrk3d(i,3,k) = dvdx(i,j,k)**3
        wrk3d(i,4,k) = dvdz(i,j,k)**3
        wrk3d(i,5,k) = dvdx(i,j,k)**4
        wrk3d(i,6,k) = dvdz(i,j,k)**4
     ENDDO; ENDDO
     var_vx(j)  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     var_vz(j)  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     skew_vx(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     skew_vz(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     flat_vx(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     flat_vz(j) = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area)
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

! Lateral terms W
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) =  dwdx(i,j,k)         **2
        wrk3d(i,2,k) = (dwdy(i,j,k)-rW_y(j))**2
        wrk3d(i,3,k) =  dwdx(i,j,k)         **3
        wrk3d(i,4,k) = (dwdy(i,j,k)-rW_y(j))**3
        wrk3d(i,5,k) =  dwdx(i,j,k)         **4
        wrk3d(i,6,k) = (dwdy(i,j,k)-rW_y(j))**4
     ENDDO; ENDDO
     var_wx(j)  = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
     var_wy(j)  = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)
     skew_wx(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
     skew_wy(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)
     flat_wx(j) = AVG_IK(imax,jmax,kmax, i5, wrk3d, dx,dz, area)
     flat_wy(j) = AVG_IK(imax,jmax,kmax, i6, wrk3d, dx,dz, area)
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

! Dilatation fluctuation
     DO k = 1,kmax; DO i = 1,imax
        wrk3d(i,1,k) = (dudx(i,j,k)+dvdy(i,j,k)-rV_y(j)+dwdz(i,j,k))**2
     ENDDO; ENDDO
     var_dil(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)

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
        re2(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)
        fe2(j) = AVG_IK(imax,jmax,kmax, i2, wrk3d, dx,dz, area)/rR(j)
        rh2(j) = AVG_IK(imax,jmax,kmax, i3, wrk3d, dx,dz, area)
        fh2(j) = AVG_IK(imax,jmax,kmax, i4, wrk3d, dx,dz, area)/rR(j)

! Correlations
        DO k=1, kmax
           DO i=1, imax
              wrk3d(i,1,k) = (rho(i,j,k)-rR(j))*(p(i,j,k)-rP(j))
           ENDDO
        ENDDO        
        rRP(j) = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        
        IF ( rR2(j) .GT. C_0_R .AND. rP2(j) .GT. C_0_R ) THEN
           rRP(j) = rRP(j)/sqrt(rR2(j)*rP2(j))
        ELSE
           rRP(j) = C_2_R
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
              rho_ac(j) = p_prime/c2(j)
              rho_en(j) = r_prime - rho_ac(j)
              
              wrk3d(i,1,k) = rho_ac(j)*rho_ac(j)
              wrk3d(i,2,k) = rho_en(j)*rho_en(j)
           ENDDO
        ENDDO

        rho_ac(j) = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        rho_en(j) = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)

! Turbulent Mach number        
        M_t(j) = SQRT((Rxx(j)+Ryy(j)+Rzz(j))/c2(j))

     ELSE
        re2(j) = C_0_R
        fe2(j) = C_0_R
        rh2(j) = C_0_R
        fh2(j) = C_0_R
        rRP(j) = C_0_R
        rho_ac(j) = C_0_R
        rho_en(j) = C_0_R
        M_t(j) = C_0_R

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
        
        rey_dil1(j)   = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        rey_dil2(j)   = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)
        rey_flux_x(j) = AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)
        rey_flux_y(j) = AVG_IK(imax, jmax, kmax, i4, wrk3d, dx, dz, area)
        rey_flux_z(j) = AVG_IK(imax, jmax, kmax, i5, wrk3d, dx, dz, area)
        rey_trp(j)    = AVG_IK(imax, jmax, kmax, i6, wrk3d, dx, dz, area)
        
        rey_prod(j) =-C_2_R*(rey_flux_y(j)*rR_y(j)+rR2(j)*rV_y(j))
        rey_conv(j) =-rV(j)*rR2_y(j)
        rey_dil1(j) = C_2_R*rR(j)*rey_dil1(j)
        
        IF( rR_y(j) .NE. C_0_R ) THEN
           eddy_diff(j) =-rey_flux_y(j)/rR_y(j)
        ELSE
           eddy_diff(j) = C_BIG_R
        ENDIF
     
     ELSE
        IF( rB_y(j) .NE. C_0_R ) THEN
           eddy_diff(j) = C_05_R*Byy(j)/rB_y(j)
        ELSE
           eddy_diff(j) = C_BIG_R
        ENDIF
     
     ENDIF

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

     Exx(j) = C_2_R*AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)/rR(j)
     Eyy(j) = C_2_R*AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)/rR(j)
     Ezz(j) = C_2_R*AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)/rR(j)
     Exy(j) =       AVG_IK(imax, jmax, kmax, i4, wrk3d, dx, dz, area)/rR(j)
     Exz(j) =       AVG_IK(imax, jmax, kmax, i7, wrk3d, dx, dz, area)/rR(j) ! new
     Eyz(j) =       AVG_IK(imax, jmax, kmax, i8, wrk3d, dx, dz, area)/rR(j) ! new

     Phi(j) = C_2_R*AVG_IK(imax, jmax, kmax, i6, wrk3d, dx, dz, area)

! -------------------------------------------------------------------
! Convective Terms 
! -------------------------------------------------------------------
     Cxx(j) =-fV(j)*Rxx_y(j)
     Cyy(j) =-fV(j)*Ryy_y(j)
     Czz(j) =-fV(j)*Rzz_y(j)
     Cxy(j) =-fV(j)*Rxy_y(j)
     Cxz(j) =-fV(j)*Rxz_y(j) ! new
     Cyz(j) =-fV(j)*Ryz_y(j) ! new 

! -------------------------------------------------------------------
! Production Terms
! -------------------------------------------------------------------
     Pxx(j) =-C_2_R*Rxy(j)*fU_y(j)
     Pyy(j) =-C_2_R*Ryy(j)*fV_y(j)
     Pzz(j) =-C_2_R*Ryz(j)*fW_y(j)                ! was set to zero before 
     Pxy(j) =-( Rxy(j)*fV_y(j) + Ryy(j)*fU_y(j) ) 
     Pxz(j) =-( Rxy(j)*fW_y(j) + Ryz(j)*fU_y(j) ) ! new
     Pyz(j) =-( Ryy(j)*fW_y(j) + Ryz(j)*fV_y(j) ) ! new

! -------------------------------------------------------------------
! Pressure Variable-Density  Terms
! -------------------------------------------------------------------
     Gxx(j) = C_0_R
     Gyy(j) = C_2_R*rVf(j)*rP_y(j)
     Gzz(j) = C_0_R
     Gxy(j) =       rUf(j)*rP_y(j)
     Gxz(j) = C_0_R
     Gyz(j) =       rWf(j)*rP_y(j)

! -------------------------------------------------------------------
! Viscous Variable-Density  Terms
! -------------------------------------------------------------------
     Dxx(j) = C_2_R*rUf(j)*Tau_xy_y(j)
     Dyy(j) = C_2_R*rVf(j)*Tau_yy_y(j)
     Dzz(j) = C_2_R*rWf(j)*Tau_yz_y(j)
     Dxy(j) = rUf(j)*Tau_yy_y(j) + rVf(j)*Tau_xy_y(j)
     Dxz(j) = rUf(j)*Tau_yz_y(j) + rWf(j)*Tau_xy_y(j) ! new
     Dyz(j) = rVf(j)*Tau_yz_y(j) + rWf(j)*Tau_yy_y(j) ! new

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
     Rxx_t(j) = -Fxx(j) + Bxx(j) + Cxx(j) + Pxx(j) - Exx(j) + ( PIxx(j) - Txxy_y(j) - Gxx(j) + Dxx(j) ) /rR(j)
     Ryy_t(j) = -Fyy(j) + Byy(j) + Cyy(j) + Pyy(j) - Eyy(j) + ( PIyy(j) - Tyyy_y(j) - Gyy(j) + Dyy(j) ) /rR(j)
     Rzz_t(j) = -Fzz(j) + Bzz(j) + Czz(j) + Pzz(j) - Ezz(j) + ( PIzz(j) - Tzzy_y(j) - Gzz(j) + Dzz(j) ) /rR(j)
     Rxy_t(j) = -Fxy(j) + Bxy(j) + Cxy(j) + Pxy(j) - Exy(j) + ( PIxy(j) - Txyy_y(j) - Gxy(j) + Dxy(j) ) /rR(j)
     Rxz_t(j) = -Fxz(j) + Bxz(j) + Cxz(j) + Pxz(j) - Exz(j) + ( PIxz(j) - Txzy_y(j) - Gxz(j) + Dxz(j) ) /rR(j)
     Ryz_t(j) = -Fyz(j) + Byz(j) + Cyz(j) + Pyz(j) - Eyz(j) + ( PIyz(j) - Tyzy_y(j) - Gyz(j) + Dyz(j) ) /rR(j)

! -------------------------------------------------------------------
! Kinetic energy equation
! -------------------------------------------------------------------
     Tke(j)  = C_05_R*(Rxx(j)    + Ryy(j)    + Rzz(j)   )

     Buo(j)  = C_05_R*(Bxx(j)    + Byy(j)    + Bzz(j)   )
     Con(j)  = C_05_R*(Cxx(j)       + Cyy(j)       + Czz(j)      )
     Prd(j)  = C_05_R*(Pxx(j)       + Pyy(j)       + Pzz(j)      )
     Pi(j)   = C_05_R*(PIxx(j)   + PIyy(j)   + PIzz(j)  )
     Eps(j)  = C_05_R*(Exx(j)       + Eyy(j)       + Ezz(j)      )
     Ty_y(j) = C_05_R*(Txxy_y(j) + Tyyy_y(j) + Tzzy_y(j))
     Gkin(j) = C_05_R*(Gxx(j)       + Gyy(j)       + Gzz(j)      )
     Dkin(j) = C_05_R*(Dxx(j)       + Dyy(j)       + Dzz(j)      )

     Tke_t(j)= Buo(j) + Con(j) + Prd(j) - Eps(j) + ( - Ty_y(j) + Pi(j) - Gkin(j) + Dkin(j) ) / rR(j)

! -------------------------------------------------------------------
! Potential energy equation
! -------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        Pot(j)       = -rB(j)*(y(j) - y(1) - scaley*ycoor_i(inb_scal))
        SourcePot(j) =-rSb(j)*(y(j) - y(1) - scaley*ycoor_i(inb_scal))

     ELSE
        Pot(j)       =-rR(j)*(y(j) - y(1) - scaley*ycoor_rho)*body_vector(2)
        SourcePot(j) = C_0_R
        
     ENDIF

     IF ( Prd(j) .NE. C_0_R ) THEN
        ri_f(j) =-Buo(j) / Prd(j) ! BuoyancyDestruction / ShearProduction
     ELSE
        ri_f(j) = C_BIG_R
     ENDIF

! -------------------------------------------------------------------
! Kolmogorov microscale and Taylor Reynolds number
! -------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC )THEN
        eta(j) = visc

     ELSE
        IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN; wrk3d(:,1,:) = visc*vis(:,j,:)/rho(:,j,:)
        ELSE;                          wrk3d(:,1,:) = visc/rho(:,j,:)            ;ENDIF
        eta(j) = AVG_IK(imax,jmax,kmax, i1, wrk3d, dx,dz, area)

     ENDIF

     IF ( eta(j) .GT. C_0_R ) THEN; re_x(j) = SQRT(Rxx(j))*lxx(j)/eta(j)
     ELSE;                       re_x(j) = C_BIG_R; ENDIF

     IF ( eta(j) .GT. C_0_R ) THEN; re_y(j) = SQRT(Ryy(j))*lyy(j)/eta(j)
     ELSE;                       re_y(j) = C_BIG_R; ENDIF

     IF ( eta(j) .GT. C_0_R ) THEN; re_z(j) = SQRT(Rzz(j))*lzz(j)/eta(j)
     ELSE;                       re_z(j) = C_BIG_R; ENDIF

     IF ( eta(j) .GT. C_0_R .AND. Eps(j) .GT. C_0_R ) THEN
        re_iso(j) = ((Rxx(j)+Ryy(j)+Rzz(j))/C_3_R)* SQRT(C_15_R/(eta(j)*Eps(j)))
        eta(j) = (eta(j)**3/Eps(j))**C_025_R
     ELSE
        re_iso(j) = C_BIG_R
        eta(j) = C_BIG_R
     ENDIF

! Independent variables
     IF      ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        VAUXPRE1 =  y(j)
        VAUXPRE2 = (y(j)-scaley*ycoor_u  -y(1))/delta_m
        VAUXPRE3 = (y(j)-scaley*ycoor_u  -y(1))/delta_w
        VAUXPRE4 = (y(j)-scaley*ycoor_rho-y(1))/delta_h01
        
     ELSE IF ( imode_flow .EQ. DNS_FLOW_JET   ) THEN
! Not developed yet; for TkStat compatibility with previous files
        VAUXPRE1 = y(j)

     ENDIF
     
  ENDDO

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

