#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# Assumes statistical homogeneity in xOz, so that the corresponding
!# partial derivative terms are assumed to be zero.
!#
!# In the incompressible case, the array p has been
!# pointed to dudz and the pressure field is stored there; do not
!# use array dudz until pressure block
!#
!# Reynolds and Favre averages
!#
!########################################################################

subroutine AVG_FLOW_XZ(q, s, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, mean2d)
    use TLAB_CONSTANTS, only: MAX_AVG_TEMPORAL
    use TLAB_CONSTANTS, only: efile, lfile, wp, wi
    use TLAB_VARS
    use TLAB_PROCS
    use TLAB_ARRAYS, only: wrk1d
    use TLAB_POINTERS_3D, only: p_wrk3d
    use Thermodynamics, only: imixture, CRATIO_INV, RRATIO
    use Thermodynamics, only: rd_ov_rv, Cd, Rv, Cvl, Lvl, Ldl, Rd, PREF_1000
    use Thermodynamics, only: Thermo_Psat_Polynomial
    use THERMO_ANELASTIC
    use THERMO_CALORIC
    use IBM_VARS, only: gamma_0, gamma_1
    use AVGS, only: AVG_IK_V
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
#endif
#ifdef USE_MPI
    use TLAB_MPI_VARS
#endif
    use FI_SOURCES, only: bbackground, FI_BUOYANCY, FI_BUOYANCY_SOURCE
    use OPR_PARTIAL

    implicit none

    real(wp), intent(IN) :: q(imax, jmax, kmax, inb_flow_array)
!    real(wp), intent(IN) :: s(imax, jmax, kmax, inb_scal_array)
    real(wp), intent(INOUT) :: s(imax, jmax, kmax, inb_scal_array) ! caluclates equi composition in airwater
    real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
    real(wp), intent(INOUT) :: mean2d(jmax, MAX_AVG_TEMPORAL)

    target q, dudz

    ! -------------------------------------------------------------------
    integer, parameter :: MAX_VARS_GROUPS = 20
    integer j, bcs(2, 2)
    real(wp) dummy
    real(wp) c23

    integer ig(MAX_VARS_GROUPS), sg(MAX_VARS_GROUPS), ng, nv

    character*32 name, groupname(MAX_VARS_GROUPS)
    character*250 line1, varname(MAX_VARS_GROUPS)

    ! Pointers to existing allocated space
    real(wp), dimension(:, :, :), pointer :: u, v, w, p, e, rho, vis

    ! ###################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    ! Define pointers
    u => q(:, :, :, 1)
    v => q(:, :, :, 2)
    w => q(:, :, :, 3)
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) then
        e => q(:, :, :, 4)
        rho => q(:, :, :, 5)
        p => q(:, :, :, 6)
        if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) vis => q(:, :, :, 8)
    else
        p => dudz
    end if

    c23 = 2.0_wp/3.0_wp

    ! Variable definition and memory management
    ! -----------------------------------------------------------------------
    ng = 1; ig(ng) = 1
#define rR(j)     mean2d(j,ig(1)  )
#define rU(j)     mean2d(j,ig(1)+1)
#define rV(j)     mean2d(j,ig(1)+2)
#define rW(j)     mean2d(j,ig(1)+3)
#define rP(j)     mean2d(j,ig(1)+4)
#define rT(j)     mean2d(j,ig(1)+5)
#define re(j)     mean2d(j,ig(1)+6)
#define rh(j)     mean2d(j,ig(1)+7)
#define rs(j)     mean2d(j,ig(1)+8)
#define rB(j)     mean2d(j,ig(1)+9)
#define fU(j)     mean2d(j,ig(1)+10)
#define fV(j)     mean2d(j,ig(1)+11)
#define fW(j)     mean2d(j,ig(1)+12)
#define fT(j)     mean2d(j,ig(1)+13)
#define fe(j)     mean2d(j,ig(1)+14)
#define fh(j)     mean2d(j,ig(1)+15)
#define fs(j)     mean2d(j,ig(1)+16)
    sg(ng) = 17

    groupname(ng) = 'Mean'
    varname(ng) = 'rR rU rV rW rP rT re rh rs rB fU fV fW fT fe fh fs'

    if (imode_ibm == 1) then
        varname(ng) = trim(adjustl(varname(ng)))//' eps_0 eps_1'
        sg(ng) = sg(ng) + 2
#define ep_0(j)   mean2d(j,ig(1)+17)
#define ep_1(j)   mean2d(j,ig(1)+18)
    end if

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Tke(j)    mean2d(j,ig(2)  )
#define Rxx(j)    mean2d(j,ig(2)+1)
#define Ryy(j)    mean2d(j,ig(2)+2)
#define Rzz(j)    mean2d(j,ig(2)+3)
#define Rxy(j)    mean2d(j,ig(2)+4)
#define Rxz(j)    mean2d(j,ig(2)+5)
#define Ryz(j)    mean2d(j,ig(2)+6)
#define rP2(j)    mean2d(j,ig(2)+7)
#define rR2(j)    mean2d(j,ig(2)+8)
#define rT2(j)    mean2d(j,ig(2)+9)
#define fT2(j)    mean2d(j,ig(2)+10)
#define re2(j)    mean2d(j,ig(2)+11)
#define fe2(j)    mean2d(j,ig(2)+12)
#define rh2(j)    mean2d(j,ig(2)+13)
#define fh2(j)    mean2d(j,ig(2)+14)
#define rs2(j)    mean2d(j,ig(2)+15)
#define fs2(j)    mean2d(j,ig(2)+16)
    sg(ng) = 17

    groupname(ng) = 'Fluctuations'
    varname(ng) = 'Tke Rxx Ryy Rzz Rxy Rxz Ryz rP2 rR2 rT2 fT2 re2 fe2 rh2 fh2 rs2 fs2'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define vortx(j)  mean2d(j,ig(3)  )
#define vorty(j)  mean2d(j,ig(3)+1)
#define vortz(j)  mean2d(j,ig(3)+2)
#define vortx2(j) mean2d(j,ig(3)+3)
#define vorty2(j) mean2d(j,ig(3)+4)
#define vortz2(j) mean2d(j,ig(3)+5)
    sg(ng) = 6

    groupname(ng) = 'Vorticity'
    varname(ng) = 'Wx Wy Wz Wx2 Wy2 Wz2'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxx_t(j)  mean2d(j,ig(4)  )
#define Bxx(j)    mean2d(j,ig(4)+1)
#define Cxx(j)    mean2d(j,ig(4)+2)
#define Pxx(j)    mean2d(j,ig(4)+3)
#define Exx(j)    mean2d(j,ig(4)+4)
#define PIxx(j)   mean2d(j,ig(4)+5)
#define Fxx(j)    mean2d(j,ig(4)+6)
#define Txxy_y(j) mean2d(j,ig(4)+7)
#define Txxy(j)   mean2d(j,ig(4)+8)
#define Gxx(j)    mean2d(j,ig(4)+9)
#define Dxx(j)    mean2d(j,ig(4)+10)
    sg(ng) = 11

    groupname(ng) = 'RxxBudget'
    varname(ng) = 'Rxx_t Bxx Cxx Pxx Exx PIxx Fxx Txxy_y Txxy Gxx Dxx'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Ryy_t(j)  mean2d(j,ig(5)  )
#define Byy(j)    mean2d(j,ig(5)+1)
#define Cyy(j)    mean2d(j,ig(5)+2)
#define Pyy(j)    mean2d(j,ig(5)+3)
#define Eyy(j)    mean2d(j,ig(5)+4)
#define PIyy(j)   mean2d(j,ig(5)+5)
#define Fyy(j)    mean2d(j,ig(5)+6)
#define Tyyy_y(j) mean2d(j,ig(5)+7)
#define Tyyy(j)   mean2d(j,ig(5)+8)
#define Gyy(j)    mean2d(j,ig(5)+9)
#define Dyy(j)    mean2d(j,ig(5)+10)
    sg(ng) = 11

    groupname(ng) = 'RyyBudget'
    varname(ng) = 'Ryy_t Byy Cyy Pyy Eyy PIyy Fyy Tyyy_y Tyyy Gyy Dyy'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rzz_t(j)  mean2d(j,ig(6)  )
#define Bzz(j)    mean2d(j,ig(6)+1)
#define Czz(j)    mean2d(j,ig(6)+2)
#define Pzz(j)    mean2d(j,ig(6)+3)
#define Ezz(j)    mean2d(j,ig(6)+4)
#define PIzz(j)   mean2d(j,ig(6)+5)
#define Fzz(j)    mean2d(j,ig(6)+6)
#define Tzzy_y(j) mean2d(j,ig(6)+7)
#define Tzzy(j)   mean2d(j,ig(6)+8)
#define Gzz(j)    mean2d(j,ig(6)+9)
#define Dzz(j)    mean2d(j,ig(6)+10)
    sg(ng) = 11

    groupname(ng) = 'RzzBudget'
    varname(ng) = 'Rzz_t Bzz Czz Pzz Ezz PIzz Fzz Tzzy_y Tzzy Gzz Dzz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxy_t(j)  mean2d(j,ig(7)  )
#define Bxy(j)    mean2d(j,ig(7)+1)
#define Cxy(j)    mean2d(j,ig(7)+2)
#define Pxy(j)    mean2d(j,ig(7)+3)
#define Exy(j)    mean2d(j,ig(7)+4)
#define PIxy(j)   mean2d(j,ig(7)+5)
#define Fxy(j)    mean2d(j,ig(7)+6)
#define Txyy_y(j) mean2d(j,ig(7)+7)
#define Txyy(j)   mean2d(j,ig(7)+8)
#define Gxy(j)    mean2d(j,ig(7)+9)
#define Dxy(j)    mean2d(j,ig(7)+10)
    sg(ng) = 11

    groupname(ng) = 'RxyBudget'
    varname(ng) = 'Rxy_t Bxy Cxy Pxy Exy PIxy Fxy Txyy_y Txyy Gxy Dxy'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Rxz_t(j)  mean2d(j,ig(8)  )
#define Bxz(j)    mean2d(j,ig(8)+1)
#define Cxz(j)    mean2d(j,ig(8)+2)
#define Pxz(j)    mean2d(j,ig(8)+3)
#define Exz(j)    mean2d(j,ig(8)+4)
#define PIxz(j)   mean2d(j,ig(8)+5)
#define Fxz(j)    mean2d(j,ig(8)+6)
#define Txzy_y(j) mean2d(j,ig(8)+7)
#define Txzy(j)   mean2d(j,ig(8)+8)
#define Gxz(j)    mean2d(j,ig(8)+9)
#define Dxz(j)    mean2d(j,ig(8)+10)
    sg(ng) = 11

    groupname(ng) = 'RxzBudget'
    varname(ng) = 'Rxz_t Bxz Cxz Pxz Exz PIxz Fxz Txzy_y Txzy Gxz Dxz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Ryz_t(j)  mean2d(j,ig(9)  )
#define Byz(j)    mean2d(j,ig(9)+1)
#define Cyz(j)    mean2d(j,ig(9)+2)
#define Pyz(j)    mean2d(j,ig(9)+3)
#define Eyz(j)    mean2d(j,ig(9)+4)
#define PIyz(j)   mean2d(j,ig(9)+5)
#define Fyz(j)    mean2d(j,ig(9)+6)
#define Tyzy_y(j) mean2d(j,ig(9)+7)
#define Tyzy(j)   mean2d(j,ig(9)+8)
#define Gyz(j)    mean2d(j,ig(9)+9)
#define Dyz(j)    mean2d(j,ig(9)+10)
    sg(ng) = 11

    groupname(ng) = 'RyzBudget'
    varname(ng) = 'Ryz_t Byz Cyz Pyz Eyz PIyz Fyz Tyzy_y Tyzy Gyz Dyz'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Tke_t(j)  mean2d(j,ig(10)  )
#define Buo(j)    mean2d(j,ig(10)+1)
#define Con(j)    mean2d(j,ig(10)+2)
#define Prd(j)    mean2d(j,ig(10)+3)
#define Eps(j)    mean2d(j,ig(10)+4)
#define Pi(j)     mean2d(j,ig(10)+5)
#define Ty_y(j)   mean2d(j,ig(10)+6)
#define Ty1(j)    mean2d(j,ig(10)+7)
#define Ty2(j)    mean2d(j,ig(10)+8)
#define Ty3(j)    mean2d(j,ig(10)+9)
#define Ty1_y(j)  mean2d(j,ig(10)+10)
#define Ty2_y(j)  mean2d(j,ig(10)+11)
#define Ty3_y(j)  mean2d(j,ig(10)+12)
#define Gkin(j)   mean2d(j,ig(10)+13)
#define Dkin(j)   mean2d(j,ig(10)+14)
#define Phi(j)    mean2d(j,ig(10)+15)
#define ugradp(j) mean2d(j,ig(10)+16)
    sg(ng) = 17

    groupname(ng) = 'TkeBudget'
    varname(ng) = 'Tke_t Buo Con Prd Eps Pi Trp Trp1 Trp2 Trp3 Trp1_y Trp2_y Trp3_y G D Phi UgradP'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rU3(j)   mean2d(j,ig(11)  )
#define rU4(j)   mean2d(j,ig(11)+1)
#define rV3(j)   mean2d(j,ig(11)+2)
#define rV4(j)   mean2d(j,ig(11)+3)
#define rW3(j)   mean2d(j,ig(11)+4)
#define rW4(j)   mean2d(j,ig(11)+5)
    sg(ng) = 6

    groupname(ng) = 'HigherOrder'
    varname(ng) = 'rU3 rU4 rV3 rV4 rW3 rW4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define U_y1(j)  mean2d(j,ig(12)  )
#define V_y1(j)  mean2d(j,ig(12)+1)
#define W_y1(j)  mean2d(j,ig(12)+2)
#define U_ii2(j) mean2d(j,ig(12)+3)
#define U_x2(j)  mean2d(j,ig(12)+4)
#define U_y2(j)  mean2d(j,ig(12)+5)
#define U_z2(j)  mean2d(j,ig(12)+6)
#define V_x2(j)  mean2d(j,ig(12)+7)
#define V_y2(j)  mean2d(j,ig(12)+8)
#define V_z2(j)  mean2d(j,ig(12)+9)
#define W_x2(j)  mean2d(j,ig(12)+10)
#define W_y2(j)  mean2d(j,ig(12)+11)
#define W_z2(j)  mean2d(j,ig(12)+12)
#define U_x3(j)  mean2d(j,ig(12)+13)
#define U_y3(j)  mean2d(j,ig(12)+14)
#define U_z3(j)  mean2d(j,ig(12)+15)
#define V_x3(j)  mean2d(j,ig(12)+16)
#define V_y3(j)  mean2d(j,ig(12)+17)
#define V_z3(j)  mean2d(j,ig(12)+18)
#define W_x3(j)  mean2d(j,ig(12)+19)
#define W_y3(j)  mean2d(j,ig(12)+20)
#define W_z3(j)  mean2d(j,ig(12)+21)
#define U_x4(j)  mean2d(j,ig(12)+22)
#define U_y4(j)  mean2d(j,ig(12)+23)
#define U_z4(j)  mean2d(j,ig(12)+24)
#define V_x4(j)  mean2d(j,ig(12)+25)
#define V_y4(j)  mean2d(j,ig(12)+26)
#define V_z4(j)  mean2d(j,ig(12)+27)
#define W_x4(j)  mean2d(j,ig(12)+28)
#define W_y4(j)  mean2d(j,ig(12)+29)
#define W_z4(j)  mean2d(j,ig(12)+30)
    sg(ng) = 31

    groupname(ng) = 'DerivativeFluctuations'
    varname(ng) = 'U_y1 V_y1 W_y1 U_ii2 ' &
                  //'U_x2 U_y2 U_z2 V_x2 V_y2 V_z2 W_x2 W_y2 W_z2 ' &
                  //'U_x3 U_y3 U_z3 V_x3 V_y3 V_z3 W_x3 W_y3 W_z3 ' &
                  //'U_x4 U_y4 U_z4 V_x4 V_y4 V_z4 W_x4 W_y4 W_z4'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rGamma(j)  mean2d(j,ig(13)  )
#define c2(j)      mean2d(j,ig(13)+1)
#define rho_ac(j)  mean2d(j,ig(13)+2)
#define rho_en(j)  mean2d(j,ig(13)+3)
#define T_ac(j)    mean2d(j,ig(13)+4)
#define T_en(j)    mean2d(j,ig(13)+5)
#define M_t(j)     mean2d(j,ig(13)+6)
#define rRP(j)     mean2d(j,ig(13)+7)
#define rRT(j)     mean2d(j,ig(13)+8)
    sg(ng) = 9

    groupname(ng) = 'Acoustics'
    varname(ng) = 'gamma C2 Rho_ac Rho_en T_ac T_en M_t rRP rRT'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rR2_flux_x(j) mean2d(j,ig(14)  )
#define rR2_flux_y(j) mean2d(j,ig(14)+1)
#define rR2_flux_z(j) mean2d(j,ig(14)+2)
#define rR2_dil1(j)   mean2d(j,ig(14)+3)
#define rR2_dil2(j)   mean2d(j,ig(14)+4)
#define rR2_trp(j)    mean2d(j,ig(14)+5)
#define rR2_prod(j)   mean2d(j,ig(14)+6)
#define rR2_conv(j)   mean2d(j,ig(14)+7)
    sg(ng) = 8

    groupname(ng) = 'RhoBudget'
    varname(ng) = 'RhoFluxX RhoFluxY RhoFluxZ RhoDil1 RhoDil2 RhoTrp RhoProd RhoConv'

    ! -----------------------------------------------------------------------
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define Pot(j)       mean2d(j,ig(15)  )
#define rref(j)      mean2d(j,ig(15)+1)
#define tref(j)      mean2d(j,ig(15)+2)
#define bfreq_fr(j)  mean2d(j,ig(15)+3)
#define bfreq_eq(j)  mean2d(j,ig(15)+4)
#define lapse_fr(j)  mean2d(j,ig(15)+5)
#define lapse_eq(j)  mean2d(j,ig(15)+6)
#define potem_fr(j)  mean2d(j,ig(15)+7)
#define potem_eq(j)  mean2d(j,ig(15)+8)
#define psat(j)      mean2d(j,ig(15)+9)
#define pref(j)      mean2d(j,ig(15)+10)
#define relhum(j)    mean2d(j,ig(15)+11)
#define dewpoint(j)  mean2d(j,ig(15)+12)
#define lapse_dew(j) mean2d(j,ig(15)+13)
    sg(ng) = 14

    groupname(ng) = 'Stratification'
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        varname(ng) = 'Pot rRref rTref BuoyFreq_fr BuoyFreq_eq LapseRate_fr LapseRate_eq ' &
                      //'PotTemp PotTemp_v SaturationPressure rPref RelativeHumidity Dewpoint LapseRate_dew'
    else
        varname(ng) = 'Pot rRref rTref BuoyFreq_fr BuoyFreq_eq LapseRate_fr LapseRate_eq ' &
                      //'PotTemp_fr PotTemp_eq SaturationPressure rPref RelativeHumidity Dewpoint LapseRate_dew'
    end if

    ! -----------------------------------------------------------------------
    ! Auxiliary variables depending on y and t; this last group is not written
    ng = ng + 1; ig(ng) = ig(ng - 1) + sg(ng - 1)
#define rUf(j)    mean2d(j,ig(16))
#define rVf(j)    mean2d(j,ig(16)+1)
#define rWf(j)    mean2d(j,ig(16)+2)

#define rU_y(j)   mean2d(j,ig(16)+3)
#define rV_y(j)   mean2d(j,ig(16)+4)
#define rW_y(j)   mean2d(j,ig(16)+5)
#define fU_y(j)   mean2d(j,ig(16)+6)
#define fV_y(j)   mean2d(j,ig(16)+7)
#define fW_y(j)   mean2d(j,ig(16)+8)
#define rP_y(j)   mean2d(j,ig(16)+9)
#define rR_y(j)   mean2d(j,ig(16)+10)
#define rT_y(j)   mean2d(j,ig(16)+11)
#define rB_y(j)   mean2d(j,ig(16)+12)

#define Rxx_y(j)  mean2d(j,ig(16)+13)
#define Ryy_y(j)  mean2d(j,ig(16)+14)
#define Rzz_y(j)  mean2d(j,ig(16)+15)
#define Rxy_y(j)  mean2d(j,ig(16)+16)
#define Rxz_y(j)  mean2d(j,ig(16)+17)
#define Ryz_y(j)  mean2d(j,ig(16)+18)
#define rR2_y(j)  mean2d(j,ig(16)+19)

#define Tau_yy(j)   mean2d(j,ig(16)+20)
#define Tau_xy(j)   mean2d(j,ig(16)+21)
#define Tau_yz(j)   mean2d(j,ig(16)+22)
#define Tau_xy_y(j) mean2d(j,ig(16)+23)
#define Tau_yy_y(j) mean2d(j,ig(16)+24)
#define Tau_yz_y(j) mean2d(j,ig(16)+25)
    sg(ng) = 26

    ! -----------------------------------------------------------------------
    nv = ig(ng) + sg(ng) - 1
    if (MAX_AVG_TEMPORAL < nv) then
        call TLAB_WRITE_ASCII(efile, 'AVG_FLOW_XZ. Not enough space in local arrays.')
        call TLAB_STOP(DNS_ERROR_AVGTMP)
    end if
    mean2d(:, 1:nv) = 0.0_wp

    ng = ng - 1
    nv = ig(ng) + sg(ng) - 1 ! the last group is not written out

    ! #######################################################################
    write (line1, *) itime; line1 = 'Calculating flow statistics at It'//trim(adjustl(line1))//'...'
    call TLAB_WRITE_ASCII(lfile, line1)

    ! #######################################################################
    ! Preliminary for IBM usage
    ! #######################################################################
    ! Asign gammas for conditional averages (c.f. Pope, p.170 [5.305])
    if (imode_ibm == 1) then
        ep_0(:) = gamma_0; ep_1(:) = gamma_1
    end if

    ! ###################################################################
    ! Averages (do not overwrite dudz; it contains p for incompressible case)
    ! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 2')
#endif

    ! Velocity
    call AVG_IK_V(imax, jmax, kmax, jmax, u, g(1)%jac, g(3)%jac, rU(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, v, g(1)%jac, g(3)%jac, rV(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, w, g(1)%jac, g(3)%jac, rW(1), wrk1d, area)

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rU(1), rU_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rV(1), rV_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rW(1), rW_y(1))

    U_y1(:) = rU_y(:)
    V_y1(:) = rV_y(:)
    W_y1(:) = rW_y(:)

    ! Density and Favre avrages
    if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE) then
        rR(:) = rbackground(:)

        fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

    else if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_DENSITY(imax, jmax, kmax, s, dwdx)
        call AVG_IK_V(imax, jmax, kmax, jmax, dwdx, g(1)%jac, g(3)%jac, rR(1), wrk1d, area)

        fU(:) = rU(:); fV(:) = rV(:); fW(:) = rW(:)

    else
        call AVG_IK_V(imax, jmax, kmax, jmax, rho, g(1)%jac, g(3)%jac, rR(1), wrk1d, area)

        dwdx = rho*u
        dwdy = rho*v
        dwdz = rho*w
        call AVG_IK_V(imax, jmax, kmax, jmax, dwdx, g(1)%jac, g(3)%jac, fU(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dwdy, g(1)%jac, g(3)%jac, fV(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dwdz, g(1)%jac, g(3)%jac, fW(1), wrk1d, area)
        fU(:) = fU(:)/rR(:)
        fV(:) = fV(:)/rR(:)
        fW(:) = fW(:)/rR(:)

    end if

    rUf(:) = rU(:) - fU(:)
    rVf(:) = rV(:) - fV(:)
    rWf(:) = rW(:) - fW(:)

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rR(1), rR_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fU(1), fU_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fV(1), fV_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), fW(1), fW_y(1))

    ! Pressure
    call AVG_IK_V(imax, jmax, kmax, jmax, p, g(1)%jac, g(3)%jac, rP(1), wrk1d, area)
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rP(1), rP_y(1))

    ! #######################################################################
    ! Main covariances (do not overwrite dudz; it contains p for incompressible case)
    ! #######################################################################
    do j = 1, jmax
        dwdx(:, j, :) = u(:, j, :) - fU(j)
        dwdy(:, j, :) = v(:, j, :) - fV(j)
        dwdz(:, j, :) = w(:, j, :) - fW(j)
    end do

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        dvdx = dwdx*dwdx
        dvdy = dwdy*dwdy
        dvdz = dwdz*dwdz
    else
        dvdx = dwdx*dwdx*rho
        dvdy = dwdy*dwdy*rho
        dvdz = dwdz*dwdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, Rxx(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, Ryy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, Rzz(1), wrk1d, area)
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == imode_eqns)) then
        Rxx(:) = Rxx(:)/rR(:)
        Ryy(:) = Ryy(:)/rR(:)
        Rzz(:) = Rzz(:)/rR(:)
    end if

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        dvdx = dwdx*dwdy
        dvdy = dwdx*dwdz
        dvdz = dwdy*dwdz
    else
        dvdx = dwdx*dwdy*rho
        dvdy = dwdx*dwdz*rho
        dvdz = dwdy*dwdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, Rxy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, Rxz(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, Ryz(1), wrk1d, area)
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == imode_eqns)) then
        Rxy(:) = Rxy(:)/rR(:)
        Rxz(:) = Rxz(:)/rR(:)
        Ryz(:) = Ryz(:)/rR(:)
    end if

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rxx(1), Rxx_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ryy(1), Ryy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rzz(1), Rzz_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rxy(1), Rxy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Rxz(1), Rxz_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ryz(1), Ryz_y(1))

    ! Density
    if (.not. (imode_eqns == DNS_EQNS_INCOMPRESSIBLE)) then
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_DENSITY(imax, jmax, kmax, s, p_wrk3d)
            do j = 1, jmax
                p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - rR(j)
            end do
        else
            do j = 1, jmax
                p_wrk3d(:, j, :) = rho(:, j, :) - rR(j)
            end do
        end if
        dvdx = p_wrk3d*p_wrk3d
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, rR2(1), wrk1d, area)

        call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rR2(1), rR2_y(1))

        ! Density Fluctuations Budget
        do j = 1, jmax
            dvdx(:, j, :) = u(:, j, :) - rU(j)
            dvdy(:, j, :) = v(:, j, :) - rV(j)
            dvdz(:, j, :) = w(:, j, :) - rW(j)
        end do
        dvdx = dvdx*p_wrk3d
        dvdy = dvdy*p_wrk3d
        dvdz = dvdz*p_wrk3d
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, rR2_flux_x(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, rR2_flux_y(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, rR2_flux_z(1), wrk1d, area)
        dvdy = dvdy*p_wrk3d
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, rR2_trp(1), wrk1d, area)

    end if

    ! higher-order moments
    p_wrk3d = dwdx*dwdx*dwdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rU3(1), wrk1d, area)
    p_wrk3d = dwdx*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rU4(1), wrk1d, area)

    p_wrk3d = dwdy*dwdy*dwdy
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rV3(1), wrk1d, area)
    p_wrk3d = dwdy*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rV4(1), wrk1d, area)

    p_wrk3d = dwdz*dwdz*dwdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rW3(1), wrk1d, area)
    p_wrk3d = dwdz*p_wrk3d
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rW4(1), wrk1d, area)

    ! Triple-velocity correlations
    dvdx = dwdx*dwdx*dwdy
    dvdy = dwdy*dwdy*dwdy
    dvdz = dwdz*dwdz*dwdy
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == imode_eqns)) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, Txxy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, Tyyy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, Tzzy(1), wrk1d, area)
    Ty1(:) = (Txxy(:) + Tyyy(:) + Tzzy(:))*0.5_wp

    dvdx = dwdx*dwdy*dwdy
    dvdy = dwdx*dwdy*dwdz
    dvdz = dwdy*dwdy*dwdz
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == imode_eqns)) then
        dvdx = dvdx*rho
        dvdy = dvdy*rho
        dvdz = dvdz*rho
    end if
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, Txyy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, Txzy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, Tyzy(1), wrk1d, area)

    ! Pressure
    do j = 1, jmax
        dvdz(:, j, :) = p(:, j, :) - rP(j)
    end do
    p_wrk3d = dvdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rP2(1), wrk1d, area)

    ! Pressure-velocity correlation in TKE transport terms
    dwdx = dwdx*dvdz
    dwdy = dwdy*dvdz
    dwdz = dwdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, jmax, dwdx, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Txyy(:) = Txyy(:) + wrk1d(1:jmax, 2)
    call AVG_IK_V(imax, jmax, kmax, jmax, dwdy, g(1)%jac, g(3)%jac, Ty2(1), wrk1d, area)
    Tyyy(:) = Tyyy(:) + Ty2(:)*2.0_wp
    call AVG_IK_V(imax, jmax, kmax, jmax, dwdz, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Tyzy(:) = Tyzy(:) + wrk1d(1:jmax, 2)

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
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), p, dwdx)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), p, dwdy)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), p, dwdz)
    p_wrk3d = u*dwdx + v*dwdy + w*dwdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, ugradp(1), wrk1d, area)

    ! Pressure Strain Terms
    ! 9 derivatives are here recomputed; ok, this routine is not called that often
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, dudx)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, dvdy)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, dwdz)
    dudx = dvdz*dudx ! dvdz contains the pressure fluctuation
    dvdy = dvdz*dvdy ! no need to substract rV_y
    dwdz = dvdz*dwdz
    call AVG_IK_V(imax, jmax, kmax, jmax, dudx, g(1)%jac, g(3)%jac, PIxx(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, PIyy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dwdz, g(1)%jac, g(3)%jac, PIzz(1), wrk1d, area)
    PIxx(:) = PIxx(:)*2.0_wp
    PIyy(:) = PIyy(:)*2.0_wp
    PIzz(:) = PIzz(:)*2.0_wp

    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, dudy)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, dvdx)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, dwdz) !dudz not free
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, dwdx)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, dudx) !dvdz not free
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, dvdy)
    dudy = dvdz*(dudy + dvdx) ! no need to substract rU_y
    dwdz = dvdz*(dwdz + dwdx)
    dudx = dvdz*(dudx + dvdy) ! no need to substract rW_y
    call AVG_IK_V(imax, jmax, kmax, jmax, dudy, g(1)%jac, g(3)%jac, PIxy(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dwdz, g(1)%jac, g(3)%jac, PIxz(1), wrk1d, area)
    call AVG_IK_V(imax, jmax, kmax, jmax, dudx, g(1)%jac, g(3)%jac, PIyz(1), wrk1d, area)

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

    if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE) then
        rT(:) = tbackground(:)

    else if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_TEMPERATURE(imax, jmax, kmax, s, T_LOC(1, 1, 1))
        call AVG_IK_V(imax, jmax, kmax, jmax, T_LOC(1, 1, 1), g(1)%jac, g(3)%jac, rT(1), wrk1d, area)

        do j = 1, jmax
            dvdx(:, j, :) = (T_LOC(:, j, :) - rT(j))**2
        end do
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, rT2(1), wrk1d, area)

        call Thermo_Psat_Polynomial(imax*jmax*kmax, T_LOC(1, 1, 1), dvdz)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, psat(1), wrk1d, area)

        call THERMO_ANELASTIC_RELATIVEHUMIDITY(imax, jmax, kmax, s, T_LOC(1, 1, 1), p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, relhum(1), wrk1d, area)

        call THERMO_ANELASTIC_THETA(imax, jmax, kmax, s, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, potem_fr(1), wrk1d, area)
        call THERMO_ANELASTIC_THETA_V(imax, jmax, kmax, s, p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, potem_eq(1), wrk1d, area)

        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), T_LOC(1, 1, 1), dudz)
        if (imixture == MIXT_TYPE_AIRWATER) &
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, 1, 1, 3), dudy)

 call THERMO_ANELASTIC_LAPSE_EQU(imax, jmax, kmax, s, dudz, dudy, GAMMA_LOC(1, 1, 1), p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, jmax, GAMMA_LOC(1, 1, 1), g(1)%jac, g(3)%jac, lapse_eq(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, bfreq_eq(1), wrk1d, area)
        bfreq_eq(:) = bfreq_eq(:)*buoyancy%vector(2)

        call THERMO_ANELASTIC_LAPSE_FR(imax, jmax, kmax, s, dudz, GAMMA_LOC(1, 1, 1), p_wrk3d)
        call AVG_IK_V(imax, jmax, kmax, jmax, GAMMA_LOC(1, 1, 1), g(1)%jac, g(3)%jac, lapse_fr(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, bfreq_fr(1), wrk1d, area)
        ! dummy = 1.0_wp /( scaleheight *gama0 )
        ! bfreq_fr(:) = -rR_y(:) /rbackground(:) -dummy *rR(:) /pbackground(:)
        bfreq_fr(:) = bfreq_fr(:)*buoyancy%vector(2)

        ! GAMMA_LOC(1,1,1) should contains lapse_fr, since lapse_dew = lapse_fr when saturated
        call THERMO_ANELASTIC_DEWPOINT(imax, jmax, kmax, s, p_wrk3d, GAMMA_LOC(1, 1, 1))
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, dewpoint(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, GAMMA_LOC(1, 1, 1), g(1)%jac, g(3)%jac, lapse_dew(1), wrk1d, area)

    else
        ! -------------------------------------------------------------------
        ! Main fields
        ! -------------------------------------------------------------------
        call THERMO_CALORIC_TEMPERATURE(imax*jmax*kmax, s, e, rho, T_LOC(1, 1, 1), p_wrk3d)
        call THERMO_GAMMA(imax*jmax*kmax, s, T_LOC(:, :, :), GAMMA_LOC(:, :, :))
        call THERMO_ENTROPY(imax*jmax*kmax, s, T_LOC(1, 1, 1), p, S_LOC(1, 1, 1))

        call AVG_IK_V(imax, jmax, kmax, jmax, T_LOC(1, 1, 1), g(1)%jac, g(3)%jac, rT(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, e, g(1)%jac, g(3)%jac, re(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, S_LOC(1, 1, 1), g(1)%jac, g(3)%jac, rs(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, GAMMA_LOC(1, 1, 1), g(1)%jac, g(3)%jac, rGamma(1), wrk1d, area)

        ! Means
        dudy = rho*e
        dudz = e + CRATIO_INV*p/rho
        dvdx = rho*dudz
        dvdy = rho*dwdz    ! rho *S_LOC
        dvdz = rho*dwdx    ! rho *T_LOC
        p_wrk3d = dudx*p/rho ! GAMMA_LOC *p /rho = speed of sound
        call AVG_IK_V(imax, jmax, kmax, jmax, dudy, g(1)%jac, g(3)%jac, fe(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dudz, g(1)%jac, g(3)%jac, rh(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, fh(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, fs(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, fT(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, c2(1), wrk1d, area)
        fe(:) = fe(:)/rR(:)
        fh(:) = fh(:)/rR(:)
        fs(:) = fs(:)/rR(:)
        fT(:) = fT(:)/rR(:)

        call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rT(1), rT_y(1))

        ! Turbulent Mach number
        M_t(:) = sqrt((Rxx(:) + Ryy(:) + Rzz(:))/c2(:))

        ! Covariances
        do j = 1, jmax
            dudy(:, j, :) = (S_LOC(:, j, :) - rs(j))**2
            dudz(:, j, :) = rho(:, j, :)*(S_LOC(:, j, :) - fs(j))**2
            dvdx(:, j, :) = (T_LOC(:, j, :) - rT(j))**2
            dvdy(:, j, :) = rho(:, j, :)*(T_LOC(:, j, :) - fT(j))**2
            dvdz(:, j, :) = (rho(:, j, :) - rR(j))*(T_LOC(:, j, :) - fT(j))
            p_wrk3d(:, j, :) = (rho(:, j, :) - rR(j))*(p(:, j, :) - rP(j))
        end do
        call AVG_IK_V(imax, jmax, kmax, jmax, dudy, g(1)%jac, g(3)%jac, rs2(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dudz, g(1)%jac, g(3)%jac, fs2(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, rT2(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, fT2(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, rRT(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rRP(1), wrk1d, area)
        fs2(:) = fs2(:)/rR(:)
        fT2(:) = fT2(:)/rR(:)

        do j = 1, jmax
            dudy(:, j, :) = (e(:, j, :) - re(j))**2
            dudz(:, j, :) = rho(:, j, :)*(e(:, j, :) - fe(j))**2
        end do
        call AVG_IK_V(imax, jmax, kmax, jmax, dudy, g(1)%jac, g(3)%jac, re2(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dudz, g(1)%jac, g(3)%jac, fe2(1), wrk1d, area)
        fe2(:) = fe2(:)/rR(:)

        p_wrk3d = e + CRATIO_INV*p/rho
        do j = 1, jmax
            dudy(:, j, :) = (p_wrk3d(:, j, :) - rh(j))**2
            dudz(:, j, :) = rho(:, j, :)*(p_wrk3d(:, j, :) - fh(j))**2
        end do
        call AVG_IK_V(imax, jmax, kmax, jmax, dudy, g(1)%jac, g(3)%jac, rh2(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dudz, g(1)%jac, g(3)%jac, fh2(1), wrk1d, area)
        fh2(:) = fh2(:)/rR(:)

        ! Acoustic and entropic density and temperature fluctuations
        do j = 1, jmax
            dudy(:, j, :) = p(:, j, :) - rP(j)                 ! pprime
            dudz(:, j, :) = dudy(:, j, :)/c2(j)               ! rho_ac
            dvdx(:, j, :) = rho(:, j, :) - rR(j) - dudz(:, j, :) ! rho_en = rprime - rho_ac
            dvdy(:, j, :) = (dudy(:, j, :)/rP(j) - dudz(:, j, :)/rR(j))*fT(j) ! T_ac
            dvdz(:, j, :) = T_LOC(:, j, :) - fT(j) - dvdy(:, j, :)               ! T_en = Tprime - T_ac
        end do
        dudz = dudz*dudz
        dvdx = dvdx*dvdx
        dvdy = dvdy*dvdy
        dvdz = dvdz*dvdz
        call AVG_IK_V(imax, jmax, kmax, jmax, dudz, g(1)%jac, g(3)%jac, rho_ac(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, rho_en(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, T_ac(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, T_en(1), wrk1d, area)

        ! -------------------------------------------------------------------
        ! Buoyancy frequency & saturation pressure
        ! -------------------------------------------------------------------
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), rho, dvdy)

        call Thermo_Psat_Polynomial(imax*jmax*kmax, T_LOC(1, 1, 1), dvdz)
        call THERMO_CP(imax*jmax*kmax, s, GAMMA_LOC(:, :, :), dvdx)

        do j = 1, jmax
            dudy(:, j, :) = dwdy(:, j, :)/p(:, j, :)/GAMMA_LOC(:, j, :) - dvdy(:, j, :)/rho(:, j, :)
            dvdx(:, j, :) = 1.0_wp/dvdx(:, j, :)
            dvdy(:, j, :) = T_LOC(:, j, :)*((p(:, j, :)/PREF_1000)**(1.0_wp/GAMMA_LOC(:, j, :) - 1.0_wp))
        end do
        call AVG_IK_V(imax, jmax, kmax, jmax, dudy, g(1)%jac, g(3)%jac, bfreq_fr(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, lapse_fr(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, potem_fr(1), wrk1d, area)
        call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, psat(1), wrk1d, area)
        bfreq_fr(:) = -bfreq_fr(:)*buoyancy%vector(2)
        lapse_fr(:) = -lapse_fr(:)*buoyancy%vector(2)*CRATIO_INV

#undef S_LOC

#define L_RATIO   dvdx
#define Q_RATIO   dvdy
#define RMEAN     dwdy
#define C_RATIO   dwdz

        if (imixture == MIXT_TYPE_AIRWATER) then
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), T_LOC(1, 1, 1), dudz)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, 1, 1, 2), dudy)

            dummy = Cvl + CRATIO_INV*Rv
            L_RATIO = -Lvl - dummy*dwdx ! dwdx is T_LOC
            L_RATIO = L_RATIO/(CRATIO_INV*Rv*dwdx)
            Q_RATIO = 1.0_wp/(p/dvdz - 1.0_wp)                ! dvdz is psat
            RMEAN = (Q_RATIO + 1.0_wp)*(1.0_wp - s(:, :, :, 1))*Rd

            p_wrk3d = (1.0_wp + Q_RATIO*L_RATIO)/RMEAN/ &
                    (dudx/(dudx - 1.0_wp) + Q_RATIO*L_RATIO*L_RATIO)  ! dudx is GAMMA_LOC
            call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, lapse_eq(1), wrk1d, area)
            lapse_eq(:) = -lapse_eq(:)*buoyancy%vector(2)/RRATIO

            p_wrk3d = (dudz - buoyancy%vector(2)/RRATIO*p_wrk3d)/dwdx &
                    *(1.0_wp + L_RATIO/rd_ov_rv/(1.0_wp - s(:, :, :, 1)))
            p_wrk3d = p_wrk3d - Rd/RMEAN*dudy
            call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, bfreq_eq(1), wrk1d, area)
            bfreq_eq(:) = -bfreq_eq(:)*buoyancy%vector(2)

            C_RATIO = Cd + s(:, :, :, 1)*Ldl
            C_RATIO = (1.0_wp - s(:, :, :, 1))*CRATIO_INV*Rv/C_RATIO
            p_wrk3d = dwdx/((p/PREF_1000)**C_RATIO)*exp(Q_RATIO*C_RATIO*L_RATIO)
            p_wrk3d = p_wrk3d*(1.0_wp + Q_RATIO)**C_RATIO/((p/dvdz)**(Q_RATIO*C_RATIO))
            call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, potem_eq(1), wrk1d, area)

        end if

#undef L_RATIO
#undef Q_RATIO
#undef RMEAN
#undef C_RATIO

    end if

#undef GAMMA_LOC
#undef T_LOC

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        pref(:) = pbackground(:)
        tref(:) = tbackground(:)
        rref(:) = rbackground(:)
    else
        pref(:) = rP(:)
        tref(:) = rT(:)
        rref(:) = rR(:)
    end if

    ! ###################################################################
    ! Potential energy
    !
    ! dudx = buoyancy
    !
    ! ###################################################################
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then

        if (buoyancy%type /= EQNS_NONE) then
            if (buoyancy%type == EQNS_EXPLICIT) then
                call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, dudx)
            else
                call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, dudx, bbackground)
            end if

            call AVG_IK_V(imax, jmax, kmax, jmax, dudx, g(1)%jac, g(3)%jac, rB(1), wrk1d, area)
            do j = 1, jmax
                dvdx(:, j, :) = (u(:, j, :) - rU(j))*(dudx(:, j, :) - rB(j))
                dvdy(:, j, :) = (v(:, j, :) - rV(j))*(dudx(:, j, :) - rB(j))
                dvdz(:, j, :) = (w(:, j, :) - rW(j))*(dudx(:, j, :) - rB(j))
            end do
            call AVG_IK_V(imax, jmax, kmax, jmax, dvdx, g(1)%jac, g(3)%jac, Bxx(1), wrk1d, area)
            call AVG_IK_V(imax, jmax, kmax, jmax, dvdy, g(1)%jac, g(3)%jac, Byy(1), wrk1d, area)
            call AVG_IK_V(imax, jmax, kmax, jmax, dvdz, g(1)%jac, g(3)%jac, Bzz(1), wrk1d, area)
            Bxy(:) = Bxx(:)*buoyancy%vector(2) + Byy(:)*buoyancy%vector(1) ! buoyancy%vector includes the Froude
            Bxz(:) = Bxx(:)*buoyancy%vector(3) + Bzz(:)*buoyancy%vector(1)
            Byz(:) = Byy(:)*buoyancy%vector(3) + Bzz(:)*buoyancy%vector(2)

            Bxx(:) = 2.0_wp*Bxx(:)*buoyancy%vector(1)
            Byy(:) = 2.0_wp*Byy(:)*buoyancy%vector(2)
            Bzz(:) = 2.0_wp*Bzz(:)*buoyancy%vector(3)

            dummy = 1.0_wp/froude
            rB(:) = rB(:)*dummy

            !        pmod(:) =-rP_y(:) + SIGN(rB(:),buoyancy%vector(2))

            call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), rB(1), rB_y(1))

        end if

    else ! Compressible case is not yet finished
        Bxx(:) = -rR(:)*rUf(:)*buoyancy%vector(1)
        Byy(:) = -rR(:)*rVf(:)*buoyancy%vector(2)
        Bzz(:) = -rR(:)*rWf(:)*buoyancy%vector(3)

        !     pmod(:) =-rP_y(:) +buoyancy%vector(2) *rR(:)

    end if

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
    call TLAB_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 3')
#endif

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, dudx)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), u, dudy)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, dudz)

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, dvdx)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, dvdy)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, dvdz)

    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, dwdx)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), w, dwdy)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, dwdz)

    ! ###################################################################
    ! Vorticity
    ! ###################################################################
    p_wrk3d = dwdy - dvdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, vortx(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - vortx(j))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, vortx2(1), wrk1d, area)

    p_wrk3d = dudz - dwdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, vorty(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - vorty(j))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, vorty2(1), wrk1d, area)

    p_wrk3d = dvdx - dudy
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, vortz(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - vortz(j))**2
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, vortz2(1), wrk1d, area)

    ! ###################################################################
    ! Derivatives Fluctuations
    ! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'AVG_FLOW_TEMPORAL_LAYER: Section 11')
#endif

    ! -------------------------------------------------------------------
    ! Longitudinal terms
    p_wrk3d = dudx*dudx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_x2(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dudx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_x3(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dudx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_x4(1), wrk1d, area)

    do j = 1, jmax
        p_wrk3d(:, j, :) = (dvdy(:, j, :) - rV_y(j))*(dvdy(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_y2(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dvdy(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_y3(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dvdy(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_y4(1), wrk1d, area)

    p_wrk3d = dwdz*dwdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_z2(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dwdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_z3(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dwdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_z4(1), wrk1d, area)

    ! -------------------------------------------------------------------
    ! Lateral terms U
    do j = 1, jmax
        p_wrk3d(:, j, :) = (dudy(:, j, :) - rU_y(j))*(dudy(:, j, :) - rU_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_y2(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dudy(:, j, :) - rU_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_y3(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dudy(:, j, :) - rU_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_y4(1), wrk1d, area)

    p_wrk3d = dudz*dudz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_z2(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dudz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_z3(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dudz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_z4(1), wrk1d, area)

    ! -------------------------------------------------------------------
    ! Lateral terms V
    p_wrk3d = dvdx*dvdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_x2(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dvdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_x3(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dvdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_x4(1), wrk1d, area)

    p_wrk3d = dvdz*dvdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_z2(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dvdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_z3(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dvdz
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, V_z4(1), wrk1d, area)

    ! -------------------------------------------------------------------
    ! Lateral terms W
    p_wrk3d = dwdx*dwdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_x2(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dwdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_x3(1), wrk1d, area)
    p_wrk3d = p_wrk3d*dwdx
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_x4(1), wrk1d, area)

    do j = 1, jmax
        p_wrk3d(:, j, :) = (dwdy(:, j, :) - rW_y(j))*(dwdy(:, j, :) - rW_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_y2(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dwdy(:, j, :) - rW_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_y3(1), wrk1d, area)
    do j = 1, jmax
        p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(dwdy(:, j, :) - rW_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, W_y4(1), wrk1d, area)

    ! -------------------------------------------------------------------
    ! Dilatation fluctuation
    p_wrk3d = dudx + dvdy + dwdz
    do j = 1, jmax
        p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - rV_y(j))*(p_wrk3d(:, j, :) - rV_y(j))
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, U_ii2(1), wrk1d, area)

    ! ###################################################################
    ! Density Fluctuations Budget
    ! ###################################################################
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) then
        p_wrk3d = dudx + dvdy + dwdz
        do j = 1, jmax
            p_wrk3d(:, j, :) = (p_wrk3d(:, j, :) - rV_y(j))*(rho(:, j, :) - rR(j))
        end do
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rR2_dil1(1), wrk1d, area)

        do j = 1, jmax
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :)*(rho(:, j, :) - rR(j))
        end do
        call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, rR2_dil2(1), wrk1d, area)

    end if

    ! ##################################################################
    ! Mean viscous dissipation rate
    ! ##################################################################
    p_wrk3d = dudx**2 + dvdy**2 + dwdz**2 + 0.5_wp*((dudy + dvdx)**2 + (dudz + dwdx)**2 + (dvdz + dwdy)**2) &
            - (dudx + dvdy + dwdz)**2/3.0_wp

    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Phi(1), wrk1d, area)
    Phi(:) = Phi(:)*visc*2.0_wp

    ! ###################################################################
    ! Dissipation Terms; final operation after viscous terms below
    ! ###################################################################
    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dudx + (dudy + dvdx)*dudy + (dudz + dwdx)*dudz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Exx(1), wrk1d, area)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dvdy*2.0_wp - p_wrk3d)*dvdy + (dudy + dvdx)*dvdx + (dvdz + dwdy)*dvdz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Eyy(1), wrk1d, area)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dwdz*2.0_wp - p_wrk3d)*dwdz + (dwdy + dvdz)*dwdy + (dwdx + dudz)*dwdx
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Ezz(1), wrk1d, area)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dvdx + (dudy + dvdx)*dvdy + (dudz + dwdx)*dvdz &
            + (dvdy*2.0_wp - p_wrk3d)*dudy + (dudy + dvdx)*dudx + (dvdz + dwdy)*dudz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Exy(1), wrk1d, area)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dudx*2.0_wp - p_wrk3d)*dwdx + (dudy + dvdx)*dwdy + (dudz + dwdx)*dwdz &
            + (dwdz*2.0_wp - p_wrk3d)*dudz + (dudz + dwdx)*dudx + (dvdz + dwdy)*dudy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Exz(1), wrk1d, area)

    p_wrk3d = (dudx + dvdy + dwdz)*c23
    p_wrk3d = (dvdy*2.0_wp - p_wrk3d)*dwdy + (dudy + dvdx)*dwdx + (dvdz + dwdy)*dwdz &
            + (dwdz*2.0_wp - p_wrk3d)*dvdz + (dudz + dwdx)*dvdx + (dvdz + dwdy)*dvdy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Eyz(1), wrk1d, area)

    ! ##################################################################
    ! Viscous shear-stress tensor
    ! ##################################################################
    p_wrk3d = dvdy*2.0_wp - dudx - dwdz
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tau_yy(1), wrk1d, area)
    do j = 1, jmax ! fluctuation tau22'
        dvdy(:, j, :) = (p_wrk3d(:, j, :) - Tau_yy(j))*c23
    end do
    Tau_yy(:) = Tau_yy(:)*visc*c23

    p_wrk3d = dudy + dvdx
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tau_xy(1), wrk1d, area)
    do j = 1, jmax ! fluctuation tau12'
        dudy(:, j, :) = p_wrk3d(:, j, :) - Tau_xy(j)
    end do
    Tau_xy(:) = Tau_xy(:)*visc

    p_wrk3d = dvdz + dwdy
    if (itransport == EQNS_TRANS_POWERLAW) p_wrk3d = p_wrk3d*vis
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, Tau_yz(1), wrk1d, area)
    do j = 1, jmax ! fluctuation tau23'
        dwdy(:, j, :) = p_wrk3d(:, j, :) - Tau_yz(j)
    end do
    Tau_yz(:) = Tau_yz(:)*visc

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_xy(1), Tau_xy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yy(1), Tau_yy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tau_yz(1), Tau_yz_y(1))

    ! -------------------------------------------------------------------
    ! Contribution to turbulent transport terms
    do j = 1, jmax
        p_wrk3d(:, j, :) = dudy(:, j, :)*(u(:, j, :) - fU(j)) ! -2*u'*tau12'
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Txxy(:) = Txxy(:) - wrk1d(1:jmax, 2)*visc*2.0_wp
    Ty3(:) = Ty3(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dvdy(:, j, :)*(v(:, j, :) - fV(j)) ! -2*v'*tau22'
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Tyyy(:) = Tyyy(:) - wrk1d(1:jmax, 2)*visc*2.0_wp
    Ty3(:) = Ty3(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dwdy(:, j, :)*(w(:, j, :) - fW(j)) ! -2*w'*tau23'
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Tzzy(:) = Tzzy(:) - wrk1d(1:jmax, 2)*visc*2.0_wp
    Ty3(:) = Ty3(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dvdy(:, j, :)*(u(:, j, :) - fU(j)) + dudy(:, j, :)*(v(:, j, :) - fV(j))! -u'*tau22' -v'*tau12'
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Txyy(:) = Txyy(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dwdy(:, j, :)*(u(:, j, :) - fU(j)) + dudy(:, j, :)*(w(:, j, :) - fW(j))! -u'*tau23' -w'*tau12'
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Txzy(:) = Txzy(:) - wrk1d(1:jmax, 2)*visc

    do j = 1, jmax
        p_wrk3d(:, j, :) = dwdy(:, j, :)*(v(:, j, :) - fV(j)) + dvdy(:, j, :)*(w(:, j, :) - fW(j))! -v'*tau23' -w'*tau22'
    end do
    call AVG_IK_V(imax, jmax, kmax, jmax, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 2), wrk1d, area)
    Tyzy(:) = Tyzy(:) - wrk1d(1:jmax, 2)*visc

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Txxy(1), Txxy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tyyy(1), Tyyy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tzzy(1), Tzzy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Txyy(1), Txyy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Txzy(1), Txzy_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Tyzy(1), Tyzy_y(1))

    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ty1(1), Ty1_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ty2(1), Ty2_y(1))
    call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), Ty3(1), Ty3_y(1))

    ! -------------------------------------------------------------------
    ! Contribution to dissipation
    Exx(:) = (Exx(:)*visc - Tau_xy(:)*rU_y(:))*2.0_wp
    Eyy(:) = (Eyy(:)*visc - Tau_yy(:)*rV_y(:))*2.0_wp
    Ezz(:) = (Ezz(:)*visc - Tau_yz(:)*rW_y(:))*2.0_wp
    Exy(:) = Exy(:)*visc - Tau_xy(:)*rV_y(:) - Tau_yy(:)*rU_y(:)
    Exz(:) = Exz(:)*visc - Tau_xy(:)*rW_y(:) - Tau_yz(:)*rU_y(:)
    Eyz(:) = Eyz(:)*visc - Tau_yy(:)*rW_y(:) - Tau_yz(:)*rV_y(:)

    ! ###################################################################
    ! Complete budget equations
    ! ###################################################################
    ! Density fluctuations budget equation
    if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) then
        rR2_prod(:) = -2.0_wp*(rR2_flux_y(:)*rR_y(:) + rR2(:)*rV_y(:))
        rR2_conv(:) = -rV(:)*rR2_y(:)
        rR2_dil1(:) = 2.0_wp*rR(:)*rR2_dil1(:)
    end if

    ! Rij Convective Terms
    Cxx(:) = -fV(:)*Rxx_y(:)
    Cyy(:) = -fV(:)*Ryy_y(:)
    Czz(:) = -fV(:)*Rzz_y(:)
    Cxy(:) = -fV(:)*Rxy_y(:)
    Cxz(:) = -fV(:)*Rxz_y(:)
    Cyz(:) = -fV(:)*Ryz_y(:)

    ! Rij Production Terms
    Pxx(:) = -2.0_wp*Rxy(:)*fU_y(:)
    Pyy(:) = -2.0_wp*Ryy(:)*fV_y(:)
    Pzz(:) = -2.0_wp*Ryz(:)*fW_y(:)
    Pxy(:) = -(Rxy(:)*fV_y(:) + Ryy(:)*fU_y(:))
    Pxz(:) = -(Rxy(:)*fW_y(:) + Ryz(:)*fU_y(:))
    Pyz(:) = -(Ryy(:)*fW_y(:) + Ryz(:)*fV_y(:))

    ! Rij Pressure Variable-Density Terms
    Gxx(:) = 0.0_wp
    Gyy(:) = 2.0_wp*rVf(:)*rP_y(:)
    Gzz(:) = 0.0_wp
    Gxy(:) = rUf(:)*rP_y(:)
    Gxz(:) = 0.0_wp
    Gyz(:) = rWf(:)*rP_y(:)

    ! Rij Viscous Variable-Density Terms
    Dxx(:) = 2.0_wp*rUf(:)*Tau_xy_y(:)
    Dyy(:) = 2.0_wp*rVf(:)*Tau_yy_y(:)
    Dzz(:) = 2.0_wp*rWf(:)*Tau_yz_y(:)
    Dxy(:) = rUf(:)*Tau_yy_y(:) + rVf(:)*Tau_xy_y(:)
    Dxz(:) = rUf(:)*Tau_yz_y(:) + rWf(:)*Tau_xy_y(:)
    Dyz(:) = rVf(:)*Tau_yz_y(:) + rWf(:)*Tau_yy_y(:)

    ! Rij Coriolis Terms
    if (coriolis%active(1) .and. coriolis%active(3)) then ! contribution from angular velocity Oy
        dummy = coriolis%vector(2)
        Fxx(:) = dummy*2.0_wp*Rxz(:)
        Fyy(:) = 0.0_wp
        Fzz(:) = -dummy*2.0_wp*Rxz(:)
        Fxy(:) = dummy*Ryz(:)
        Fxz(:) = dummy*(Rzz(:) - Rxx(:))
        Fyz(:) = -dummy*Rxy(:)
    end if

    ! Rij Buoyancy Terms; Calculated in Section Potential Energy

    ! Rij Transient terms
    Rxx_t(:) = -Fxx(:) + Bxx(:) + Cxx(:) + Pxx(:) - Exx(:) + (PIxx(:) - Txxy_y(:) - Gxx(:) + Dxx(:))/rR(:)
    Ryy_t(:) = -Fyy(:) + Byy(:) + Cyy(:) + Pyy(:) - Eyy(:) + (PIyy(:) - Tyyy_y(:) - Gyy(:) + Dyy(:))/rR(:)
    Rzz_t(:) = -Fzz(:) + Bzz(:) + Czz(:) + Pzz(:) - Ezz(:) + (PIzz(:) - Tzzy_y(:) - Gzz(:) + Dzz(:))/rR(:)
    Rxy_t(:) = -Fxy(:) + Bxy(:) + Cxy(:) + Pxy(:) - Exy(:) + (PIxy(:) - Txyy_y(:) - Gxy(:) + Dxy(:))/rR(:)
    Rxz_t(:) = -Fxz(:) + Bxz(:) + Cxz(:) + Pxz(:) - Exz(:) + (PIxz(:) - Txzy_y(:) - Gxz(:) + Dxz(:))/rR(:)
    Ryz_t(:) = -Fyz(:) + Byz(:) + Cyz(:) + Pyz(:) - Eyz(:) + (PIyz(:) - Tyzy_y(:) - Gyz(:) + Dyz(:))/rR(:)

    ! Kinetic energy equation
    Tke(:) = 0.5_wp*(Rxx(:) + Ryy(:) + Rzz(:))

    Buo(:) = 0.5_wp*(Bxx(:) + Byy(:) + Bzz(:))
    Con(:) = 0.5_wp*(Cxx(:) + Cyy(:) + Czz(:))
    Prd(:) = 0.5_wp*(Pxx(:) + Pyy(:) + Pzz(:))
    Pi(:) = 0.5_wp*(PIxx(:) + PIyy(:) + PIzz(:))
    Eps(:) = 0.5_wp*(Exx(:) + Eyy(:) + Ezz(:))
    Ty_y(:) = 0.5_wp*(Txxy_y(:) + Tyyy_y(:) + Tzzy_y(:))
    Gkin(:) = 0.5_wp*(Gxx(:) + Gyy(:) + Gzz(:))
    Dkin(:) = 0.5_wp*(Dxx(:) + Dyy(:) + Dzz(:))

    Tke_t(:) = Buo(:) + Con(:) + Prd(:) - Eps(:) + (-Ty_y(:) + Pi(:) - Gkin(:) + Dkin(:))/rR(:)

    ! Potential energy equation
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        Pot(:) = -rB(:)*(g(2)%nodes(:) - sbg(inb_scal)%ymean)

    else
        Pot(:) = -rR(:)*(g(2)%nodes(:) - rbg%ymean)*buoyancy%vector(2)

    end if

    ! ###################################################################
    ! Output
    ! ###################################################################
    ! 14 t-dependent variables, for consistency with old format
    ! ng = ng +1
    ! groupname(ng) = ''
    ! varname(ng)   = 'dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy dummy'
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''
    ! ng = ng +1; groupname(ng) = ''; varname(ng) = ''

    write (name, *) itime; name = 'avg'//trim(adjustl(name))
    call IO_WRITE_AVERAGES(name, itime, rtime, jmax, nv, ng, g(2)%nodes, varname, groupname, mean2d)

    return
end subroutine AVG_FLOW_XZ
