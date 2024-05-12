#include "dns_const.h"
#include "dns_error.h"
#include "avgij_map.h"

!########################################################################
!#
!# Post-processing statistical data accumulated in mean1d. Based on
!# mappings define in file avgij_map.h
!#
!########################################################################
subroutine AVG_FLOW_SPATIAL_LAYER(itxc, jmin_loc, jmax_loc, mean1d, stat)
    use TLAB_CONSTANTS, only: efile, tfile, wp, wi, big_wp
    use TLAB_VARS
    use TLAB_PROCS
    use TLAB_ARRAYS, only: wrk1d, wrk2d
    use Thermodynamics, only: MRATIO
    use OPR_PARTIAL
    implicit none

    integer(wi), intent(in) :: itxc                     ! size of array stat containing postprocess data
    integer(wi), intent(in) :: jmin_loc, jmax_loc       ! interval in the y direction
    real(wp), intent(in) :: mean1d(nstatavg, jmax, *)   ! array with raw mean data
    real(wp), intent(inout) :: stat(nstatavg, jmax, *)     ! local array with postprocess (final) mean data

! -------------------------------------------------------------------
#define rU(A,B)     stat(A,B,1)
#define rV(A,B)     stat(A,B,2)
#define rW(A,B)     stat(A,B,3)
#define rP(A,B)     stat(A,B,4)
#define rR(A,B)     stat(A,B,5)
#define rT(A,B)     stat(A,B,6)
#define rUf2(A,B)   stat(A,B,7)
#define rVf2(A,B)   stat(A,B,8)
#define rWf2(A,B)   stat(A,B,9)
#define rPf2(A,B)   stat(A,B,10)
#define rRf2(A,B)   stat(A,B,11)
#define rTf2(A,B)   stat(A,B,12)
#define rUfVf(A,B)  stat(A,B,13)
#define rUfWf(A,B)  stat(A,B,14)
#define rVfWf(A,B)  stat(A,B,15)

#define fU(A,B)     stat(A,B,16)
#define fV(A,B)     stat(A,B,17)
#define fW(A,B)     stat(A,B,18)
#define fT(A,B)     stat(A,B,19)
#define fRxx(A,B)   stat(A,B,20)
#define fRyy(A,B)   stat(A,B,21)
#define fRzz(A,B)   stat(A,B,22)
#define fTf2(A,B)   stat(A,B,23)
#define fRxy(A,B)   stat(A,B,24)
#define fRxz(A,B)   stat(A,B,25)
#define fRyz(A,B)   stat(A,B,26)

#define rTKE(A,B)   stat(A,B,27)
#define fTKE(A,B)   stat(A,B,28)

#define rbxx(A,B)   stat(A,B,29)
#define rbyy(A,B)   stat(A,B,30)
#define rbzz(A,B)   stat(A,B,31)
#define rbxy(A,B)   stat(A,B,32)
#define rbxz(A,B)   stat(A,B,33)
#define rbyz(A,B)   stat(A,B,34)

#define fbxx(A,B)   stat(A,B,35)
#define fbyy(A,B)   stat(A,B,36)
#define fbzz(A,B)   stat(A,B,37)
#define fbxy(A,B)   stat(A,B,38)
#define fbxz(A,B)   stat(A,B,39)
#define fbyz(A,B)   stat(A,B,40)

#define eta(A,B)       stat(A,B,41)
#define tau(A,B)       stat(A,B,42)
#define lambda(A,B)    stat(A,B,43)
#define lambda_x(A,B)  stat(A,B,44)
#define lambda_y(A,B)  stat(A,B,45)
#define lambda_z(A,B)  stat(A,B,46)
#define equi(A,B)      stat(A,B,47)

#define rdUdx(A,B)     stat(A,B,48)
#define rdUdy(A,B)     stat(A,B,49)
#define rdUdz(A,B)     stat(A,B,50)
#define rdVdx(A,B)     stat(A,B,51)
#define rdVdy(A,B)     stat(A,B,52)
#define rdVdz(A,B)     stat(A,B,53)
#define rdWdx(A,B)     stat(A,B,54)
#define rdWdy(A,B)     stat(A,B,55)
#define rdWdz(A,B)     stat(A,B,56)

#define rdUdxf2(A,B)   stat(A,B,57)
#define rdUdyf2(A,B)   stat(A,B,58)
#define rdUdzf2(A,B)   stat(A,B,59)
#define rdVdxf2(A,B)   stat(A,B,60)
#define rdVdyf2(A,B)   stat(A,B,61)
#define rdVdzf2(A,B)   stat(A,B,62)
#define rdWdxf2(A,B)   stat(A,B,63)
#define rdWdyf2(A,B)   stat(A,B,64)
#define rdWdzf2(A,B)   stat(A,B,65)

#define rdVdxfdUdyf(A,B) stat(A,B,66)
#define rdWdxfdUdzf(A,B) stat(A,B,67)
#define rdWdyfdVdzf(A,B) stat(A,B,68)
#define rdUdxfdVdyf(A,B) stat(A,B,69)
#define rdUdxfdWdzf(A,B) stat(A,B,70)
#define rdVdyfdWdzf(A,B) stat(A,B,71)

#define fdUdx(A,B)     stat(A,B,72)
#define fdUdy(A,B)     stat(A,B,73)
#define fdUdz(A,B)     stat(A,B,74)
#define fdVdx(A,B)     stat(A,B,75)
#define fdVdy(A,B)     stat(A,B,76)
#define fdVdz(A,B)     stat(A,B,77)
#define fdWdx(A,B)     stat(A,B,78)
#define fdWdy(A,B)     stat(A,B,79)
#define fdWdz(A,B)     stat(A,B,80)

#define dRxxdx(A,B)    stat(A,B,81)
#define dRxxdy(A,B)    stat(A,B,82)
#define dRxxdz(A,B)    stat(A,B,83)
#define dRyydx(A,B)    stat(A,B,84)
#define dRyydy(A,B)    stat(A,B,85)
#define dRyydz(A,B)    stat(A,B,86)
#define dRzzdx(A,B)    stat(A,B,87)
#define dRzzdy(A,B)    stat(A,B,88)
#define dRzzdz(A,B)    stat(A,B,89)
#define dRxydx(A,B)    stat(A,B,90)
#define dRxydy(A,B)    stat(A,B,91)
#define dRxydz(A,B)    stat(A,B,92)
#define dRxzdx(A,B)    stat(A,B,93)
#define dRxzdy(A,B)    stat(A,B,94)
#define dRxzdz(A,B)    stat(A,B,95)
#define dRyzdx(A,B)    stat(A,B,96)
#define dRyzdy(A,B)    stat(A,B,97)
#define dRyzdz(A,B)    stat(A,B,98)

#define Vortx(A,B)       stat(A,B,99)
#define Vorty(A,B)       stat(A,B,100)
#define Vortz(A,B)       stat(A,B,101)
#define Dil(A,B)         stat(A,B,102)
#define Vortxf2(A,B)     stat(A,B,103)
#define Vortyf2(A,B)     stat(A,B,104)
#define Vortzf2(A,B)     stat(A,B,105)
#define Dilf2(A,B)       stat(A,B,106)

#define rVis(A,B)        stat(A,B,107)

#define Conv_tt(A,B)     stat(A,B,108)
#define Prod_tt(A,B)     stat(A,B,109)
#define Diss_tt(A,B)     stat(A,B,110)
#define Tran_tt(A,B)     stat(A,B,111)
#define Pres_tt(A,B)     stat(A,B,112)
#define MnFl_tt(A,B)     stat(A,B,113)
#define Resi_tt(A,B)     stat(A,B,114)

#define fTKE_nf(A,B)     stat(A,B,115)
#define eps_f(A,B)       stat(A,B,116)

#define dRdx(A,B)        stat(A,B,117)
#define dRdy(A,B)        stat(A,B,118)
#define dRdz(A,B)        stat(A,B,119)

#define T1xx(A,B)        stat(A,B,120)
#define T1yy(A,B)        stat(A,B,121)
#define T1zz(A,B)        stat(A,B,122)
#define T1xy(A,B)        stat(A,B,123)
#define T1xz(A,B)        stat(A,B,124)
#define T1yz(A,B)        stat(A,B,125)

#define tau_xx(A,B)      stat(A,B,126)
#define tau_yy(A,B)      stat(A,B,127)
#define tau_zz(A,B)      stat(A,B,128)
#define tau_xy(A,B)      stat(A,B,129)
#define tau_xz(A,B)      stat(A,B,130)
#define tau_yz(A,B)      stat(A,B,131)

#define T4xx(A,B)        stat(A,B,132)
#define T4yy(A,B)        stat(A,B,133)
#define T4zz(A,B)        stat(A,B,134)
#define T4xy(A,B)        stat(A,B,135)
#define T4xz(A,B)        stat(A,B,136)
#define T4yz(A,B)        stat(A,B,137)
#define T4yx(A,B)        stat(A,B,138)
#define T4zx(A,B)        stat(A,B,139)
#define T4zy(A,B)        stat(A,B,140)

#define phi(A,B)         stat(A,B,141)

#define dPdx(A,B)        stat(A,B,142)
#define dPdy(A,B)        stat(A,B,143)
#define dPdz(A,B)        stat(A,B,144)

#define Conv_xx(A,B)     stat(A,B,145)
#define Prod_xx(A,B)     stat(A,B,146)
#define Diss_xx(A,B)     stat(A,B,147)
#define Tran_xx(A,B)     stat(A,B,148)
#define Pres_xx(A,B)     stat(A,B,149)
#define MnFl_xx(A,B)     stat(A,B,150)
#define Resi_xx(A,B)     stat(A,B,151)

#define Conv_yy(A,B)     stat(A,B,152)
#define Prod_yy(A,B)     stat(A,B,153)
#define Diss_yy(A,B)     stat(A,B,154)
#define Tran_yy(A,B)     stat(A,B,155)
#define Pres_yy(A,B)     stat(A,B,156)
#define MnFl_yy(A,B)     stat(A,B,157)
#define Resi_yy(A,B)     stat(A,B,158)

#define Conv_zz(A,B)     stat(A,B,159)
#define Prod_zz(A,B)     stat(A,B,160)
#define Diss_zz(A,B)     stat(A,B,161)
#define Tran_zz(A,B)     stat(A,B,162)
#define Pres_zz(A,B)     stat(A,B,163)
#define MnFl_zz(A,B)     stat(A,B,164)
#define Resi_zz(A,B)     stat(A,B,165)

#define Conv_xy(A,B)     stat(A,B,166)
#define Prod_xy(A,B)     stat(A,B,167)
#define Diss_xy(A,B)     stat(A,B,168)
#define Tran_xy(A,B)     stat(A,B,169)
#define Pres_xy(A,B)     stat(A,B,170)
#define MnFl_xy(A,B)     stat(A,B,171)
#define Resi_xy(A,B)     stat(A,B,172)

#define Conv(A,B)        stat(A,B,173)
#define Prod(A,B)        stat(A,B,174)
#define Diss(A,B)        stat(A,B,175)
#define Tran(A,B)        stat(A,B,176)
#define Pres(A,B)        stat(A,B,177)
#define MnFl(A,B)        stat(A,B,178)
#define Resi(A,B)        stat(A,B,179)

#define Conv_u(A,B)      stat(A,B,180)
#define Tran_u(A,B)      stat(A,B,181)
#define Reyn_u(A,B)      stat(A,B,182)
#define Resi_u(A,B)      stat(A,B,183)

#define Conv_v(A,B)      stat(A,B,184)
#define Tran_v(A,B)      stat(A,B,185)
#define Reyn_v(A,B)      stat(A,B,186)
#define Resi_v(A,B)      stat(A,B,187)

#define Conv_w(A,B)      stat(A,B,188)
#define Tran_w(A,B)      stat(A,B,189)
#define Reyn_w(A,B)      stat(A,B,190)
#define Resi_w(A,B)      stat(A,B,191)

#define Conv_p(A,B)      stat(A,B,192)
#define Reve_p(A,B)      stat(A,B,193)
#define Diss_p(A,B)      stat(A,B,194)
#define Tran_p(A,B)      stat(A,B,195)
#define Reyn_p(A,B)      stat(A,B,196)
#define Resi_p(A,B)      stat(A,B,197)

#define rho_ac(A,B)      stat(A,B,198)
#define rho_en(A,B)      stat(A,B,199)
#define T_ac(A,B)        stat(A,B,200)
#define T_en(A,B)        stat(A,B,201)
#define rho_p(A,B)       stat(A,B,202)
#define rho_T(A,B)       stat(A,B,203)

#define S_rho(A,B)       stat(A,B,204)
#define F_rho(A,B)       stat(A,B,205)
#define S_u(A,B)         stat(A,B,206)
#define F_u(A,B)         stat(A,B,207)
#define S_v(A,B)         stat(A,B,208)
#define F_v(A,B)         stat(A,B,209)
#define S_w(A,B)         stat(A,B,210)
#define F_w(A,B)         stat(A,B,211)
#define S_p(A,B)         stat(A,B,212)
#define F_p(A,B)         stat(A,B,213)
#define S_T(A,B)         stat(A,B,214)
#define F_T(A,B)         stat(A,B,215)

#define rRuT(A,B)        stat(A,B,216)
#define rRvT(A,B)        stat(A,B,217)
#define rRwT(A,B)        stat(A,B,218)
#define fRuT(A,B)        stat(A,B,219)
#define fRvT(A,B)        stat(A,B,220)
#define fRwT(A,B)        stat(A,B,221)

#define Conv_T(A,B)      stat(A,B,222)
#define Reve_T(A,B)      stat(A,B,223)
#define Diss_T(A,B)      stat(A,B,224)
#define Tran_T(A,B)      stat(A,B,225)
#define Reyn_T(A,B)      stat(A,B,226)
#define Resi_T(A,B)      stat(A,B,227)

#define LAST_INDEX 228

#define delta_m_u(A)      stat(A,1,LAST_INDEX)
#define delta_m_d(A)      stat(A,2,LAST_INDEX)
#define delta_w_u(A)      stat(A,3,LAST_INDEX)
#define delta_w_d(A)      stat(A,4,LAST_INDEX)
#define delta_u_u(A)      stat(A,5,LAST_INDEX)
#define delta_u_d(A)      stat(A,6,LAST_INDEX)
#define delta_01_u(A)     stat(A,7,LAST_INDEX)
#define delta_01_d(A)     stat(A,8,LAST_INDEX)
#define delta_t_u(A)      stat(A,9,LAST_INDEX)
#define delta_t_d(A)      stat(A,10,LAST_INDEX)
#define delta_r_u(A)      stat(A,11,LAST_INDEX)
#define delta_r_d(A)      stat(A,12,LAST_INDEX)
#define delta_u_center(A) stat(A,13,LAST_INDEX)
#define Reynolds_d(A)     stat(A,14,LAST_INDEX)
#define Reynolds_l(A)     stat(A,15,LAST_INDEX)
#define Reynolds_i(A)     stat(A,16,LAST_INDEX)
#define IntExcMomU(A)     stat(A,17,LAST_INDEX)
#define IntExcMomP(A)     stat(A,18,LAST_INDEX)
#define IntExcMomRxx(A)   stat(A,19,LAST_INDEX)
#define IntFluxT(A)       stat(A,20,LAST_INDEX)
#define IntMassU(A)       stat(A,21,LAST_INDEX)
#define IntMassV(A)       stat(A,22,LAST_INDEX)
#define IntTkeK(A)        stat(A,23,LAST_INDEX)
#define IntTkePi(A)       stat(A,24,LAST_INDEX)
#define IntTkeP(A)        stat(A,25,LAST_INDEX)
#define IntTkeF(A)        stat(A,26,LAST_INDEX)

#define VARMX1D 26

    integer(wi) i, j, k, n, bcs(2, 2)
    real(wp) pts, c13, zero
    real(wp) dum1, dum2, dum3, dum4, dum5
    real(wp) SIMPSON_NU
    real(wp) U2, DU, UC, r05, r005, r09, T2, DH, R2
    real(wp) y_center, dt_mean
    real(wp) delta_05, delta_w, delta_t
    real(wp) dfTdx, dfTdy, dRTTdx, dRTTdy, dfTf2dx, dfTf2dy
    real(wp) tranttx, trantty, dRUTdx, dRVTdy
    real(wp) fdTdx, fdTdy, fdTdz
    integer(wi) nj, jloc_max(1)

    real(wp) VAUXPRE(4), VAUXPOS(28)
    integer(wi) ivauxpre, ivauxpos, ivauxdum
    character*32 name
    character*400 line1
    character*2750 line2

    integer, parameter :: i23 = 23
    
! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING AVG_FLOW_SPATIAL_LAYER')
#endif

    bcs = 0

    r05 = 0.5_wp
    r005 = 5.0_wp*1.0e-2_wp
    r09 = 9.0_wp/10.0_wp
    c13 = 1.0_wp/3.0_wp
    zero = 1.0e-6_wp

    if (nstatavg_points == 0) then
        call TLAB_WRITE_ASCII(efile, 'AVG_FLOW_SPATIAL_LAYER: Zero number of points')
        call TLAB_STOP(DNS_ERROR_STATZERO)
    else
        pts = 1.0_wp/real(nstatavg_points, wp)
    end if

    dt_mean = (rtime - rstattimeorg)/real(itime - istattimeorg, wp)
    U2 = qbg(1)%mean - 0.5_wp*qbg(1)%delta
    T2 = tbg%mean - 0.5_wp*tbg%delta
    R2 = rbg%mean - 0.5_wp*rbg%delta

    if (itxc < nstatavg*jmax*LAST_INDEX) then
        call TLAB_WRITE_ASCII(efile, 'AVG_FLOW_SPATIAL_LAYER: Not enough space in stat')
        call TLAB_STOP(DNS_ERROR_WRKSIZE)
    end if

    nj = jmax_loc - jmin_loc + 1

! ###################################################################
! Main loop
! ###################################################################
    do j = 1, jmax*nstatavg

! ###################################################################
! Reynolds Averages
! ###################################################################
        rU(j, 1) = MA_U(j)*pts
        rV(j, 1) = MA_V(j)*pts
        rW(j, 1) = MA_W(j)*pts
        rP(j, 1) = MA_P(j)*pts
        rR(j, 1) = MA_R(j)*pts
        rT(j, 1) = MA_T(j)*pts

        rUf2(j, 1) = MA_UU(j)*pts - rU(j, 1)*rU(j, 1)
        rVf2(j, 1) = MA_VV(j)*pts - rV(j, 1)*rV(j, 1)
        rWf2(j, 1) = MA_WW(j)*pts - rW(j, 1)*rW(j, 1)
        rUfVf(j, 1) = MA_UV(j)*pts - rU(j, 1)*rV(j, 1)
        rUfWf(j, 1) = MA_UW(j)*pts - rU(j, 1)*rW(j, 1)
        rVfWf(j, 1) = MA_VW(j)*pts - rV(j, 1)*rW(j, 1)
        rTKE(j, 1) = 0.5_wp*(rUf2(j, 1) + rVf2(j, 1) + rWf2(j, 1))
        rbxx(j, 1) = 0.5_wp*rUf2(j, 1)/rTKE(j, 1) - c13
        rbyy(j, 1) = 0.5_wp*rVf2(j, 1)/rTKE(j, 1) - c13
        rbzz(j, 1) = 0.5_wp*rWf2(j, 1)/rTKE(j, 1) - c13
        rbxy(j, 1) = 0.5_wp*rUfVf(j, 1)/rTKE(j, 1)
        rbxz(j, 1) = 0.5_wp*rUfWf(j, 1)/rTKE(j, 1)
        rbyz(j, 1) = 0.5_wp*rVfWf(j, 1)/rTKE(j, 1)

        rPf2(j, 1) = MA_PP(j)*pts - rP(j, 1)*rP(j, 1)
        rRf2(j, 1) = MA_RR(j)*pts - rR(j, 1)*rR(j, 1)
        rTf2(j, 1) = MA_TT(j)*pts - rT(j, 1)*rT(j, 1)

        rRuT(j, 1) = MA_TU(j)*pts - rT(j, 1)*rU(j, 1)
        rRvT(j, 1) = MA_TV(j)*pts - rT(j, 1)*rV(j, 1)
        rRwT(j, 1) = MA_TW(j)*pts - rT(j, 1)*rW(j, 1)

! ###################################################################
! Favre Averages
! ###################################################################
        fU(j, 1) = MA_RU(j)*pts/rR(j, 1)
        fV(j, 1) = MA_RV(j)*pts/rR(j, 1)
        fW(j, 1) = MA_RW(j)*pts/rR(j, 1)
        fT(j, 1) = MA_RT(j)*pts/rR(j, 1)

        fRxx(j, 1) = MA_RUU(j)*pts/rR(j, 1) - fU(j, 1)*fU(j, 1)
        fRyy(j, 1) = MA_RVV(j)*pts/rR(j, 1) - fV(j, 1)*fV(j, 1)
        fRzz(j, 1) = MA_RWW(j)*pts/rR(j, 1) - fW(j, 1)*fW(j, 1)
        fRxy(j, 1) = MA_RUV(j)*pts/rR(j, 1) - fU(j, 1)*fV(j, 1)
        fRxz(j, 1) = MA_RUW(j)*pts/rR(j, 1) - fU(j, 1)*fW(j, 1)
        fRyz(j, 1) = MA_RVW(j)*pts/rR(j, 1) - fV(j, 1)*fW(j, 1)
        fTKE(j, 1) = 0.5_wp*(fRxx(j, 1) + fRyy(j, 1) + fRzz(j, 1))
        fbxx(j, 1) = 0.5_wp*fRxx(j, 1)/fTKE(j, 1) - c13
        fbyy(j, 1) = 0.5_wp*fRyy(j, 1)/fTKE(j, 1) - c13
        fbzz(j, 1) = 0.5_wp*fRzz(j, 1)/fTKE(j, 1) - c13
        fbxy(j, 1) = 0.5_wp*fRxy(j, 1)/fTKE(j, 1)
        fbxz(j, 1) = 0.5_wp*fRxz(j, 1)/fTKE(j, 1)
        fbyz(j, 1) = 0.5_wp*fRyz(j, 1)/fTKE(j, 1)

        fTf2(j, 1) = MA_RTT(j)*pts/rR(j, 1) - fT(j, 1)*fT(j, 1)

        fRuT(j, 1) = MRATIO*MA_PU(j)*pts/rR(j, 1) - fU(j, 1)*fT(j, 1)
        fRvT(j, 1) = MRATIO*MA_PV(j)*pts/rR(j, 1) - fV(j, 1)*fT(j, 1)
        fRwT(j, 1) = MRATIO*MA_PW(j)*pts/rR(j, 1) - fW(j, 1)*fT(j, 1)

! the TKE before filtering is stored every iteration
        dum1 = 1.0_wp/real((itime - istattimeorg)*g(3)%size, wp)
        fTKE_nf(j, 1) = 0.5_wp*(MA_FLT_RUU(j) + MA_FLT_RVV(j) + MA_FLT_RWW(j) - &
                                (MA_FLT_RU(j)*MA_FLT_RU(j) &
                                 + MA_FLT_RV(j)*MA_FLT_RV(j) &
                                 + MA_FLT_RW(j)*MA_FLT_RW(j))*dum1/rR(j, 1))*dum1/rR(j, 1)

!     eps_f(j,1) = (fTKE_nf(j,1)-fTKE(j,1))/dt_mean/real(ifilt_step, wp) ! to be rechecked

! ###################################################################
! First derivative terms (Reynolds)
! ###################################################################
        dRdx(j, 1) = MA_Rx(j)*pts
        dRdy(j, 1) = MA_Ry(j)*pts
        dRdz(j, 1) = MA_Rz(j)*pts
        dPdx(j, 1) = MA_Px(j)*pts
        dPdy(j, 1) = MA_Py(j)*pts
        dPdz(j, 1) = MA_Pz(j)*pts

! velocities
        rdUdx(j, 1) = MA_Ux(j)*pts
        rdUdy(j, 1) = MA_Uy(j)*pts
        rdUdz(j, 1) = MA_Uz(j)*pts
        rdVdx(j, 1) = MA_Vx(j)*pts
        rdVdy(j, 1) = MA_Vy(j)*pts
        rdVdz(j, 1) = MA_Vz(j)*pts
        rdWdx(j, 1) = MA_Wx(j)*pts
        rdWdy(j, 1) = MA_Wy(j)*pts
        rdWdz(j, 1) = MA_Wz(j)*pts

        rdUdxf2(j, 1) = MA_Ux2(j)*pts - rdUdx(j, 1)*rdUdx(j, 1)
        rdUdyf2(j, 1) = MA_Uy2(j)*pts - rdUdy(j, 1)*rdUdy(j, 1)
        rdUdzf2(j, 1) = MA_Uz2(j)*pts - rdUdz(j, 1)*rdUdz(j, 1)
        rdVdxf2(j, 1) = MA_Vx2(j)*pts - rdVdx(j, 1)*rdVdx(j, 1)
        rdVdyf2(j, 1) = MA_Vy2(j)*pts - rdVdy(j, 1)*rdVdy(j, 1)
        rdVdzf2(j, 1) = MA_Vz2(j)*pts - rdVdz(j, 1)*rdVdz(j, 1)
        rdWdxf2(j, 1) = MA_Wx2(j)*pts - rdWdx(j, 1)*rdWdx(j, 1)
        rdWdyf2(j, 1) = MA_Wy2(j)*pts - rdWdy(j, 1)*rdWdy(j, 1)
        rdWdzf2(j, 1) = MA_Wz2(j)*pts - rdWdz(j, 1)*rdWdz(j, 1)

        rdVdxfdUdyf(j, 1) = MA_VxUy(j)*pts - rdVdx(j, 1)*rdUdy(j, 1)
        rdWdxfdUdzf(j, 1) = MA_WxUz(j)*pts - rdWdx(j, 1)*rdUdz(j, 1)
        rdWdyfdVdzf(j, 1) = MA_WyVz(j)*pts - rdWdy(j, 1)*rdVdz(j, 1)
        rdUdxfdVdyf(j, 1) = MA_UXVY(j)*pts - rdUdx(j, 1)*rdVdy(j, 1)
        rdUdxfdWdzf(j, 1) = MA_UxWz(j)*pts - rdUdx(j, 1)*rdWdz(j, 1)
        rdVdyfdWdzf(j, 1) = MA_VyWz(j)*pts - rdVdy(j, 1)*rdWdz(j, 1)

! vorticity and dilatation
        Vortx(j, 1) = rdWdy(j, 1) - rdVdz(j, 1)
        Vorty(j, 1) = rdUdz(j, 1) - rdWdx(j, 1)
        Vortz(j, 1) = rdVdx(j, 1) - rdUdy(j, 1)
        Dil(j, 1) = rdUdx(j, 1) + rdVdy(j, 1) + rdWdz(j, 1)

        Vortxf2(j, 1) = rdWdyf2(j, 1) + rdVdzf2(j, 1) - 2.0_wp*rdWdyfdVdzf(j, 1)
        Vortyf2(j, 1) = rdUdzf2(j, 1) + rdWdxf2(j, 1) - 2.0_wp*rdWdxfdUdzf(j, 1)
        Vortzf2(j, 1) = rdVdxf2(j, 1) + rdUdyf2(j, 1) - 2.0_wp*rdVdxfdUdyf(j, 1)

        Dilf2(j, 1) = rdUdxf2(j, 1) + rdVdyf2(j, 1) + rdWdzf2(j, 1) + &
                      2.0_wp*(rdUdxfdVdyf(j, 1) + rdUdxfdWdzf(j, 1) + rdVdyfdWdzf(j, 1))

! ###################################################################
! First derivative terms (Favre)
! ###################################################################
        fdUdx(j, 1) = (MA_RUx(j) + MA_URx(j) - fU(j, 1)*MA_Rx(j))*pts/rR(j, 1)
        fdUdy(j, 1) = (MA_RUy(j) + MA_URy(j) - fU(j, 1)*MA_Ry(j))*pts/rR(j, 1)
        fdUdz(j, 1) = (MA_RUz(j) + MA_URz(j) - fU(j, 1)*MA_Rz(j))*pts/rR(j, 1)
        fdVdx(j, 1) = (MA_RVx(j) + MA_VRx(j) - fV(j, 1)*MA_Rx(j))*pts/rR(j, 1)
        fdVdy(j, 1) = (MA_RVy(j) + MA_VRy(j) - fV(j, 1)*MA_Ry(j))*pts/rR(j, 1)
        fdVdz(j, 1) = (MA_RVz(j) + MA_VRz(j) - fV(j, 1)*MA_Rz(j))*pts/rR(j, 1)
        fdWdx(j, 1) = (MA_RWx(j) + MA_WRx(j) - fW(j, 1)*MA_Rx(j))*pts/rR(j, 1)
        fdWdy(j, 1) = (MA_RWy(j) + MA_WRy(j) - fW(j, 1)*MA_Ry(j))*pts/rR(j, 1)
        fdWdz(j, 1) = (MA_RWz(j) + MA_WRz(j) - fW(j, 1)*MA_Rz(j))*pts/rR(j, 1)

! ###################################################################
! Derivatives of the Reynolds stresses
! ###################################################################
        dRxxdx(j, 1) = (MA_RUUx(j) - MA_RUU(j)/rR(j, 1)*dRdx(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdUdx(j, 1) - fU(j, 1)*fdUdx(j, 1)
        dRxxdy(j, 1) = (MA_RUUy(j) - MA_RUU(j)/rR(j, 1)*dRdy(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdUdy(j, 1) - fU(j, 1)*fdUdy(j, 1)
        dRxxdz(j, 1) = (MA_RUUz(j) - MA_RUU(j)/rR(j, 1)*dRdz(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdUdz(j, 1) - fU(j, 1)*fdUdz(j, 1)

        dRyydx(j, 1) = (MA_RVVx(j) - MA_RVV(j)/rR(j, 1)*dRdx(j, 1))*pts/rR(j, 1) - &
                       fV(j, 1)*fdVdx(j, 1) - fV(j, 1)*fdVdx(j, 1)
        dRyydy(j, 1) = (MA_RVVy(j) - MA_RVV(j)/rR(j, 1)*dRdy(j, 1))*pts/rR(j, 1) - &
                       fV(j, 1)*fdVdy(j, 1) - fV(j, 1)*fdVdy(j, 1)
        dRyydz(j, 1) = (MA_RVVz(j) - MA_RVV(j)/rR(j, 1)*dRdz(j, 1))*pts/rR(j, 1) - &
                       fV(j, 1)*fdVdz(j, 1) - fV(j, 1)*fdVdz(j, 1)

        dRzzdx(j, 1) = (MA_RWWx(j) - MA_RWW(j)/rR(j, 1)*dRdx(j, 1))*pts/rR(j, 1) - &
                       fW(j, 1)*fdWdx(j, 1) - fW(j, 1)*fdWdx(j, 1)
        dRzzdy(j, 1) = (MA_RWWy(j) - MA_RWW(j)/rR(j, 1)*dRdy(j, 1))*pts/rR(j, 1) - &
                       fW(j, 1)*fdWdy(j, 1) - fW(j, 1)*fdWdy(j, 1)
        dRzzdz(j, 1) = (MA_RWWz(j) - MA_RWW(j)/rR(j, 1)*dRdz(j, 1))*pts/rR(j, 1) - &
                       fW(j, 1)*fdWdz(j, 1) - fW(j, 1)*fdWdz(j, 1)

        dRxydx(j, 1) = (MA_RUVx(j) - MA_RUV(j)/rR(j, 1)*dRdx(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdVdx(j, 1) - fV(j, 1)*fdUdx(j, 1)
        dRxydy(j, 1) = (MA_RUVy(j) - MA_RUV(j)/rR(j, 1)*dRdy(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdVdy(j, 1) - fV(j, 1)*fdUdy(j, 1)
        dRxydz(j, 1) = (MA_RUVz(j) - MA_RUV(j)/rR(j, 1)*dRdz(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdVdz(j, 1) - fV(j, 1)*fdUdz(j, 1)

        dRxzdx(j, 1) = (MA_RUWx(j) - MA_RUW(j)/rR(j, 1)*dRdx(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdWdx(j, 1) - fW(j, 1)*fdUdx(j, 1)
        dRxzdy(j, 1) = (MA_RUWy(j) - MA_RUW(j)/rR(j, 1)*dRdy(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdWdy(j, 1) - fW(j, 1)*fdUdy(j, 1)
        dRxzdz(j, 1) = (MA_RUWz(j) - MA_RUW(j)/rR(j, 1)*dRdz(j, 1))*pts/rR(j, 1) - &
                       fU(j, 1)*fdWdz(j, 1) - fW(j, 1)*fdUdz(j, 1)

        dRyzdx(j, 1) = (MA_RVWx(j) - MA_RVW(j)/rR(j, 1)*dRdx(j, 1))*pts/rR(j, 1) - &
                       fV(j, 1)*fdWdx(j, 1) - fW(j, 1)*fdVdx(j, 1)
        dRyzdy(j, 1) = (MA_RVWy(j) - MA_RVW(j)/rR(j, 1)*dRdy(j, 1))*pts/rR(j, 1) - &
                       fV(j, 1)*fdWdy(j, 1) - fW(j, 1)*fdVdy(j, 1)
        dRyzdz(j, 1) = (MA_RVWz(j) - MA_RVW(j)/rR(j, 1)*dRdz(j, 1))*pts/rR(j, 1) - &
                       fV(j, 1)*fdWdz(j, 1) - fW(j, 1)*fdVdz(j, 1)

! ###################################################################
! Viscous shear stress tensor
! ###################################################################
        rVis(j, 1) = MA_VIS(j)*pts

        tau_xx(j, 1) = MA_TAUxx(j)*pts
        tau_yy(j, 1) = MA_TAUyy(j)*pts
        tau_zz(j, 1) = MA_TAUzz(j)*pts
        tau_xy(j, 1) = MA_TAUxy(j)*pts
        tau_xz(j, 1) = MA_TAUxz(j)*pts
        tau_yz(j, 1) = MA_TAUyz(j)*pts

        phi(j, 1) = (MA_TAUXkUk(j) + MA_TAUYkVk(j) + MA_TAUZkWk(j))*pts

! ###################################################################
! Transport equations
! ###################################################################
! -------------------------------------------------------------------
! auxiliar quantities
! -------------------------------------------------------------------
        dum1 = fU(j, 1)*dRdx(j, 1) + fV(j, 1)*dRdy(j, 1) + fW(j, 1)*dRdz(j, 1)

        dum2 = fU(j, 1)*fdUdx(j, 1) + fV(j, 1)*fdUdy(j, 1) + fW(j, 1)*fdUdz(j, 1)
        dum3 = fU(j, 1)*fdVdx(j, 1) + fV(j, 1)*fdVdy(j, 1) + fW(j, 1)*fdVdz(j, 1)
        dum4 = fU(j, 1)*fdWdx(j, 1) + fV(j, 1)*fdWdy(j, 1) + fW(j, 1)*fdWdz(j, 1)

        dum5 = fdUdx(j, 1) + fdVdy(j, 1) + fdWdz(j, 1)

! -------------------------------------------------------------------
! X-, Y-, and Z-Momentum equation
! -------------------------------------------------------------------
        Conv_u(j, 1) = -dum2
        Tran_u(j, 1) = (-dPdx(j, 1) + MA_TAUXkk(j)*pts)/rR(j, 1)
        Reyn_u(j, 1) = -dRxxdx(j, 1) - dRxydy(j, 1) - dRxzdz(j, 1) - &
                       (fRxx(j, 1)*dRdx(j, 1) + fRxy(j, 1)*dRdy(j, 1) + fRxz(j, 1)*dRdz(j, 1))/rR(j, 1)
        Resi_u(j, 1) = Conv_u(j, 1) + Tran_u(j, 1) + Reyn_u(j, 1)

        Conv_v(j, 1) = -dum3
        Tran_v(j, 1) = (-dPdy(j, 1) + MA_TAUYkk(j)*pts)/rR(j, 1)
        Reyn_v(j, 1) = -dRxydx(j, 1) - dRyydy(j, 1) - dRyzdz(j, 1) - &
                       (fRxy(j, 1)*dRdx(j, 1) + fRyy(j, 1)*dRdy(j, 1) + fRyz(j, 1)*dRdz(j, 1))/rR(j, 1)
        Resi_v(j, 1) = Conv_v(j, 1) + Tran_v(j, 1) + Reyn_v(j, 1)

        Conv_w(j, 1) = -dum4
        Tran_w(j, 1) = (-dPdz(j, 1) + MA_TAUZkk(j)*pts)/rR(j, 1)
        Reyn_w(j, 1) = -dRxzdx(j, 1) - dRyzdy(j, 1) - dRzzdz(j, 1) - &
                       (fRxz(j, 1)*dRdx(j, 1) + fRyz(j, 1)*dRdy(j, 1) + fRzz(j, 1)*dRdz(j, 1))/rR(j, 1)
        Resi_w(j, 1) = Conv_w(j, 1) + Tran_w(j, 1) + Reyn_w(j, 1)

! -------------------------------------------------------------------
! Convective element of transport term of Reynolds equations
! -------------------------------------------------------------------
        T1xx(j, 1) = (MA_RUUUkk(j) - MA_RUU(j)*dum5 - &
                      MA_RUU(j)*fdUdx(j, 1) - MA_RUV(j)*fdUdy(j, 1) - MA_RUW(j)*fdUdz(j, 1) - &
                      MA_RUU(j)*fdUdx(j, 1) - MA_RUV(j)*fdUdy(j, 1) - MA_RUW(j)*fdUdz(j, 1) - &
                      (MA_RUUx(j) + MA_RUVy(j) + MA_RUWz(j))*fU(j, 1) - &
                      (MA_RUUx(j) + MA_RUVy(j) + MA_RUWz(j))*fU(j, 1) - &
                      MA_RUUx(j)*fU(j, 1) - MA_RUUy(j)*fV(j, 1) - MA_RUUz(j)*fW(j, 1))*pts + &
                     2.0_wp*(fU(j, 1)*fU(j, 1)*dum1 + &
                             rR(j, 1)*(fU(j, 1)*fU(j, 1)*dum5 + fU(j, 1)*dum2 + fU(j, 1)*dum2))

        T1yy(j, 1) = (MA_RVVUkk(j) - MA_RVV(j)*dum5 - &
                      MA_RUV(j)*fdVdx(j, 1) - MA_RVV(j)*fdVdy(j, 1) - MA_RVW(j)*fdVdz(j, 1) - &
                      MA_RUV(j)*fdVdx(j, 1) - MA_RVV(j)*fdVdy(j, 1) - MA_RVW(j)*fdVdz(j, 1) - &
                      (MA_RUVx(j) + MA_RVVy(j) + MA_RVWz(j))*fV(j, 1) - &
                      (MA_RUVx(j) + MA_RVVy(j) + MA_RVWz(j))*fV(j, 1) - &
                      MA_RVVx(j)*fU(j, 1) - MA_RVVy(j)*fV(j, 1) - MA_RVVz(j)*fW(j, 1))*pts + &
                     2.0_wp*(fV(j, 1)*fV(j, 1)*dum1 + &
                             rR(j, 1)*(fV(j, 1)*fV(j, 1)*dum5 + fV(j, 1)*dum3 + fV(j, 1)*dum3))

        T1zz(j, 1) = (MA_RWWUkk(j) - MA_RWW(j)*dum5 - &
                      MA_RUW(j)*fdWdx(j, 1) - MA_RVW(j)*fdWdy(j, 1) - MA_RWW(j)*fdWdz(j, 1) - &
                      MA_RUW(j)*fdWdx(j, 1) - MA_RVW(j)*fdWdy(j, 1) - MA_RWW(j)*fdWdz(j, 1) - &
                      (MA_RUWx(j) + MA_RVWy(j) + MA_RWWz(j))*fW(j, 1) - &
                      (MA_RUWx(j) + MA_RVWy(j) + MA_RWWz(j))*fW(j, 1) - &
                      MA_RWWx(j)*fU(j, 1) - MA_RWWy(j)*fV(j, 1) - MA_RWWz(j)*fW(j, 1))*pts + &
                     2.0_wp*(fW(j, 1)*fW(j, 1)*dum1 + &
                             rR(j, 1)*(fW(j, 1)*fW(j, 1)*dum5 + fW(j, 1)*dum4 + fW(j, 1)*dum4))

        T1xy(j, 1) = (MA_RUVUkk(j) - MA_RUV(j)*dum5 - &
                      MA_RUU(j)*fdVdx(j, 1) - MA_RUV(j)*fdVdy(j, 1) - MA_RUW(j)*fdVdz(j, 1) - &
                      MA_RUV(j)*fdUdx(j, 1) - MA_RVV(j)*fdUdy(j, 1) - MA_RVW(j)*fdUdz(j, 1) - &
                      (MA_RUVx(j) + MA_RVVy(j) + MA_RVWz(j))*fU(j, 1) - &
                      (MA_RUUx(j) + MA_RUVy(j) + MA_RUWz(j))*fV(j, 1) - &
                      MA_RUVx(j)*fU(j, 1) - MA_RUVy(j)*fV(j, 1) - MA_RUVz(j)*fW(j, 1))*pts + &
                     2.0_wp*(fU(j, 1)*fV(j, 1)*dum1 + &
                             rR(j, 1)*(fU(j, 1)*fV(j, 1)*dum5 + fU(j, 1)*dum3 + fV(j, 1)*dum2))

        T1xz(j, 1) = (MA_RUWUkk(j) - MA_RUW(j)*dum5 - &
                      MA_RUU(j)*fdWdx(j, 1) - MA_RUV(j)*fdWdy(j, 1) - MA_RUW(j)*fdWdz(j, 1) - &
                      MA_RUW(j)*fdUdx(j, 1) - MA_RVW(j)*fdUdy(j, 1) - MA_RWW(j)*fdUdz(j, 1) - &
                      (MA_RUWx(j) + MA_RVWy(j) + MA_RWWz(j))*fU(j, 1) - &
                      (MA_RUUx(j) + MA_RUVy(j) + MA_RUWz(j))*fW(j, 1) - &
                      MA_RUWx(j)*fU(j, 1) - MA_RUWy(j)*fV(j, 1) - MA_RUWz(j)*fW(j, 1))*pts + &
                     2.0_wp*(fU(j, 1)*fW(j, 1)*dum1 + &
                             rR(j, 1)*(fU(j, 1)*fW(j, 1)*dum5 + fU(j, 1)*dum4 + fW(j, 1)*dum2))

        T1yz(j, 1) = (MA_RVWUkk(j) - MA_RVW(j)*dum5 - &
                      MA_RUV(j)*fdWdx(j, 1) - MA_RVV(j)*fdWdy(j, 1) - MA_RVW(j)*fdWdz(j, 1) - &
                      MA_RUW(j)*fdVdx(j, 1) - MA_RVW(j)*fdVdy(j, 1) - MA_RWW(j)*fdVdz(j, 1) - &
                      (MA_RUWx(j) + MA_RVWy(j) + MA_RWWz(j))*fV(j, 1) - &
                      (MA_RUVx(j) + MA_RVVy(j) + MA_RVWz(j))*fW(j, 1) - &
                      MA_RVWx(j)*fU(j, 1) - MA_RVWy(j)*fV(j, 1) - MA_RVWz(j)*fW(j, 1))*pts + &
                     2.0_wp*(fV(j, 1)*fW(j, 1)*dum1 + &
                             rR(j, 1)*(fV(j, 1)*fW(j, 1)*dum5 + fV(j, 1)*dum4 + fW(j, 1)*dum3))

! -------------------------------------------------------------------
! Viscous element of transport term of Reynolds equations
! -------------------------------------------------------------------
        T4xx(j, 1) = (MA_TAUXkUk(j) + MA_UTAUXkk(j) - rU(j, 1)*MA_TAUXkk(j))*pts - &
                     tau_xx(j, 1)*fdUdx(j, 1) - tau_xy(j, 1)*fdUdy(j, 1) - tau_xz(j, 1)*fdUdz(j, 1)

        T4xy(j, 1) = (MA_TAUYkUk(j) + MA_UTAUYkk(j) - rU(j, 1)*MA_TAUYkk(j))*pts - &
                     tau_xy(j, 1)*fdUdx(j, 1) - tau_yy(j, 1)*fdUdy(j, 1) - tau_yz(j, 1)*fdUdz(j, 1)

        T4xz(j, 1) = (MA_TAUZkUk(j) + MA_UTAUZkk(j) - rU(j, 1)*MA_TAUZkk(j))*pts - &
                     tau_xz(j, 1)*fdUdx(j, 1) - tau_yz(j, 1)*fdUdy(j, 1) - tau_zz(j, 1)*fdUdz(j, 1)

        T4yx(j, 1) = (MA_TAUXkVk(j) + MA_VTAUXkk(j) - rV(j, 1)*MA_TAUXkk(j))*pts - &
                     tau_xx(j, 1)*fdVdx(j, 1) - tau_xy(j, 1)*fdVdy(j, 1) - tau_xz(j, 1)*fdVdz(j, 1)

        T4yy(j, 1) = (MA_TAUYkVk(j) + MA_VTAUYkk(j) - rV(j, 1)*MA_TAUYkk(j))*pts - &
                     tau_xy(j, 1)*fdVdx(j, 1) - tau_yy(j, 1)*fdVdy(j, 1) - tau_yz(j, 1)*fdVdz(j, 1)

        T4yz(j, 1) = (MA_TAUZkVk(j) + MA_VTAUZkk(j) - rV(j, 1)*MA_TAUZkk(j))*pts - &
                     tau_xz(j, 1)*fdVdx(j, 1) - tau_yz(j, 1)*fdVdy(j, 1) - tau_zz(j, 1)*fdVdz(j, 1)

        T4zx(j, 1) = (MA_TAUXkWk(j) + MA_WTAUXkk(j) - rW(j, 1)*MA_TAUXkk(j))*pts - &
                     tau_xx(j, 1)*fdWdx(j, 1) - tau_xy(j, 1)*fdWdy(j, 1) - tau_xz(j, 1)*fdWdz(j, 1)

        T4zy(j, 1) = (MA_TAUYkWk(j) + MA_WTAUYkk(j) - rW(j, 1)*MA_TAUYkk(j))*pts - &
                     tau_xy(j, 1)*fdWdx(j, 1) - tau_yy(j, 1)*fdWdy(j, 1) - tau_yz(j, 1)*fdWdz(j, 1)

        T4zz(j, 1) = (MA_TAUZkWk(j) + MA_WTAUZkk(j) - rW(j, 1)*MA_TAUZkk(j))*pts - &
                     tau_xz(j, 1)*fdWdx(j, 1) - tau_yz(j, 1)*fdWdy(j, 1) - tau_zz(j, 1)*fdWdz(j, 1)

! -------------------------------------------------------------------
! Rxx Reynolds stress equation
! -------------------------------------------------------------------
        Conv_xx(j, 1) = -fU(j, 1)*dRxxdx(j, 1) - fV(j, 1)*dRxxdy(j, 1) - fW(j, 1)*dRxxdz(j, 1)
        Prod_xx(j, 1) = -2.0_wp*(fRxx(j, 1)*fdUdx(j, 1) + fRxy(j, 1)*fdUdy(j, 1) + fRxz(j, 1)*fdUdz(j, 1))
        Diss_xx(j, 1) = -2.0_wp*(MA_TAUXkUk(j)*pts - &
                                 tau_xx(j, 1)*rdUdx(j, 1) - tau_xy(j, 1)*rdUdy(j, 1) - tau_xz(j, 1)*rdUdz(j, 1))/rR(j, 1)

! pressure terms with Reynolds average
        Tran_xx(j, 1) = -(T1xx(j, 1) + 2.0_wp*(-T4xx(j, 1) + (MA_PUx(j) + MA_UPx(j))*pts - &
                                               rP(j, 1)*rdUdx(j, 1) - rU(j, 1)*dPdx(j, 1)))/rR(j, 1)
        Pres_xx(j, 1) = 2.0_wp*(MA_PUx(j)*pts - rP(j, 1)*rdUdx(j, 1))/rR(j, 1)
        MnFl_xx(j, 1) = 2.0_wp*(rU(j, 1) - fU(j, 1))*(MA_TAUXkk(j)*pts - dPdx(j, 1))/rR(j, 1)

        Resi_xx(j, 1) = Conv_xx(j, 1) + Prod_xx(j, 1) + Diss_xx(j, 1) + Tran_xx(j, 1) + &
                        Pres_xx(j, 1) + MnFl_xx(j, 1)

! -------------------------------------------------------------------
! Ryy Reynolds stress equation
! -------------------------------------------------------------------
        Conv_yy(j, 1) = -fU(j, 1)*dRyydx(j, 1) - fV(j, 1)*dRyydy(j, 1) - fW(j, 1)*dRyydz(j, 1)
        Prod_yy(j, 1) = -2.0_wp*(fRxy(j, 1)*fdVdx(j, 1) + fRyy(j, 1)*fdVdy(j, 1) + fRyz(j, 1)*fdVdz(j, 1))
        Diss_yy(j, 1) = -2.0_wp*(MA_TAUYkVk(j)*pts - &
                                 tau_xy(j, 1)*rdVdx(j, 1) - tau_yy(j, 1)*rdVdy(j, 1) - tau_yz(j, 1)*rdVdz(j, 1))/rR(j, 1)

! pressure terms with Reynolds average
        Tran_yy(j, 1) = -(T1yy(j, 1) + 2.0_wp*(-T4yy(j, 1) + (MA_PVY(j) + MA_VPy(j))*pts - &
                                               rP(j, 1)*rdVdy(j, 1) - rV(j, 1)*dPdy(j, 1)))/rR(j, 1)
        Pres_yy(j, 1) = 2.0_wp*(MA_PVY(j)*pts - rP(j, 1)*rdVdy(j, 1))/rR(j, 1)
        MnFl_yy(j, 1) = 2.0_wp*(rV(j, 1) - fV(j, 1))*(MA_TAUYkk(j)*pts - dPdy(j, 1))/rR(j, 1)

        Resi_yy(j, 1) = Conv_yy(j, 1) + Prod_yy(j, 1) + Diss_yy(j, 1) + Tran_yy(j, 1) + &
                        Pres_yy(j, 1) + MnFl_yy(j, 1)

! -------------------------------------------------------------------
! Rzz Reynolds stress equation
! -------------------------------------------------------------------
        Conv_zz(j, 1) = -fU(j, 1)*dRzzdx(j, 1) - fV(j, 1)*dRzzdy(j, 1) - fW(j, 1)*dRzzdz(j, 1)
        Prod_zz(j, 1) = -2.0_wp*(fRxz(j, 1)*fdWdx(j, 1) + fRyz(j, 1)*fdWdy(j, 1) + fRzz(j, 1)*fdWdz(j, 1))
        Diss_zz(j, 1) = -2.0_wp*(MA_TAUZkWk(j)*pts - &
                                 tau_xz(j, 1)*rdWdx(j, 1) - tau_yz(j, 1)*rdWdy(j, 1) - tau_zz(j, 1)*rdWdz(j, 1))/rR(j, 1)

! pressure terms with Reynolds average
        Tran_zz(j, 1) = -(T1zz(j, 1) + 2.0_wp*(-T4zz(j, 1) + (MA_PWz(j) + MA_WPz(j))*pts - &
                                               rP(j, 1)*rdWdz(j, 1) - rW(j, 1)*dPdz(j, 1)))/rR(j, 1)
        Pres_zz(j, 1) = 2.0_wp*(MA_PWz(j)*pts - rP(j, 1)*rdWdz(j, 1))/rR(j, 1)
        MnFl_zz(j, 1) = 2.0_wp*(rW(j, 1) - fW(j, 1))*(MA_TAUZkk(j)*pts - dPdz(j, 1))/rR(j, 1)

        Resi_zz(j, 1) = Conv_zz(j, 1) + Prod_zz(j, 1) + Diss_zz(j, 1) + Tran_zz(j, 1) + &
                        Pres_zz(j, 1) + MnFl_zz(j, 1)

! -------------------------------------------------------------------
! Rxy Reynolds stress equation
! -------------------------------------------------------------------
        Conv_xy(j, 1) = -fU(j, 1)*dRxydx(j, 1) - fV(j, 1)*dRxydy(j, 1) - fW(j, 1)*dRxydz(j, 1)
        Prod_xy(j, 1) = &
            -fRxx(j, 1)*fdVdx(j, 1) - fRxy(j, 1)*fdVdy(j, 1) - fRxz(j, 1)*fdVdz(j, 1) &
            - fRxy(j, 1)*fdUdx(j, 1) - fRyy(j, 1)*fdUdy(j, 1) - fRyz(j, 1)*fdUdz(j, 1)
        Diss_xy(j, 1) = -( &
                        MA_TAUXkVk(j)*pts - &
                        tau_xx(j, 1)*rdVdx(j, 1) - tau_xy(j, 1)*rdVdy(j, 1) - tau_xz(j, 1)*rdVdz(j, 1) + &
                        MA_TAUYkUk(j)*pts - &
                        tau_xy(j, 1)*rdUdx(j, 1) - tau_yy(j, 1)*rdUdy(j, 1) - tau_yz(j, 1)*rdUdz(j, 1))/rR(j, 1)

! pressure terms with Reynolds average
        Tran_xy(j, 1) = -(T1xy(j, 1) - T4xy(j, 1) - T4yx(j, 1) + &
                          (MA_PUy(j) + MA_UPy(j))*pts - rP(j, 1)*rdUdy(j, 1) - rU(j, 1)*dPdy(j, 1) + &
                          (MA_PVX(j) + MA_VPx(j))*pts - rP(j, 1)*rdVdx(j, 1) - rV(j, 1)*dPdx(j, 1))/rR(j, 1)
        Pres_xy(j, 1) = (MA_PUy(j)*pts - rP(j, 1)*rdUdy(j, 1) + &
                         MA_PVX(j)*pts - rP(j, 1)*rdVdx(j, 1))/rR(j, 1)
        MnFl_xy(j, 1) = ((rU(j, 1) - fU(j, 1))*(MA_TAUYkk(j)*pts - dPdy(j, 1)) + &
                         (rV(j, 1) - fV(j, 1))*(MA_TAUXkk(j)*pts - dPdx(j, 1)))/rR(j, 1)

        Resi_xy(j, 1) = Conv_xy(j, 1) + Prod_xy(j, 1) + Diss_xy(j, 1) + Tran_xy(j, 1) + &
                        Pres_xy(j, 1) + MnFl_xy(j, 1)

! -------------------------------------------------------------------
! Turbulent kinetic energy equation
! -------------------------------------------------------------------
        Conv(j, 1) = 0.5_wp*(Conv_xx(j, 1) + Conv_yy(j, 1) + Conv_zz(j, 1))
        Prod(j, 1) = 0.5_wp*(Prod_xx(j, 1) + Prod_yy(j, 1) + Prod_zz(j, 1))
        Diss(j, 1) = 0.5_wp*(Diss_xx(j, 1) + Diss_yy(j, 1) + Diss_zz(j, 1))
        Pres(j, 1) = 0.5_wp*(Pres_xx(j, 1) + Pres_yy(j, 1) + Pres_zz(j, 1))
        Tran(j, 1) = 0.5_wp*(Tran_xx(j, 1) + Tran_yy(j, 1) + Tran_zz(j, 1))
        MnFl(j, 1) = 0.5_wp*(MnFl_xx(j, 1) + MnFl_yy(j, 1) + MnFl_zz(j, 1))

        Resi(j, 1) = 0.5_wp*(Resi_xx(j, 1) + Resi_yy(j, 1) + Resi_zz(j, 1))

! -------------------------------------------------------------------
! Energy equation in terms of p
! -------------------------------------------------------------------
! pressure terms with Reynolds average
        Conv_p(j, 1) = -(fU(j, 1)*dPdx(j, 1) + fV(j, 1)*dPdy(j, 1) + fW(j, 1)*dPdz(j, 1))
        Reve_p(j, 1) = -gama0*rP(j, 1)*Dil(j, 1)
        Diss_p(j, 1) = (gama0 - 1.0_wp)*phi(j, 1)
        Tran_p(j, 1) = MA_Tkk(j)*pts*gama0*visc/prandtl
        Reyn_p(j, 1) = -((MA_UkPk(j) + MA_PUx(j) + MA_PVY(j) + MA_PWz(j))*pts - &
                         rP(j, 1)*Dil(j, 1) + Conv_p(j, 1))

        Resi_p(j, 1) = Conv_p(j, 1) + Reve_p(j, 1) + Diss_p(j, 1) + Tran_p(j, 1) + &
                       Reyn_p(j, 1) - (gama0 - 1)*rR(j, 1)*Pres(j, 1)

! -------------------------------------------------------------------
! Energy equation in terms of T
! -------------------------------------------------------------------
! using MRATIO*p/rho=T
        fdTdx = (MRATIO*dPdx(j, 1) - fT(j, 1)*dRdx(j, 1))/rR(j, 1)
        fdTdy = (MRATIO*dPdy(j, 1) - fT(j, 1)*dRdy(j, 1))/rR(j, 1)
        fdTdz = (MRATIO*dPdz(j, 1) - fT(j, 1)*dRdz(j, 1))/rR(j, 1)

        Conv_T(j, 1) = -(fU(j, 1)*fdTdx + fV(j, 1)*fdTdy + fW(j, 1)*fdTdz)
! dilatation-pressure terms with Reynolds average
        Reve_T(j, 1) = -MRATIO*(gama0 - 1)*rP(j, 1)*Dil(j, 1)/rR(j, 1)
        Diss_T(j, 1) = gama0*phi(j, 1)/rR(j, 1)
        Tran_T(j, 1) = MA_Tkk(j)*pts*gama0*visc/prandtl/rR(j, 1)
        Reyn_T(j, 1) = -(MRATIO*(MA_UkPk(j) + MA_PUx(j) + MA_PVY(j) + MA_PWz(j))*pts/rR(j, 1) + &
                         Conv_T(j, 1))

        Resi_T(j, 1) = Conv_T(j, 1) + Reve_T(j, 1) + Diss_T(j, 1) + Tran_T(j, 1) + &
                       Reyn_T(j, 1) - MRATIO*(gama0 - 1)*Pres(j, 1)

! -------------------------------------------------------------------
! Turbulent temperature equation
! -------------------------------------------------------------------
! !!! Not complete
        dfTdx = (MRATIO*dPdx(j, 1) - fT(j, 1)*dRdx(j, 1))/rR(j, 1)
        dfTdy = (MRATIO*dPdy(j, 1) - fT(j, 1)*dRdy(j, 1))/rR(j, 1)

        dRTTdx = MRATIO*(MA_PTx(j) + MA_TPx(j))*pts
        dRTTdy = MRATIO*(MA_PTy(j) + MA_TPy(j))*pts

        dfTf2dx = (dRTTdx - (fT(j, 1)*fT(j, 1) + fTf2(j, 1))*dRdx(j, 1))/rR(j, 1) - 2.0_wp*fT(j, 1)*dfTdx
        dfTf2dy = (dRTTdy - (fT(j, 1)*fT(j, 1) + fTf2(j, 1))*dRdy(j, 1))/rR(j, 1) - 2.0_wp*fT(j, 1)*dfTdy

        Conv_tt(j, 1) = -fU(j, 1)*dfTf2dx - fV(j, 1)*dfTf2dy
        Prod_tt(j, 1) = -2.0_wp*(fRuT(j, 1)*dfTdx + fRvT(j, 1)*dfTdy)

        dRUTdx = MRATIO*(MA_PUx(j) + MA_UPx(j))*pts
        dRVTdy = MRATIO*(MA_PVY(j) + MA_VPy(j))*pts

        tranttx = MA_RUTTx(j)*pts - &
                  fU(j, 1)*dRTTdx - &
                  rR(j, 1)*(fT(j, 1)**2 + fTf2(j, 1))* &
                  fdUdx(j, 1) - &
                  2.0_wp*fT(j, 1)*dRUTdx - &
                  2.0_wp*rR(j, 1)*(fU(j, 1)*fT(j, 1) + &
                                   fRuT(j, 1))*dfTdx + &
                  2.0_wp*fU(j, 1)*fT(j, 1)**2*dRdx(j, 1) + &
                  2.0_wp*rR(j, 1)*fT(j, 1)**2*fdUdx(j, 1) + &
                  4.0_wp*rR(j, 1)*fU(j, 1)*fT(j, 1)*dfTdx
        trantty = MA_RVTTy(j)*pts - &
                  fV(j, 1)*dRTTdy - &
                  rR(j, 1)*(fT(j, 1)**2 + fTf2(j, 1))* &
                  fdVdy(j, 1) - &
                  2.0_wp*fT(j, 1)*dRVTdy - &
                  2.0_wp*rR(j, 1)*(fV(j, 1)*fT(j, 1) + &
                                   fRvT(j, 1))*dfTdy + &
                  2.0_wp*fV(j, 1)*fT(j, 1)**2*dRdy(j, 1) + &
                  2.0_wp*rR(j, 1)*fT(j, 1)**2*fdVdy(j, 1) + &
                  4.0_wp*rR(j, 1)*fV(j, 1)*fT(j, 1)*dfTdy

        Tran_tt(j, 1) = -2.0_wp*(tranttx + trantty)

        Diss_tt(j, 1) = 0.0_wp
        Pres_tt(j, 1) = 0.0_wp
        MnFl_tt(j, 1) = 0.0_wp

        Resi_tt(j, 1) = Conv_tt(j, 1) + Prod_tt(j, 1) + Tran_tt(j, 1) + Diss_tt(j, 1) + MnFl_tt(j, 1)

! ###################################################################
! Variable density quantities
! ###################################################################
! speed of sound. Using Reynolds, not Favre, for T
        dum1 = rT(j, 1)/(mach*mach)
        dum2 = rT(j, 1)*(1.0_wp/rP(j, 1) - 1.0_wp/(rR(j, 1)*dum1))

        rho_p(j, 1) = MA_RP(j)*pts - rR(j, 1)*rP(j, 1)
        rho_T(j, 1) = MA_RT(j)*pts - rR(j, 1)*rT(j, 1)
! T-p correlation
        dum3 = MA_RTT(j)*pts/MRATIO - rT(j, 1)*rP(j, 1)

        rho_ac(j, 1) = rPf2(j, 1)/(dum1*dum1)
        rho_en(j, 1) = rRf2(j, 1) + rho_ac(j, 1) - 2.0_wp*rho_p(j, 1)/dum1

        T_ac(j, 1) = rPf2(j, 1)*dum2*dum2
        T_en(j, 1) = rTf2(j, 1) + T_ac(j, 1) - 2.0_wp*dum3*dum2

! ###################################################################
! Scales
! ###################################################################
        if (Diss(j, 1) == 0.0_wp) then
            eta(j, 1) = big_wp
            tau(j, 1) = big_wp
            lambda(j, 1) = big_wp
        else
            eta(j, 1) = ((visc/rR(j, 1))**3.0_wp/abs(Diss(j, 1)))**0.25_wp
            tau(j, 1) = sqrt(visc/(rR(j, 1)*abs(Diss(j, 1))))
            lambda(j, 1) = sqrt(10.0_wp*rTKE(j, 1)/(rR(j, 1)*abs(Diss(j, 1))/visc))
        end if

        if (rdUdxf2(j, 1) == 0.0_wp) then
            lambda_x(j, 1) = big_wp
        else
            lambda_x(j, 1) = sqrt(rUf2(j, 1)/rdUdxf2(j, 1))
        end if

        if (rdVdyf2(j, 1) == 0.0_wp) then
            lambda_y(j, 1) = big_wp
        else
            lambda_y(j, 1) = sqrt(rVf2(j, 1)/rdVdyf2(j, 1))
        end if

        if (rdWdzf2(j, 1) == 0.0_wp) then
            lambda_z(j, 1) = big_wp
        else
            lambda_z(j, 1) = sqrt(rWf2(j, 1)/rdWdzf2(j, 1))
        end if

! ###################################################################
! Skewness and flatness
! ###################################################################
        S_rho(j, 1) = MA_R3(j)*pts - rR(j, 1)**3.0_wp - 3.0_wp*rR(j, 1)*rRf2(j, 1)
        S_u(j, 1) = MA_U3(j)*pts - rU(j, 1)**3.0_wp - 3.0_wp*rU(j, 1)*rUf2(j, 1)
        S_v(j, 1) = MA_V3(j)*pts - rV(j, 1)**3.0_wp - 3.0_wp*rV(j, 1)*rVf2(j, 1)
        S_w(j, 1) = MA_W3(j)*pts - rW(j, 1)**3.0_wp - 3.0_wp*rW(j, 1)*rWf2(j, 1)
        S_p(j, 1) = MA_P3(j)*pts - rP(j, 1)**3.0_wp - 3.0_wp*rP(j, 1)*rPf2(j, 1)
        S_T(j, 1) = MA_T3(j)*pts - rT(j, 1)**3.0_wp - 3.0_wp*rT(j, 1)*rTf2(j, 1)

        F_rho(j, 1) = MA_R4(j)*pts - rR(j, 1)**4.0_wp - 4.0_wp*rR(j, 1)*S_rho(j, 1) - &
                      6.0_wp*rR(j, 1)**2.0_wp*rRf2(j, 1)
        F_u(j, 1) = MA_U4(j)*pts - rU(j, 1)**4.0_wp - 4.0_wp*rU(j, 1)*S_u(j, 1) - &
                    6.0_wp*rU(j, 1)**2.0_wp*rUf2(j, 1)
        F_v(j, 1) = MA_V4(j)*pts - rV(j, 1)**4.0_wp - 4.0_wp*rV(j, 1)*S_v(j, 1) - &
                    6.0_wp*rV(j, 1)**2.0_wp*rVf2(j, 1)
        F_w(j, 1) = MA_W4(j)*pts - rW(j, 1)**4.0_wp - 4.0_wp*rW(j, 1)*S_w(j, 1) - &
                    6.0_wp*rW(j, 1)**2.0_wp*rWf2(j, 1)
        F_p(j, 1) = MA_P4(j)*pts - rP(j, 1)**4.0_wp - 4.0_wp*rP(j, 1)*S_p(j, 1) - &
                    6.0_wp*rP(j, 1)**2.0_wp*rPf2(j, 1)
        F_T(j, 1) = MA_T4(j)*pts - rT(j, 1)**4.0_wp - 4.0_wp*rT(j, 1)*S_T(j, 1) - &
                    6.0_wp*rT(j, 1)**2.0_wp*rTf2(j, 1)

! Normalization
        if (rRf2(j, 1) == 0.0_wp) then
            S_rho(j, 1) = big_wp
            F_rho(j, 1) = big_wp
        else
            S_rho(j, 1) = S_rho(j, 1)/rRf2(j, 1)**(3.0_wp/2.0_wp)
            F_rho(j, 1) = F_rho(j, 1)/rRf2(j, 1)**2.0_wp
        end if
        if (rUf2(j, 1) == 0.0_wp) then
            S_u(j, 1) = big_wp
            F_u(j, 1) = big_wp
        else
            S_u(j, 1) = S_u(j, 1)/rUf2(j, 1)**(3.0_wp/2.0_wp)
            F_u(j, 1) = F_u(j, 1)/rUf2(j, 1)**2.0_wp
        end if
        if (rVf2(j, 1) == 0.0_wp) then
            S_v(j, 1) = big_wp
            F_v(j, 1) = big_wp
        else
            S_v(j, 1) = S_v(j, 1)/rVf2(j, 1)**(3.0_wp/2.0_wp)
            F_v(j, 1) = F_v(j, 1)/rVf2(j, 1)**2.0_wp
        end if
        if (rWf2(j, 1) == 0.0_wp) then
            S_w(j, 1) = big_wp
            F_w(j, 1) = big_wp
        else
            S_w(j, 1) = S_w(j, 1)/rWf2(j, 1)**(3.0_wp/2.0_wp)
            F_w(j, 1) = F_w(j, 1)/rWf2(j, 1)**2.0_wp
        end if
        if (rPf2(j, 1) == 0.0_wp) then
            S_p(j, 1) = big_wp
            F_p(j, 1) = big_wp
        else
            S_p(j, 1) = S_p(j, 1)/rPf2(j, 1)**(3.0_wp/2.0_wp)
            F_p(j, 1) = F_p(j, 1)/rPf2(j, 1)**2.0_wp
        end if
        if (rTf2(j, 1) == 0.0_wp) then
            S_T(j, 1) = big_wp
            F_T(j, 1) = big_wp
        else
            S_T(j, 1) = S_T(j, 1)/rTf2(j, 1)**(3.0_wp/2.0_wp)
            F_T(j, 1) = F_T(j, 1)/rTf2(j, 1)**2.0_wp
        end if

    end do

! ###################################################################
! Integral quantities shear layer
! ###################################################################
!   IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
!
! ! Vorticity Thickness
!      DO n = 1,nstatavg
!         DO j = 1,jmax
!            wrk1d(j,1) = fU(n,j)
!         ENDDO
!         CALL OPR_PARTIAL_Y(OPR_P1, i1,jmax,i1, bcs, g(2), wrk1d(1,1), wrk1d(1,2))
!         delta_w_u(n) = (fU(n,jmax)-fU(n,1)) / MINVAL(wrk1d(1:jmax,2))
!      ENDDO
!
! ! Momentum thickness
!      DO n = 1,nstatavg
!         UC = qbg(1)%mean
!         DU = qbg(1)%delta
!         DO j = jmin_loc, jmax_loc
!            wrk1d(j,1) = rR(n,j)*( 0.25_wp - ((fU(n,j)-UC)/DU)**2 )
!         ENDDO
!         delta_m_u(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
!      ENDDO
!
! ! Mixing layer limit (U=0.1dU and U=0.9dU)
!      y_center = g(2)%nodes(1) + qbg(1)%ymean_rel*g(2)%scale
!      DO n = 1,nstatavg
!         fU_05 = U2 + C_01_R*qbg(1)%delta
!         DO j = 1,jmax
!            IF ( fU(n,j) .GT. fU_05 .AND. fU(n,j+1) .LE. fU_05 ) THEN
!               delta_01_u(n) = g(2)%nodes(j) + (fU_05-fU(n,j))*(g(2)%nodes(j+1) -g(2)%nodes(j))/(fU(n,j+1)-fU(n,j))
!            ENDIF
!         ENDDO
!         delta_01_u(n) = delta_01_u(n) - y_center
!
!         fU_05 = U2 + r09*qbg(1)%delta
!         DO j = 1,jmax
!            IF ( fU(n,j) .GT. fU_05 .AND. fU(n,j+1) .LE. fU_05 ) THEN
!               delta_01_d(n) = g(2)%nodes(j) + (fU_05-fU(n,j))*(g(2)%nodes(j+1) -g(2)%nodes(j))/(fU(n,j+1)-fU(n,j))
!            ENDIF
!         ENDDO
!         delta_01_d(n) = delta_01_d(n) - y_center
!
!         delta_u_u(n) = 0.5_wp*( delta_01_u(n) + delta_01_d(n) )
!      ENDDO
!
! ! ###################################################################
! ! 1D quantities of the jet
! ! ###################################################################
!   ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
! -------------------------------------------------------------------
! Integral balance of mass
! -------------------------------------------------------------------
    do n = 1, nstatavg
! axial
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = rR(n, j)*fU(n, j)
        end do
        IntMassU(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
! lateral
        do k = 1, n
            i = statavg(k)
            wrk1d(k, 1) = rR(k, jmin_loc)*fV(k, jmin_loc) - rR(k, jmax_loc)*fV(k, jmax_loc)
            wrk1d(k, 2) = g(1)%nodes(i)
        end do
        if (n == 1) then
            IntMassV(n) = 0.0_wp
        else if (n == 2) then
            IntMassV(n) = 0.5_wp*(wrk1d(1, 1) + wrk1d(2, 1))*(wrk1d(2, 2) - wrk1d(1, 2))
        else
            IntMassV(n) = SIMPSON_NU(n, wrk1d(1, 1), wrk1d(1, 2))
        end if
    end do

! -------------------------------------------------------------------
! Integral balance of axial momentum (conserved)
! -------------------------------------------------------------------
    do n = 1, nstatavg
! mean velocity part
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = rR(n, j)*fU(n, j)*(fU(n, j) - U2)
        end do
        IntExcMomU(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
! pressure part
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = (rP(n, j) - pbg%mean)
        end do
        IntExcMomP(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
! Reynolds stress part
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = rR(n, j)*fRxx(n, j)
        end do
        IntExcMomRxx(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
    end do

! -------------------------------------------------------------------
! Integral balance of turbulent kinetic energy
! -------------------------------------------------------------------
    do n = 1, nstatavg
! TKE flux
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = rR(n, j)*fU(n, j)*fTKE(n, j)
        end do
        IntTkeK(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
! Integral of production term
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = rR(n, j)*Prod(n, j)
        end do
        IntTkeP(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
! Integral of numerical dissipation
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = -rR(n, j)*eps_f(n, j)
        end do
        IntTkeF(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
! Integral of production term
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = Pres(n, j)
        end do
        IntTkePi(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
    end do

! -------------------------------------------------------------------
! Axial flux of heat
! -------------------------------------------------------------------
    do n = 1, nstatavg
        do j = jmin_loc, jmax_loc
            wrk1d(j, 1) = rR(n, j)*fU(n, j)*(fT(n, j) - T2)
        end do
        IntFluxT(n) = SIMPSON_NU(nj, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))
    end do

! -------------------------------------------------------------------
! Jet thickness
! -------------------------------------------------------------------
! Vorticity thickness
    do n = 1, nstatavg
        do j = 1, jmax
            wrk1d(j, 1) = fU(n, j)
        end do
        call OPR_PARTIAL_Y(OPR_P1, 1, jmax, 1, bcs, g(2), wrk1d(1, 1), wrk1d(1, 2))
        delta_w_u(n) = (fU(n, jmax/2 + 1) - U2)/abs(minval(wrk1d(1:jmax, 2)))
        delta_w_d(n) = (fU(n, jmax/2) - U2)/abs(maxval(wrk1d(1:jmax, 2)))
    end do

! Momentum thickness
    do n = 1, nstatavg
        UC = 0.5_wp*(U2 + fU(n, jmax/2))
        DU = fU(n, jmax/2) - U2
        do j = jmin_loc, jmax/2
            wrk1d(j, 1) = rR(n, j)*(0.25_wp - ((fU(n, j) - UC)/DU)**2)
        end do
        delta_m_d(n) = SIMPSON_NU(jmax/2 - jmin_loc + 1, wrk1d(jmin_loc, 1), g(2)%nodes(jmin_loc))

        UC = 0.5_wp*(U2 + fU(n, jmax/2 + 1))
        DU = fU(n, jmax/2 + 1) - U2
        do j = jmax/2 + 1, jmax_loc
            wrk1d(j, 1) = rR(n, j)*(0.25_wp - ((fU(n, j) - UC)/DU)**2)
        end do
        delta_m_u(n) = SIMPSON_NU(jmax_loc - jmax/2, wrk1d(jmax/2 + 1, 1), g(2)%nodes(jmax/2 + 1))
    end do

! Jet half-width based on velocity
    call DELTA_X(nstatavg, jmax, g(2)%nodes, fU(1, 1), wrk1d(1, 1), &
                 delta_u_d(1), delta_u_u(1), U2, r05)

! Jet Limit (U=0.05Uc)
    call DELTA_X(nstatavg, jmax, g(2)%nodes, fU(1, 1), wrk1d(1, 1), &
                 delta_01_d(1), delta_01_u(1), U2, r005)

! Jet half-width based on temperature/density
    if (rbg%delta /= 0.0_wp) then
        do j = 1, jmax*nstatavg
            wrk2d(j, 1) = abs(fT(j, 1) - T2) + T2 ! we can have hot or cold jet
        end do
        call DELTA_X(nstatavg, jmax, g(2)%nodes, wrk2d(1, 1), wrk1d(1, 1), &
                     delta_t_d(1), delta_t_u(1), T2, r05)

        do j = 1, jmax*nstatavg
            wrk2d(j, 1) = abs(rR(j, 1) - R2) + R2 ! we can have hot or cold jet
        end do
        call DELTA_X(nstatavg, jmax, g(2)%nodes, wrk2d(1, 1), wrk1d(1, 1), &
                     delta_r_d(1), delta_r_u(1), R2, r05)
    else
        do n = 1, nstatavg
            delta_t_d(n) = 1.0_wp
            delta_t_u(n) = 1.0_wp
            delta_r_d(n) = 1.0_wp
            delta_r_u(n) = 1.0_wp
        end do
    end if

! Jet center line based on velocity
    y_center = g(2)%nodes(1) + qbg(1)%ymean_rel*g(2)%scale
    do n = 1, nstatavg
        do j = 1, jmax
            wrk1d(j, 1) = fU(n, j)
        end do
        jloc_max = maxloc(wrk1d(1:jmax, 1)); j = jloc_max(1)
        if (wrk1d(j - 1, 1) > wrk1d(j + 1, 1)) then
            delta_u_center(n) = 0.5_wp*(g(2)%nodes(j) + g(2)%nodes(j - 1))
        else
            delta_u_center(n) = 0.5_wp*(g(2)%nodes(j) + g(2)%nodes(j + 1))
        end if
        delta_u_center(n) = delta_u_center(n) - y_center
    end do

    ! ENDIF

! ###################################################################
! Scaling of the quatities
! ###################################################################
#define simuc(A) wrk1d(A,2)
#define simtc(A) wrk1d(A,3)
#define simrc(A) wrk1d(A,4)

    do n = 1, nstatavg

        ! IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        !    DU = qbg(1)%delta
        !    delta_05 = delta_01_u(n) - delta_01_d(n)
        ! ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
        delta_05 = 0.5_wp*(delta_u_u(n) + delta_u_d(n))

        simuc(n) = 0.5_wp*(fU(n, jmax/2) + fU(n, jmax/2 + 1)) - U2
        if (rbg%delta /= 0.0_wp) then
            simtc(n) = 0.5_wp*(fT(n, jmax/2) + fT(n, jmax/2 + 1)) - T2
            simrc(n) = 0.5_wp*(rR(n, jmax/2) + rR(n, jmax/2 + 1)) - R2
        else
            simtc(n) = 1.0_wp
            simrc(n) = 1.0_wp
        end if

        DU = simuc(n)
        DH = abs(simtc(n))
        ! ENDIF

! reynolds based on half-width
        Reynolds_d(n) = rR(n, jmax/2)*2.0_wp*delta_05*DU/visc
! reynolds based on isotropic lambda
        Reynolds_i(n) = rR(n, jmax/2)*lambda(n, jmax/2)*sqrt(2.0_wp*fTKE(n, jmax/2)/3.0_wp)/visc
! reynolds based on longitudinal lambda
        Reynolds_l(n) = rR(n, jmax/2)*lambda_x(n, jmax/2)*sqrt(fRxx(n, jmax/2))/visc

        do j = 1, jmax
            Vortx(n, j) = Vortx(n, j)/DU*delta_05
            Vorty(n, j) = Vorty(n, j)/DU*delta_05
            Vortz(n, j) = Vortz(n, j)/DU*delta_05
            Dil(n, j) = Dil(n, j)/DU*delta_05
            Vortxf2(n, j) = sqrt(Vortxf2(n, j))/DU*delta_05
            Vortyf2(n, j) = sqrt(Vortyf2(n, j))/DU*delta_05
            Vortzf2(n, j) = sqrt(Vortzf2(n, j))/DU*delta_05
            Dilf2(n, j) = Dilf2(n, j)/DU/DU*delta_05*delta_05

            Conv_xx(n, j) = Conv_xx(n, j)/(DU*DU*DU)*delta_05
            Prod_xx(n, j) = Prod_xx(n, j)/(DU*DU*DU)*delta_05
            Diss_xx(n, j) = Diss_xx(n, j)/(DU*DU*DU)*delta_05
            Tran_xx(n, j) = Tran_xx(n, j)/(DU*DU*DU)*delta_05
            Pres_xx(n, j) = Pres_xx(n, j)/(DU*DU*DU)*delta_05
            MnFl_xx(n, j) = MnFl_xx(n, j)/(DU*DU*DU)*delta_05
            Resi_xx(n, j) = Resi_xx(n, j)/(DU*DU*DU)*delta_05

            Conv_yy(n, j) = Conv_yy(n, j)/(DU*DU*DU)*delta_05
            Prod_yy(n, j) = Prod_yy(n, j)/(DU*DU*DU)*delta_05
            Diss_yy(n, j) = Diss_yy(n, j)/(DU*DU*DU)*delta_05
            Tran_yy(n, j) = Tran_yy(n, j)/(DU*DU*DU)*delta_05
            Pres_yy(n, j) = Pres_yy(n, j)/(DU*DU*DU)*delta_05
            MnFl_yy(n, j) = MnFl_yy(n, j)/(DU*DU*DU)*delta_05
            Resi_yy(n, j) = Resi_yy(n, j)/(DU*DU*DU)*delta_05

            Conv_zz(n, j) = Conv_zz(n, j)/(DU*DU*DU)*delta_05
            Prod_zz(n, j) = Prod_zz(n, j)/(DU*DU*DU)*delta_05
            Diss_zz(n, j) = Diss_zz(n, j)/(DU*DU*DU)*delta_05
            Tran_zz(n, j) = Tran_zz(n, j)/(DU*DU*DU)*delta_05
            Pres_zz(n, j) = Pres_zz(n, j)/(DU*DU*DU)*delta_05
            MnFl_zz(n, j) = MnFl_zz(n, j)/(DU*DU*DU)*delta_05
            Resi_zz(n, j) = Resi_zz(n, j)/(DU*DU*DU)*delta_05

            Conv_xy(n, j) = Conv_xy(n, j)/(DU*DU*DU)*delta_05
            Prod_xy(n, j) = Prod_xy(n, j)/(DU*DU*DU)*delta_05
            Diss_xy(n, j) = Diss_xy(n, j)/(DU*DU*DU)*delta_05
            Tran_xy(n, j) = Tran_xy(n, j)/(DU*DU*DU)*delta_05
            Pres_xy(n, j) = Pres_xy(n, j)/(DU*DU*DU)*delta_05
            MnFl_xy(n, j) = MnFl_xy(n, j)/(DU*DU*DU)*delta_05
            Resi_xy(n, j) = Resi_xy(n, j)/(DU*DU*DU)*delta_05

            fTKE(n, j) = fTKE(n, j)/(DU*DU)
            Conv(n, j) = Conv(n, j)/(DU*DU*DU)*delta_05
            Prod(n, j) = Prod(n, j)/(DU*DU*DU)*delta_05
            Diss(n, j) = Diss(n, j)/(DU*DU*DU)*delta_05
            Tran(n, j) = Tran(n, j)/(DU*DU*DU)*delta_05
            Pres(n, j) = Pres(n, j)/(DU*DU*DU)*delta_05
            MnFl(n, j) = MnFl(n, j)/(DU*DU*DU)*delta_05
            Resi(n, j) = Resi(n, j)/(DU*DU*DU)*delta_05

            equi(n, j) = fTKE(n, j)/abs(Diss(n, j))

            eps_f(n, j) = -eps_f(n, j)/(DU*DU*DU)*delta_05

            Conv_u(n, j) = Conv_u(n, j)/(DU*DU)*delta_05
            Tran_u(n, j) = Tran_u(n, j)/(DU*DU)*delta_05
            Reyn_u(n, j) = Reyn_u(n, j)/(DU*DU)*delta_05
            Resi_u(n, j) = Resi_u(n, j)/(DU*DU)*delta_05

            Conv_v(n, j) = Conv_v(n, j)/(DU*DU)*delta_05
            Tran_v(n, j) = Tran_v(n, j)/(DU*DU)*delta_05
            Reyn_v(n, j) = Reyn_v(n, j)/(DU*DU)*delta_05
            Resi_v(n, j) = Resi_v(n, j)/(DU*DU)*delta_05

            Conv_w(n, j) = Conv_w(n, j)/(DU*DU)*delta_05
            Tran_w(n, j) = Tran_w(n, j)/(DU*DU)*delta_05
            Reyn_w(n, j) = Reyn_w(n, j)/(DU*DU)*delta_05
            Resi_w(n, j) = Resi_w(n, j)/(DU*DU)*delta_05

            Conv_p(n, j) = Conv_p(n, j)/(DU*DU*DU)*delta_05
            Reve_p(n, j) = Reve_p(n, j)/(DU*DU*DU)*delta_05
            Diss_p(n, j) = Diss_p(n, j)/(DU*DU*DU)*delta_05
            Tran_p(n, j) = Tran_p(n, j)/(DU*DU*DU)*delta_05
            Reyn_p(n, j) = Reyn_p(n, j)/(DU*DU*DU)*delta_05
            Resi_p(n, j) = Resi_p(n, j)/(DU*DU*DU)*delta_05

            Conv_T(n, j) = Conv_T(n, j)/(DH*DU)*delta_05
            Reve_T(n, j) = Reve_T(n, j)/(DH*DU)*delta_05
            Diss_T(n, j) = Diss_T(n, j)/(DH*DU)*delta_05
            Tran_T(n, j) = Tran_T(n, j)/(DH*DU)*delta_05
            Reyn_T(n, j) = Reyn_T(n, j)/(DH*DU)*delta_05
            Resi_T(n, j) = Resi_T(n, j)/(DH*DU)*delta_05

            Conv_tt(n, j) = Conv_tt(n, j)/(DH*DH*DU)*delta_05
            Prod_tt(n, j) = Prod_tt(n, j)/(DH*DH*DU)*delta_05
            Diss_tt(n, j) = Diss_tt(n, j)/(DH*DH*DU)*delta_05
            Tran_tt(n, j) = Tran_tt(n, j)/(DH*DH*DU)*delta_05
            Pres_tt(n, j) = Pres_tt(n, j)/(DH*DH*DU)*delta_05
            MnFl_tt(n, j) = MnFl_tt(n, j)/(DH*DH*DU)*delta_05
            Resi_tt(n, j) = Resi_tt(n, j)/(DH*DH*DU)*delta_05

        end do

    end do

! ###################################################################
! Saving the data in TkStat format
! ###################################################################
    write (name, *) itime
    ! IF      ( imode_flow .EQ. DNS_FLOW_JET   ) THEN; name = 'jetavg'//TRIM(ADJUSTL(name))
    ! ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN; name = 'shravg'//TRIM(ADJUSTL(name)); ENDIF
    name = 'avg'//trim(adjustl(name))

#ifdef USE_RECLEN
    open (UNIT=i23, RECL=3260, FILE=name, STATUS='unknown')
#else
    open (UNIT=i23, FILE=name, STATUS='unknown')
#endif

! -------------------------------------------------------------------
! Header
! -------------------------------------------------------------------
    write (i23, '(A8,E14.7E3)') 'RTIME = ', rtime

! Independent variables
    line2 = 'I J X Y SU ST'

! Dependent variables depending on y and x
    line1 = 'Xg Yg'
    write (i23, 1010) 'GROUP = Grid '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'rU rV rW rP rR rT rUf2 rVf2 rWf2 rPf2 rRf2 rTf2 rUfVf rUfWf rVfWf rTKE ' &
            //'rbxx rbyy rbzz rbxy rbxz rbyz rRuT rRvT rRwT'
    write (i23, 1010) 'GROUP = Reynolds_Avgs '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'fU fV fW fT fTf2 fRxy fRxz fRyz fRxx fRyy fRzz ' &
            //'fbxx fbyy fbzz fbxy fbxz fbyz fRuT fRvT fRwT'
    write (i23, 1010) 'GROUP = Favre_Avgs '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'rdUdx rdUdy rdUdz rdVdx rdVdy rdVdz rdWdx rdWdy rdWdz ' &
            //'rdUdxf2 rdUdyf2 rdUdzf2 rdVdxf2 rdVdyf2 rdVdzf2 ' &
            //'rdWdxf2 rdWdyf2 rdWdzf2 ' &
            //'rdVdxfdUdyf rdWdxfdUdzf rdWdyfdVdzf ' &
            //'rdUdxfdVdyf rdUdxfdWdzf rdVdyfdWdzf ' &
            //'dPdx dPdy dPdz dRdx dRdy dRdz'
    write (i23, 1010) 'GROUP = Derivatives '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'Vortx Vorty Vortz Dil fDil Vortxf2 Vortyf2 Vortzf2 Dilf2'
    write (i23, 1010) 'GROUP = Vort_Dil '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'eta tau lambda lambda_x lambda_y lambda_z equi'
    write (i23, 1010) 'GROUP = Scales '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'Rxx Conv_xx Prod_xx Diss_xx Tran_xx Pres_xx MnFl_xx Resi_xx'
    write (i23, 1010) 'GROUP = Rxx_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'Ryy Conv_yy Prod_yy Diss_yy Tran_yy Pres_yy MnFl_yy Resi_yy'
    write (i23, 1010) 'GROUP = Ryy_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'Rzz Conv_zz Prod_zz Diss_zz Tran_zz Pres_zz MnFl_zz Resi_zz '
    write (i23, 1010) 'GROUP = Rzz_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'Rxy Conv_xy Prod_xy Diss_xy Tran_xy Pres_xy MnFl_xy Resi_xy'
    write (i23, 1010) 'GROUP = Rxy_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'TKE Conv Prod Diss Tran Pres MnFl Resi'
    write (i23, 1010) 'GROUP = TKE_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'Rtt Conv_tt Prod_tt Diss_tt Tran_tt Pres_tt MnFl_tt Resi_tt'
    write (i23, 1010) 'GROUP = Rtt_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'U Conv_u Tran_u Reyn_u Resi_u'
    write (i23, 1010) 'GROUP = U_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'V Conv_v Tran_v Reyn_v Resi_v'
    write (i23, 1010) 'GROUP = V_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'W Conv_w Tran_w Reyn_w Resi_w'
    write (i23, 1010) 'GROUP = W_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'p Conv_p Reve_p Diss_p Tran_p Reyn_p Pres_p Resi_p'
    write (i23, 1010) 'GROUP = p_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'T Conv_T Reve_T Diss_T Tran_T Reyn_T Pres_T Resi_T'
    write (i23, 1010) 'GROUP = T_Eqn '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'fTKE_nf eps_f'
    write (i23, 1010) 'GROUP = Filter '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'tau_xx tau_yy tau_zz tau_xy tau_xz tau_yz phi rVis'
    write (i23, 1010) 'GROUP = Mean_Stresses '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'Corr_RP Corr_RT R_ac R_en T_ac T_en RuT RvT RwT Rur Rvr Rwr'
    write (i23, 1010) 'GROUP = VarDensity '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    line1 = 'S_R S_U S_V S_W S_P S_T F_R F_U F_V F_W F_P F_T'
    write (i23, 1010) 'GROUP = Skewness_Flatness '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

! dependent variables dependent on t only
    ! IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
    line1 = 'Del_mom_u Del_mom_d Del_vor_u Del_vor_d ' &
            //'Del_half_u Del_half_d Del_lim_u Del_lim_d ' &
            //'Del_tem_u Del_tem_d Del_rho_u Del_rho_d Del_Umax ' &
            //'Sim_U Sim_T ' &
            //'Re_half Re_lambda_iso Re_lambda_lon ' &
            //'Int_mom_U Int_mom_P Int_mom_Rxx ' &
            //'Int_mass_U Int_mass_V Int_flux_T ' &
            //'Int_tke_K Int_tke_Pi Int_tke_P Int_tke_F'
    write (i23, 1010) 'GROUP = 1D_Quantities '//trim(adjustl(line1))
    line2 = trim(adjustl(line2))//' '//trim(adjustl(line1))

    ! ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
    !    line1 = 'Delta_m Delta_w y_01 y_09 y_05 Re Re_l'
    !    WRITE(i23,1010) 'GROUP = 1D_Quantities '//TRIM(ADJUSTL(line1))
    !    line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))
    ! ENDIF

    write (i23, 1010) trim(adjustl(line2))

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
    do n = 1, nstatavg
        i = statavg(n)

        ! IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        !    delta_05 = delta_01_u(n) - delta_01_d(n)
        !    delta_w  = delta_w_u(n)
        ! ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
        delta_05 = 0.5_wp*(delta_u_u(n) + delta_u_d(n))
        delta_w = 0.5_wp*(delta_w_u(n) + delta_w_d(n))
        delta_t = 0.5_wp*(delta_t_u(n) + delta_t_d(n))
        ! ENDIF

        ! IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
        ivauxpos = VARMX1D + 2
        VAUXPOS(1) = delta_m_u(n)
        VAUXPOS(2) = delta_m_d(n)
        VAUXPOS(3) = delta_w_u(n)
        VAUXPOS(4) = delta_w_d(n)
        VAUXPOS(5) = delta_u_u(n)
        VAUXPOS(6) = delta_u_d(n)
        VAUXPOS(7) = delta_01_u(n)
        VAUXPOS(8) = delta_01_d(n)
        VAUXPOS(9) = delta_t_u(n)
        VAUXPOS(10) = delta_t_d(n)
        VAUXPOS(11) = delta_r_u(n)
        VAUXPOS(12) = delta_r_d(n)
        VAUXPOS(13) = delta_u_center(n)
        VAUXPOS(14) = (simuc(1)/simuc(n))**2.0_wp
        VAUXPOS(15) = (simtc(1)/simtc(n))**2.0_wp
        VAUXPOS(16) = Reynolds_d(n)
        VAUXPOS(17) = Reynolds_i(n)
        VAUXPOS(18) = Reynolds_l(n)
        VAUXPOS(19) = IntExcMomU(n)
        VAUXPOS(20) = IntExcMomP(n)
        VAUXPOS(21) = IntExcMomRxx(n)
        VAUXPOS(22) = IntMassU(n)
        VAUXPOS(23) = IntMassV(n)
        VAUXPOS(24) = IntFluxT(n)
        VAUXPOS(25) = IntTkeK(n)
        VAUXPOS(26) = IntTkePi(n)
        VAUXPOS(27) = IntTkeP(n)
        VAUXPOS(28) = IntTkeF(n)
        ! ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        !    ivauxpos = 7
        !    VAUXPOS(1) = delta_m_u(n)
        !    VAUXPOS(2) = delta_w_u(n)
        !    VAUXPOS(3) = delta_01_u(n)
        !    VAUXPOS(4) = delta_01_d(n)
        !    VAUXPOS(5) = delta_u_u(n)
        !    VAUXPOS(6) = Reynolds_d(n)
        !    VAUXPOS(7) = Reynolds_l(n)
        ! ENDIF

        do j = 1, jmax
            ivauxpre = 4
            VAUXPRE(1) = g(1)%nodes(i)/qbg(1)%diam
            VAUXPRE(2) = g(2)%nodes(j)/qbg(1)%diam
            VAUXPRE(3) = (g(2)%nodes(j) - qbg(1)%ymean)/delta_05
            VAUXPRE(4) = (g(2)%nodes(j) - tbg%ymean)/delta_t

            if (j == jmax/2) then
                ivauxdum = ivauxpos
            else
                ivauxdum = 0
            end if

            write (i23, 1100) i, j, (VAUXPRE(k), k=1, ivauxpre), &
                ! Grid&
                g(1)%nodes(i), g(2)%nodes(j), &
                ! Reynolds averages&
                rU(n, j), rV(n, j), rW(n, j), &
                rP(n, j), rR(n, j), rT(n, j), &
                rUf2(n, j), rVf2(n, j), rWf2(n, j), &
                rPf2(n, j), rRf2(n, j), rTf2(n, j), &
                rUfVf(n, j), rUfWf(n, j), rVfWf(n, j), &
                rTKE(n, j), rbxx(n, j), rbyy(n, j), &
                rbzz(n, j), rbxy(n, j), rbxz(n, j), &
                rbyz(n, j), rRuT(n, j), rRvT(n, j), rRwT(n, j), &
                ! Favre averages&
                fU(n, j), fV(n, j), fW(n, j), fT(n, j), fTf2(n, j), &
                fRxy(n, j), fRxz(n, j), fRyz(n, j), &
                fRxx(n, j), fRyy(n, j), fRzz(n, j), &
                fbxx(n, j), fbyy(n, j), fbzz(n, j), &
                fbxy(n, j), fbxz(n, j), fbyz(n, j), &
                fRuT(n, j), fRvT(n, j), fRwT(n, j), &
                ! Derivatives&
                fdUdx(n, j), fdUdy(n, j), fdUdz(n, j), &
                fdVdx(n, j), fdVdy(n, j), fdVdz(n, j), &
                fdWdx(n, j), fdWdy(n, j), fdWdz(n, j), &
                rdUdxf2(n, j), rdUdyf2(n, j), rdUdzf2(n, j), &
                rdVdxf2(n, j), rdVdyf2(n, j), rdVdzf2(n, j), &
                rdWdxf2(n, j), rdWdyf2(n, j), rdWdzf2(n, j), &
                rdVdxfdUdyf(n, j), rdWdxfdUdzf(n, j), rdWdyfdVdzf(n, j), &
                rdUdxfdVdyf(n, j), rdUdxfdWdzf(n, j), rdVdyfdWdzf(n, j), &
                dPdx(n, j), dPdy(n, j), dPdz(n, j), &
                dRdx(n, j), dRdy(n, j), dRdz(n, j), &
                ! Vorticity $ Dilation&
                Vortx(n, j), Vorty(n, j), Vortz(n, j), Dil(n, j), &
                fdUdx(n, j) + fdVdy(n, j) + fdWdz(n, j), &
                Vortxf2(n, j), Vortyf2(n, j), Vortzf2(n, j), Dilf2(n, j), &
                ! Scales&
                eta(n, j), tau(n, j), lambda(n, j), &
                lambda_x(n, j), lambda_y(n, j), lambda_z(n, j), &
                equi(n, j), &
                ! Rxx equation&
                sqrt(fRxx(n, j))/simuc(n), &
                Conv_xx(n, j), Prod_xx(n, j), Diss_xx(n, j), &
                Tran_xx(n, j), Pres_xx(n, j), MnFl_xx(n, j), Resi_xx(n, j), &
                ! Ryy equation&
                sqrt(fRyy(n, j))/simuc(n), &
                Conv_yy(n, j), Prod_yy(n, j), Diss_yy(n, j), &
                Tran_yy(n, j), Pres_yy(n, j), MnFl_yy(n, j), Resi_yy(n, j), &
                ! Rzz equation&
                sqrt(fRzz(n, j))/simuc(n), &
                Conv_zz(n, j), Prod_zz(n, j), Diss_zz(n, j), &
                Tran_zz(n, j), Pres_zz(n, j), MnFl_zz(n, j), Resi_zz(n, j), &
                ! Rxy equation&
                fRxy(n, j)/simuc(n)/simuc(n), &
                Conv_xy(n, j), Prod_xy(n, j), Diss_xy(n, j), &
                Tran_xy(n, j), Pres_xy(n, j), MnFl_xy(n, j), Resi_xy(n, j), &
                ! TKE equation&
                fTKE(n, j), Conv(n, j), Prod(n, j), Diss(n, j), &
                Tran(n, j), Pres(n, j), MnFl(n, j), Resi(n, j), &
                ! temperate equation&
                sqrt(fTf2(n, j))/abs(simtc(n)), &
                Conv_tt(n, j), Prod_tt(n, j), &
                Diss_tt(n, j), Tran_tt(n, j), Pres_tt(n, j), &
                MnFl_tt(n, j), Resi_tt(n, j), &
                ! X-momentum equation&
                (fU(n, j) - U2)/simuc(n), Conv_u(n, j), Tran_u(n, j), &
                Reyn_u(n, j), Resi_u(n, j), &
                ! Y-momentum equation&
                fV(n, j)/simuc(n), Conv_v(n, j), Tran_v(n, j), &
                Reyn_v(n, j), Resi_v(n, j), &
                ! Z-momentum equation&
                fW(n, j)/simuc(n), Conv_w(n, j), Tran_w(n, j), &
                Reyn_w(n, j), Resi_w(n, j), &
                ! energy equation in p&
                (rP(n, j) - rP(n, 1))/ &
                (rP(n, jmax/2) - rP(n, 1)), &
                Conv_p(n, j), Reve_p(n, j), Diss_p(n, j), &
                Tran_p(n, j), Reyn_p(n, j), &
                -(gama0 - 1)*rR(n, j)*Pres(n, j), Resi_p(n, j), &
                ! energy equation in T&
                (fT(n, j) - T2)/abs(simtc(n)), &
                Conv_T(n, j), Reve_T(n, j), Diss_T(n, j), &
                Tran_T(n, j), Reyn_T(n, j), &
                -MRATIO*(gama0 - 1)*Pres(n, j)* &
                simuc(n)*simuc(n)/abs(simtc(n)), Resi_T(n, j), &
                ! Filtering&
                fTKE_nf(n, j), eps_f(n, j), &
                ! Stress tensor&
                tau_xx(n, j), tau_yy(n, j), tau_zz(n, j), &
                tau_xy(n, j), tau_xz(n, j), tau_yz(n, j), &
                phi(n, j), rVis(n, j), &
                ! Variable density quantities&
                rho_p(n, j), rho_T(n, j), rho_ac(n, j), rho_en(n, j), &
                T_ac(n, j), T_en(n, j), &
                fRuT(n, j)/abs(simtc(n)*simuc(n)), &
                fRvT(n, j)/abs(simtc(n)*simuc(n)), &
                fRwT(n, j)/abs(simtc(n)*simuc(n)), &
                (fU(n, j) - rU(n, j))*rR(n, j)/ &
                abs(simrc(n)*simuc(n)), &
                (fV(n, j) - rV(n, j))*rR(n, j)/ &
                abs(simrc(n)*simuc(n)), &
                (fW(n, j) - rW(n, j))*rR(n, j)/ &
                abs(simrc(n)*simuc(n)), &
                ! Skewness&Flatness&
                S_rho(n, j), S_u(n, j), S_v(n, j), S_w(n, j), &
                S_p(n, j), S_T(n, j), F_rho(n, j), F_u(n, j), &
                F_v(n, j), F_w(n, j), F_p(n, j), F_T(n, j), &
                ! 1D quantities&
                (VAUXPOS(k), k=1, ivauxdum)

        end do
    end do

    close (i23)

#undef simuc
#undef simtc
#undef simrc

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING AVG_FLOW_SPATIAL_LAYER')
#endif

    return

1010 format(A)
1100 format(I3, 1x, I3, 4(1x, G_FORMAT_R), 206(1x, G_FORMAT_R), VARMX1D(1x, G_FORMAT_R), 2(1x, G_FORMAT_R))

end subroutine AVG_FLOW_SPATIAL_LAYER
