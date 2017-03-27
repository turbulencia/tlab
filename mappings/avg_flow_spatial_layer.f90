#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "avgij_map.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2000/09/25 - J.P. Mellado
!#              Created
!# 2008/01/11 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Post-processing statistical data accumulated in mean1d. Based on 
!# mappings define in file avgij_map.h
!#
!########################################################################
!# ARGUMENTS 
!#
!# itxc    In     size of array stat, containing postprocess data
!# mean1d  In     array with raw mean data
!# stat           aux local array with postprocess (final) mean data
!#
!########################################################################
SUBROUTINE AVG_FLOW_SPATIAL_LAYER(itxc, jmin_loc,jmax_loc, mean1d, stat, wrk1d,wrk2d)

  USE DNS_CONSTANTS, ONLY : efile, tfile
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0, MRATIO
  IMPLICIT NONE

#include "integers.h"

  TINTEGER itxc, jmin_loc, jmax_loc
  TREAL, DIMENSION(nstatavg,jmax,*)  :: mean1d, stat

  TREAL wrk1d(isize_wrk1d,*)
  TREAL wrk2d(isize_wrk2d,*)

  TREAL dy(1) ! To use old wrappers to calculate derivatives

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

  TINTEGER i, j, k, n
  TREAL pts, c13, zero
  TREAL dum1, dum2, dum3, dum4, dum5
  TREAL SIMPSON_NU
  TREAL U2, DU, UC, fU_05, r05, r005, r09, T2, DH, R2
  TREAL y_center, dt_mean
  TREAL delta_05, delta_w, delta_t
  TREAL dfTdx, dfTdy, dRTTdx, dRTTdy, dfTf2dx, dfTf2dy
  TREAL tranttx, trantty, dRUTdx, dRVTdy
  TREAL fdTdx, fdTdy, fdTdz
  TINTEGER nj, jloc_max(1)

  TREAL  VAUXPRE(4), VAUXPOS(28)
  TINTEGER ivauxpre, ivauxpos, ivauxdum
  CHARACTER*32 name
  CHARACTER*400 line1
  CHARACTER*2750 line2

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING AVG_FLOW_SPATIAL_LAYER' )
#endif

  r05  = C_05_R
  r005 = C_5_R*C_1EM2_R
  r09  = C_9_R/C_10_R
  c13  = C_1_R/C_3_R
  zero = C_1EM6_R

  if ( nstatavg_points .EQ. 0 ) then
     CALL IO_WRITE_ASCII(efile,'AVG_FLOW_SPATIAL_LAYER: Zero number of points')
     CALL DNS_STOP(DNS_ERROR_STATZERO)
  ELSE
     pts = C_1_R/M_REAL(nstatavg_points)
  endif

  dt_mean = (rtime-rstattimeorg)/M_REAL(itime-istattimeorg)
  U2 = mean_u   - C_05_R*delta_u
  T2 = mean_tem - C_05_R*delta_tem
  R2 = rbg%mean - C_05_R*rbg%delta

  IF ( itxc .LT. nstatavg*jmax*LAST_INDEX ) THEN
     CALL IO_WRITE_ASCII(efile,'AVG_FLOW_SPATIAL_LAYER: Not enough space in stat')
     CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

  nj = jmax_loc-jmin_loc+1

! ###################################################################
! Main loop
! ###################################################################
  DO j = 1,jmax*nstatavg

! ###################################################################
! Reynolds Averages
! ###################################################################
     rU(j,1) = MA_U(j)*pts
     rV(j,1) = MA_V(j)*pts
     rW(j,1) = MA_W(j)*pts
     rP(j,1) = MA_P(j)*pts
     rR(j,1) = MA_R(j)*pts
     rT(j,1) = MA_T(j)*pts

     rUf2(j,1) = MA_UU(j)*pts - rU(j,1)*rU(j,1)
     rVf2(j,1) = MA_VV(j)*pts - rV(j,1)*rV(j,1)
     rWf2(j,1) = MA_WW(j)*pts - rW(j,1)*rW(j,1)
     rUfVf(j,1)= MA_UV(j)*pts - rU(j,1)*rV(j,1)
     rUfWf(j,1)= MA_UW(j)*pts - rU(j,1)*rW(j,1)
     rVfWf(j,1)= MA_VW(j)*pts - rV(j,1)*rW(j,1)
     rTKE(j,1) = C_05_R*( rUf2(j,1) + rVf2(j,1) + rWf2(j,1) )
     rbxx(j,1) = C_05_R*rUf2(j,1)/rTKE(j,1) - c13  
     rbyy(j,1) = C_05_R*rVf2(j,1)/rTKE(j,1) - c13 
     rbzz(j,1) = C_05_R*rWf2(j,1)/rTKE(j,1) - c13
     rbxy(j,1) = C_05_R*rUfVf(j,1)/rTKE(j,1)
     rbxz(j,1) = C_05_R*rUfWf(j,1)/rTKE(j,1)
     rbyz(j,1) = C_05_R*rVfWf(j,1)/rTKE(j,1)

     rPf2(j,1) = MA_PP(j)*pts - rP(j,1)*rP(j,1)
     rRf2(j,1) = MA_RR(j)*pts - rR(j,1)*rR(j,1)
     rTf2(j,1) = MA_TT(j)*pts - rT(j,1)*rT(j,1)

     rRuT(j,1) = MA_TU(j)*pts - rT(j,1)*rU(j,1)
     rRvT(j,1) = MA_TV(j)*pts - rT(j,1)*rV(j,1)
     rRwT(j,1) = MA_TW(j)*pts - rT(j,1)*rW(j,1)

! ###################################################################
! Favre Averages
! ###################################################################
     fU(j,1) = MA_RU(j)*pts/rR(j,1)
     fV(j,1) = MA_RV(j)*pts/rR(j,1)
     fW(j,1) = MA_RW(j)*pts/rR(j,1)
     fT(j,1) = MA_RT(j)*pts/rR(j,1)

     fRxx(j,1) = MA_RUU(j)*pts/rR(j,1) - fU(j,1)*fU(j,1)
     fRyy(j,1) = MA_RVV(j)*pts/rR(j,1) - fV(j,1)*fV(j,1)
     fRzz(j,1) = MA_RWW(j)*pts/rR(j,1) - fW(j,1)*fW(j,1)
     fRxy(j,1) = MA_RUV(j)*pts/rR(j,1) - fU(j,1)*fV(j,1)
     fRxz(j,1) = MA_RUW(j)*pts/rR(j,1) - fU(j,1)*fW(j,1)
     fRyz(j,1) = MA_RVW(j)*pts/rR(j,1) - fV(j,1)*fW(j,1)
     fTKE(j,1) = C_05_R*( fRxx(j,1) + fRyy(j,1) + fRzz(j,1) )
     fbxx(j,1) = C_05_R*fRxx(j,1)/fTKE(j,1) - c13  
     fbyy(j,1) = C_05_R*fRyy(j,1)/fTKE(j,1) - c13 
     fbzz(j,1) = C_05_R*fRzz(j,1)/fTKE(j,1) - c13
     fbxy(j,1) = C_05_R*fRxy(j,1)/fTKE(j,1)
     fbxz(j,1) = C_05_R*fRxz(j,1)/fTKE(j,1)
     fbyz(j,1) = C_05_R*fRyz(j,1)/fTKE(j,1)

     fTf2(j,1) = MA_RTT(j)*pts/rR(j,1) - fT(j,1)*fT(j,1)

     fRuT(j,1)= MRATIO*MA_PU(j)*pts/rR(j,1) - fU(j,1)*fT(j,1)
     fRvT(j,1)= MRATIO*MA_PV(j)*pts/rR(j,1) - fV(j,1)*fT(j,1)
     fRwT(j,1)= MRATIO*MA_PW(j)*pts/rR(j,1) - fW(j,1)*fT(j,1)

! the TKE before filtering is stored every iteration
     dum1 = C_1_R/M_REAL((itime-istattimeorg)*kmax_total)
     fTKE_nf(j,1) = C_05_R*( MA_FLT_RUU(j) + MA_FLT_RVV(j) + MA_FLT_RWW(j) -&
          ( MA_FLT_RU(j)*MA_FLT_RU(j) &
          + MA_FLT_RV(j)*MA_FLT_RV(j)&
          + MA_FLT_RW(j)*MA_FLT_RW(j) )*dum1/rR(j,1) )*dum1/rR(j,1)   

!     eps_f(j,1) = (fTKE_nf(j,1)-fTKE(j,1))/dt_mean/M_REAL(ifilt_step) ! to be rechecked

! ###################################################################
! First derivative terms (Reynolds)
! ###################################################################
     dRdx(j,1) = MA_Rx(j)*pts
     dRdy(j,1) = MA_Ry(j)*pts
     dRdz(j,1) = MA_Rz(j)*pts
     dPdx(j,1) = MA_Px(j)*pts
     dPdy(j,1) = MA_Py(j)*pts
     dPdz(j,1) = MA_Pz(j)*pts

! velocities
     rdUdx(j,1) = MA_Ux(j)*pts
     rdUdy(j,1) = MA_Uy(j)*pts
     rdUdz(j,1) = MA_Uz(j)*pts
     rdVdx(j,1) = MA_Vx(j)*pts
     rdVdy(j,1) = MA_Vy(j)*pts
     rdVdz(j,1) = MA_Vz(j)*pts
     rdWdx(j,1) = MA_Wx(j)*pts
     rdWdy(j,1) = MA_Wy(j)*pts
     rdWdz(j,1) = MA_Wz(j)*pts

     rdUdxf2(j,1) = MA_Ux2(j)*pts - rdUdx(j,1)*rdUdx(j,1)
     rdUdyf2(j,1) = MA_Uy2(j)*pts - rdUdy(j,1)*rdUdy(j,1)
     rdUdzf2(j,1) = MA_Uz2(j)*pts - rdUdz(j,1)*rdUdz(j,1)
     rdVdxf2(j,1) = MA_Vx2(j)*pts - rdVdx(j,1)*rdVdx(j,1)
     rdVdyf2(j,1) = MA_Vy2(j)*pts - rdVdy(j,1)*rdVdy(j,1)
     rdVdzf2(j,1) = MA_Vz2(j)*pts - rdVdz(j,1)*rdVdz(j,1)
     rdWdxf2(j,1) = MA_Wx2(j)*pts - rdWdx(j,1)*rdWdx(j,1)
     rdWdyf2(j,1) = MA_Wy2(j)*pts - rdWdy(j,1)*rdWdy(j,1)
     rdWdzf2(j,1) = MA_Wz2(j)*pts - rdWdz(j,1)*rdWdz(j,1)

     rdVdxfdUdyf(j,1) = MA_VxUy(j)*pts - rdVdx(j,1)*rdUdy(j,1)
     rdWdxfdUdzf(j,1) = MA_WxUz(j)*pts - rdWdx(j,1)*rdUdz(j,1)
     rdWdyfdVdzf(j,1) = MA_WyVz(j)*pts - rdWdy(j,1)*rdVdz(j,1)
     rdUdxfdVdyf(j,1) = MA_UXVY(j)*pts - rdUdx(j,1)*rdVdy(j,1)
     rdUdxfdWdzf(j,1) = MA_UxWz(j)*pts - rdUdx(j,1)*rdWdz(j,1)
     rdVdyfdWdzf(j,1) = MA_VyWz(j)*pts - rdVdy(j,1)*rdWdz(j,1)

! vorticity and dilatation
     Vortx(j,1) = rdWdy(j,1) - rdVdz(j,1)
     Vorty(j,1) = rdUdz(j,1) - rdWdx(j,1)
     Vortz(j,1) = rdVdx(j,1) - rdUdy(j,1)
     Dil(j,1)   = rdUdx(j,1) + rdVdy(j,1) + rdWdz(j,1)

     Vortxf2(j,1) = rdWdyf2(j,1) + rdVdzf2(j,1) - C_2_R*rdWdyfdVdzf(j,1)
     Vortyf2(j,1) = rdUdzf2(j,1) + rdWdxf2(j,1) - C_2_R*rdWdxfdUdzf(j,1)
     Vortzf2(j,1) = rdVdxf2(j,1) + rdUdyf2(j,1) - C_2_R*rdVdxfdUdyf(j,1)

     Dilf2(j,1) =  rdUdxf2(j,1) + rdVdyf2(j,1) + rdWdzf2(j,1) + &
          C_2_R*(rdUdxfdVdyf(j,1) + rdUdxfdWdzf(j,1) + rdVdyfdWdzf(j,1))

! ###################################################################
! First derivative terms (Favre)
! ###################################################################
     fdUdx(j,1) = (MA_RUx(j) + MA_URx(j) - fU(j,1)*MA_Rx(j))*pts/rR(j,1)
     fdUdy(j,1) = (MA_RUy(j) + MA_URy(j) - fU(j,1)*MA_Ry(j))*pts/rR(j,1)
     fdUdz(j,1) = (MA_RUz(j) + MA_URz(j) - fU(j,1)*MA_Rz(j))*pts/rR(j,1)
     fdVdx(j,1) = (MA_RVx(j) + MA_VRx(j) - fV(j,1)*MA_Rx(j))*pts/rR(j,1)
     fdVdy(j,1) = (MA_RVy(j) + MA_VRy(j) - fV(j,1)*MA_Ry(j))*pts/rR(j,1)
     fdVdz(j,1) = (MA_RVz(j) + MA_VRz(j) - fV(j,1)*MA_Rz(j))*pts/rR(j,1)
     fdWdx(j,1) = (MA_RWx(j) + MA_WRx(j) - fW(j,1)*MA_Rx(j))*pts/rR(j,1)
     fdWdy(j,1) = (MA_RWy(j) + MA_WRy(j) - fW(j,1)*MA_Ry(j))*pts/rR(j,1)
     fdWdz(j,1) = (MA_RWz(j) + MA_WRz(j) - fW(j,1)*MA_Rz(j))*pts/rR(j,1)

! ###################################################################
! Derivatives of the Reynolds stresses
! ###################################################################
     dRxxdx(j,1) = (MA_RUUx(j) - MA_RUU(j)/rR(j,1)*dRdx(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdUdx(j,1) - fU(j,1)*fdUdx(j,1)
     dRxxdy(j,1) = (MA_RUUy(j) - MA_RUU(j)/rR(j,1)*dRdy(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdUdy(j,1) - fU(j,1)*fdUdy(j,1)
     dRxxdz(j,1) = (MA_RUUz(j) - MA_RUU(j)/rR(j,1)*dRdz(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdUdz(j,1) - fU(j,1)*fdUdz(j,1)

     dRyydx(j,1) = (MA_RVVx(j) - MA_RVV(j)/rR(j,1)*dRdx(j,1))*pts/rR(j,1) -&
          fV(j,1)*fdVdx(j,1) - fV(j,1)*fdVdx(j,1)
     dRyydy(j,1) = (MA_RVVy(j) - MA_RVV(j)/rR(j,1)*dRdy(j,1))*pts/rR(j,1) -&
          fV(j,1)*fdVdy(j,1) - fV(j,1)*fdVdy(j,1)
     dRyydz(j,1) = (MA_RVVz(j) - MA_RVV(j)/rR(j,1)*dRdz(j,1))*pts/rR(j,1) -&
          fV(j,1)*fdVdz(j,1) - fV(j,1)*fdVdz(j,1)

     dRzzdx(j,1) = (MA_RWWx(j) - MA_RWW(j)/rR(j,1)*dRdx(j,1))*pts/rR(j,1) -&
          fW(j,1)*fdWdx(j,1) - fW(j,1)*fdWdx(j,1)
     dRzzdy(j,1) = (MA_RWWy(j) - MA_RWW(j)/rR(j,1)*dRdy(j,1))*pts/rR(j,1) -&
          fW(j,1)*fdWdy(j,1) - fW(j,1)*fdWdy(j,1)
     dRzzdz(j,1) = (MA_RWWz(j) - MA_RWW(j)/rR(j,1)*dRdz(j,1))*pts/rR(j,1) -&
          fW(j,1)*fdWdz(j,1) - fW(j,1)*fdWdz(j,1)

     dRxydx(j,1) = (MA_RUVx(j) - MA_RUV(j)/rR(j,1)*dRdx(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdVdx(j,1) - fV(j,1)*fdUdx(j,1)
     dRxydy(j,1) = (MA_RUVy(j) - MA_RUV(j)/rR(j,1)*dRdy(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdVdy(j,1) - fV(j,1)*fdUdy(j,1)
     dRxydz(j,1) = (MA_RUVz(j) - MA_RUV(j)/rR(j,1)*dRdz(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdVdz(j,1) - fV(j,1)*fdUdz(j,1)

     dRxzdx(j,1) = (MA_RUWx(j) - MA_RUW(j)/rR(j,1)*dRdx(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdWdx(j,1) - fW(j,1)*fdUdx(j,1)
     dRxzdy(j,1) = (MA_RUWy(j) - MA_RUW(j)/rR(j,1)*dRdy(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdWdy(j,1) - fW(j,1)*fdUdy(j,1)
     dRxzdz(j,1) = (MA_RUWz(j) - MA_RUW(j)/rR(j,1)*dRdz(j,1))*pts/rR(j,1) -&
          fU(j,1)*fdWdz(j,1) - fW(j,1)*fdUdz(j,1)

     dRyzdx(j,1) = (MA_RVWx(j) - MA_RVW(j)/rR(j,1)*dRdx(j,1))*pts/rR(j,1) -&
          fV(j,1)*fdWdx(j,1) - fW(j,1)*fdVdx(j,1)
     dRyzdy(j,1) = (MA_RVWy(j) - MA_RVW(j)/rR(j,1)*dRdy(j,1))*pts/rR(j,1) -&
          fV(j,1)*fdWdy(j,1) - fW(j,1)*fdVdy(j,1)
     dRyzdz(j,1) = (MA_RVWz(j) - MA_RVW(j)/rR(j,1)*dRdz(j,1))*pts/rR(j,1) -&
          fV(j,1)*fdWdz(j,1) - fW(j,1)*fdVdz(j,1)

! ###################################################################
! Viscous shear stress tensor 
! ###################################################################
     rVis(j,1) = MA_VIS(j)*pts

     tau_xx(j,1) = MA_TAUxx(j)*pts
     tau_yy(j,1) = MA_TAUyy(j)*pts
     tau_zz(j,1) = MA_TAUzz(j)*pts
     tau_xy(j,1) = MA_TAUxy(j)*pts
     tau_xz(j,1) = MA_TAUxz(j)*pts
     tau_yz(j,1) = MA_TAUyz(j)*pts

     phi(j,1) = (MA_TAUXkUk(j) + MA_TAUYkVk(j) + MA_TAUZkWk(j))*pts

! ###################################################################
! Transport equations
! ###################################################################
! -------------------------------------------------------------------
! auxiliar quantities
! -------------------------------------------------------------------
     dum1 = fU(j,1)*dRdx(j,1) + fV(j,1)*dRdy(j,1) + fW(j,1)*dRdz(j,1)

     dum2 = fU(j,1)*fdUdx(j,1) + fV(j,1)*fdUdy(j,1) + fW(j,1)*fdUdz(j,1)
     dum3 = fU(j,1)*fdVdx(j,1) + fV(j,1)*fdVdy(j,1) + fW(j,1)*fdVdz(j,1)
     dum4 = fU(j,1)*fdWdx(j,1) + fV(j,1)*fdWdy(j,1) + fW(j,1)*fdWdz(j,1)

     dum5 = fdUdx(j,1) + fdVdy(j,1) + fdWdz(j,1)

! -------------------------------------------------------------------
! X-, Y-, and Z-Momentum equation
! -------------------------------------------------------------------
     Conv_u(j,1) =-dum2 
     Tran_u(j,1) = (-dPdx(j,1)+MA_TAUXkk(j)*pts)/rR(j,1)
     Reyn_u(j,1) =-dRxxdx(j,1)-dRxydy(j,1)-dRxzdz(j,1)-&
          (fRxx(j,1)*dRdx(j,1) + fRxy(j,1)*dRdy(j,1) + fRxz(j,1)*dRdz(j,1))/rR(j,1) 
     Resi_u(j,1) = Conv_u(j,1) + Tran_u(j,1) + Reyn_u(j,1)

     Conv_v(j,1) =-dum3
     Tran_v(j,1) = (-dPdy(j,1)+MA_TAUYkk(j)*pts)/rR(j,1)
     Reyn_v(j,1) =-dRxydx(j,1)-dRyydy(j,1)-dRyzdz(j,1)-&
          (fRxy(j,1)*dRdx(j,1) + fRyy(j,1)*dRdy(j,1) + fRyz(j,1)*dRdz(j,1))/rR(j,1) 
     Resi_v(j,1) = Conv_v(j,1) + Tran_v(j,1) + Reyn_v(j,1)

     Conv_w(j,1) =-dum4
     Tran_w(j,1) = (-dPdz(j,1)+MA_TAUZkk(j)*pts)/rR(j,1)
     Reyn_w(j,1) =-dRxzdx(j,1)-dRyzdy(j,1)-dRzzdz(j,1)-&
          (fRxz(j,1)*dRdx(j,1) + fRyz(j,1)*dRdy(j,1) + fRzz(j,1)*dRdz(j,1))/rR(j,1) 
     Resi_w(j,1) = Conv_w(j,1) + Tran_w(j,1) + Reyn_w(j,1)

! -------------------------------------------------------------------
! Convective element of transport term of Reynolds equations
! -------------------------------------------------------------------
     T1xx(j,1) = ( MA_RUUUkk(j) - MA_RUU(j)*dum5   -&
          MA_RUU(j)*fdUdx(j,1) - MA_RUV(j)*fdUdy(j,1) - MA_RUW(j)*fdUdz(j,1) -&
          MA_RUU(j)*fdUdx(j,1) - MA_RUV(j)*fdUdy(j,1) - MA_RUW(j)*fdUdz(j,1) -&
          (MA_RUUx(j)+MA_RUVy(j)+MA_RUWz(j))*fU(j,1) -&
          (MA_RUUx(j)+MA_RUVy(j)+MA_RUWz(j))*fU(j,1) -&
          MA_RUUx(j)*fU(j,1) - MA_RUUy(j)*fV(j,1) - MA_RUUz(j)*fW(j,1) )*pts +&
          C_2_R*( fU(j,1)*fU(j,1)*dum1 + &
          rR(j,1)*(fU(j,1)*fU(j,1)*dum5 + fU(j,1)*dum2 + fU(j,1)*dum2) )

     T1yy(j,1) = (MA_RVVUkk(j) - MA_RVV(j)*dum5   -&
          MA_RUV(j)*fdVdx(j,1) - MA_RVV(j)*fdVdy(j,1) - MA_RVW(j)*fdVdz(j,1) - &
          MA_RUV(j)*fdVdx(j,1) - MA_RVV(j)*fdVdy(j,1) - MA_RVW(j)*fdVdz(j,1) -&
          (MA_RUVx(j)+MA_RVVy(j)+MA_RVWz(j))*fV(j,1)-&
          (MA_RUVx(j)+MA_RVVy(j)+MA_RVWz(j))*fV(j,1)-&
          MA_RVVx(j)*fU(j,1) - MA_RVVy(j)*fV(j,1) - MA_RVVz(j)*fW(j,1) )*pts +&
          C_2_R*( fV(j,1)*fV(j,1)*dum1 +&
          rR(j,1)*(fV(j,1)*fV(j,1)*dum5 + fV(j,1)*dum3 + fV(j,1)*dum3) ) 

     T1zz(j,1) = (MA_RWWUkk(j) - MA_RWW(j)*dum5   -&
          MA_RUW(j)*fdWdx(j,1) - MA_RVW(j)*fdWdy(j,1) - MA_RWW(j)*fdWdz(j,1) -&
          MA_RUW(j)*fdWdx(j,1) - MA_RVW(j)*fdWdy(j,1) - MA_RWW(j)*fdWdz(j,1) -&
          (MA_RUWx(j)+MA_RVWy(j)+MA_RWWz(j))*fW(j,1)-&
          (MA_RUWx(j)+MA_RVWy(j)+MA_RWWz(j))*fW(j,1)-&
          MA_RWWx(j)*fU(j,1) - MA_RWWy(j)*fV(j,1) - MA_RWWz(j)*fW(j,1) )*pts +&
          C_2_R*( fW(j,1)*fW(j,1)*dum1 +&
          rR(j,1)*(fW(j,1)*fW(j,1)*dum5 + fW(j,1)*dum4 + fW(j,1)*dum4) )   

     T1xy(j,1) = (MA_RUVUkk(j) - MA_RUV(j)*dum5   -&
          MA_RUU(j)*fdVdx(j,1) - MA_RUV(j)*fdVdy(j,1) - MA_RUW(j)*fdVdz(j,1) -&
          MA_RUV(j)*fdUdx(j,1) - MA_RVV(j)*fdUdy(j,1) - MA_RVW(j)*fdUdz(j,1) -&
          (MA_RUVx(j)+MA_RVVy(j)+MA_RVWz(j))*fU(j,1)-&
          (MA_RUUx(j)+MA_RUVy(j)+MA_RUWz(j))*fV(j,1)-&
          MA_RUVx(j)*fU(j,1) - MA_RUVy(j)*fV(j,1) - MA_RUVz(j)*fW(j,1) )*pts +&
          C_2_R*( fU(j,1)*fV(j,1)*dum1 +&
          rR(j,1)*(fU(j,1)*fV(j,1)*dum5 + fU(j,1)*dum3 + fV(j,1)*dum2) )   

     T1xz(j,1) = (MA_RUWUkk(j) - MA_RUW(j)*dum5   -&
          MA_RUU(j)*fdWdx(j,1) - MA_RUV(j)*fdWdy(j,1) - MA_RUW(j)*fdWdz(j,1) -&
          MA_RUW(j)*fdUdx(j,1) - MA_RVW(j)*fdUdy(j,1) - MA_RWW(j)*fdUdz(j,1) -&
          (MA_RUWx(j)+MA_RVWy(j)+MA_RWWz(j))*fU(j,1)-&
          (MA_RUUx(j)+MA_RUVy(j)+MA_RUWz(j))*fW(j,1) -&
          MA_RUWx(j)*fU(j,1) - MA_RUWy(j)*fV(j,1) - MA_RUWz(j)*fW(j,1) )*pts +&
          C_2_R*( fU(j,1)*fW(j,1)*dum1 +&
          rR(j,1)*(fU(j,1)*fW(j,1)*dum5 + fU(j,1)*dum4 + fW(j,1)*dum2) )   

     T1yz(j,1) = (MA_RVWUkk(j) - MA_RVW(j)*dum5   -&
          MA_RUV(j)*fdWdx(j,1) - MA_RVV(j)*fdWdy(j,1) - MA_RVW(j)*fdWdz(j,1) -&
          MA_RUW(j)*fdVdx(j,1) - MA_RVW(j)*fdVdy(j,1) - MA_RWW(j)*fdVdz(j,1) -&
          (MA_RUWx(j)+MA_RVWy(j)+MA_RWWz(j))*fV(j,1)-&
          (MA_RUVx(j)+MA_RVVy(j)+MA_RVWz(j))*fW(j,1)-&
          MA_RVWx(j)*fU(j,1) - MA_RVWy(j)*fV(j,1) - MA_RVWz(j)*fW(j,1) )*pts +&
          C_2_R*( fV(j,1)*fW(j,1)*dum1 + &
          rR(j,1)*(fV(j,1)*fW(j,1)*dum5 + fV(j,1)*dum4 + fW(j,1)*dum3) )   

! -------------------------------------------------------------------
! Viscous element of transport term of Reynolds equations
! -------------------------------------------------------------------
     T4xx(j,1) = (MA_TAUXkUk(j) + MA_UTAUXkk(j) - rU(j,1)*MA_TAUXkk(j))*pts -&
          tau_xx(j,1)*fdUdx(j,1) - tau_xy(j,1)*fdUdy(j,1) - tau_xz(j,1)*fdUdz(j,1)

     T4xy(j,1) = (MA_TAUYkUk(j) + MA_UTAUYkk(j) - rU(j,1)*MA_TAUYkk(j))*pts -&
          tau_xy(j,1)*fdUdx(j,1) - tau_yy(j,1)*fdUdy(j,1) - tau_yz(j,1)*fdUdz(j,1)

     T4xz(j,1) = (MA_TAUZkUk(j) + MA_UTAUZkk(j) - rU(j,1)*MA_TAUZkk(j))*pts -&
          tau_xz(j,1)*fdUdx(j,1) - tau_yz(j,1)*fdUdy(j,1) - tau_zz(j,1)*fdUdz(j,1)

     T4yx(j,1) = (MA_TAUXkVk(j) + MA_VTAUXkk(j) - rV(j,1)*MA_TAUXkk(j))*pts -&
          tau_xx(j,1)*fdVdx(j,1) - tau_xy(j,1)*fdVdy(j,1) - tau_xz(j,1)*fdVdz(j,1)

     T4yy(j,1) = (MA_TAUYkVk(j) + MA_VTAUYkk(j) - rV(j,1)*MA_TAUYkk(j))*pts -&
          tau_xy(j,1)*fdVdx(j,1) - tau_yy(j,1)*fdVdy(j,1) - tau_yz(j,1)*fdVdz(j,1)

     T4yz(j,1) = (MA_TAUZkVk(j) + MA_VTAUZkk(j) - rV(j,1)*MA_TAUZkk(j))*pts -&
          tau_xz(j,1)*fdVdx(j,1) - tau_yz(j,1)*fdVdy(j,1) - tau_zz(j,1)*fdVdz(j,1)

     T4zx(j,1) = (MA_TAUXkWk(j) + MA_WTAUXkk(j) - rW(j,1)*MA_TAUXkk(j))*pts -&
          tau_xx(j,1)*fdWdx(j,1) - tau_xy(j,1)*fdWdy(j,1) - tau_xz(j,1)*fdWdz(j,1)

     T4zy(j,1) = (MA_TAUYkWk(j) + MA_WTAUYkk(j) - rW(j,1)*MA_TAUYkk(j))*pts -&
          tau_xy(j,1)*fdWdx(j,1) - tau_yy(j,1)*fdWdy(j,1) - tau_yz(j,1)*fdWdz(j,1)

     T4zz(j,1) = (MA_TAUZkWk(j) + MA_WTAUZkk(j) - rW(j,1)*MA_TAUZkk(j))*pts -&
          tau_xz(j,1)*fdWdx(j,1) - tau_yz(j,1)*fdWdy(j,1) - tau_zz(j,1)*fdWdz(j,1)

! -------------------------------------------------------------------
! Rxx Reynolds stress equation
! -------------------------------------------------------------------
     Conv_xx(j,1) =-fU(j,1)*dRxxdx(j,1) - fV(j,1)*dRxxdy(j,1) - fW(j,1)*dRxxdz(j,1)
     Prod_xx(j,1) =-C_2_R*( fRxx(j,1)*fdUdx(j,1)+fRxy(j,1)*fdUdy(j,1)+fRxz(j,1)*fdUdz(j,1) )
     Diss_xx(j,1) =-C_2_R*( MA_TAUXkUk(j)*pts -&
          tau_xx(j,1)*rdUdx(j,1) - tau_xy(j,1)*rdUdy(j,1) - tau_xz(j,1)*rdUdz(j,1) )/rR(j,1)

! pressure terms with Reynolds average
     Tran_xx(j,1) =-(T1xx(j,1) + C_2_R*(-T4xx(j,1) + (MA_PUx(j)+MA_UPx(j))*pts -&
          rP(j,1)*rdUdx(j,1) - rU(j,1)*dPdx(j,1)) )/rR(j,1)
     Pres_xx(j,1) = C_2_R*( MA_PUx(j)*pts - rP(j,1)*rdUdx(j,1) )/rR(j,1)
     MnFl_xx(j,1) = C_2_R*(rU(j,1)-fU(j,1))*(MA_TAUXkk(j)*pts-dPdx(j,1))/rR(j,1)

     Resi_xx(j,1) = Conv_xx(j,1) + Prod_xx(j,1) + Diss_xx(j,1) + Tran_xx(j,1) +&
          Pres_xx(j,1) + MnFl_xx(j,1) 

! -------------------------------------------------------------------
! Ryy Reynolds stress equation
! -------------------------------------------------------------------
     Conv_yy(j,1) =-fU(j,1)*dRyydx(j,1) - fV(j,1)*dRyydy(j,1) - fW(j,1)*dRyydz(j,1)
     Prod_yy(j,1) =-C_2_R*( fRxy(j,1)*fdVdx(j,1)+fRyy(j,1)*fdVdy(j,1)+fRyz(j,1)*fdVdz(j,1) )
     Diss_yy(j,1) =-C_2_R*( MA_TAUYkVk(j)*pts -&
          tau_xy(j,1)*rdVdx(j,1) - tau_yy(j,1)*rdVdy(j,1) - tau_yz(j,1)*rdVdz(j,1) )/rR(j,1)

! pressure terms with Reynolds average
     Tran_yy(j,1) =-(T1yy(j,1) + C_2_R*(-T4yy(j,1) + (MA_PVY(j)+MA_VPy(j))*pts -&
          rP(j,1)*rdVdy(j,1) - rV(j,1)*dPdy(j,1) ) )/rR(j,1)
     Pres_yy(j,1) = C_2_R*( MA_PVY(j)*pts - rP(j,1)*rdVdy(j,1) )/rR(j,1)
     MnFl_yy(j,1) = C_2_R*(rV(j,1)-fV(j,1))*(MA_TAUYkk(j)*pts-dPdy(j,1))/rR(j,1)        

     Resi_yy(j,1) = Conv_yy(j,1) + Prod_yy(j,1) + Diss_yy(j,1) + Tran_yy(j,1) +&
          Pres_yy(j,1) + MnFl_yy(j,1) 

! -------------------------------------------------------------------
! Rzz Reynolds stress equation
! -------------------------------------------------------------------
     Conv_zz(j,1) =-fU(j,1)*dRzzdx(j,1) - fV(j,1)*dRzzdy(j,1) - fW(j,1)*dRzzdz(j,1)
     Prod_zz(j,1) =-C_2_R*( fRxz(j,1)*fdWdx(j,1)+fRyz(j,1)*fdWdy(j,1)+fRzz(j,1)*fdWdz(j,1) )
     Diss_zz(j,1) =-C_2_R*( MA_TAUZkWk(j)*pts -&
          tau_xz(j,1)*rdWdx(j,1) - tau_yz(j,1)*rdWdy(j,1) - tau_zz(j,1)*rdWdz(j,1) )/rR(j,1)

! pressure terms with Reynolds average
     Tran_zz(j,1) =-(T1zz(j,1) + C_2_R*(-T4zz(j,1) + (MA_PWz(j)+MA_WPz(j))*pts -&
          rP(j,1)*rdWdz(j,1) - rW(j,1)*dPdz(j,1) ) )/rR(j,1)
     Pres_zz(j,1) = C_2_R*( MA_PWz(j)*pts - rP(j,1)*rdWdz(j,1) )/rR(j,1)
     MnFl_zz(j,1) = C_2_R*(rW(j,1)-fW(j,1))*(MA_TAUZkk(j)*pts-dPdz(j,1))/rR(j,1)        

     Resi_zz(j,1) = Conv_zz(j,1) + Prod_zz(j,1) + Diss_zz(j,1) + Tran_zz(j,1) +&
          Pres_zz(j,1) + MnFl_zz(j,1) 

! -------------------------------------------------------------------
! Rxy Reynolds stress equation
! -------------------------------------------------------------------
     Conv_xy(j,1) =-fU(j,1)*dRxydx(j,1) - fV(j,1)*dRxydy(j,1) - fW(j,1)*dRxydz(j,1)
     Prod_xy(j,1) =&
          - fRxx(j,1)*fdVdx(j,1) - fRxy(j,1)*fdVdy(j,1) - fRxz(j,1)*fdVdz(j,1)&
          - fRxy(j,1)*fdUdx(j,1) - fRyy(j,1)*fdUdy(j,1) - fRyz(j,1)*fdUdz(j,1)
     Diss_xy(j,1) =-( &
          MA_TAUXkVk(j)*pts  -&
          tau_xx(j,1)*rdVdx(j,1) - tau_xy(j,1)*rdVdy(j,1) - tau_xz(j,1)*rdVdz(j,1) +&
          MA_TAUYkUk(j)*pts  -&
          tau_xy(j,1)*rdUdx(j,1) - tau_yy(j,1)*rdUdy(j,1) - tau_yz(j,1)*rdUdz(j,1) )/rR(j,1)

! pressure terms with Reynolds average
     Tran_xy(j,1) =-( T1xy(j,1) - T4xy(j,1) - T4yx(j,1) + &
          (MA_PUy(j)+MA_UPy(j))*pts - rP(j,1)*rdUdy(j,1) - rU(j,1)*dPdy(j,1) + &
          (MA_PVX(j)+MA_VPx(j))*pts - rP(j,1)*rdVdx(j,1) - rV(j,1)*dPdx(j,1) )/rR(j,1)
     Pres_xy(j,1) =( MA_PUy(j)*pts - rP(j,1)*rdUdy(j,1) + &
          MA_PVX(j)*pts - rP(j,1)*rdVdx(j,1) )/rR(j,1)
     MnFl_xy(j,1) = ( (rU(j,1)-fU(j,1))*(MA_TAUYkk(j)*pts-dPdy(j,1)) +&
          (rV(j,1)-fV(j,1))*(MA_TAUXkk(j)*pts-dPdx(j,1)) )/rR(j,1)     

     Resi_xy(j,1) = Conv_xy(j,1) + Prod_xy(j,1) + Diss_xy(j,1) + Tran_xy(j,1) +&
          Pres_xy(j,1) + MnFl_xy(j,1) 

! -------------------------------------------------------------------
! Turbulent kinetic energy equation
! -------------------------------------------------------------------
     Conv(j,1) = C_05_R*(Conv_xx(j,1) + Conv_yy(j,1) + Conv_zz(j,1))
     Prod(j,1) = C_05_R*(Prod_xx(j,1) + Prod_yy(j,1) + Prod_zz(j,1))
     Diss(j,1) = C_05_R*(Diss_xx(j,1) + Diss_yy(j,1) + Diss_zz(j,1))
     Pres(j,1) = C_05_R*(Pres_xx(j,1) + Pres_yy(j,1) + Pres_zz(j,1))
     Tran(j,1) = C_05_R*(Tran_xx(j,1) + Tran_yy(j,1) + Tran_zz(j,1))
     MnFl(j,1) = C_05_R*(MnFl_xx(j,1) + MnFl_yy(j,1) + MnFl_zz(j,1))

     Resi(j,1) = C_05_R*(Resi_xx(j,1) + Resi_yy(j,1) + Resi_zz(j,1))

! -------------------------------------------------------------------
! Energy equation in terms of p
! -------------------------------------------------------------------
! pressure terms with Reynolds average
     Conv_p(j,1) =-(fU(j,1)*dPdx(j,1) + fV(j,1)*dPdy(j,1) + fW(j,1)*dPdz(j,1))
     Reve_p(j,1) =-gama0*rP(j,1)*Dil(j,1)
     Diss_p(j,1) = (gama0-C_1_R)*phi(j,1)
     Tran_p(j,1) = MA_Tkk(j)*pts*gama0/reynolds/prandtl
     Reyn_p(j,1) =-( (MA_UkPk(j)+MA_PUx(j)+MA_PVY(j)+MA_PWz(j))*pts-&
          rP(j,1)*Dil(j,1)+ Conv_p(j,1) )

     Resi_p(j,1) = Conv_p(j,1) + Reve_p(j,1) + Diss_p(j,1) + Tran_p(j,1) +&
          Reyn_p(j,1) - (gama0-1)*rR(j,1)*Pres(j,1)  

! -------------------------------------------------------------------
! Energy equation in terms of T
! -------------------------------------------------------------------
! using MRATIO*p/rho=T
     fdTdx = (MRATIO*dPdx(j,1)-fT(j,1)*dRdx(j,1))/rR(j,1)
     fdTdy = (MRATIO*dPdy(j,1)-fT(j,1)*dRdy(j,1))/rR(j,1)
     fdTdz = (MRATIO*dPdz(j,1)-fT(j,1)*dRdz(j,1))/rR(j,1)

     Conv_T(j,1) =-(fU(j,1)*fdTdx + fV(j,1)*fdTdy + fW(j,1)*fdTdz )
! dilatation-pressure terms with Reynolds average
     Reve_T(j,1) =-MRATIO*(gama0-1)*rP(j,1)*Dil(j,1)/rR(j,1)
     Diss_T(j,1) = gama0*phi(j,1)/rR(j,1)
     Tran_T(j,1) = MA_Tkk(j)*pts*gama0/reynolds/prandtl/rR(j,1)
     Reyn_T(j,1) = -( MRATIO*(MA_UkPk(j) + MA_PUx(j) + MA_PVY(j) + MA_PWz(j))*pts/rR(j,1)+&
          Conv_T(j,1) )

     Resi_T(j,1) = Conv_T(j,1) + Reve_T(j,1) + Diss_T(j,1) + Tran_T(j,1) +&
          Reyn_T(j,1) - MRATIO*(gama0-1)*Pres(j,1)      

! -------------------------------------------------------------------
! Turbulent temperature equation 
! -------------------------------------------------------------------
! !!! Not complete 
     dfTdx = (MRATIO*dPdx(j,1)-fT(j,1)*dRdx(j,1))/rR(j,1)
     dfTdy = (MRATIO*dPdy(j,1)-fT(j,1)*dRdy(j,1))/rR(j,1)

     dRTTdx = MRATIO*(MA_PTx(j)+MA_TPx(j))*pts
     dRTTdy = MRATIO*(MA_PTy(j)+MA_TPy(j))*pts

     dfTf2dx = (dRTTdx-(fT(j,1)*fT(j,1)+fTf2(j,1))*dRdx(j,1))/rR(j,1)-C_2_R*fT(j,1)*dfTdx
     dfTf2dy = (dRTTdy-(fT(j,1)*fT(j,1)+fTf2(j,1))*dRdy(j,1))/rR(j,1)-C_2_R*fT(j,1)*dfTdy

     Conv_tt(j,1) = -fU(j,1)*dfTf2dx-fV(j,1)*dfTf2dy
     Prod_tt(j,1) = -C_2_R*(fRuT(j,1)*dfTdx+fRvT(j,1)*dfTdy)

     dRUTdx = MRATIO*(MA_PUx(j)+MA_UPx(j))*pts
     dRVTdy = MRATIO*(MA_PVY(j)+MA_VPy(j))*pts

     tranttx = MA_RUTTx(j)*pts - &
          fU(j,1)*dRTTdx -&
          rR(j,1)*(fT(j,1)**2+fTf2(j,1))*&
          fdUdx(j,1) - &
          C_2_R*fT(j,1)*dRUTdx -&
          C_2_R*rR(j,1)*(fU(j,1)*fT(j,1)+&
          fRuT(j,1))*dfTdx +&
          C_2_R*fU(j,1)*fT(j,1)**2*dRdx(j,1) +&
          C_2_R*rR(j,1)*fT(j,1)**2*fdUdx(j,1) +&
          C_4_R*rR(j,1)*fU(j,1)*fT(j,1)*dfTdx
     trantty = MA_RVTTy(j)*pts - &
          fV(j,1)*dRTTdy -&
          rR(j,1)*(fT(j,1)**2+fTf2(j,1))*&
          fdVdy(j,1) - &
          C_2_R*fT(j,1)*dRVTdy -&
          C_2_R*rR(j,1)*(fV(j,1)*fT(j,1)+&
          fRvT(j,1))*dfTdy +&
          C_2_R*fV(j,1)*fT(j,1)**2*dRdy(j,1) +&
          C_2_R*rR(j,1)*fT(j,1)**2*fdVdy(j,1) +&
          C_4_R*rR(j,1)*fV(j,1)*fT(j,1)*dfTdy

     Tran_tt(j,1) =-C_2_R*(tranttx + trantty)

     Diss_tt(j,1) = C_0_R
     Pres_tt(j,1) = C_0_R
     MnFl_tt(j,1) = C_0_R

     Resi_tt(j,1) = Conv_tt(j,1) + Prod_tt(j,1) + Tran_tt(j,1) + Diss_tt(j,1) + MnFl_tt(j,1)

! ###################################################################
! Variable density quantities
! ###################################################################
! speed of sound. Using Reynolds, not Favre, for T
     dum1 = rT(j,1)/(mach*mach)
     dum2 = rT(j,1)*(C_1_R/rP(j,1)-C_1_R/(rR(j,1)*dum1))

     rho_p(j,1) = MA_RP(j)*pts - rR(j,1)*rP(j,1)
     rho_T(j,1) = MA_RT(j)*pts - rR(j,1)*rT(j,1)
! T-p correlation
     dum3 = MA_RTT(j)*pts/MRATIO- rT(j,1)*rP(j,1)

     rho_ac(j,1) = rPf2(j,1)/(dum1*dum1)
     rho_en(j,1) = rRf2(j,1)+rho_ac(j,1)-C_2_R*rho_p(j,1)/dum1

     T_ac(j,1) = rPf2(j,1)*dum2*dum2
     T_en(j,1) = rTf2(j,1)+T_ac(j,1)-C_2_R*dum3*dum2

! ###################################################################
! Scales
! ###################################################################
     IF ( Diss(j,1) .EQ. C_0_R ) THEN
        eta(j,1)    = C_BIG_R
        tau(j,1)    = C_BIG_R
        lambda(j,1) = C_BIG_R
     ELSE
        eta(j,1) = C_1_R/( ABS(Diss(j,1))*(reynolds*rR(j,1))**C_3_R )**C_025_R
        tau(j,1) = C_1_R/SQRT( reynolds*rR(j,1)*ABS(Diss(j,1)) )
        lambda(j,1) = SQRT( C_10_R*rTKE(j,1) / (rR(j,1)*ABS(Diss(j,1))*reynolds) )
     ENDIF

     IF (rdUdxf2(j,1) .EQ. C_0_R ) THEN
        lambda_x(j,1) = C_BIG_R
     ELSE
        lambda_x(j,1) = SQRT( rUf2(j,1) / rdUdxf2(j,1) )
     ENDIF

     IF (rdVdyf2(j,1) .EQ. C_0_R ) THEN
        lambda_y(j,1) = C_BIG_R
     ELSE
        lambda_y(j,1) = SQRT( rVf2(j,1) / rdVdyf2(j,1) )
     ENDIF

     IF (rdWdzf2(j,1) .EQ. C_0_R ) THEN
        lambda_z(j,1) = C_BIG_R
     ELSE
        lambda_z(j,1) = SQRT( rWf2(j,1) / rdWdzf2(j,1) )
     ENDIF

! ###################################################################
! Skewness and flatness
! ###################################################################
     S_rho(j,1) = MA_R3(j)*pts - rR(j,1)**C_3_R - C_3_R*rR(j,1)*rRf2(j,1)
     S_u(j,1) =   MA_U3(j)*pts - rU(j,1)**C_3_R - C_3_R*rU(j,1)*rUf2(j,1)
     S_v(j,1) =   MA_V3(j)*pts - rV(j,1)**C_3_R - C_3_R*rV(j,1)*rVf2(j,1)
     S_w(j,1) =   MA_W3(j)*pts - rW(j,1)**C_3_R - C_3_R*rW(j,1)*rWf2(j,1)
     S_p(j,1) =   MA_P3(j)*pts - rP(j,1)**C_3_R - C_3_R*rP(j,1)*rPf2(j,1)
     S_T(j,1) =   MA_T3(j)*pts - rT(j,1)**C_3_R - C_3_R*rT(j,1)*rTf2(j,1)

     F_rho(j,1) = MA_R4(j)*pts - rR(j,1)**C_4_R - C_4_R*rR(j,1)*S_rho(j,1)-&
          C_6_R*rR(j,1)**C_2_R*rRf2(j,1)
     F_u(j,1) =   MA_U4(j)*pts - rU(j,1)**C_4_R - C_4_R*rU(j,1)*S_u(j,1)-&
          C_6_R*rU(j,1)**C_2_R*rUf2(j,1)
     F_v(j,1) =   MA_V4(j)*pts - rV(j,1)**C_4_R - C_4_R*rV(j,1)*S_v(j,1)-&
          C_6_R*rV(j,1)**C_2_R*rVf2(j,1)
     F_w(j,1) =   MA_W4(j)*pts - rW(j,1)**C_4_R - C_4_R*rW(j,1)*S_w(j,1)-&
          C_6_R*rW(j,1)**C_2_R*rWf2(j,1)
     F_p(j,1) =   MA_P4(j)*pts - rP(j,1)**C_4_R - C_4_R*rP(j,1)*S_p(j,1)-&
          C_6_R*rP(j,1)**C_2_R*rPf2(j,1)
     F_T(j,1) =   MA_T4(j)*pts - rT(j,1)**C_4_R - C_4_R*rT(j,1)*S_T(j,1)-&
          C_6_R*rT(j,1)**C_2_R*rTf2(j,1)

! Normalization
     IF ( rRf2(j,1) .EQ. C_0_R) THEN
        S_rho(j,1) = C_BIG_R
        F_rho(j,1) = C_BIG_R
     ELSE
        S_rho(j,1) = S_rho(j,1)/rRf2(j,1)**(C_3_R/C_2_R)
        F_rho(j,1) = F_rho(j,1)/rRf2(j,1)**C_2_R
     ENDIF
     IF ( rUf2(j,1) .EQ. C_0_R) THEN
        S_u(j,1) = C_BIG_R
        F_u(j,1) = C_BIG_R
     ELSE
        S_u(j,1) = S_u(j,1)/rUf2(j,1)**(C_3_R/C_2_R)
        F_u(j,1) = F_u(j,1)/rUf2(j,1)**C_2_R
     ENDIF
     IF ( rVf2(j,1) .EQ. C_0_R) THEN
        S_v(j,1) = C_BIG_R
        F_v(j,1) = C_BIG_R
     ELSE
        S_v(j,1) = S_v(j,1)/rVf2(j,1)**(C_3_R/C_2_R)
        F_v(j,1) = F_v(j,1)/rVf2(j,1)**C_2_R
     ENDIF
     IF ( rWf2(j,1) .EQ. C_0_R) THEN
        S_w(j,1) = C_BIG_R
        F_w(j,1) = C_BIG_R
     ELSE
        S_w(j,1) = S_w(j,1)/rWf2(j,1)**(C_3_R/C_2_R)
        F_w(j,1) = F_w(j,1)/rWf2(j,1)**C_2_R
     ENDIF
     IF ( rPf2(j,1) .EQ. C_0_R) THEN
        S_p(j,1) = C_BIG_R
        F_p(j,1) = C_BIG_R
     ELSE
        S_p(j,1) = S_p(j,1)/rPf2(j,1)**(C_3_R/C_2_R)
        F_p(j,1) = F_p(j,1)/rPf2(j,1)**C_2_R
     ENDIF
     IF ( rTf2(j,1) .EQ. C_0_R) THEN
        S_T(j,1) = C_BIG_R
        F_T(j,1) = C_BIG_R
     ELSE
        S_T(j,1) = S_T(j,1)/rTf2(j,1)**(C_3_R/C_2_R)
        F_T(j,1) = F_T(j,1)/rTf2(j,1)**C_2_R
     ENDIF

  ENDDO

! ###################################################################
! Integral quantities shear layer
! ###################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN

! Vorticity Thickness
     DO n = 1,nstatavg
        DO j = 1,jmax
           wrk1d(j,1) = fU(n,j)
        ENDDO
        CALL PARTIAL_Y(imode_fdm, i1, jmax, i1, j1bc, dy, wrk1d(1,1), &
             wrk1d(1,2), i0, i0, wrk1d(1,3), wrk2d, wrk2d(1,2))
        delta_w_u(n) = (fU(n,jmax)-fU(n,1)) / MINVAL(wrk1d(1:jmax,2))
     ENDDO

! Momentum thickness
     DO n = 1,nstatavg
        UC = mean_u
        DU = delta_u
        DO j = jmin_loc, jmax_loc
           wrk1d(j,1) = rR(n,j)*( C_025_R - ((fU(n,j)-UC)/DU)**2 )
        ENDDO
        delta_m_u(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
     ENDDO

! Mixing layer limit (U=0.1dU and U=0.9dU)
     y_center = g(2)%nodes(1) + ycoor_u*scaley
     DO n = 1,nstatavg
        fU_05 = U2 + C_01_R*delta_u
        DO j = 1,jmax
           IF ( fU(n,j) .GT. fU_05 .AND. fU(n,j+1) .LE. fU_05 ) THEN
              delta_01_u(n) = g(2)%nodes(j) + (fU_05-fU(n,j))*(g(2)%nodes(j+1) -g(2)%nodes(j))/(fU(n,j+1)-fU(n,j))
           ENDIF
        ENDDO
        delta_01_u(n) = delta_01_u(n) - y_center

        fU_05 = U2 + r09*delta_u
        DO j = 1,jmax
           IF ( fU(n,j) .GT. fU_05 .AND. fU(n,j+1) .LE. fU_05 ) THEN
              delta_01_d(n) = g(2)%nodes(j) + (fU_05-fU(n,j))*(g(2)%nodes(j+1) -g(2)%nodes(j))/(fU(n,j+1)-fU(n,j))
           ENDIF
        ENDDO
        delta_01_d(n) = delta_01_d(n) - y_center

        delta_u_u(n) = C_05_R*( delta_01_u(n) + delta_01_d(n) )
     ENDDO

! ###################################################################
! 1D quantities of the jet
! ###################################################################
  ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
! -------------------------------------------------------------------
! Integral balance of mass
! -------------------------------------------------------------------
     DO n = 1,nstatavg
! axial
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = rR(n,j)*fU(n,j)
        ENDDO
        IntMassU(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
! lateral
        DO k = 1,n
           i = statavg(k)
           wrk1d(k,1) = rR(k,jmin_loc)*fV(k,jmin_loc) - rR(k,jmax_loc)*fV(k,jmax_loc)
           wrk1d(k,2) = g(1)%nodes(i)
        ENDDO
        IF ( n .EQ. 1 ) THEN
           IntMassV(n) = C_0_R
        ELSE IF ( n .EQ. 2 ) THEN
           IntMassV(n) = C_05_R*(wrk1d(1,1)+wrk1d(2,1))*(wrk1d(2,2)-wrk1d(1,2))
        ELSE 
           IntMassV(n) = SIMPSON_NU(n,wrk1d(1,1),wrk1d(1,2))
        ENDIF
     ENDDO

! -------------------------------------------------------------------
! Integral balance of axial momentum (conserved)
! -------------------------------------------------------------------
     DO n = 1,nstatavg
! mean velocity part
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = rR(n,j)*fU(n,j)*(fU(n,j)-U2)
        ENDDO
        IntExcMomU(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
! pressure part
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = (rP(n,j)-pbg%mean)
        ENDDO
        IntExcMomP(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
! Reynolds stress part
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = rR(n,j)*fRxx(n,j)
        ENDDO
        IntExcMomRxx(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
     ENDDO

! -------------------------------------------------------------------
! Integral balance of turbulent kinetic energy
! -------------------------------------------------------------------
     DO n = 1,nstatavg
! TKE flux
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = rR(n,j)*fU(n,j)*fTKE(n,j)
        ENDDO
        IntTkeK(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
! Integral of production term
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = rR(n,j)*Prod(n,j)
        ENDDO
        IntTkeP(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
! Integral of numerical dissipation
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) =-rR(n,j)*eps_f(n,j)
        ENDDO
        IntTkeF(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
! Integral of production term
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = Pres(n,j)
        ENDDO
        IntTkePi(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
     ENDDO

! -------------------------------------------------------------------
! Axial flux of heat
! -------------------------------------------------------------------
     DO n = 1,nstatavg
        DO j = jmin_loc,jmax_loc
           wrk1d(j,1) = rR(n,j)*fU(n,j)*(fT(n,j)-T2)
        ENDDO
        IntFluxT(n) = SIMPSON_NU(nj,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))
     ENDDO

! -------------------------------------------------------------------
! Jet thickness
! -------------------------------------------------------------------
! Vorticity thickness
     DO n = 1,nstatavg
        DO j = 1,jmax
           wrk2d(j,1) = fU(n,j)
        ENDDO
        CALL PARTIAL_Y( imode_fdm, i1, jmax, i1, j1bc, dy, wrk2d(1,1), &
             wrk2d(1,2), i0, i0, wrk1d, wrk2d(1,3), wrk2d(1,4) )
        delta_w_u(n) = (fU(n,jmax/2+1)-U2)/ ABS(MINVAL(wrk2d(1:jmax,2)))
        delta_w_d(n) = (fU(n,jmax/2)  -U2)/ ABS(MAXVAL(wrk2d(1:jmax,2)))
     ENDDO

! Momentum thickness
     DO n = 1,nstatavg
        UC = C_05_R*(U2+fU(n,jmax/2))
        DU = fU(n,jmax/2)-U2
        DO j = jmin_loc,jmax/2
           wrk1d(j,1) = rR(n,j)*( C_025_R - ((fU(n,j)-UC)/DU)**2 )
        ENDDO
        delta_m_d(n) = SIMPSON_NU(jmax/2-jmin_loc+1,wrk1d(jmin_loc,1), g(2)%nodes(jmin_loc))

        UC = C_05_R*(U2+fU(n,jmax/2+1))
        DU = fU(n,jmax/2+1)-U2
        DO j = jmax/2+1,jmax_loc
           wrk1d(j,1) = rR(n,j)*( C_025_R - ((fU(n,j)-UC)/DU)**2 )
        ENDDO
        delta_m_u(n) = SIMPSON_NU(jmax_loc-jmax/2, wrk1d(jmax/2+1,1), g(2)%nodes(jmax/2+1))
     ENDDO

! Jet half-width based on velocity
     CALL DELTA_X(nstatavg, jmax, g(2)%nodes, fU(1,1), wrk1d(1,1),&
          delta_u_d(1), delta_u_u(1), U2, r05)

! Jet Limit (U=0.05Uc)
     CALL DELTA_X(nstatavg, jmax, g(2)%nodes, fU(1,1), wrk1d(1,1),&
          delta_01_d(1), delta_01_u(1), U2, r005)

! Jet half-width based on temperature/density
     IF ( rbg%delta .NE. C_0_R ) THEN
        DO j = 1,jmax*nstatavg
           wrk2d(j,1) = ABS(fT(j,1)-T2) + T2 ! we can have hot or cold jet
        ENDDO
        CALL DELTA_X(nstatavg, jmax, g(2)%nodes, wrk2d(1,1), wrk1d(1,1),&
             delta_t_d(1), delta_t_u(1), T2, r05)

        DO j = 1,jmax*nstatavg
           wrk2d(j,1) = ABS(rR(j,1)-R2) + R2 ! we can have hot or cold jet
        ENDDO
        CALL DELTA_X(nstatavg, jmax, g(2)%nodes, wrk2d(1,1), wrk1d(1,1),&
             delta_r_d(1), delta_r_u(1), R2, r05)
     ELSE
        DO n = 1,nstatavg
           delta_t_d(n) = C_1_R
           delta_t_u(n) = C_1_R
           delta_r_d(n) = C_1_R
           delta_r_u(n) = C_1_R
        ENDDO
     ENDIF

! Jet center line based on velocity
     y_center = g(2)%nodes(1) + ycoor_u*scaley
     DO n = 1,nstatavg
        DO j = 1,jmax
           wrk1d(j,1) = fU(n,j)
        ENDDO
        jloc_max = MAXLOC(wrk1d(1:jmax,1)); j = jloc_max(1)
        IF ( wrk1d(j-1,1) .GT. wrk1d(j+1,1) ) THEN
           delta_u_center(n) = C_05_R*(g(2)%nodes(j) +g(2)%nodes(j-1))
        ELSE
           delta_u_center(n) = C_05_R*(g(2)%nodes(j) +g(2)%nodes(j+1))
        ENDIF
        delta_u_center(n) = delta_u_center(n) - y_center
     ENDDO

  ENDIF

! ###################################################################
! Scaling of the quatities
! ###################################################################
#define simuc(A) wrk1d(A,2)
#define simtc(A) wrk1d(A,3)
#define simrc(A) wrk1d(A,4)

  DO n = 1,nstatavg

     IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        DU = delta_u
        delta_05 = delta_01_u(n) - delta_01_d(n)
     ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
        delta_05 = C_05_R*(delta_u_u(n)+delta_u_d(n))

        simuc(n) = C_05_R*(fU(n,jmax/2)+fU(n,jmax/2+1)) - U2
        IF ( rbg%delta .NE. C_0_R ) THEN
           simtc(n) = C_05_R*(fT(n,jmax/2)+fT(n,jmax/2+1)) - T2
           simrc(n) = C_05_R*(rR(n,jmax/2)+rR(n,jmax/2+1)) - R2
        ELSE
           simtc(n) = C_1_R
           simrc(n) = C_1_R
        ENDIF

        DU = simuc(n)
        DH = ABS(simtc(n))
     ENDIF

! reynolds based on half-width
     Reynolds_d(n) = reynolds*rR(n,jmax/2)* C_2_R*delta_05*DU
! reynolds based on isotropic lambda
     Reynolds_i(n) = reynolds*rR(n,jmax/2)* lambda(n,jmax/2)*SQRT(C_2_R*fTKE(n,jmax/2)/C_3_R)
! reynolds based on longitudinal lambda
     Reynolds_l(n) = reynolds*rR(n,jmax/2)* lambda_x(n,jmax/2)*SQRT(fRxx(n,jmax/2))

     DO j = 1,jmax
        Vortx(n,j)   = Vortx(n,j)/DU*delta_05
        Vorty(n,j)   = Vorty(n,j)/DU*delta_05
        Vortz(n,j)   = Vortz(n,j)/DU*delta_05
        Dil(n,j)     = Dil(n,j)/DU*delta_05
        Vortxf2(n,j) = SQRT(Vortxf2(n,j))/DU*delta_05
        Vortyf2(n,j) = SQRT(Vortyf2(n,j))/DU*delta_05
        Vortzf2(n,j) = SQRT(Vortzf2(n,j))/DU*delta_05
        Dilf2(n,j)   = Dilf2(n,j)/DU/DU*delta_05*delta_05

        Conv_xx(n,j) = Conv_xx(n,j)/(DU*DU*DU)*delta_05 
        Prod_xx(n,j) = Prod_xx(n,j)/(DU*DU*DU)*delta_05
        Diss_xx(n,j) = Diss_xx(n,j)/(DU*DU*DU)*delta_05
        Tran_xx(n,j) = Tran_xx(n,j)/(DU*DU*DU)*delta_05
        Pres_xx(n,j) = Pres_xx(n,j)/(DU*DU*DU)*delta_05
        MnFl_xx(n,j) = MnFl_xx(n,j)/(DU*DU*DU)*delta_05
        Resi_xx(n,j) = Resi_xx(n,j)/(DU*DU*DU)*delta_05

        Conv_yy(n,j) = Conv_yy(n,j)/(DU*DU*DU)*delta_05 
        Prod_yy(n,j) = Prod_yy(n,j)/(DU*DU*DU)*delta_05
        Diss_yy(n,j) = Diss_yy(n,j)/(DU*DU*DU)*delta_05
        Tran_yy(n,j) = Tran_yy(n,j)/(DU*DU*DU)*delta_05
        Pres_yy(n,j) = Pres_yy(n,j)/(DU*DU*DU)*delta_05
        MnFl_yy(n,j) = MnFl_yy(n,j)/(DU*DU*DU)*delta_05
        Resi_yy(n,j) = Resi_yy(n,j)/(DU*DU*DU)*delta_05

        Conv_zz(n,j) = Conv_zz(n,j)/(DU*DU*DU)*delta_05 
        Prod_zz(n,j) = Prod_zz(n,j)/(DU*DU*DU)*delta_05
        Diss_zz(n,j) = Diss_zz(n,j)/(DU*DU*DU)*delta_05
        Tran_zz(n,j) = Tran_zz(n,j)/(DU*DU*DU)*delta_05
        Pres_zz(n,j) = Pres_zz(n,j)/(DU*DU*DU)*delta_05
        MnFl_zz(n,j) = MnFl_zz(n,j)/(DU*DU*DU)*delta_05
        Resi_zz(n,j) = Resi_zz(n,j)/(DU*DU*DU)*delta_05

        Conv_xy(n,j) = Conv_xy(n,j)/(DU*DU*DU)*delta_05 
        Prod_xy(n,j) = Prod_xy(n,j)/(DU*DU*DU)*delta_05
        Diss_xy(n,j) = Diss_xy(n,j)/(DU*DU*DU)*delta_05
        Tran_xy(n,j) = Tran_xy(n,j)/(DU*DU*DU)*delta_05
        Pres_xy(n,j) = Pres_xy(n,j)/(DU*DU*DU)*delta_05
        MnFl_xy(n,j) = MnFl_xy(n,j)/(DU*DU*DU)*delta_05
        Resi_xy(n,j) = Resi_xy(n,j)/(DU*DU*DU)*delta_05

        fTKE(n,j) = fTKE(n,j)/(DU*DU)
        Conv(n,j) = Conv(n,j)/(DU*DU*DU)*delta_05 
        Prod(n,j) = Prod(n,j)/(DU*DU*DU)*delta_05
        Diss(n,j) = Diss(n,j)/(DU*DU*DU)*delta_05
        Tran(n,j) = Tran(n,j)/(DU*DU*DU)*delta_05
        Pres(n,j) = Pres(n,j)/(DU*DU*DU)*delta_05
        MnFl(n,j) = MnFl(n,j)/(DU*DU*DU)*delta_05
        Resi(n,j) = Resi(n,j)/(DU*DU*DU)*delta_05

        equi(n,j) = fTKE(n,j)/ABS(Diss(n,j))

        eps_f(n,j) =-eps_f(n,j)/(DU*DU*DU)*delta_05

        Conv_u(n,j) = Conv_u(n,j)/(DU*DU)*delta_05
        Tran_u(n,j) = Tran_u(n,j)/(DU*DU)*delta_05
        Reyn_u(n,j) = Reyn_u(n,j)/(DU*DU)*delta_05
        Resi_u(n,j) = Resi_u(n,j)/(DU*DU)*delta_05

        Conv_v(n,j) = Conv_v(n,j)/(DU*DU)*delta_05
        Tran_v(n,j) = Tran_v(n,j)/(DU*DU)*delta_05
        Reyn_v(n,j) = Reyn_v(n,j)/(DU*DU)*delta_05
        Resi_v(n,j) = Resi_v(n,j)/(DU*DU)*delta_05

        Conv_w(n,j) = Conv_w(n,j)/(DU*DU)*delta_05
        Tran_w(n,j) = Tran_w(n,j)/(DU*DU)*delta_05
        Reyn_w(n,j) = Reyn_w(n,j)/(DU*DU)*delta_05
        Resi_w(n,j) = Resi_w(n,j)/(DU*DU)*delta_05

        Conv_p(n,j) = Conv_p(n,j)/(DU*DU*DU)*delta_05 
        Reve_p(n,j) = Reve_p(n,j)/(DU*DU*DU)*delta_05
        Diss_p(n,j) = Diss_p(n,j)/(DU*DU*DU)*delta_05
        Tran_p(n,j) = Tran_p(n,j)/(DU*DU*DU)*delta_05
        Reyn_p(n,j) = Reyn_p(n,j)/(DU*DU*DU)*delta_05
        Resi_p(n,j) = Resi_p(n,j)/(DU*DU*DU)*delta_05

        Conv_T(n,j) = Conv_T(n,j)/(DH*DU)*delta_05 
        Reve_T(n,j) = Reve_T(n,j)/(DH*DU)*delta_05
        Diss_T(n,j) = Diss_T(n,j)/(DH*DU)*delta_05
        Tran_T(n,j) = Tran_T(n,j)/(DH*DU)*delta_05
        Reyn_T(n,j) = Reyn_T(n,j)/(DH*DU)*delta_05
        Resi_T(n,j) = Resi_T(n,j)/(DH*DU)*delta_05

        Conv_tt(n,j) = Conv_tt(n,j)/(DH*DH*DU)*delta_05 
        Prod_tt(n,j) = Prod_tt(n,j)/(DH*DH*DU)*delta_05
        Diss_tt(n,j) = Diss_tt(n,j)/(DH*DH*DU)*delta_05
        Tran_tt(n,j) = Tran_tt(n,j)/(DH*DH*DU)*delta_05
        Pres_tt(n,j) = Pres_tt(n,j)/(DH*DH*DU)*delta_05
        MnFl_tt(n,j) = MnFl_tt(n,j)/(DH*DH*DU)*delta_05
        Resi_tt(n,j) = Resi_tt(n,j)/(DH*DH*DU)*delta_05

     ENDDO

  ENDDO

! ###################################################################
! Saving the data in TkStat format
! ###################################################################
  WRITE(name,*) itime
  IF      ( imode_flow .EQ. DNS_FLOW_JET   ) THEN; name = 'jetavg'//TRIM(ADJUSTL(name))
  ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN; name = 'shravg'//TRIM(ADJUSTL(name)); ENDIF

#ifdef USE_RECLEN
  OPEN(UNIT=i23,RECL=3260,FILE=name,STATUS='unknown')
#else
  OPEN(UNIT=i23,FILE=name,STATUS='unknown')
#endif     

! -------------------------------------------------------------------
! Header
! -------------------------------------------------------------------
  WRITE(i23,'(A8,E14.7E3)') 'RTIME = ', rtime

! Independent variables
  line2 = 'I J X Y SU ST'

! Dependent variables depending on y and x
  line1 = 'Xg Yg'
  WRITE(i23,1010) 'GROUP = Grid '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'rU rV rW rP rR rT rUf2 rVf2 rWf2 rPf2 rRf2 rTf2 rUfVf rUfWf rVfWf rTKE ' &
       //'rbxx rbyy rbzz rbxy rbxz rbyz rRuT rRvT rRwT'
  WRITE(i23,1010) 'GROUP = Reynolds_Avgs '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'fU fV fW fT fTf2 fRxy fRxz fRyz fRxx fRyy fRzz '&
       //'fbxx fbyy fbzz fbxy fbxz fbyz fRuT fRvT fRwT'
  WRITE(i23,1010) 'GROUP = Favre_Avgs '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'rdUdx rdUdy rdUdz rdVdx rdVdy rdVdz rdWdx rdWdy rdWdz '&
       //'rdUdxf2 rdUdyf2 rdUdzf2 rdVdxf2 rdVdyf2 rdVdzf2 '&
       //'rdWdxf2 rdWdyf2 rdWdzf2 '&
       //'rdVdxfdUdyf rdWdxfdUdzf rdWdyfdVdzf '&
       //'rdUdxfdVdyf rdUdxfdWdzf rdVdyfdWdzf '&
       //'dPdx dPdy dPdz dRdx dRdy dRdz'
  WRITE(i23,1010) 'GROUP = Derivatives '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'Vortx Vorty Vortz Dil fDil Vortxf2 Vortyf2 Vortzf2 Dilf2'
  WRITE(i23,1010) 'GROUP = Vort_Dil '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'eta tau lambda lambda_x lambda_y lambda_z equi'
  WRITE(i23,1010) 'GROUP = Scales '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'Rxx Conv_xx Prod_xx Diss_xx Tran_xx Pres_xx MnFl_xx Resi_xx'
  WRITE(i23,1010) 'GROUP = Rxx_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'Ryy Conv_yy Prod_yy Diss_yy Tran_yy Pres_yy MnFl_yy Resi_yy'
  WRITE(i23,1010) 'GROUP = Ryy_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'Rzz Conv_zz Prod_zz Diss_zz Tran_zz Pres_zz MnFl_zz Resi_zz '
  WRITE(i23,1010) 'GROUP = Rzz_Eqn '//TRIM(ADJUSTL(line1)) 
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'Rxy Conv_xy Prod_xy Diss_xy Tran_xy Pres_xy MnFl_xy Resi_xy'
  WRITE(i23,1010) 'GROUP = Rxy_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'TKE Conv Prod Diss Tran Pres MnFl Resi'
  WRITE(i23,1010) 'GROUP = TKE_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'Rtt Conv_tt Prod_tt Diss_tt Tran_tt Pres_tt MnFl_tt Resi_tt'
  WRITE(i23,1010) 'GROUP = Rtt_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'U Conv_u Tran_u Reyn_u Resi_u'
  WRITE(i23,1010) 'GROUP = U_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'V Conv_v Tran_v Reyn_v Resi_v'
  WRITE(i23,1010) 'GROUP = V_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'W Conv_w Tran_w Reyn_w Resi_w'
  WRITE(i23,1010) 'GROUP = W_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'p Conv_p Reve_p Diss_p Tran_p Reyn_p Pres_p Resi_p'
  WRITE(i23,1010) 'GROUP = p_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'T Conv_T Reve_T Diss_T Tran_T Reyn_T Pres_T Resi_T'
  WRITE(i23,1010) 'GROUP = T_Eqn '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'fTKE_nf eps_f'
  WRITE(i23,1010) 'GROUP = Filter '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'tau_xx tau_yy tau_zz tau_xy tau_xz tau_yz phi rVis'
  WRITE(i23,1010) 'GROUP = Mean_Stresses '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'Corr_RP Corr_RT R_ac R_en T_ac T_en RuT RvT RwT Rur Rvr Rwr'
  WRITE(i23,1010) 'GROUP = VarDensity '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  line1 = 'S_R S_U S_V S_W S_P S_T F_R F_U F_V F_W F_P F_T'
  WRITE(i23,1010) 'GROUP = Skewness_Flatness '//TRIM(ADJUSTL(line1))
  line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

! dependent variables dependent on t only
  IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
     line1 ='Del_mom_u Del_mom_d Del_vor_u Del_vor_d '&
          //'Del_half_u Del_half_d Del_lim_u Del_lim_d '&
          //'Del_tem_u Del_tem_d Del_rho_u Del_rho_d Del_Umax '&
          //'Sim_U Sim_T '&
          //'Re_half Re_lambda_iso Re_lambda_lon '&
          //'Int_mom_U Int_mom_P Int_mom_Rxx '&
          //'Int_mass_U Int_mass_V Int_flux_T '&
          //'Int_tke_K Int_tke_Pi Int_tke_P Int_tke_F'
     WRITE(i23,1010) 'GROUP = 1D_Quantities '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

  ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
     line1 = 'Delta_m Delta_w y_01 y_09 y_05 Re Re_l'
     WRITE(i23,1010) 'GROUP = 1D_Quantities '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))
  ENDIF

  WRITE(i23,1010) TRIM(ADJUSTL(line2))

! -------------------------------------------------------------------
! Output
! -------------------------------------------------------------------
  DO n = 1,nstatavg
     i = statavg(n)

     IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        delta_05 = delta_01_u(n) - delta_01_d(n)
        delta_w  = delta_w_u(n)
     ELSE IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
        delta_05 = C_05_R*(delta_u_u(n)+delta_u_d(n))
        delta_w  = C_05_R*(delta_w_u(n)+delta_w_d(n))
        delta_t  = C_05_R*(delta_t_u(n)+delta_t_d(n))
     ENDIF

     IF ( imode_flow .EQ. DNS_FLOW_JET ) THEN
        ivauxpos = VARMX1D+2
        VAUXPOS(1)  = delta_m_u(n)
        VAUXPOS(2)  = delta_m_d(n)
        VAUXPOS(3)  = delta_w_u(n)
        VAUXPOS(4)  = delta_w_d(n)
        VAUXPOS(5)  = delta_u_u(n)
        VAUXPOS(6)  = delta_u_d(n)
        VAUXPOS(7)  = delta_01_u(n)
        VAUXPOS(8)  = delta_01_d(n)
        VAUXPOS(9)  = delta_t_u(n)
        VAUXPOS(10) = delta_t_d(n)
        VAUXPOS(11) = delta_r_u(n)
        VAUXPOS(12) = delta_r_d(n)
        VAUXPOS(13) = delta_u_center(n)
        VAUXPOS(14) = (simuc(1)/simuc(n))**C_2_R
        VAUXPOS(15) = (simtc(1)/simtc(n))**C_2_R
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
     ELSE IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
        ivauxpos = 7
        VAUXPOS(1) = delta_m_u(n)
        VAUXPOS(2) = delta_w_u(n)
        VAUXPOS(3) = delta_01_u(n)
        VAUXPOS(4) = delta_01_d(n)
        VAUXPOS(5) = delta_u_u(n)
        VAUXPOS(6) = Reynolds_d(n)
        VAUXPOS(7) = Reynolds_l(n)
     ENDIF

     DO j = 1,jmax
        ivauxpre = 4
        VAUXPRE(1) = g(1)%nodes(i)/diam_u
        VAUXPRE(2) = g(2)%nodes(j)/diam_u
        VAUXPRE(3) = (g(2)%nodes(j)- g(2)%nodes(1) - ycoor_u*  scaley)/delta_05
        VAUXPRE(4) = (g(2)%nodes(j)- g(2)%nodes(1) - ycoor_tem*scaley)/delta_t

        IF ( j .EQ. jmax/2 ) THEN
           ivauxdum = ivauxpos
        ELSE
           ivauxdum = 0
        ENDIF

        WRITE(i23,1100) i,j,(VAUXPRE(k), k=1,ivauxpre), &
! Grid&
             g(1)%nodes(i),g(2)%nodes(j),&
! Reynolds averages&
             rU(n,j),rV(n,j),rW(n,j),&
             rP(n,j),rR(n,j),rT(n,j), &
             rUf2(n,j),rVf2(n,j),rWf2(n,j),   &
             rPf2(n,j),rRf2(n,j),rTf2(n,j), &
             rUfVf(n,j),rUfWf(n,j),rVfWf(n,j), &
             rTKE(n,j),rbxx(n,j), rbyy(n,j), &
             rbzz(n,j),rbxy(n,j), rbxz(n,j), &
             rbyz(n,j),rRuT(n,j),rRvT(n,j),rRwT(n,j),&
! Favre averages&
             fU(n,j),fV(n,j),fW(n,j),fT(n,j),fTf2(n,j),   &
             fRxy(n,j),fRxz(n,j),fRyz(n,j), &
             fRxx(n,j),fRyy(n,j),fRzz(n,j), &
             fbxx(n,j),fbyy(n,j),fbzz(n,j),&
             fbxy(n,j),fbxz(n,j),fbyz(n,j),&
             fRuT(n,j),fRvT(n,j),fRwT(n,j),&
! Derivatives&
             fdUdx(n,j),fdUdy(n,j),fdUdz(n,j),&
             fdVdx(n,j),fdVdy(n,j),fdVdz(n,j),&
             fdWdx(n,j),fdWdy(n,j),fdWdz(n,j),&
             rdUdxf2(n,j),rdUdyf2(n,j),rdUdzf2(n,j),&
             rdVdxf2(n,j),rdVdyf2(n,j),rdVdzf2(n,j),&
             rdWdxf2(n,j),rdWdyf2(n,j),rdWdzf2(n,j),&
             rdVdxfdUdyf(n,j),rdWdxfdUdzf(n,j),rdWdyfdVdzf(n,j),&
             rdUdxfdVdyf(n,j),rdUdxfdWdzf(n,j),rdVdyfdWdzf(n,j),&
             dPdx(n,j),dPdy(n,j),dPdz(n,j),&
             dRdx(n,j),dRdy(n,j),dRdz(n,j),&
! Vorticity $ Dilation&
             Vortx(n,j),Vorty(n,j),Vortz(n,j),Dil(n,j),&
             fdUdx(n,j)+fdVdy(n,j)+fdWdz(n,j),&
             Vortxf2(n,j),Vortyf2(n,j),Vortzf2(n,j),Dilf2(n,j),&
! Scales&
             eta(n,j),tau(n,j),lambda(n,j), &
             lambda_x(n,j),lambda_y(n,j),lambda_z(n,j), &
             equi(n,j),&
! Rxx equation&
             SQRT(fRxx(n,j))/simuc(n),&
             Conv_xx(n,j),Prod_xx(n,j),Diss_xx(n,j),&
             Tran_xx(n,j),Pres_xx(n,j),MnFl_xx(n,j),Resi_xx(n,j),&
! Ryy equation&
             SQRT(fRyy(n,j))/simuc(n),&
             Conv_yy(n,j),Prod_yy(n,j),Diss_yy(n,j),&
             Tran_yy(n,j),Pres_yy(n,j),MnFl_yy(n,j),Resi_yy(n,j),&
! Rzz equation&
             SQRT(fRzz(n,j))/simuc(n),&
             Conv_zz(n,j),Prod_zz(n,j),Diss_zz(n,j),&
             Tran_zz(n,j),Pres_zz(n,j),MnFl_zz(n,j),Resi_zz(n,j),&
! Rxy equation&
             fRxy(n,j)/simuc(n)/simuc(n),&
             Conv_xy(n,j),Prod_xy(n,j),Diss_xy(n,j),&
             Tran_xy(n,j),Pres_xy(n,j),MnFl_xy(n,j),Resi_xy(n,j),&
! TKE equation&
             fTKE(n,j),Conv(n,j),Prod(n,j),Diss(n,j),&
             Tran(n,j),Pres(n,j),MnFl(n,j),Resi(n,j),&
! temperate equation&
             SQRT(fTf2(n,j))/ABS(simtc(n)),&
             Conv_tt(n,j),Prod_tt(n,j),&
             Diss_tt(n,j),Tran_tt(n,j),Pres_tt(n,j),&
             MnFl_tt(n,j),Resi_tt(n,j),&
! X-momentum equation&
             (fU(n,j)-U2)/simuc(n),Conv_u(n,j),Tran_u(n,j),&
             Reyn_u(n,j),Resi_u(n,j),&
! Y-momentum equation&
             fV(n,j)/simuc(n),Conv_v(n,j),Tran_v(n,j),&
             Reyn_v(n,j),Resi_v(n,j),&
! Z-momentum equation&
             fW(n,j)/simuc(n),Conv_w(n,j),Tran_w(n,j),&
             Reyn_w(n,j),Resi_w(n,j),&
! energy equation in p&
             (rP(n,j)-rP(n,1))/&
             (rP(n,jmax/2)-rP(n,1)),&
             Conv_p(n,j),Reve_p(n,j),Diss_p(n,j),&
             Tran_p(n,j),Reyn_p(n,j),&
             -(gama0-1)*rR(n,j)*Pres(n,j),Resi_p(n,j),&
! energy equation in T&
             (fT(n,j)-T2)/ABS(simtc(n)),&
             Conv_T(n,j),Reve_T(n,j),Diss_T(n,j),&
             Tran_T(n,j),Reyn_T(n,j),&
             -MRATIO*(gama0-1)*Pres(n,j)*&
             simuc(n)*simuc(n)/ABS(simtc(n)),Resi_T(n,j),&
! Filtering&
             fTKE_nf(n,j),eps_f(n,j),&
! Stress tensor&
             tau_xx(n,j),tau_yy(n,j),tau_zz(n,j),&
             tau_xy(n,j),tau_xz(n,j),tau_yz(n,j),&
             phi(n,j), rVis(n,j),&
! Variable density quantities&
             rho_p(n,j), rho_T(n,j),rho_ac(n,j), rho_en(n,j),&
             T_ac(n,j), T_en(n,j),&
             fRuT(n,j)/ABS(simtc(n)*simuc(n)),&
             fRvT(n,j)/ABS(simtc(n)*simuc(n)),&
             fRwT(n,j)/ABS(simtc(n)*simuc(n)),&
             (fU(n,j)-rU(n,j))*rR(n,j)/&
             ABS(simrc(n)*simuc(n)),&
             (fV(n,j)-rV(n,j))*rR(n,j)/&
             ABS(simrc(n)*simuc(n)),&
             (fW(n,j)-rW(n,j))*rR(n,j)/&
             ABS(simrc(n)*simuc(n)),&
! Skewness&Flatness&
             S_rho(n,j),S_u(n,j),S_v(n,j),S_w(n,j),&
             S_p(n,j),S_T(n,j),F_rho(n,j),F_u(n,j),&
             F_v(n,j),F_w(n,j),F_p(n,j),F_T(n,j),&
! 1D quantities&
             (VAUXPOS(k),k=1,ivauxdum)

     ENDDO
  ENDDO

  CLOSE(i23)

#undef simuc
#undef simtc
#undef simrc

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING AVG_FLOW_SPATIAL_LAYER' )
#endif

  RETURN

1010 FORMAT(A)
1100 FORMAT(I3,1X,I3,4(1X,E12.5E3),206(1X,E12.5E3),VARMX1D(1X,E12.5E3),2(1X,E12.5E3))

END SUBROUTINE AVG_FLOW_SPATIAL_LAYER
