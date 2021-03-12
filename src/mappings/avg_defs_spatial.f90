SUBROUTINE AVG_DEFS_SPATIAL

  IMPLICIT NONE

#include "types.h"

#ifdef USE_RECLEN
  OPEN(UNIT=22, RECL=1024, FILE='dns.def', STATUS='unknown')
#else
  OPEN(UNIT=22, FILE='dns.def', STATUS='unknown')
#endif

  WRITE(22,1010) '# Statistics File Definition & Description'
  WRITE(22,1010) '# Date = 10/03/2000'
  WRITE(22,1010) '# Written by Juan Pedro Mellado'
  WRITE(22,1010) '# Format <TEXT>:<XMGR>:<LATEX>:<FORMULA>'
  WRITE(22,1010) 'I:i:i:'
  WRITE(22,1010) 'J:j:j:'
  WRITE(22,1010) 'X:x/h:x/h:'
  WRITE(22,1010) 'Y:y/h:y/h:'
  WRITE(22,1010) 'SU:x\s2\N/\xd\s0.5:x_2/\delta_{0.5}:'
  WRITE(22,1010) 'SW:x\s2\N/\xd\sw\N\0:x_2/\delta_\omega:'
  WRITE(22,1010) 'SS:x\s2\N/\xd\4\ss\N:x_2/\delta_s:'
!  WRITE(22,1010) 'SL1:x\s2\N/\xd\s0.5:x_2/\delta_{0.5}:'
!  WRITE(22,1010) 'SL2:x\s2\N/\xd\sw\N\0:x_2/\delta_s:'

! Scales

  WRITE(22,1010) 'lambda:\xl\4:\lambda:'&
       //'\bigg(\frac{5\nu q^2}{\epsilon}\bigg)^{1/2}'&
       //'\qquad q^2=\overline{u^{''2}_1}+\overline{u^{''2}_2}'&
       //'+\overline{u^{''2}_3}'
  WRITE(22,1010) 'lambda_x:\xl\4\sx\N:\lambda_x:'&
       //'\sqrt{\overline{u^{''2}_1}\bigg/\overline'&
       //'{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_1}\bigg)^2}}'
  WRITE(22,1010) 'lambda_y:\xl\4\sy\N:\lambda_y:'&
       //'\sqrt{\overline{u^{''2}_2}\bigg/\overline'&
       //'{\bigg(\frac{\partial '&
       //'u^{''}_2}{\partial x_2}\bigg)^2}}'
  WRITE(22,1010) 'lambda_z:\xl\4\sz\N:\lambda_z:'&
       //'\sqrt{\overline{u^{''2}_3}\bigg/\overline'&
       //'{\bigg(\frac{\partial '&
       //'u^{''}_3}{\partial x_3}\bigg)^2}}'
  WRITE(22,1010) 'eta:\xh\4:\eta:'&
       //'\bigg(\frac{\nu^3}{\epsilon}\bigg)^{1/4}'
  WRITE(22,1010) 'tau:\xt\4:\tau:'&
       //'\bigg(\frac{\nu}{\epsilon}\bigg)^{1/2}'
  WRITE(22,1010) 'equi:\xk/e\4:\k/\e:'

  WRITE(22,1010) 'Re:Re\sx\N:Re_x:Re_h\frac{2\delta}{h}'&
       //'\frac{\Delta U_c}{\Delta U_0}\frac{\bar{\rho}}{\rho_0}'
  WRITE(22,1010) 'Re_l:Re\s\xl\4\N:Re_\lambda:Re_h\frac{\lambda}{h}'&
       //'\frac{q}{\Delta U_0}\frac{\bar{\rho}}{\rho_0}'


! Reynolds averages

  WRITE(22,1010) 'rR:\xr\4:\bar{\rho}:'
  WRITE(22,1010) 'rU:rU:\bar{u}_1:'
  WRITE(22,1010) 'rV:rV:\bar{u}_2:'
  WRITE(22,1010) 'rW:rW:\bar{u}_3:'
  WRITE(22,1010) 'rP:P:\bar{p}:'
  WRITE(22,1010) 'rT:rT:\bar{T}:'
  WRITE(22,1010) 'rS:r\xx\4:\bar{\xi}:'

  WRITE(22,1010) 'rUf2:u\S2\N\srms\N:\overline{u^{''2}}:'
  WRITE(22,1010) 'rVf2:v\S2\N\srms\N:\overline{v^{''2}}:'
  WRITE(22,1010) 'rWf2:w\S2\N\srms\N:\overline{w^{''2}}:'
  WRITE(22,1010) 'rPf2:p\S2\N\srms\N:\overline{p^{''2}}:'
  WRITE(22,1010) 'rRf2:\xr\4\S2\N\srms\N:\overline{\rho^{''2}}:'
  WRITE(22,1010) 'rTf2:T\S2\N\srms\N:\overline{T^{''2}}:'
  WRITE(22,1010) 'rSf2:r\xx\4\S2\N\srms\N:\overline{\xi^{''2}}:'
  WRITE(22,1010) 'rUfVf:rR\sxy\N:\overline{u^{''}v^{''}}:'
  WRITE(22,1010) 'rUfWf:rR\sxz\N:\overline{u^{''}w^{''}}:'
  WRITE(22,1010) 'rVfWf:rR\syz\N:\overline{v^{''}w^{''}}:'
  WRITE(22,1010) 'rUfSf:rR\su\xx\4\N:\overline{u^{''}\xi^{''}}:'
  WRITE(22,1010) 'rVfSf:rR\sv\xx\4\N:\overline{v^{''}\xi^{''}}:'
  WRITE(22,1010) 'rWfSf:rR\sw\xx\4\N:\overline{w^{''}\xi^{''}}:'
  WRITE(22,1010) 'rTKE:k:k:'

  WRITE(22,1010) 'rbxx:rb\sxx\N:rb_{xx}:'&
       //'\frac{\overline{u^{''2}}}{q^2}-\frac{1}{3}'
  WRITE(22,1010) 'rbyy:rb\syy\N:rb_{yy}:'&
       //'\frac{\overline{v^{''2}}}{q^2}-\frac{1}{3}'
  WRITE(22,1010) 'rbzz:rb\szz\N:rb_{zz}:'&
       //'\frac{\overline{w^{''2}}}{q^2}-\frac{1}{3}'
  WRITE(22,1010) 'rbxy:rb\sxy\N:rb_{xy}:'&
       //'\frac{\overline{u^{''}v^{''}}}{q^2}'
  WRITE(22,1010) 'rbxz:rb\sxz\N:rb_{xz}:'&
       //'\frac{\overline{u^{''}w^{''}}}{q^2}'
  WRITE(22,1010) 'rbyz:rb\syz\N:rb_{yz}:'&
       //'\frac{\overline{v^{''}w^{''}}}{q^2}'

! Favre averages

  WRITE(22,1010) 'fU:fU:\tilde{u}_1:'
  WRITE(22,1010) 'fV:fV:\tilde{u}_2:'
  WRITE(22,1010) 'fW:fW:\tilde{u}_3:'
  WRITE(22,1010) 'fT:fT:\tilde{T}:'
  WRITE(22,1010) 'fS:f\xx\0:\tilde{\xi}:'
  WRITE(22,1010) 'fTf2:fT\S2\N\srms\N:\overline{T^{''''2}}:'&
       //'\frac{\overline{\rho T^{''''2}}}{\bar{\rho}}'
  WRITE(22,1010) 'fRxz:R\sxz\N:R_{13}:'&
       //'\frac{\overline{\rho u^{''''}_1 u^{''''}_3}}{\bar{\rho}}'
  WRITE(22,1010) 'fRyz:R\syz\N:R_{23}:'&
       //'\frac{\overline{\rho u^{''''}_2 u^{''''}_3}}{\bar{\rho}}'

  WRITE(22,1010) 'fbxx:b\sxx\N:b_{xx}:'&
       //'\frac{R_{xx}}{2k}-\frac{1}{3}'
  WRITE(22,1010) 'fbyy:b\syy\N:b_{yy}:'&
       //'\frac{R_{yy}}{2k}-\frac{1}{3}'
  WRITE(22,1010) 'fbzz:b\szz\N:b_{zz}:'&
       //'\frac{R_{zz}}{2k}-\frac{1}{3}'
  WRITE(22,1010) 'fbxy:b\sxy\N:b_{xy}:\frac{R_{xy}}{2k}'
  WRITE(22,1010) 'fbxz:b\sxz\N:b_{xz}:\frac{R_{xz}}{2k}'
  WRITE(22,1010) 'fbyz:b\syz\N:b_{yz}:\frac{R_{yz}}{2k}'

! R11 Budget Equation

  WRITE(22,1010) 'sRxx:R\sxx\N\S1/2\N:R_{11}:'&
       //'\frac{\overline{\rho u^{''''2}_1}}{\bar{\rho}}'
  WRITE(22,1010) 'Resi_xx:R\sxx,t\N:'&
       //'\frac{\partial R_{11}}{\partial t}:'&
       //'C_{11} + P_{11} - \epsilon_{11} + (-T_{11k,k}+\Pi_{11}+'&
       //'\Sigma_{11})/\bar{\rho}'
  WRITE(22,1010) 'Conv_xx:C\sxx\N:C_{11}:'&
       //'-\tilde{u}_k \frac{\partial R_{11}}{\partial x_k}'
  WRITE(22,1010) 'Prod_xx:P\sxx\N:P_{11}:'&
       //'-2 R_{1k} \frac{\partial \tilde{u}_1}{\partial x_k}'
  WRITE(22,1010) 'Diss_xx:-\xe\4\sxx\N:\epsilon_{11}:'&
       //'\frac{2}{\bar{\rho}} \overline{ \tau^{''}_{1k}'&
       //'\frac{\partial u^{''}_1}{\partial x_k}}'
  WRITE(22,1010) 'Tran_xx:-T\sxxk,k\N/\xr\4:-T_{11k,k}/\rho:'&
       //'T_{11k}=\overline{\rho u^{''''2}_1 u^{''''}_k} '&
       //'- 2 \overline{\tau^{''}_{1k} u^{''''}_1}'&
       //'+ 2 \overline{p^{''}u^{''''}_1}\delta_{1k}'
  WRITE(22,1010) 'Pres_xx:\xP\4\sxx\N/\xr\4:\Pi_{11}/\rho:'&
       //'2 \overline{p^{''}\frac{\partial u^{''}_1}{\partial x_1}}'
  WRITE(22,1010) 'MnFl_xx:\xS\4\sxx\N/\xr\4:\Sigma_{11}/\rho:'&
       //'2 \overline{u^{''''}_1}\frac{\partial \bar{\tau}_{1k}}{'&
       //'\partial x_k}-\overline{u^{''''}_1}\frac{\partial '&
       //'\bar{p}}{\partial x_1}'

! R22 Budget Equation

  WRITE(22,1010) 'sRyy:R\syy\N\S1/2\N:R_{22}:'&
       //'\frac{\overline{\rho u^{''''2}_2}}{\bar{\rho}}'
  WRITE(22,1010) 'Resi_yy:R\syy,t\N:'&
       //'\frac{\partial R_{22}}{\partial t}:'&
       //'C_{22} + P_{22} - \epsilon_{22} + (-T_{22k,k}+\Pi_{22}+'&
       //'\Sigma_{22})/\bar{\rho}'
  WRITE(22,1010) 'Conv_yy:C\syy\N:C_{22}:'&
       //'-\tilde{u}_k \frac{\partial R_{22}}{\partial x_k}'
  WRITE(22,1010) 'Prod_yy:P\syy\N:P_{22}:'&
       //'-2 R_{2k} \frac{\partial \tilde{u}_2}{\partial x_k}'
  WRITE(22,1010) 'Diss_yy:-\xe\4\syy\N:\epsilon_{22}:'&
       //'\frac{2}{\bar{\rho}} \overline{ \tau^{''}_{2k}'&
       //'\frac{\partial u^{''}_2}{\partial x_k}}'
  WRITE(22,1010) 'Tran_yy:-T\syyk,k\N/\xr\4:-T_{22k,k}/\rho:'&
       //'T_{22k}=\overline{\rho u^{''''2}_1 u^{''''}_k} '&
       //'- 2 \overline{\tau^{''}_{2k} u^{''''}_2}'&
       //'+ 2 \overline{p^{''}u^{''''}_2}\delta_{2k}'
  WRITE(22,1010) 'Pres_yy:\xP\4\syy\N/\xr\4:\Pi_{22}/\rho:'&
       //'2 \overline{p^{''}\frac{\partial u^{''}_2}{\partial x_2}}'
  WRITE(22,1010) 'MnFl_yy:\xS\4\syy\N/\xr\4:\Sigma_{22}/\rho:'&
       //'2 \overline{u^{''''}_2}\frac{\partial \bar{\tau}_{2k}}{'&
       //'\partial x_k}-\overline{u^{''''}_2}\frac{\partial '&
       //'\bar{p}}{\partial x_2}'

! R33 Budget Equation

  WRITE(22,1010) 'sRzz:R\szz\N\S1/2\N:R_{33}:'&
       //'\frac{\overline{\rho u^{''''2}_2}}{\bar{\rho}}'
  WRITE(22,1010) 'Resi_zz:R\szz,t\N:'&
       //'\frac{\partial R_{22}}{\partial t}:'&
       //'C_{33} + P_{33} - \epsilon_{33} + (-T_{33k,k}+\Pi_{33}+'&
       //'\Sigma_{33})/\bar{\rho}'
  WRITE(22,1010) 'Conv_zz:C\szz\N:C_{33}:'&
       //'-\tilde{u}_k \frac{\partial R_{33}}{\partial x_k}'
  WRITE(22,1010) 'Prod_zz:P\szz\N:P_{33}:'&
       //'-2 R_{3k} \frac{\partial \tilde{u}_3}{\partial x_k}'
  WRITE(22,1010) 'Diss_zz:-\xe\4\szz\N:\epsilon_{33}:'&
       //'\frac{2}{\bar{\rho}} \overline{ \tau^{''}_{3k}'&
       //'\frac{\partial u^{''}_3}{\partial x_k}}'
  WRITE(22,1010) 'Tran_zz:-T\szzk,k\N/\xr\4:-T_{33k,k}/\rho:'&
       //'T_{33k}=\overline{\rho u^{''''3}_1 u^{''''}_k} '&
       //'- 2 \overline{\tau^{''}_{3k} u^{''''}_3}'&
       //'+ 2 \overline{p^{''}u^{''''}_3}\delta_{3k}'
  WRITE(22,1010) 'Pres_zz:\xP\4\szz\N/\xr\4:\Pi_{33}/\rho:'&
       //'2 \overline{p^{''}\frac{\partial u^{''}_3}{\partial x_3}}'
  WRITE(22,1010) 'MnFl_zz:\xS\4\szz\N/\xr\4:\Sigma_{33}/\rho:'&
       //'2 \overline{u^{''''}_3}\frac{\partial \bar{\tau}_{3k}}{'&
       //'\partial x_k}-\overline{u^{''''}_3}\frac{\partial '&
       //'\bar{p}}{\partial x_3}'

! R12 Budget Equation

  WRITE(22,1010) 'sRxy:R\sxy\N:R_{12}:'&
       //'\frac{\overline{\rho u^{''''}_1 u^{''''}_2}}{\bar{\rho}}'
  WRITE(22,1010) 'Resi_xy:R\sxy,t\N:'&
       //'\frac{\partial R_{12}}{\partial t}:'&
       //'C_{12} + P_{12} - \epsilon_{12} + (-T_{12k,k}+\Pi_{12}+'&
       //'\Sigma_{12})/\bar{\rho}'
  WRITE(22,1010) 'Conv_xy:C\sxy\N:C_{12}:'&
       //'-\tilde{u}_k \frac{\partial R_{12}}{\partial x_k}'
  WRITE(22,1010) 'Prod_xy:P\sxy\N:P_{12}:'&
       //'-\left( R_{1k} \frac{\partial \tilde{u}_2}{\partial x_k}+'&
       //' R_{2k} \frac{\partial \tilde{u}_1}{\partial x_k} \right)'
  WRITE(22,1010) 'Diss_xy:-\xe\4\sxy\N:\epsilon_{12}:'&
       //'\frac{1}{\bar{\rho}} \bigg(\overline{ \tau^{''}_{1k}'&
       //'\frac{\partial u^{''}_2}{\partial x_k}} + \overline{ '&
       //'\tau^{''}_{2k}'&
       //'\frac{\partial u^{''}_1}{\partial x_k}}\bigg)'
  WRITE(22,1010) 'Tran_xy:-T\sxyk,k\N/\xr\4:-T_{12k,k}/\rho:'&
       //'T_{12k}=\overline{\rho u^{''''}_1 u^{''''}_2 '&
       //'u^{''''}_k} + '&
       //'\overline{p^{''} u^{''''}_1}\delta_{2k} +'&
       //'\overline{p^{''} u^{''''}_2}\delta_{1k} -'&
       //'(\overline{\tau^{''}_{1k} u^{''''}_2}+'&
       //' \overline{\tau^{''}_{2k} u^{''''}_1})'
  WRITE(22,1010) 'Pres_xy:\xP\4\sxy\N/\xr\4:\Pi_{12}/\rho:'&
       //'\overline{p^{''}\bigg(\frac{\partial u^{''}_1}{'&
       //'\partial x_2}+ \frac{\partial u^{''}_2}{'&
       //'\partial x_1} \bigg)}'
  WRITE(22,1010) 'MnFl_xy:\xS\4\sxy\N/\xr\4:\Sigma_{12}/\rho:'&
       //'\overline{u^{''''}_1}\frac{\partial \bar{\tau}_{2k}}{'&
       //'\partial x_k}+ \overline{u^{''''}_2}\frac{\partial '&
       //'\bar{\tau}_{1k}}{\partial x_k} -'&
       //'\bigg(\overline{u^{''''}_1}\frac{\partial \bar{p}}'&
       //'{\partial x_2} + \overline{u^{''''}_2}'&
       //'\frac{\partial \bar{p}}{\partial x_1} \bigg)'

! Kinetic Energy

  WRITE(22,1010) 'fTKE:k:k:'&
       //'(R_{11} + R_{22} + R_{33})/2'
  WRITE(22,1010) 'Resi:k\s,t\N:'&
       //'\frac{\partial k}{\partial t}:'&
       //'C + P - \epsilon + (-T_{k,k}+\Pi+\Sigma)/\bar{\rho}'
  WRITE(22,1010) 'Conv:C:C:'&
       //'(C_{11}+C_{22}+C_{33})/2'
  WRITE(22,1010) 'Prod:P:P:'&
       //'(P_{11}+P_{22}+P_{33})/2'
  WRITE(22,1010) 'Diss:-\xe\4:\epsilon:'&
       //'(\epsilon_{11}+\epsilon_{22}+\epsilon_{33})/2'
  WRITE(22,1010) 'Pres:\xP\4/\xr\4:\Pi/\rho:'&
       //'(\Pi_{11}+\Pi_{22}+\Pi_{33})/2'
  WRITE(22,1010) 'Tran:-T\sk,k\N/\xr\4:T_{k,k}/\rho:'&
       //'T_k=(T_{11k}+T_{22k}+T_{33k})/2'
  WRITE(22,1010) 'MnFl:\xS\4/\xr\4:\Sigma/\rho:'&
       //'(\Sigma_{11}+\Sigma_{22}+\Sigma_{33})/2'

! Filtering

  WRITE(22,1010) 'fTKE_nf:k\sbf\N:k_{bf}:'&
       //'(R_{11,bf} + R_{22,bf} + R_{33,bf})/2'
  WRITE(22,1010) 'eps_f:\xe\4\sf\N/\xe\4:\epsilon_f/\epsilon:'&
       //'\frac{k-k_{bf}}{\epsilon \Delta t}' 

! Shear stresses

  WRITE(22,1010) 'tau_xx:\xt\4\sxx\N:\overline{\tau_{xx}}:'
  WRITE(22,1010) 'tau_yy:\xt\4\syy\N:\overline{\tau_{yy}}:'
  WRITE(22,1010) 'tau_zz:\xt\4\szz\N:\overline{\tau_{zz}}:'
  WRITE(22,1010) 'tau_xy:\xt\4\sxy\N:\overline{\tau_{xy}}:'
  WRITE(22,1010) 'tau_xz:\xt\4\sxz\N:\overline{\tau_{xz}}:'
  WRITE(22,1010) 'tau_yz:\xt\4\syz\N:\overline{\tau_{yz}}:'

  WRITE(22,1010) 'Phi:\xf\0:\bar{\phi}:'&
       //'\overline{\tau_{lm}\frac{\partial u_l}{\partial x_m}}'

! Similarity variables

  WRITE(22,1010) 'sU:\xD\4U/\xD\4U\sc\N:\Delta U/\Delta U_c:'
  WRITE(22,1010) 'sV:V/\xD\4U\sc\N:\Delta U/\Delta U_c:'
  WRITE(22,1010) 'sW:W/\xD\4U\sc\N:\Delta U/\Delta U_c:'
  WRITE(22,1010) 'sS:\xDx\4/\xDx\4\sc\N:\Delta \xi/\Delta \xi_c:'

  WRITE(22,1010) 'SimSC:(\xx\4\s0\N/\xx\4\sc\N)\S2\N:'&
       //'\Delta \xi/\Delta \xi_c:'

! Scalar Budget

  WRITE(22,1010) 'sRss:R\s\xxx\4\N\S1/2\N:R_{\xi\xi}:'&
       //'\frac{\overline{\rho \xi^{''''2}}}{\bar{\rho}}'
  WRITE(22,1010) 'Rsu:R\s\xx\4u\N:R_{\xi 1}:'&
       //'\frac{\overline{\rho \xi^{''''} u^{''''}_1}}{\bar{\rho}}'
  WRITE(22,1010) 'Rsv:R\s\xx\4v\N:R_{\xi 2}:'&
       //'\frac{\overline{\rho \xi^{''''} u^{''''}_2}}{\bar{\rho}}'
  WRITE(22,1010) 'Rsw:R\s\xx\4w\N:R_{\xi 3}:'&
       //'\frac{\overline{\rho \xi^{''''} u^{''''}_3}}{\bar{\rho}}'

! Dilatation

  WRITE(22,1010) 'Vortx:\xw\4\sx\N/'&
       //'(\xD\4U\sc\N/\xd\4\s0.5\N):\omega_{x}/'&
       //'(\Delta U_c/\delta_{0.5}):'
  WRITE(22,1010) 'Vorty:\xw\4\sy\N/'&
       //'(\xD\4U\sc\N/\xd\4\s0.5\N):\omega_{y}/'&
       //'(\Delta U_c/\delta_{0.5}):'
  WRITE(22,1010) 'Vortz:\xw\4\sz\N/'&
       //'(\xD\4U\sc\N/\xd\4\s0.5\N):\omega_{z}/'&
       //'(\Delta U_c/\delta_{0.5}):'
  WRITE(22,1010) 'Dil:\9V\4.\5u\4:\nabla . \overline{u}:'
  WRITE(22,1010) 'Vortxf2:\xw\4\sx,rms\N/'&
       //'(\xD\4U\sc\N/\xd\4\s0.5\N):\omega_{x,rms}/'&
       //'(\Delta U_c/\delta_{0.5}):'
  WRITE(22,1010) 'Vortyf2:\xw\4\sy,rms\N/'&
       //'(\xD\4U\sc\N/\xd\4\s0.5\N):\omega_{x,rms}/'&
       //'(\Delta U_c/\delta_{0.5}):'
  WRITE(22,1010) 'Vortzf2:\xw\4\sz,rms\N/'&
       //'(\xD\4U\sc\N/\xd\4\s0.5\N):\omega_{x,rms}/'&
       //'(\Delta U_c/\delta_{0.5}):'
  WRITE(22,1010) 'Dilf2:(\9V\4.u'')\S2\N:'&
       //'\overline{(\nabla . u^{''})^2}:'

! Others

  WRITE(22,1010) 'F_y:F\sy\N:F_2:'&
       //'\frac{1}{Re Sc}\overline{\Gamma^{*} \frac{\partial \xi}{'&
       //'\partial x_2}}'
  WRITE(22,1010) 'Tp_y:Tp\sy\N:T_{p,2}:'&
       //'\frac{\partial}{\partial x_2}\bigg( '&
       //'\overline{p^{''}u^{''}_2} - \frac{\gamma}{Re Pr}'&
       //'\frac{\partial }{\partial x_2}\bigg( '&
       //'\overline{\frac{p}{\rho}}\bigg)\bigg)'
  WRITE(22,1010) 'rRP:rRP:\overline{\rho^{''}p^{''}}:'
  WRITE(22,1010) 'rRT:rRT:\overline{\rho^{''}T^{''}}:'
  WRITE(22,1010) 'Pp2:P\sp\N:P_p:'&
       //'- 2 \bigg( \overline{p^{''}u^{''}_2}' &
       //'\frac{\partial \bar{p}}{\partial x_2} + \gamma'&
       //'\overline{p^{''2}}\frac{\partial \bar{u}_2 '&
       //'}{\partial x_2} \bigg)'
  WRITE(22,1010) 'Ep2:\xe\0\sp\N:\epsilon_p:'&
       //'\frac{2 \gamma}{Re Pr} \overline{\frac{1}{\rho}'&
       //'\frac{\partial p^{''}}{\partial x_k}\frac{\partial'&
       //' p^{''} }{\partial x_k}}'
  WRITE(22,1010) 'Tp2_y:T\spy,2\N:T_{p2,2}:'&
       //'\frac{\partial}{\partial x_2} \bigg( '&
       //'\overline{u^{''}_2 p^{''2}} - \frac{2 \gamma}{Re Pr}' &
       //'\overline{p^{''}\frac{\partial}{\partial x_2}'&
       //'\bigg(\frac{p}{\rho}\bigg)}\bigg)'
  WRITE(22,1010) 'Pip2:\xP\0\sp\N:\Pi_p:'&
       //'\overline{p^{''2}\bigg(\frac{\partial u^{''}_1}{'&
       //'\partial x_1}+ \frac{\partial u^{''}_2}{'&
       //'\partial x_2} + \frac{\partial u^{''}_3}{'&
       //'\partial x_3}\bigg)}'
  WRITE(22,1010) 'PPhi:p\xf\0:\overline{p^{''}\phi}:'
  WRITE(22,1010) 'Rho_ac:\xr\0\sac\N:\rho^{''}_{ac}:'
  WRITE(22,1010) 'Rho_en:\xr\0\sen\N:\rho^{''}_{en}:'
  WRITE(22,1010) 'T_ac:T\sac\N:T^{''}_{ac}:'
  WRITE(22,1010) 'T_en:T\sen\N:T^{''}_{en}:'
  WRITE(22,1010) 'RhoFluxX:\xr\0\S''\Nu\S''\Nu:'&
       //'\overline{\rho^{''}u^{''}}:'
  WRITE(22,1010) 'RhoFluxY:\xr\0\S''\Nv\S''\Nu:'&
       //'\overline{\rho^{''}v^{''}}:'
  WRITE(22,1010) 'RhoFluxZ:\xr\0\S''\Nw\S''\Nu:'&
       //'\overline{\rho^{''}w^{''}}:'
  WRITE(22,1010) 'Sp2:\xS\0\sp\N:\Sigma_p:'&
       //'\frac{2 \gamma}{Re Pr} \bigg(\overline{\frac{p}'&
       //'{\rho^2}\frac{\partial \rho}{\partial x_k}'&
       //'\frac{\partial p^{''}}{\partial x_k}} - '&
       //'\overline{\frac{1}{\rho}\frac{\partial p^{''}}{'&
       //'\partial x_k}} \frac{\partial \bar{p}}{\partial x_k}'&
       //'\bigg)'

  WRITE(22,1010) 'Skew_x:Skew\sx\N:Sk_1:'&
       //'\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_1}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_1}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'Skew_y:Skew\sy\N:Sk_2:'&
       //'\overline{\bigg(\frac{\partial u^{''}_2}{\partial x_2}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_2}{\partial x_2}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'Skew_z:Skew\sz\N:Sk_3:'&
       //'\overline{\bigg(\frac{\partial u^{''}_3}{\partial x_3}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_3}{\partial x_3}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'Entropy:s:\bar{s}:'&
       //'\overline{\log{p/p_o}-\gamma \log{\rho/\rho_o}}'
  WRITE(22,1010) 'Entropy2:s\S2\N\srms\N:\overline{s^{''2}}:'

! Growth Rates and Thickness

  WRITE(22,1010) 'Delta_m_u:\xd\sq\N\4\Su\N:\delta_\theta:'&
       //'\frac{1}{\rho_o}\int\bar{\rho}(1/4-\tilde{u}^2_1/'&
       //'\Delta u^2) d x_2'
  WRITE(22,1010) 'Delta_m_d:\xd\sq\N\4\Sd\N:\delta_\theta:'&
       //'\frac{1}{\rho_o}\int\bar{\rho}(1/4-\tilde{u}^2_1/'&
       //'\Delta u^2) d x_2'
  WRITE(22,1010) 'Delta_w_u:\xd\sw\N\4\Su\N:\delta_\omega:'&
       //'\frac{\Delta U_c}{|\frac{d\tilde{u}_1}{dx_2}|_{max}}'
  WRITE(22,1010) 'Delta_w_d:\xd\sw\N\4\Sd\N:\delta_\omega:'&
       //'\frac{\Delta U_c}{|\frac{d\tilde{u}_1}{dx_2}|_{max}}'
  WRITE(22,1010) 'Delta_05_u:\xd\s0.5\N\4\Su\N:\delta_{0.5}:'&
       //'\tilde{u}_1(\delta_{0.5})=\frac{U_2+U_c}{2}'
  WRITE(22,1010) 'Delta_05_d:\xd\s0.5\N\4\Sd\N:\delta_{0.5}:'&
       //'\tilde{u}_1(\delta_{0.5})=\frac{U_2+U_c}{2}'
  WRITE(22,1010) 'Delta_s_u:\xd\sx\N\4\Su\N:\delta_\xi:'&
       //'\tilde{\xi}_1(\delta_{\xi})=\frac{\xi_2+\xi_c}{2}'
  WRITE(22,1010) 'Delta_s_d:\xd\sx\N\4\Sd\N:\delta_\xi:'&
       //'\tilde{\xi}_1(\delta_{\xi})=\frac{\xi_2+\xi_c}{2}'
  WRITE(22,1010) 'ExcMom:J:J:'&
       //'\int\bar{\rho}\tilde{u}_1(\tilde{u}_1-U_2)'
  WRITE(22,1010) 'SimUc:(\xD\4U\s0\N/\xD\4U\sc\N\4)\S2\N:'&
       //'\left(\frac{\Delta U_{c,0}}{\Delta U_c}\right)^2'
  WRITE(22,1010) 'Center:y\sc\N:y_c:'

1010 FORMAT(A)
  CLOSE(22)

  RETURN
END SUBROUTINE AVG_DEFS_SPATIAL

