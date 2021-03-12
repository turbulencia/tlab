#include "types.h"
#include "dns_const.h"

SUBROUTINE AVG_DEFS_TEMPORAL

  IMPLICIT NONE

#ifdef USE_RECLEN
  OPEN(UNIT=22, RECL=1024, FILE='dns.def', STATUS='unknown')
#else
  OPEN(UNIT=22, FILE='dns.def', STATUS='unknown')
#endif

  WRITE(22,1010) '# Statistics File Definition & Description'
  WRITE(22,1010) '# Format <TEXT>:<XMGR>:<LATEX>:<FORMULA>'
  WRITE(22,1010) 'I:i:i:'
  WRITE(22,1010) 'J:j:j:'
  WRITE(22,1010) 'Y:y:x_2:'
  WRITE(22,1010) 'SM:(y-y\sc\N)/\xd\sq\N\0:(x_2-x_{2,c})/\delta_\theta:'
  WRITE(22,1010) 'SW:(y-y\sc\N)/\xd\sw\N\0:(x_2-x_{2,c})/\delta_\omega:'
  WRITE(22,1010) 'SR:(y-y\sc\N)/h         :(x_2-x_{2,c})/h:'
  WRITE(22,1010) 'SS:(y-y\sc\N)/\xd\sz\N\0:(x_2-x_{2,c})/\delta_\zeta:'

  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'rR:\x\ca\Cr\cq\C\0    :\bar{\rho}:'
  WRITE(22,1010) 'rU:\x\ca\C\1u\x\cq\C\0:\bar{u}_1:'
  WRITE(22,1010) 'rV:\x\ca\C\1v\x\cq\C\0:\bar{u}_2:'
  WRITE(22,1010) 'rW:\x\ca\C\1w\x\cq\C\0:\bar{u}_3:'
  WRITE(22,1010) 'rP:\x\ca\C\1p\x\cq\C\0:\bar{p}:'
  WRITE(22,1010) 'rT:\x\ca\C\1T\x\cq\C\0:\bar{T}:'
  WRITE(22,1010) 're:\x\ca\C\1e\x\cq\C\0:\bar{e}:'
  WRITE(22,1010) 'rh:\x\ca\C\1h\x\cq\C\0:\bar{h}:'
  WRITE(22,1010) 'rs:\x\ca\C\1s\x\cq\C\0:\bar{s}:'
  WRITE(22,1010) 'rB:\x\ca\C\1b\x\cq\C\0:\bar{b}:'
  WRITE(22,1010) 'fU:\x\ca\C\1u\x\cq\C\0:\tilde{u}_1:'
  WRITE(22,1010) 'fV:\x\ca\C\1v\x\cq\C\0:\tilde{u}_2:'
  WRITE(22,1010) 'fW:\x\ca\C\1w\x\cq\C\0:\tilde{u}_3:'
  WRITE(22,1010) 'fT:\x\ca\C\1T\x\cq\C\0:\tilde{T}:'
  WRITE(22,1010) 'fe:\x\ca\C\1e\x\cq\C\0:\tilde{e}:'
  WRITE(22,1010) 'fh:\x\ca\C\1h\x\cq\C\0:\tilde{h}:'
  WRITE(22,1010) 'fs:\x\ca\C\1s\x\cq\C\0:\tilde{s}:'

  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Tke:\1TKE\0    :K     :(R_{11}+R_{22}+R_{33})/2'
  WRITE(22,1010) 'Rxx:\1R\0\sxx\N:R_{11}:'
  WRITE(22,1010) 'Ryy:\1R\0\syy\N:R_{22}:'
  WRITE(22,1010) 'Rzz:\1R\0\szz\N:R_{33}:'
  WRITE(22,1010) 'Rxy:\1R\0\sxy\N:R_{12}:'
  WRITE(22,1010) 'Rxz:\1R\0\sxz\N:R_{13}:'
  WRITE(22,1010) 'Ryz:\1R\0\syz\N:R_{23}:'

  WRITE(22,1010) 'rP2:\1p\0\S2\N\srms\N:\overline{p^{''2}}:'
  WRITE(22,1010) 'rR2:\xr\0\S2\N\srms\N:\overline{\rho^{''2}}:'
  WRITE(22,1010) 'rT2:\1T\0\S2\N\srms\N:\overline{T^{''2}}:'
  WRITE(22,1010) 'fT2:\1T\0\S2\N\srms\N:\widetilde{T^{''''2}}:'
  WRITE(22,1010) 're2:\1e\0\S2\N\srms\N:\overline{e^{''2}}:'
  WRITE(22,1010) 'fe2:\1e\0\S2\N\srms\N:\widetilde{e^{''''2}}:'
  WRITE(22,1010) 'rh2:\1h\0\S2\N\srms\N:\overline{h^{''2}}:'
  WRITE(22,1010) 'fh2:\1h\0\S2\N\srms\N:\widetilde{h^{''''2}}:'
  WRITE(22,1010) 'rs2:\1s\0\S2\N\srms\N:\overline{s^{''2}}:'
  WRITE(22,1010) 'fs2:\1s\0\S2\N\srms\N:\widetilde{s^{''''2}}:'

! Vorticity
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Wx :\x\ca\Cw\0\sx\N\x\cq\C\0:\overline{\omega}_1:'&
       //'\partial \bar{u}_3/\partial x_2 -\partial \bar{u}_2/\partial x_3'
  WRITE(22,1010) 'Wy :\x\ca\Cw\0\sy\N\x\cq\C\0:\overline{\omega}_2:'&
       //'\partial \bar{u}_1/\partial x_3 -\partial \bar{u}_3/\partial x_1'
  WRITE(22,1010) 'Wz :\x\ca\Cw\0\sz\N\x\cq\C\0:\overline{\omega}_3:'&
       //'\partial \bar{u}_2/\partial x_1 -\partial \bar{u}_1/\partial x_2'
  WRITE(22,1010) 'Wx2:(\xw\0\sx\N)\S2\N\srms\N:\overline{\omega_1^{''2}}:'
  WRITE(22,1010) 'Wy2:(\xw\0\sy\N)\S2\N\srms\N:\overline{\omega_2^{''2}}:'
  WRITE(22,1010) 'Wz2:(\xw\0\sz\N)\S2\N\srms\N:\overline{\omega_3^{''2}}:'

! R11 Budget Equation
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Rxx_t:R\sxx,t\N:\partial R_{11} / \partial t:'&
       //'C_{11} -F_{11} +P_{11} +B_{11} -\varepsilon_{11} + (\Pi_{11}-T_{112,2}+D_{11}-G_{11})/\bar{\rho}'
  WRITE(22,1010) 'Cxx:C\sxx\N:C_{11}:'
  WRITE(22,1010) 'Fxx:F\sxx\N:F_{11}:'
  WRITE(22,1010) 'Pxx:P\sxx\N:P_{11}:'
  WRITE(22,1010) 'Bxx:B\sxx\N:B_{11}:'
  WRITE(22,1010) 'Exx:\xe\0\sxx\N:\varepsilon_{11}:'
  WRITE(22,1010) 'PIxx:\xP\0\sxx\N:\Pi_{11}:'
  WRITE(22,1010) 'Txxy_y:T\sxxy,y\N:\partial T_{112}/\partial x_2:'
  WRITE(22,1010) 'Txxy:T\sxxy\N:T_{112}:'
  WRITE(22,1010) 'Gxx:G\sxx\N:G_{11}:'
  WRITE(22,1010) 'Dxx:D\sxx\N:D_{11}:'

! R22 Budget Equation
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Ryy_t:R\syy,t\N:\partial R_{22} / \partial t:'&
       //'C_{22} -F_{22} +P_{22} +B_{22} -\varepsilon_{22} + (\Pi_{22}-T_{222,2}+D_{22}-G_{22})/\bar{\rho}'
  WRITE(22,1010) 'Cyy:C\syy\N:C_{22}:'
  WRITE(22,1010) 'Fyy:F\syy\N:F_{22}:'
  WRITE(22,1010) 'Pyy:P\syy\N:P_{22}:'
  WRITE(22,1010) 'Byy:B\syy\N:B_{22}:'
  WRITE(22,1010) 'Eyy:\xe\0\syy\N:\varepsilon_{22}:'
  WRITE(22,1010) 'PIyy:\xP\0\syy\N:\Pi_{22}:'
  WRITE(22,1010) 'Tyyy_y:T\syyy,y\N:\partial T_{222}/\partial x_2:'
  WRITE(22,1010) 'Tyyy:T\syyy\N:T_{222}:'
  WRITE(22,1010) 'Gyy:G\syy\N:G_{22}:'
  WRITE(22,1010) 'Dyy:D\syy\N:D_{22}:'

! R33 Budget Equation
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Rzz_t:R\szz,t\N:\partial R_{33} / \partial t:'&
       //'C_{33} -F_{33} +P_{33} +B_{33} -\varepsilon_{33} + (\Pi_{33}-T_{332,2}+D_{33}-G_{33})/\bar{\rho}'
  WRITE(22,1010) 'Czz:C\szz\N:C_{33}:'
  WRITE(22,1010) 'Fzz:F\szz\N:F_{33}:'
  WRITE(22,1010) 'Pzz:P\szz\N:P_{33}:'
  WRITE(22,1010) 'Bzz:B\szz\N:B_{33}:'
  WRITE(22,1010) 'Ezz:\xe\0\szz\N:\varepsilon_{33}:'
  WRITE(22,1010) 'PIzz:\xP\0\szz\N:\Pi_{33}:'
  WRITE(22,1010) 'Tzzy_y:T\szzy,y\N:\partial T_{332}/\partial x_2:'
  WRITE(22,1010) 'Tzzy:T\szzy\N:T_{332}:'
  WRITE(22,1010) 'Gzz:G\szz\N:G_{33}:'
  WRITE(22,1010) 'Dzz:D\szz\N:D_{33}:'

! R12 Budget Equation
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Rxy_t:R\sxy,t\N:\partial R_{12} / \partial t:'&
       //'C_{12} -F_{12} +P_{12} +B_{12} -\varepsilon_{12} + (\Pi_{12}-T_{122,2}+D_{12}-G_{12})/\bar{\rho}'
  WRITE(22,1010) 'Cxy:C\sxy\N:C_{12}:'
  WRITE(22,1010) 'Fxy:F\sxy\N:F_{12}:'
  WRITE(22,1010) 'Pxy:P\sxy\N:P_{12}:'
  WRITE(22,1010) 'Bxy:B\sxy\N:B_{12}:'
  WRITE(22,1010) 'Exy:\xe\0\sxy\N:\varepsilon_{12}:'
  WRITE(22,1010) 'PIxy:\xP\0\sxy\N:\Pi_{12}:'
  WRITE(22,1010) 'Txyy_y:T\sxyy,y\N:\partial T_{122}/\partial x_2:'
  WRITE(22,1010) 'Txyy:T\sxyy\N:T_{122}:'
  WRITE(22,1010) 'Gxy:G\sxy\N:G_{12}:'
  WRITE(22,1010) 'Dxy:D\sxy\N:D_{12}:'

! R13 Budget Equation

! R23 Budget Equation

! TKE budget
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Tke_t:K\s,t\N:\partial K / \partial t:'&
       //'C +P +B -\varepsilon + (\Pi-T_{,2}+D-G)/\bar{\rho}'
  WRITE(22,1010) 'Prd:P:P:'
  WRITE(22,1010) 'Buo:B:B:'
  WRITE(22,1010) 'Eps:\xe\0:\varepsilon:'
  WRITE(22,1010) 'Pi:\xP\0:\Pi:'
  WRITE(22,1010) 'Trp:T\sy,y\N:\partial T / \partial x_2:'
  WRITE(22,1010) 'Phi:\xf\0:\bar{\phi}:'&
       //'\overline{\tau_{lm}\frac{\partial u_l}{\partial x_m}}'

! Scales
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Eta        :\xh\0         :\eta:'&
       //'\left(\nu^3/\varepsilon\right)^{1/4}'
  WRITE(22,1010) 'LambdaUx   :\xl\0\sx\N:\lambda_1:'&
       //'\sqrt{R_{11}\bigg/\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_1}\bigg)^2}}'
  WRITE(22,1010) 'LambdaVy   :\xl\0\sy\N:\lambda_2:'&
       //'\sqrt{R_{22}\bigg/\overline{\bigg(\frac{\partial u^{''}_2}{\partial x_2}\bigg)^2}}'
  WRITE(22,1010) 'LambdaWz   :\xl\0\sz\N:\lambda_3:'&
       //'\sqrt{R_{33}\bigg/\overline{\bigg(\frac{\partial u^{''}_3}{\partial x_3}\bigg)^2}}'
  WRITE(22,1010) 'ReLambdaUx :\1Re\s\xl\0x\N:Re_{\lambda,x}:\sqrt{R_{11}}\lambda_1/\nu'
  WRITE(22,1010) 'ReLambdaVy :\1Re\s\xl\0y\N:Re_{\lambda,y}:\sqrt{R_{22}}\lambda_2/\nu'
  WRITE(22,1010) 'ReLambdaWz :\1Re\s\xl\0z\N:Re_{\lambda,z}:\sqrt{R_{33}}\lambda_3/\nu'
  WRITE(22,1010) 'ReLambdaIso:\1Re\s\xl\N   :Re_{\lambda}  :'&
       //'[(R_{11}+R_{22}+R_{33})/3]\sqrt{\frac{15}{\nu\varepsilon}}'

! Derivatives
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'VarDilatation:\xq\0\S2\N\srms\N:\overline{\theta^{''2}}:'&
       //'\overline{\left(\partial u^{''}_i/\partial x_i\right)^2}'

  WRITE(22,1010) 'SkewUx:\1S\su,x\N:S_{11}:'&
       //'\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_1}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_1}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'SkewUy:\1S\su,y\N:S_{12}:'&
       //'\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_2}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_2}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'SkewUz:\1S\su,z\N:S_{13}:'&
       //'\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_3}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_3}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'SkewVy:\1S\sv,y\N:S_{22}:'
  WRITE(22,1010) 'SkewVx:\1S\sv,x\N:S_{21}:'
  WRITE(22,1010) 'SkewVz:\1S\sv,z\N:S_{23}:'
  WRITE(22,1010) 'SkewWz:\1S\sw,z\N:S_{33}:'
  WRITE(22,1010) 'SkewWx:\1S\sw,x\N:S_{31}:'
  WRITE(22,1010) 'SkewWy:\1S\sw,y\N:S_{32}:'

  WRITE(22,1010) 'FlatUx:\1F\su,x\N:F_{11}:'&
       //'\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_1}'&
       //'\bigg)^4}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_1}\bigg)^2}\bigg)^{2}'
  WRITE(22,1010) 'FlatUy:\1F\su,y\N:F_{12}:'&
       //'\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_2}'&
       //'\bigg)^4}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_2}\bigg)^2}\bigg)^{2}'
  WRITE(22,1010) 'FlatUz:\1F\su,z\N:F_{13}:'&
       //'\overline{\bigg(\frac{\partial u^{''}_1}{\partial x_3}'&
       //'\bigg)^4}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'u^{''}_1}{\partial x_3}\bigg)^2}\bigg)^{2}'
  WRITE(22,1010) 'FlatVy:\1F\sv,y\N:F_{22}:'
  WRITE(22,1010) 'FlatVx:\1F\sv,x\N:F_{21}:'
  WRITE(22,1010) 'FlatVz:\1F\sv,z\N:F_{23}:'
  WRITE(22,1010) 'FlatWz:\1F\sw,z\N:F_{33}:'
  WRITE(22,1010) 'FlatWx:\1F\sw,x\N:F_{31}:'
  WRITE(22,1010) 'FlatWy:\1F\sw,y\N:F_{32}:'

! Compressibility
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'gamma :\xg\0      :\gamma:'
  WRITE(22,1010) 'C2    :\1c\0\S2\N\srms\N:\overline{\gamma p/\rho}'
  WRITE(22,1010) 'Rho_ac:\xr\0\sac\N:\rho^{''}_{ac}:'
  WRITE(22,1010) 'Rho_en:\xr\0\sen\N:\rho^{''}_{en}:'
  WRITE(22,1010) 'T_ac  :\1T\sac\N:T^{''}_{ac}:'
  WRITE(22,1010) 'T_en  :\1T\sen\N:T^{''}_{en}:'
  WRITE(22,1010) 'M_t   :\xr\0\sac\N:\rho^{''}_{ac}:'
  WRITE(22,1010) 'rRP   :rRP:\overline{\rho^{''}p^{''}}:'
  WRITE(22,1010) 'rRT   :rRT:\overline{\rho^{''}T^{''}}:'

  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

! Density
  WRITE(22,1010) 'RhoFluxX:\xr\0\S''\Nu\S'':\overline{\rho^{''}u^{''}}:'
  WRITE(22,1010) 'RhoFluxY:\xr\0\S''\Nv\S'':\overline{\rho^{''}v^{''}}:'
  WRITE(22,1010) 'RhoFluxZ:\xr\0\S''\Nw\S'':\overline{\rho^{''}w^{''}}:'

! Stratification

! Scalar
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'rS:\x\ca\Cz\cq\C\0    :\bar{\zeta}:'
  WRITE(22,1010) 'fS:\x\ca\Cz\cq\C\0    :\tilde{\zeta}:'

  WRITE(22,1010) 'Rsu:\1R\0\s\xz\0x\N:R_{\zeta 1}:'
  WRITE(22,1010) 'Rsv:\1R\0\s\xz\0y\N:R_{\zeta 2}:'
  WRITE(22,1010) 'Rsw:\1R\0\s\xz\0z\N:R_{\zeta 3}:'

  WRITE(22,1010) 'rS2:\1R\0\s\xzz\0\N:R_{\zeta\zeta}:'
  WRITE(22,1010) 'fS2:\1R\0\s\xzz\0\N:R_{\zeta\zeta}:'
  WRITE(22,1010) 'rS3:\xz\0''\S3\N:\overline{\zeta^{''3}}:'
  WRITE(22,1010) 'fS3:\xz\0''\S3\N:\widetilde{\zeta^{''''3}}:'
  WRITE(22,1010) 'rS4:\xz\0''\S4\N:\overline{\zeta^{''4}}:'
  WRITE(22,1010) 'fS4:\xz\0''\S4\N:\widetilde{\zeta^{''''4}}:'
  WRITE(22,1010) 'rS5:\xz\0''\S5\N:\overline{\zeta^{''5}}:'
  WRITE(22,1010) 'fS5:\xz\0''\S5\N:\widetilde{\zeta^{''''5}}:'
  WRITE(22,1010) 'rS6:\xz\0''\S6\N:\overline{\zeta^{''6}}:'
  WRITE(22,1010) 'fS6:\xz\0''\S6\N:\widetilde{\zeta^{''''6}}:'

! Scalar Budget
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Rss_t:R\s\xzz\0,t\N:\partial R_{\zeta\zeta} / \partial t:'&
       //'C_{\zeta\zeta}+P_{\zeta\zeta} -\chi + (-T_{\zeta\zeta 2,2}+D_{\zeta\zeta}+Q_{\zeta\zeta})/\bar{\rho}'
  WRITE(22,1010) 'Css:C\s\xzz\0\N:C_{\zeta\zeta}:'
  WRITE(22,1010) 'Pss:P\s\xzz\0\N:P_{\zeta\zeta}:'
  WRITE(22,1010) 'Ess:\xc\0:\chi:'
  WRITE(22,1010) 'Tssy_y:T\s\xzz\0,y\N:T_{\zeta\zeta 2,2}:'
  WRITE(22,1010) 'Tssy:T\s\xzz\0,y\N:T_{\zeta\zeta 2}:'
  WRITE(22,1010) 'Dss:\xS\szz\N\0:D_{\zeta\zeta}:'
  WRITE(22,1010) 'Qss:Q\s\xzz\N\0:Q_{\zeta\zeta}:'

! Derivatives
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'SkewSx:\1S\s\xz\0,x\N:S_{\zeta 1}:'&
       //'\overline{\bigg(\frac{\partial \zeta^{''}}{\partial x_1}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'\zeta^{''}}{\partial x_1}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'SkewSy:\1S\s\xz\0,y\N:S_{\zeta 2}:'&
       //'\overline{\bigg(\frac{\partial \zeta^{''}}{\partial x_2}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'\zeta^{''}}{\partial x_2}\bigg)^2}\bigg)^{3/2}'
  WRITE(22,1010) 'SkewSz:\1S\s\xz\0,z\N:S_{\zeta 3}:'&
       //'\overline{\bigg(\frac{\partial \zeta^{''}}{\partial x_3}'&
       //'\bigg)^3}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'\zeta^{''}}{\partial x_3}\bigg)^2}\bigg)^{3/2}'

  WRITE(22,1010) 'FlatSx:\1F\s\xz\0,x\N:F_{\zeta 1}:'&
       //'\overline{\bigg(\frac{\partial \zeta^{''}}{\partial x_1}'&
       //'\bigg)^4}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'\zeta^{''}}{\partial x_1}\bigg)^2}\bigg)^{2}'
  WRITE(22,1010) 'FlatSy:\1F\s\xz\0,y\N:F_{\zeta 2}:'&
       //'\overline{\bigg(\frac{\partial \zeta^{''}}{\partial x_2}'&
       //'\bigg)^4}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'\zeta^{''}}{\partial x_2}\bigg)^2}\bigg)^{2}'
  WRITE(22,1010) 'FlatSz:\1F\s\xz\0,z\N:F_{\zeta 3}:'&
       //'\overline{\bigg(\frac{\partial \zeta^{''}}{\partial x_3}'&
       //'\bigg)^4}\bigg/\bigg(\overline{\bigg(\frac{\partial '&
       //'\zeta^{''}}{\partial x_3}\bigg)^2}\bigg)^{2}'

! Growth Rates and Thickness
  WRITE(22,1010) '\medskip\hrule{\hfill}\medskip::\,:'

  WRITE(22,1010) 'Delta_w:\xd\sw\N\0:\delta_\omega:'&
       //'\frac{\Delta u}{\partial\tilde{u}/\partial x_2}'
  WRITE(22,1010) 'Delta_m:\xd\sq\N\0:\delta_\theta:'&
       //'\frac{1}{\rho_o}\int\bar{\rho}(1/4-\tilde{u}^2_1/'&
       //'(\Delta u)^2) d x_2'
  WRITE(22,1010) 'Delta_m_p:d\xd\sq\N\0/dt:\dot{\delta}_\theta:'&
       //'\frac{2}{\rho_o (\Delta u)^3}\int(\tau_{12} - \bar{\rho}'&
       //'R_{12}) \frac{\partial \tilde{u}_1}{\partial x_2}d x_2'
  WRITE(22,1010) 'Delta_s:\xd\sz\N\0:\delta_\zeta:'&
       //'\frac{\Delta \zeta}{\partial\tilde{\zeta}/\partial x_2}'

1010 FORMAT(A)
  CLOSE(22)

  RETURN
END SUBROUTINE AVG_DEFS_TEMPORAL

