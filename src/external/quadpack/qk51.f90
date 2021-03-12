      subroutine qk51(f,p,a,b,result,abserr,resabs,resasc)
      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_200_L   0.2e+03
#define C_1P5_L   1.5e+00
#define C_50_L    0.5e+02


#define C_XGKC1_L    0.9992621049926098e+00
#define C_XGKC2_L    0.9955569697904981e+00
#define C_XGKC3_L    0.9880357945340772e+00
#define C_XGKC4_L    0.9766639214595175e+00
#define C_XGKC5_L    0.9616149864258425e+00
#define C_XGKC6_L    0.9429745712289743e+00
#define C_XGKC7_L    0.9207471152817016e+00
#define C_XGKC8_L    0.8949919978782754e+00
#define C_XGKC9_L    0.8658470652932756e+00
#define C_XGKC10_L   0.8334426287608340e+00
#define C_XGKC11_L   0.7978737979985001e+00
#define C_XGKC12_L   0.7592592630373576e+00
#define C_XGKC13_L   0.7177664068130844e+00
#define C_XGKC14_L   0.6735663684734684e+00
#define C_XGKC15_L   0.6268100990103174e+00
#define C_XGKC16_L   0.5776629302412230e+00
#define C_XGKC17_L   0.5263252843347192e+00
#define C_XGKC18_L   0.4730027314457150e+00
#define C_XGKC19_L   0.4178853821930377e+00
#define C_XGKC20_L   0.3611723058093878e+00
#define C_XGKC21_L   0.3030895389311078e+00
#define C_XGKC22_L   0.2438668837209884e+00
#define C_XGKC23_L   0.1837189394210489e+00
#define C_XGKC24_L   0.1228646926107104e+00
#define C_XGKC25_L   0.6154448300568508e-01   

#define C_WGKC1_L    0.1987383892330316e-02
#define C_WGKC2_L    0.5561932135356714e-02
#define C_WGKC3_L    0.9473973386174152e-02
#define C_WGKC4_L    0.1323622919557167e-01
#define C_WGKC5_L    0.1684781770912830e-01
#define C_WGKC6_L    0.2043537114588284e-01
#define C_WGKC7_L    0.2400994560695322e-01
#define C_WGKC8_L    0.2747531758785174e-01
#define C_WGKC9_L    0.3079230016738749e-01
#define C_WGKC10_L   0.3400213027432934e-01
#define C_WGKC11_L   0.3711627148341554e-01
#define C_WGKC12_L   0.4008382550403238e-01
#define C_WGKC13_L   0.4287284502017005e-01
#define C_WGKC14_L   0.4550291304992179e-01
#define C_WGKC15_L   0.4798253713883671e-01
#define C_WGKC16_L   0.5027767908071567e-01
#define C_WGKC17_L   0.5236288580640748e-01
#define C_WGKC18_L   0.5425112988854549e-01
#define C_WGKC19_L   0.5595081122041232e-01
#define C_WGKC20_L   0.5743711636156783e-01
#define C_WGKC21_L   0.5868968002239421e-01
#define C_WGKC22_L   0.5972034032417406e-01
#define C_WGKC23_L   0.6053945537604586e-01
#define C_WGKC24_L   0.6112850971705305e-01
#define C_WGKC25_L   0.6147118987142532e-01
#define C_WGKC26_L   0.6158081806783294e-01

#define C_WGC1_L     0.1139379850102629e-01
#define C_WGC2_L     0.2635498661503214e-01
#define C_WGC3_L     0.4093915670130631e-01
#define C_WGC4_L     0.5490469597583519e-01
#define C_WGC5_L     0.6803833381235692e-01
#define C_WGC6_L     0.8014070033500102e-01
#define C_WGC7_L     0.9102826198296365e-01
#define C_WGC8_L     0.1005359490670506e+00
#define C_WGC9_L     0.1085196244742637e+00
#define C_WGC10_L    0.1148582591457116e+00
#define C_WGC11_L    0.1194557635357848e+00
#define C_WGC12_L    0.1222424429903100e+00
#define C_WGC13_L    0.1231760537267155e+00


#else
#define C_200_L   0.2d+03
#define C_1P5_L   1.5d+00
#define C_50_L    0.5d+02


#define C_XGKC1_L    0.9992621049926098d+00
#define C_XGKC2_L    0.9955569697904981d+00
#define C_XGKC3_L    0.9880357945340772d+00
#define C_XGKC4_L    0.9766639214595175d+00
#define C_XGKC5_L    0.9616149864258425d+00
#define C_XGKC6_L    0.9429745712289743d+00
#define C_XGKC7_L    0.9207471152817016d+00
#define C_XGKC8_L    0.8949919978782754d+00
#define C_XGKC9_L    0.8658470652932756d+00
#define C_XGKC10_L   0.8334426287608340d+00
#define C_XGKC11_L   0.7978737979985001d+00
#define C_XGKC12_L   0.7592592630373576d+00
#define C_XGKC13_L   0.7177664068130844d+00
#define C_XGKC14_L   0.6735663684734684d+00
#define C_XGKC15_L   0.6268100990103174d+00
#define C_XGKC16_L   0.5776629302412230d+00
#define C_XGKC17_L   0.5263252843347192d+00
#define C_XGKC18_L   0.4730027314457150d+00
#define C_XGKC19_L   0.4178853821930377d+00
#define C_XGKC20_L   0.3611723058093878d+00
#define C_XGKC21_L   0.3030895389311078d+00
#define C_XGKC22_L   0.2438668837209884d+00
#define C_XGKC23_L   0.1837189394210489d+00
#define C_XGKC24_L   0.1228646926107104d+00
#define C_XGKC25_L   0.6154448300568508d-01   

#define C_WGKC1_L    0.1987383892330316d-02
#define C_WGKC2_L    0.5561932135356714d-02
#define C_WGKC3_L    0.9473973386174152d-02
#define C_WGKC4_L    0.1323622919557167d-01
#define C_WGKC5_L    0.1684781770912830d-01
#define C_WGKC6_L    0.2043537114588284d-01
#define C_WGKC7_L    0.2400994560695322d-01
#define C_WGKC8_L    0.2747531758785174d-01
#define C_WGKC9_L    0.3079230016738749d-01
#define C_WGKC10_L   0.3400213027432934d-01
#define C_WGKC11_L   0.3711627148341554d-01
#define C_WGKC12_L   0.4008382550403238d-01
#define C_WGKC13_L   0.4287284502017005d-01
#define C_WGKC14_L   0.4550291304992179d-01
#define C_WGKC15_L   0.4798253713883671d-01
#define C_WGKC16_L   0.5027767908071567d-01
#define C_WGKC17_L   0.5236288580640748d-01
#define C_WGKC18_L   0.5425112988854549d-01
#define C_WGKC19_L   0.5595081122041232d-01
#define C_WGKC20_L   0.5743711636156783d-01
#define C_WGKC21_L   0.5868968002239421d-01
#define C_WGKC22_L   0.5972034032417406d-01
#define C_WGKC23_L   0.6053945537604586d-01
#define C_WGKC24_L   0.6112850971705305d-01
#define C_WGKC25_L   0.6147118987142532d-01
#define C_WGKC26_L   0.6158081806783294d-01

#define C_WGC1_L     0.1139379850102629d-01
#define C_WGC2_L     0.2635498661503214d-01
#define C_WGC3_L     0.4093915670130631d-01
#define C_WGC4_L     0.5490469597583519d-01
#define C_WGC5_L     0.6803833381235692d-01
#define C_WGC6_L     0.8014070033500102d-01
#define C_WGC7_L     0.9102826198296365d-01
#define C_WGC8_L     0.1005359490670506d+00
#define C_WGC9_L     0.1085196244742637d+00
#define C_WGC10_L    0.1148582591457116d+00
#define C_WGC11_L    0.1194557635357848d+00
#define C_WGC12_L    0.1222424429903100d+00
#define C_WGC13_L    0.1231760537267155d+00

#endif

!***begin prologue  qk51
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  51-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!   de doncker,elise,appl. math & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                   estimate
!               j = integral of abs(f) over (a,b)
!***description
!
!   integration rules
!   standard fortran subroutine
!   real version
!
!   parameters
!    on entry
!      f      - real
!               function subroutine defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the calling program.
!
!      a      - real
!               lower limit of integration
!
!      b      - real
!               upper limit of integration
!
!    on return
!      result - real
!               approximation to the integral i
!               result is computed by applying the 51-point
!               kronrod rule (resk) obtained by optimal addition
!               of abscissae to the 25-point gauss rule (resg).
!
!      abserr - real
!               estimate of the modulus of the absolute error,
!               which should not exceed abs(i-result)
!
!      resabs - real
!               approximation to the integral j
!
!      resasc - real
!               approximation to the integral of abs(f-i/(b-a))
!               over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  qk51
!
      TREAL a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,&
        fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,M_MACH,&
        uflow,&
        wg,wgk,xgk,p
      TINTEGER j,jtw,jtwm1
      TINTEGER i1,i4
      external f
!
      dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)
!
!   the abscissae and weights are given for the interval (-1,1).
!   because of symmetry only the positive abscissae and their
!   corresponding weights are given.
!
!   xgk    - abscissae of the 51-point kronrod rule
!            xgk(2), xgk(4), ...  abscissae of the 25-point
!            gauss rule
!            xgk(1), xgk(3), ...  abscissae which are optimally
!            added to the 25-point gauss rule
!
!   wgk    - weights of the 51-point kronrod rule
!
!   wg     - weights of the 25-point gauss rule
!
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),&
        xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14)/&
           C_XGKC1_L,   C_XGKC2_L,&
           C_XGKC3_L,   C_XGKC4_L,&
           C_XGKC5_L,   C_XGKC6_L,&
           C_XGKC7_L,   C_XGKC8_L,&
           C_XGKC9_L,   C_XGKC10_L,&
           C_XGKC11_L,   C_XGKC12_L,&
           C_XGKC13_L,   C_XGKC14_L/
       data xgk(15),xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21),&
        xgk(22),xgk(23),xgk(24),xgk(25),xgk(26)/&
           C_XGKC15_L,   C_XGKC16_L,&
           C_XGKC17_L,   C_XGKC18_L,&
           C_XGKC19_L,   C_XGKC20_L,&
           C_XGKC21_L,   C_XGKC22_L,&
           C_XGKC23_L,   C_XGKC24_L,&
           C_XGKC25_L,   C_0_R               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),&
        wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14)/&
           C_WGKC1_L,   C_WGKC2_L,&
           C_WGKC3_L,   C_WGKC4_L,&
           C_WGKC5_L,   C_WGKC6_L,&
           C_WGKC7_L,   C_WGKC8_L,&
           C_WGKC9_L,   C_WGKC10_L,&
           C_WGKC11_L,   C_WGKC12_L,&
           C_WGKC13_L,   C_WGKC14_L/
       data wgk(15),wgk(16),wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)&
        ,wgk(22),wgk(23),wgk(24),wgk(25),wgk(26)/&
           C_WGKC15_L,   C_WGKC16_L,&
           C_WGKC17_L,   C_WGKC18_L,&
           C_WGKC19_L,   C_WGKC20_L,&
           C_WGKC21_L,   C_WGKC22_L,&
           C_WGKC23_L,   C_WGKC24_L,&
           C_WGKC25_L,   C_WGKC26_L/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),&
        wg(10),wg(11),wg(12),wg(13)/&
           C_WGC1_L,    C_WGC2_L,&
           C_WGC3_L,    C_WGC4_L,&
           C_WGC5_L,    C_WGC6_L,&
           C_WGC7_L,    C_WGC8_L,&
           C_WGC9_L,    C_WGC10_L,&
           C_WGC11_L,    C_WGC12_L,&
           C_WGC13_L/
!
!
!   list of major variables
!   -----------------------
!
!   centr  - mid point of the interval
!   hlgth  - half-length of the interval
!   absc   - abscissa
!   fval*  - function value
!   resg   - result of the 25-point gauss formula
!   resk   - result of the 51-point kronrod formula
!   reskh  - approximation to the mean value of f over (a,b),
!            i.e. to i/(b-a)
!
!   machine dependent constants
!   ---------------------------
!
!   epmach is the largest relative spacing.
!   uflow is the smallest positive magnitude.
!
!***first executable statement  qk51
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
      centr = C_05_R*(a+b)
      hlgth = C_05_R*(b-a)
      dhlgth = abs(hlgth)
!
!   compute the 51-point kronrod approximation to
!   the integral, and estimate the absolute error.
!
      fc = f(centr,p)
      resg = wg(13)*fc
      resk = wgk(26)*fc
      resabs = abs(resk)
      do 10 j=1,12
        jtw = j*2
        absc = hlgth*xgk(jtw)
        fval1 = f(centr-absc,p)
        fval2 = f(centr+absc,p)
        fv1(jtw) = fval1
        fv2(jtw) = fval2
        fsum = fval1+fval2
        resg = resg+wg(j)*fsum
        resk = resk+wgk(jtw)*fsum
        resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
   10 continue
      do 15 j = 1,13
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,p)
        fval2 = f(centr+absc,p)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
   15 continue
      reskh = resk*C_05_R
      resasc = wgk(26)*abs(fc-reskh)
      do 20 j=1,25
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.C_0_R.and.abserr.ne. C_0_R)&
        abserr = resasc*MIN(C_1_R,&
        (C_200_L*abserr/resasc)**C_1P5_L)
      if(resabs.gt.uflow/(C_50_L*epmach)) abserr = MAX&
        ((epmach*C_50_L)*resabs,abserr)
      return
      end
