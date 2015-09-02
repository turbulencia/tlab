      subroutine qk41(f,p,a,b,result,abserr,resabs,resasc)

      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_1_L     1.0e+00
#define C_200_L   0.2e+03
#define C_1P5_L   1.5e+00
#define C_50_L    0.5e+02


#define C_XGKC1_L   0.9988590315882777e+00
#define C_XGKC2_L   0.9931285991850949e+00
#define C_XGKC3_L   0.9815078774502503e+00
#define C_XGKC4_L   0.9639719272779138e+00
#define C_XGKC5_L   0.9408226338317548e+00
#define C_XGKC6_L   0.9122344282513259e+00
#define C_XGKC7_L   0.8782768112522820e+00
#define C_XGKC8_L   0.8391169718222188e+00
#define C_XGKC9_L   0.7950414288375512e+00
#define C_XGKC10_L  0.7463319064601508e+00
#define C_XGKC11_L  0.6932376563347514e+00
#define C_XGKC12_L  0.6360536807265150e+00
#define C_XGKC13_L  0.5751404468197103e+00
#define C_XGKC14_L  0.5108670019508271e+00
#define C_XGKC15_L  0.4435931752387251e+00
#define C_XGKC16_L  0.3737060887154196e+00
#define C_XGKC17_L  0.3016278681149130e+00
#define C_XGKC18_L  0.2277858511416451e+00
#define C_XGKC19_L  0.1526054652409227e+00
#define C_XGKC20_L  0.7652652113349733e-01

#define C_WGKC1_L   0.3073583718520532e-02
#define C_WGKC2_L   0.8600269855642942e-02
#define C_WGKC3_L   0.1462616925697125e-01
#define C_WGKC4_L   0.2038837346126652e-01
#define C_WGKC5_L   0.2588213360495116e-01
#define C_WGKC6_L   0.3128730677703280e-01
#define C_WGKC7_L   0.3660016975820080e-01
#define C_WGKC8_L   0.4166887332797369e-01
#define C_WGKC9_L   0.4643482186749767e-01
#define C_WGKC10_L  0.5094457392372869e-01
#define C_WGKC11_L  0.5519510534828599e-01
#define C_WGKC12_L  0.5911140088063957e-01
#define C_WGKC13_L  0.6265323755478117e-01
#define C_WGKC14_L  0.6583459713361842e-01
#define C_WGKC15_L  0.6864867292852162e-01
#define C_WGKC16_L  0.7105442355344407e-01
#define C_WGKC17_L  0.7303069033278667e-01
#define C_WGKC18_L  0.7458287540049919e-01
#define C_WGKC19_L  0.7570449768455667e-01
#define C_WGKC20_L  0.7637786767208074e-01
#define C_WGKC21_L  0.7660071191799966e-01

#define C_WGC1_L    0.1761400713915212e-01
#define C_WGC2_L    0.4060142980038694e-01
#define C_WGC3_L    0.6267204833410906e-01
#define C_WGC4_L    0.8327674157670475e-01
#define C_WGC5_L    0.1019301198172404e+00
#define C_WGC6_L    0.1181945319615184e+00
#define C_WGC7_L    0.1316886384491766e+00
#define C_WGC8_L    0.1420961093183821e+00
#define C_WGC9_L    0.1491729864726037e+00
#define C_WGC10_L   0.1527533871307259e+00


#else
#define C_1_L     1.0d+00
#define C_200_L   0.2d+03
#define C_1P5_L   1.5d+00
#define C_50_L    0.5d+02

#define C_XGKC1_L   0.9988590315882777d+00
#define C_XGKC2_L   0.9931285991850949d+00
#define C_XGKC3_L   0.9815078774502503d+00
#define C_XGKC4_L   0.9639719272779138d+00
#define C_XGKC5_L   0.9408226338317548d+00
#define C_XGKC6_L   0.9122344282513259d+00
#define C_XGKC7_L   0.8782768112522820d+00
#define C_XGKC8_L   0.8391169718222188d+00
#define C_XGKC9_L   0.7950414288375512d+00
#define C_XGKC10_L  0.7463319064601508d+00
#define C_XGKC11_L  0.6932376563347514d+00
#define C_XGKC12_L  0.6360536807265150d+00
#define C_XGKC13_L  0.5751404468197103d+00
#define C_XGKC14_L  0.5108670019508271d+00
#define C_XGKC15_L  0.4435931752387251d+00
#define C_XGKC16_L  0.3737060887154196d+00
#define C_XGKC17_L  0.3016278681149130d+00
#define C_XGKC18_L  0.2277858511416451d+00
#define C_XGKC19_L  0.1526054652409227d+00
#define C_XGKC20_L  0.7652652113349733d-01

#define C_WGKC1_L   0.3073583718520532d-02
#define C_WGKC2_L   0.8600269855642942d-02
#define C_WGKC3_L   0.1462616925697125d-01
#define C_WGKC4_L   0.2038837346126652d-01
#define C_WGKC5_L   0.2588213360495116d-01
#define C_WGKC6_L   0.3128730677703280d-01
#define C_WGKC7_L   0.3660016975820080d-01
#define C_WGKC8_L   0.4166887332797369d-01
#define C_WGKC9_L   0.4643482186749767d-01
#define C_WGKC10_L  0.5094457392372869d-01
#define C_WGKC11_L  0.5519510534828599d-01
#define C_WGKC12_L  0.5911140088063957d-01
#define C_WGKC13_L  0.6265323755478117d-01
#define C_WGKC14_L  0.6583459713361842d-01
#define C_WGKC15_L  0.6864867292852162d-01
#define C_WGKC16_L  0.7105442355344407d-01
#define C_WGKC17_L  0.7303069033278667d-01
#define C_WGKC18_L  0.7458287540049919d-01
#define C_WGKC19_L  0.7570449768455667d-01
#define C_WGKC20_L  0.7637786767208074d-01
#define C_WGKC21_L  0.7660071191799966d-01

#define C_WGC1_L    0.1761400713915212d-01
#define C_WGC2_L    0.4060142980038694d-01
#define C_WGC3_L    0.6267204833410906d-01
#define C_WGC4_L    0.8327674157670475d-01
#define C_WGC5_L    0.1019301198172404d+00
#define C_WGC6_L    0.1181945319615184d+00
#define C_WGC7_L    0.1316886384491766d+00
#define C_WGC8_L    0.1420961093183821d+00
#define C_WGC9_L    0.1491729864726037d+00
#define C_WGC10_L   0.1527533871307259d+00


#endif

!***begin prologue  qk41
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  41-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!   de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b), with error
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
!               function subprogram defining the integrand
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
!               result is computed by applying the 41-point
!               gauss-kronrod rule (resk) obtained by optimal
!               addition of abscissae to the 20-point gauss
!               rule (resg).
!
!      abserr - real
!               estimate of the modulus of the absolute error,
!               which should not exceed abs(i-result)
!
!      resabs - real
!               approximation to the integral j
!
!      resasc - real
!               approximation to the integal of abs(f-i/(b-a))
!               over (a,b)
!
!***references  (none)
!***routines called  d1mach
!***end prologue  qk41
!
      TREAL a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,&
        fv1,fv2,hlgth,resabs,&
        resasc,resg,resk,reskh,result,M_MACH,uflow,&
        wg,wgk,xgk,p
      TINTEGER j,jtw,jtwm1
      TINTEGER i1,i4
      external f
!
      dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)
!
!   the abscissae and weights are given for the interval (-1,1).
!   because of symmetry only the positive abscissae and their
!   corresponding weights are given.
!
!   xgk    - abscissae of the 41-point gauss-kronrod rule
!            xgk(2), xgk(4), ...  abscissae of the 20-point
!            gauss rule
!            xgk(1), xgk(3), ...  abscissae which are optimally
!            added to the 20-point gauss rule
!
!   wgk    - weights of the 41-point gauss-kronrod rule
!
!   wg     - weights of the 20-point gauss rule
!
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),&
        xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),&
        xgk(16),xgk(17),xgk(18),xgk(19),xgk(20),xgk(21)/&
           C_XGKC1_L,   C_XGKC2_L,&
           C_XGKC3_L,   C_XGKC4_L,&
           C_XGKC5_L,   C_XGKC6_L,&
           C_XGKC7_L,   C_XGKC8_L,&
           C_XGKC9_L,   C_XGKC10_L,&
           C_XGKC11_L,   C_XGKC12_L,&
           C_XGKC13_L,   C_XGKC14_L,&
           C_XGKC15_L,   C_XGKC16_L,&
           C_XGKC17_L,   C_XGKC18_L,&
           C_XGKC19_L,   C_XGKC20_L,&
           C_0_R               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),&
        wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),&
        wgk(17),wgk(18),wgk(19),wgk(20),wgk(21)/&
           C_WGKC1_L,   C_WGKC2_L,&
           C_WGKC3_L,   C_WGKC4_L,&
           C_WGKC5_L,   C_WGKC6_L,&
           C_WGKC7_L,   C_WGKC8_L,&
           C_WGKC9_L,   C_WGKC10_L,&
           C_WGKC11_L,   C_WGKC12_L,&
           C_WGKC13_L,   C_WGKC14_L,&
           C_WGKC15_L,   C_WGKC16_L,&
           C_WGKC17_L,   C_WGKC18_L,&
           C_WGKC19_L,   C_WGKC20_L,&
           C_WGKC21_L/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8),wg(9),wg(10)/&
           C_WGC1_L,    C_WGC2_L,&
           C_WGC3_L,    C_WGC4_L,&
           C_WGC5_L,    C_WGC6_L,&
           C_WGC7_L,    C_WGC8_L,&
           C_WGC9_L,    C_WGC10_L/
!
!
!   list of major variables
!   -----------------------
!
!   centr  - mid point of the interval
!   hlgth  - half-length of the interval
!   absc   - abscissa
!   fval*  - function value
!   resg   - result of the 20-point gauss formula
!   resk   - result of the 41-point kronrod formula
!   reskh  - approximation to mean value of f over (a,b), i.e.
!            to i/(b-a)
!
!   machine dependent constants
!   ---------------------------
!
!   epmach is the largest relative spacing.
!   uflow is the smallest positive magnitude.
!
!***first executable statement  qk41
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
      centr = C_05_R*(a+b)
      hlgth = C_05_R*(b-a)
      dhlgth = abs(hlgth)
!
!   compute the 41-point gauss-kronrod approximation to
!   the integral, and estimate the absolute error.
!
      resg = C_0_R
      fc = f(centr,p)
      resk = wgk(21)*fc
      resabs = abs(resk)
      do 10 j=1,10
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
      do 15 j = 1,10
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
      resasc = wgk(21)*abs(fc-reskh)
      do 20 j=1,20
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne.C_0_R.and.abserr.ne. C_0_R)&
        abserr = resasc*MIN(C_1_L,&
        (C_200_L*abserr/resasc)**C_1P5_L)
      if(resabs.gt.uflow/(C_50_L*epmach)) abserr = MAX&
        ((epmach*C_50_L)*resabs,abserr)
      return
      end
