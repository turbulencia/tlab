      subroutine qk31(f,p,a,b,result,abserr,resabs,resasc)
      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_200_L   0.2e+03
#define C_1P5_L   1.5e+00
#define C_50_L    0.5e+02

#define C_XGKC1_L    0.9980022986933971e+00
#define C_XGKC2_L    0.9879925180204854e+00
#define C_XGKC3_L    0.9677390756791391e+00
#define C_XGKC4_L    0.9372733924007059e+00
#define C_XGKC5_L    0.8972645323440819e+00
#define C_XGKC6_L    0.8482065834104272e+00
#define C_XGKC7_L    0.7904185014424659e+00
#define C_XGKC8_L    0.7244177313601700e+00
#define C_XGKC9_L    0.6509967412974170e+00
#define C_XGKC10_L   0.5709721726085388e+00
#define C_XGKC11_L   0.4850818636402397e+00
#define C_XGKC12_L   0.3941513470775634e+00
#define C_XGKC13_L   0.2991800071531688e+00
#define C_XGKC14_L   0.2011940939974345e+00
#define C_XGKC15_L   0.1011420669187175e+00
                                           
#define C_WGKC1_L    0.5377479872923349e-02
#define C_WGKC2_L    0.1500794732931612e-01
#define C_WGKC3_L    0.2546084732671532e-01
#define C_WGKC4_L    0.3534636079137585e-01
#define C_WGKC5_L    0.4458975132476488e-01
#define C_WGKC6_L    0.5348152469092809e-01
#define C_WGKC7_L    0.6200956780067064e-01
#define C_WGKC8_L    0.6985412131872826e-01
#define C_WGKC9_L    0.7684968075772038e-01
#define C_WGKC10_L   0.8308050282313302e-01
#define C_WGKC11_L   0.8856444305621177e-01
#define C_WGKC12_L   0.9312659817082532e-01
#define C_WGKC13_L   0.9664272698362368e-01
#define C_WGKC14_L   0.9917359872179196e-01
#define C_WGKC15_L   0.1007698455238756e+00
#define C_WGKC16_L   0.1013300070147915e+00
                                           
#define C_WGC1_L     0.3075324199611727e-01
#define C_WGC2_L     0.7036604748810812e-01
#define C_WGC3_L     0.1071592204671719e+00
#define C_WGC4_L     0.1395706779261543e+00
#define C_WGC5_L     0.1662692058169939e+00
#define C_WGC6_L     0.1861610000155622e+00
#define C_WGC7_L     0.1984314853271116e+00
#define C_WGC8_L     0.2025782419255613e+00


#else
#define C_200_L   0.2d+03
#define C_1P5_L   1.5d+00
#define C_50_L    0.5d+02

#define C_XGKC1_L    0.9980022986933971d+00
#define C_XGKC2_L    0.9879925180204854d+00
#define C_XGKC3_L    0.9677390756791391d+00
#define C_XGKC4_L    0.9372733924007059d+00
#define C_XGKC5_L    0.8972645323440819d+00
#define C_XGKC6_L    0.8482065834104272d+00
#define C_XGKC7_L    0.7904185014424659d+00
#define C_XGKC8_L    0.7244177313601700d+00
#define C_XGKC9_L    0.6509967412974170d+00
#define C_XGKC10_L   0.5709721726085388d+00
#define C_XGKC11_L   0.4850818636402397d+00
#define C_XGKC12_L   0.3941513470775634d+00
#define C_XGKC13_L   0.2991800071531688d+00
#define C_XGKC14_L   0.2011940939974345d+00
#define C_XGKC15_L   0.1011420669187175d+00
                                           
#define C_WGKC1_L    0.5377479872923349d-02
#define C_WGKC2_L    0.1500794732931612d-01
#define C_WGKC3_L    0.2546084732671532d-01
#define C_WGKC4_L    0.3534636079137585d-01
#define C_WGKC5_L    0.4458975132476488d-01
#define C_WGKC6_L    0.5348152469092809d-01
#define C_WGKC7_L    0.6200956780067064d-01
#define C_WGKC8_L    0.6985412131872826d-01
#define C_WGKC9_L    0.7684968075772038d-01
#define C_WGKC10_L   0.8308050282313302d-01
#define C_WGKC11_L   0.8856444305621177d-01
#define C_WGKC12_L   0.9312659817082532d-01
#define C_WGKC13_L   0.9664272698362368d-01
#define C_WGKC14_L   0.9917359872179196d-01
#define C_WGKC15_L   0.1007698455238756d+00
#define C_WGKC16_L   0.1013300070147915d+00
                                           
#define C_WGC1_L     0.3075324199611727d-01
#define C_WGC2_L     0.7036604748810812d-01
#define C_WGC3_L     0.1071592204671719d+00
#define C_WGC4_L     0.1395706779261543d+00
#define C_WGC5_L     0.1662692058169939d+00
#define C_WGC6_L     0.1861610000155622d+00
#define C_WGC7_L     0.1984314853271116d+00
#define C_WGC8_L     0.2025782419255613d+00
#endif

!***begin prologue  qk31
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  31-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!   de doncker,elise,appl. math. & progr. div. - k.u.leuven
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
!               result is computed by applying the 31-point
!               gauss-kronrod rule (resk), obtained by optimal
!               addition of abscissae to the 15-point gauss
!               rule (resg).
!
!      abserr - real
!               estimate of the modulus of the modulus,
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
!***end prologue  qk31
      TREAL a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,&
        fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,M_MACH,&
        uflow,&
        wg,wgk,xgk,p
      TINTEGER j,jtw,jtwm1
      TINTEGER i1,i4
      external f
!
      dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)
!
!   the abscissae and weights are given for the interval (-1,1).
!   because of symmetry only the positive abscissae and their
!   corresponding weights are given.
!
!   xgk    - abscissae of the 31-point kronrod rule
!            xgk(2), xgk(4), ...  abscissae of the 15-point
!            gauss rule
!            xgk(1), xgk(3), ...  abscissae which are optimally
!            added to the 15-point gauss rule
!
!   wgk    - weights of the 31-point kronrod rule
!
!   wg     - weights of the 15-point gauss rule
!
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),&
        xgk(9),xgk(10),xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),&
        xgk(16)/&
           C_XGKC1_L,   C_XGKC2_L,&
           C_XGKC3_L,   C_XGKC4_L,&
           C_XGKC5_L,   C_XGKC6_L,&
           C_XGKC7_L,   C_XGKC8_L,&
           C_XGKC9_L,   C_XGKC10_L,&
           C_XGKC11_L,   C_XGKC12_L,&
           C_XGKC13_L,   C_XGKC14_L,&
           C_XGKC15_L,   C_0_R               /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),&
        wgk(9),wgk(10),wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),&
        wgk(16)/&
           C_WGKC1_L,   C_WGKC2_L,&
           C_WGKC3_L,   C_WGKC4_L,&
           C_WGKC5_L,   C_WGKC6_L,&
           C_WGKC7_L,   C_WGKC8_L,&
           C_WGKC9_L,   C_WGKC10_L,&
           C_WGKC11_L,   C_WGKC12_L,&
           C_WGKC13_L,   C_WGKC14_L,&
           C_WGKC15_L,   C_WGKC16_L/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/&
           C_WGC1_L,   C_WGC2_L,&
           C_WGC3_L,   C_WGC4_L,&
           C_WGC5_L,   C_WGC6_L,&
           C_WGC7_L,   C_WGC8_L/
!
!
!   list of major variables
!   -----------------------
!   centr  - mid point of the interval
!   hlgth  - half-length of the interval
!   absc   - abscissa
!   fval*  - function value
!   resg   - result of the 15-point gauss formula
!   resk   - result of the 31-point kronrod formula
!   reskh  - approximation to the mean value of f over (a,b),
!            i.e. to i/(b-a)
!
!   machine dependent constants
!   ---------------------------
!   epmach is the largest relative spacing.
!   uflow is the smallest positive magnitude.
!
!***first executable statement  qk31
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
      centr = C_05_R*(a+b)
      hlgth = C_05_R*(b-a)
      dhlgth = abs(hlgth)
!
!   compute the 31-point kronrod approximation to
!   the integral, and estimate the absolute error.
!
      fc = f(centr,p)
      resg = wg(8)*fc
      resk = wgk(16)*fc
      resabs = abs(resk)
      do 10 j=1,7
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
      do 15 j = 1,8
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
      resasc = wgk(16)*abs(fc-reskh)
      do 20 j=1,15
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
