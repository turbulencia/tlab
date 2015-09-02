      subroutine qk61(f,p,a,b,result,abserr,resabs,resasc)
      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_200_L   0.2e+03
#define C_1P5_L   1.5e+00
#define C_50_L    0.5e+02


#define C_XGKC1_L     0.9994844100504906e+00
#define C_XGKC2_L     0.9968934840746495e+00
#define C_XGKC3_L     0.9916309968704046e+00
#define C_XGKC4_L     0.9836681232797472e+00
#define C_XGKC5_L     0.9731163225011263e+00
#define C_XGKC6_L     0.9600218649683075e+00
#define C_XGKC7_L     0.9443744447485600e+00
#define C_XGKC8_L     0.9262000474292743e+00
#define C_XGKC9_L     0.9055733076999078e+00
#define C_XGKC10_L    0.8825605357920527e+00
#define C_XGKC11_L    0.8572052335460611e+00
#define C_XGKC12_L    0.8295657623827684e+00
#define C_XGKC13_L    0.7997278358218391e+00
#define C_XGKC14_L    0.7677774321048262e+00
#define C_XGKC15_L    0.7337900624532268e+00
#define C_XGKC16_L    0.6978504947933158e+00
#define C_XGKC17_L    0.6600610641266270e+00
#define C_XGKC18_L    0.6205261829892429e+00
#define C_XGKC19_L    0.5793452358263617e+00
#define C_XGKC20_L    0.5366241481420199e+00
#define C_XGKC21_L    0.4924804678617786e+00
#define C_XGKC22_L    0.4470337695380892e+00
#define C_XGKC23_L    0.4004012548303944e+00
#define C_XGKC24_L    0.3527047255308781e+00
#define C_XGKC25_L    0.3040732022736251e+00
#define C_XGKC26_L    0.2546369261678898e+00
#define C_XGKC27_L    0.2045251166823099e+00
#define C_XGKC28_L    0.1538699136085835e+00
#define C_XGKC29_L    0.1028069379667370e+00
#define C_XGKC30_L    0.5147184255531770e-01
     
#define C_WGKC1_L     0.1389013698677008e-02
#define C_WGKC2_L     0.3890461127099884e-02
#define C_WGKC3_L     0.6630703915931292e-02
#define C_WGKC4_L     0.9273279659517763e-02
#define C_WGKC5_L     0.1182301525349634e-01
#define C_WGKC6_L     0.1436972950704580e-01
#define C_WGKC7_L     0.1692088918905327e-01
#define C_WGKC8_L     0.1941414119394238e-01
#define C_WGKC9_L     0.2182803582160919e-01
#define C_WGKC10_L    0.2419116207808060e-01
#define C_WGKC11_L    0.2650995488233310e-01
#define C_WGKC12_L    0.2875404876504129e-01
#define C_WGKC13_L    0.3090725756238776e-01
#define C_WGKC14_L    0.3298144705748373e-01
#define C_WGKC15_L    0.3497933802806002e-01
#define C_WGKC16_L    0.3688236465182123e-01
#define C_WGKC17_L    0.3867894562472759e-01 
#define C_WGKC18_L    0.4037453895153596e-01
#define C_WGKC19_L    0.4196981021516425e-01
#define C_WGKC20_L    0.4345253970135607e-01
#define C_WGKC21_L    0.4481480013316266e-01
#define C_WGKC22_L    0.4605923827100699e-01
#define C_WGKC23_L    0.4718554656929915e-01
#define C_WGKC24_L    0.4818586175708713e-01
#define C_WGKC25_L    0.4905543455502978e-01
#define C_WGKC26_L    0.4979568342707421e-01
#define C_WGKC27_L    0.5040592140278235e-01
#define C_WGKC28_L    0.5088179589874961e-01
#define C_WGKC29_L    0.5122154784925877e-01
#define C_WGKC30_L    0.5142612853745903e-01
#define C_WGKC31_L    0.5149472942945157e-01


#define C_WGC1_L     0.7968192496166606e-02
#define C_WGC2_L     0.1846646831109096e-01
#define C_WGC3_L     0.2878470788332337e-01
#define C_WGC4_L     0.3879919256962705e-01
#define C_WGC5_L     0.4840267283059405e-01
#define C_WGC6_L     0.5749315621761907e-01
#define C_WGC7_L     0.6597422988218050e-01
#define C_WGC8_L     0.7375597473770521e-01
#define C_WGC9_L     0.8075589522942022e-01
#define C_WGC10_L    0.8689978720108298e-01
#define C_WGC11_L    0.9212252223778613e-01
#define C_WGC12_L    0.9636873717464426e-01
#define C_WGC13_L    0.9959342058679527e-01
#define C_WGC14_L    0.1017623897484055e+00
#define C_WGC15_L    0.1028526528935588e+00

#else
#define C_200_L   0.2d+03
#define C_1P5_L   1.5d+00
#define C_50_L    0.5d+02

#define C_XGKC1_L     0.9994844100504906d+00
#define C_XGKC2_L     0.9968934840746495d+00
#define C_XGKC3_L     0.9916309968704046d+00
#define C_XGKC4_L     0.9836681232797472d+00
#define C_XGKC5_L     0.9731163225011263d+00
#define C_XGKC6_L     0.9600218649683075d+00
#define C_XGKC7_L     0.9443744447485600d+00
#define C_XGKC8_L     0.9262000474292743d+00
#define C_XGKC9_L     0.9055733076999078d+00
#define C_XGKC10_L    0.8825605357920527d+00
#define C_XGKC11_L    0.8572052335460611d+00
#define C_XGKC12_L    0.8295657623827684d+00
#define C_XGKC13_L    0.7997278358218391d+00
#define C_XGKC14_L    0.7677774321048262d+00
#define C_XGKC15_L    0.7337900624532268d+00
#define C_XGKC16_L    0.6978504947933158d+00
#define C_XGKC17_L    0.6600610641266270d+00
#define C_XGKC18_L    0.6205261829892429d+00
#define C_XGKC19_L    0.5793452358263617d+00
#define C_XGKC20_L    0.5366241481420199d+00
#define C_XGKC21_L    0.4924804678617786d+00
#define C_XGKC22_L    0.4470337695380892d+00
#define C_XGKC23_L    0.4004012548303944d+00
#define C_XGKC24_L    0.3527047255308781d+00
#define C_XGKC25_L    0.3040732022736251d+00
#define C_XGKC26_L    0.2546369261678898d+00
#define C_XGKC27_L    0.2045251166823099d+00
#define C_XGKC28_L    0.1538699136085835d+00
#define C_XGKC29_L    0.1028069379667370d+00
#define C_XGKC30_L    0.5147184255531770d-01
     
#define C_WGKC1_L     0.1389013698677008d-02
#define C_WGKC2_L     0.3890461127099884d-02
#define C_WGKC3_L     0.6630703915931292d-02
#define C_WGKC4_L     0.9273279659517763d-02
#define C_WGKC5_L     0.1182301525349634d-01
#define C_WGKC6_L     0.1436972950704580d-01
#define C_WGKC7_L     0.1692088918905327d-01
#define C_WGKC8_L     0.1941414119394238d-01
#define C_WGKC9_L     0.2182803582160919d-01
#define C_WGKC10_L    0.2419116207808060d-01
#define C_WGKC11_L    0.2650995488233310d-01
#define C_WGKC12_L    0.2875404876504129d-01
#define C_WGKC13_L    0.3090725756238776d-01
#define C_WGKC14_L    0.3298144705748373d-01
#define C_WGKC15_L    0.3497933802806002d-01
#define C_WGKC16_L    0.3688236465182123d-01
#define C_WGKC17_L    0.3867894562472759d-01 
#define C_WGKC18_L    0.4037453895153596d-01
#define C_WGKC19_L    0.4196981021516425d-01
#define C_WGKC20_L    0.4345253970135607d-01
#define C_WGKC21_L    0.4481480013316266d-01
#define C_WGKC22_L    0.4605923827100699d-01
#define C_WGKC23_L    0.4718554656929915d-01
#define C_WGKC24_L    0.4818586175708713d-01
#define C_WGKC25_L    0.4905543455502978d-01
#define C_WGKC26_L    0.4979568342707421d-01
#define C_WGKC27_L    0.5040592140278235d-01
#define C_WGKC28_L    0.5088179589874961d-01
#define C_WGKC29_L    0.5122154784925877d-01
#define C_WGKC30_L    0.5142612853745903d-01
#define C_WGKC31_L    0.5149472942945157d-01


#define C_WGC1_L     0.7968192496166606d-02
#define C_WGC2_L     0.1846646831109096d-01
#define C_WGC3_L     0.2878470788332337d-01
#define C_WGC4_L     0.3879919256962705d-01
#define C_WGC5_L     0.4840267283059405d-01
#define C_WGC6_L     0.5749315621761907d-01
#define C_WGC7_L     0.6597422988218050d-01
#define C_WGC8_L     0.7375597473770521d-01
#define C_WGC9_L     0.8075589522942022d-01
#define C_WGC10_L    0.8689978720108298d-01
#define C_WGC11_L    0.9212252223778613d-01
#define C_WGC12_L    0.9636873717464426d-01
#define C_WGC13_L    0.9959342058679527d-01
#define C_WGC14_L    0.1017623897484055d+00
#define C_WGC15_L    0.1028526528935588d+00
#endif


!***begin prologue  qk61
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  61-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!   de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  to compute i = integral of f over (a,b) with error
!                   estimate
!               j = integral of dabs(f) over (a,b)
!***description
!
!integration rule
!standard fortran subroutine
!real version
!
!
!parameters
! on entry
!   f      - real
!            function subprogram defining the integrand
!            function f(x). the actual name for f needs to be
!            declared e x t e r n a l in the calling program.
!
!   a      - real
!            lower limit of integration
!
!   b      - real
!            upper limit of integration
!
! on return
!   result - real
!            approximation to the integral i
!            result is computed by applying the 61-point
!            kronrod rule (resk) obtained by optimal addition of
!            abscissae to the 30-point gauss rule (resg).
!
!   abserr - real
!            estimate of the modulus of the absolute error,
!            which should equal or exceed dabs(i-result)
!
!   resabs - real
!            approximation to the integral j
!
!   resasc - real
!            approximation to the integral of dabs(f-i/(b-a))
!
!
!***references  (none)
!***routines called  d1mach
!***end prologue  qk61
!
      TREAL a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,&
        fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,M_MACH,&
        uflow,&
        wg,wgk,xgk,p
      TINTEGER j,jtw,jtwm1
      TINTEGER i1,i4
      external f
!
      dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)
!
!   the abscissae and weights are given for the
!   interval (-1,1). because of symmetry only the positive
!   abscissae and their corresponding weights are given.
!
!   xgk   - abscissae of the 61-point kronrod rule
!           xgk(2), xgk(4)  ... abscissae of the 30-point
!           gauss rule
!           xgk(1), xgk(3)  ... optimally added abscissae
!           to the 30-point gauss rule
!
!   wgk   - weights of the 61-point kronrod rule
!
!   wg    - weigths of the 30-point gauss rule
!
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8),&
         xgk(9),xgk(10)/&
           C_XGKC1_L,     C_XGKC2_L,&
           C_XGKC3_L,     C_XGKC4_L,&
           C_XGKC5_L,     C_XGKC6_L,&
           C_XGKC7_L,     C_XGKC8_L,&
           C_XGKC9_L,     C_XGKC10_L/
      data xgk(11),xgk(12),xgk(13),xgk(14),xgk(15),xgk(16),&
        xgk(17),xgk(18),xgk(19),xgk(20)/&
           C_XGKC11_L,     C_XGKC12_L,&
           C_XGKC13_L,     C_XGKC14_L,&
           C_XGKC15_L,     C_XGKC16_L,&
           C_XGKC17_L,     C_XGKC18_L,&
           C_XGKC19_L,     C_XGKC20_L/
      data xgk(21),xgk(22),xgk(23),xgk(24),&
        xgk(25),xgk(26),xgk(27),xgk(28),xgk(29),xgk(30),xgk(31)/&
           C_XGKC21_L,     C_XGKC22_L,&
           C_XGKC23_L,     C_XGKC24_L,&
           C_XGKC25_L,     C_XGKC26_L,&
           C_XGKC27_L,     C_XGKC28_L,&
           C_XGKC29_L,     C_XGKC30_L,&
           C_0_R                   /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8),&
        wgk(9),wgk(10)/&
           C_WGKC1_L,     C_WGKC2_L,&
           C_WGKC3_L,     C_WGKC4_L,&
           C_WGKC5_L,     C_WGKC6_L,&
           C_WGKC7_L,     C_WGKC8_L,&
           C_WGKC9_L,     C_WGKC10_L/
      data wgk(11),wgk(12),wgk(13),wgk(14),wgk(15),wgk(16),&
        wgk(17),wgk(18),wgk(19),wgk(20)/&
           C_WGKC11_L,     C_WGKC12_L,&
           C_WGKC13_L,     C_WGKC14_L,&
           C_WGKC15_L,     C_WGKC16_L,&
           C_WGKC17_L,     C_WGKC18_L,&
           C_WGKC19_L,     C_WGKC20_L/
      data wgk(21),wgk(22),wgk(23),wgk(24),&
        wgk(25),wgk(26),wgk(27),wgk(28),wgk(29),wgk(30),wgk(31)/&
           C_WGKC21_L,     C_WGKC22_L,&
           C_WGKC23_L,     C_WGKC24_L,&
           C_WGKC25_L,     C_WGKC26_L,&
           C_WGKC27_L,     C_WGKC28_L,&
           C_WGKC29_L,     C_WGKC30_L,&
           C_WGKC31_L/
      data wg(1),wg(2),wg(3),wg(4),wg(5),wg(6),wg(7),wg(8)/&
           C_WGC1_L,     C_WGC2_L,&
           C_WGC3_L,     C_WGC4_L,&
           C_WGC5_L,     C_WGC6_L,&
           C_WGC7_L,     C_WGC8_L/
      data wg(9),wg(10),wg(11),wg(12),wg(13),wg(14),wg(15)/&
           C_WGC9_L,     C_WGC10_L,&
           C_WGC11_L,     C_WGC12_L,&
           C_WGC13_L,     C_WGC14_L,&
           C_WGC15_L/
!
!   list of major variables
!   -----------------------
!
!   centr  - mid point of the interval
!   hlgth  - half-length of the interval
!   absc   - abscissa
!   fval*  - function value
!   resg   - result of the 30-point gauss rule
!   resk   - result of the 61-point kronrod rule
!   reskh  - approximation to the mean value of f
!            over (a,b), i.e. to i/(b-a)
!
!   machine dependent constants
!   ---------------------------
!
!   epmach is the largest relative spacing.
!   uflow is the smallest positive magnitude.
!
!***first executable statement  qk61
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
      centr = C_05_R*(b+a)
      hlgth = C_05_R*(b-a)
      dhlgth = abs(hlgth)
!
!   compute the 61-point kronrod approximation to the
!   integral, and estimate the absolute error.
!
      resg = C_0_R
      fc = f(centr,p)
      resk = wgk(31)*fc
      resabs = abs(resk)
      do 10 j=1,15
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
      do 15 j=1,15
        jtwm1 = j*2-1
        absc = hlgth*xgk(jtwm1)
        fval1 = f(centr-absc,p)
        fval2 = f(centr+absc,p)
        fv1(jtwm1) = fval1
        fv2(jtwm1) = fval2
        fsum = fval1+fval2
        resk = resk+wgk(jtwm1)*fsum
        resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
  15    continue
      reskh = resk*C_05_R
      resasc = wgk(31)*abs(fc-reskh)
      do 20 j=1,30
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
