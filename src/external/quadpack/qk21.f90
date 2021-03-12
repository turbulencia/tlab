      subroutine qk21(f,p,a,b,result,abserr,resabs,resasc)

      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_200_L   0.2e+03
#define C_1P5_L   1.5e+00
#define C_50_L    0.5e+02


#define C_XGKC1_L   0.9956571630258081e+00
#define C_XGKC2_L   0.9739065285171717e+00
#define C_XGKC3_L   0.9301574913557082e+00
#define C_XGKC4_L   0.8650633666889845e+00
#define C_XGKC5_L   0.7808177265864169e+00
#define C_XGKC6_L   0.6794095682990244e+00
#define C_XGKC7_L   0.5627571346686047e+00
#define C_XGKC8_L   0.4333953941292472e+00
#define C_XGKC9_L   0.2943928627014602e+00
#define C_XGKC10_L  0.1488743389816312e+00
     
#define C_WGKC1_L   0.1169463886737187e-01 
#define C_WGKC2_L   0.3255816230796473e-01
#define C_WGKC3_L   0.5475589657435200e-01
#define C_WGKC4_L   0.7503967481091995e-01
#define C_WGKC5_L   0.9312545458369761e-01
#define C_WGKC6_L   0.1093871588022976e+00
#define C_WGKC7_L   0.1234919762620659e+00
#define C_WGKC8_L   0.1347092173114733e+00
#define C_WGKC9_L   0.1427759385770601e+00
#define C_WGKC10_L  0.1477391049013385e+00
#define C_WGKC11_L  0.1494455540029169e+00

#define C_WGC1_L    0.6667134430868814e-01 
#define C_WGC2_L    0.1494513491505806e+00
#define C_WGC3_L    0.2190863625159820e+00 
#define C_WGC4_L    0.2692667193099964e+00
#define C_WGC5_L    0.2955242247147529e+00


#else
#define C_200_L   0.2d+03
#define C_1P5_L   1.5d+00
#define C_50_L    0.5d+02

#define C_XGKC1_L   0.9956571630258081d+00
#define C_XGKC2_L   0.9739065285171717d+00
#define C_XGKC3_L   0.9301574913557082d+00
#define C_XGKC4_L   0.8650633666889845d+00
#define C_XGKC5_L   0.7808177265864169d+00
#define C_XGKC6_L   0.6794095682990244d+00
#define C_XGKC7_L   0.5627571346686047d+00
#define C_XGKC8_L   0.4333953941292472d+00
#define C_XGKC9_L   0.2943928627014602d+00
#define C_XGKC10_L  0.1488743389816312d+00
     
#define C_WGKC1_L   0.1169463886737187d-01 
#define C_WGKC2_L   0.3255816230796473d-01
#define C_WGKC3_L   0.5475589657435200d-01
#define C_WGKC4_L   0.7503967481091995d-01
#define C_WGKC5_L   0.9312545458369761d-01
#define C_WGKC6_L   0.1093871588022976d+00
#define C_WGKC7_L   0.1234919762620659d+00
#define C_WGKC8_L   0.1347092173114733d+00
#define C_WGKC9_L   0.1427759385770601d+00
#define C_WGKC10_L  0.1477391049013385d+00
#define C_WGKC11_L  0.1494455540029169d+00

#define C_WGC1_L    0.6667134430868814d-01 
#define C_WGC2_L    0.1494513491505806d+00
#define C_WGC3_L    0.2190863625159820d+00 
#define C_WGC4_L    0.2692667193099964d+00
#define C_WGC5_L    0.2955242247147529d+00
#endif

!***begin prologue  qk21
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  21-point gauss-kronrod rules
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
!               declared e x t e r n a l in the driver program.
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
!               result is computed by applying the 21-point
!               kronrod rule (resk) obtained by optimal addition
!               of abscissae to the 10-point gauss rule (resg).
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
!***end prologue  qk21
!
      TREAL a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,&
        fv1,fv2,hlgth,resabs,resg,resk,reskh,result,M_MACH,&
        uflow,wg,wgk,&
        xgk,resasc,p
      TINTEGER j,jtw,jtwm1
      TINTEGER i1,i4
      external f
!
      dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)
!
!   the abscissae and weights are given for the interval (-1,1).
!   because of symmetry only the positive abscissae and their
!   corresponding weights are given.
!
!   xgk    - abscissae of the 21-point kronrod rule
!            xgk(2), xgk(4), ...  abscissae of the 10-point
!            gauss rule
!            xgk(1), xgk(3), ...  abscissae which are optimally
!            added to the 10-point gauss rule
!
!   wgk    - weights of the 21-point kronrod rule
!
!   wg     - weights of the 10-point gauss rule
!
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),&
           xgk(8),xgk(9),xgk(10),xgk(11)/&
           C_XGKC1_L,     C_XGKC2_L,&
           C_XGKC3_L,     C_XGKC4_L,&
           C_XGKC5_L,     C_XGKC6_L,&
           C_XGKC7_L,     C_XGKC8_L,&
           C_XGKC9_L,     C_XGKC10_L,&
           C_0_R/
!
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),&
        wgk(8),wgk(9),wgk(10),wgk(11)/&
           C_WGKC1_L,     C_WGKC2_L,&
           C_WGKC3_L,     C_WGKC4_L,&
           C_WGKC5_L,     C_WGKC6_L,&
           C_WGKC7_L,     C_WGKC8_L,&
           C_WGKC9_L,     C_WGKC10_L,&
           C_WGKC11_L/
!
      data wg(1),wg(2),wg(3),wg(4),wg(5)/&
           C_WGC1_L,     C_WGC2_L,&
           C_WGC3_L,     C_WGC4_L,&
           C_WGC5_L/
!
!
!   list of major variables
!   -----------------------
!
!   centr  - mid point of the interval
!   hlgth  - half-length of the interval
!   absc   - abscissa
!   fval*  - function value
!   resg   - result of the 10-point gauss formula
!   resk   - result of the 21-point kronrod formula
!   reskh  - approximation to the mean value of f over (a,b),
!            i.e. to i/(b-a)
!
!
!   machine dependent constants
!   ---------------------------
!
!   epmach is the largest relative spacing.
!   uflow is the smallest positive magnitude.
!
!***first executable statement  qk21
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
      centr = C_05_R*(a+b)
      hlgth = C_05_R*(b-a)
      dhlgth = abs(hlgth)
!
!   compute the 21-point kronrod approximation to
!   the integral, and estimate the absolute error.
!
      resg = C_0_R
      fc = f(centr,p)
      resk = wgk(11)*fc
      resabs = abs(resk)
      do 10 j=1,5
        jtw = 2*j
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
      do 15 j = 1,5
        jtwm1 = 2*j-1
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
      resasc = wgk(11)*abs(fc-reskh)
      do 20 j=1,10
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
