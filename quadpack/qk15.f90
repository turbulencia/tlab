      subroutine qk15(f,p,a,b,result,abserr,resabs,resasc)

      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_02EP3_L 0.2e+03
#define C_1P5_L   1.5e+00
#define C_05EP2_L 0.5e+02

#define C_XGK1_L 0.9914553711208126e+00
#define C_XGK2_L 0.9491079123427585e+00
#define C_XGK3_L 0.8648644233597691e+00
#define C_XGK4_L 0.7415311855993944e+00
#define C_XGK5_L 0.5860872354676911e+00
#define C_XGK6_L 0.4058451513773972e+00
#define C_XGK7_L 0.2077849550078985e+00

#define C_WGK1_L 0.2293532201052922e-01
#define C_WGK2_L 0.6309209262997855e-01
#define C_WGK3_L 0.1047900103222502e+00
#define C_WGK4_L 0.1406532597155259e+00
#define C_WGK5_L 0.1690047266392679e+00
#define C_WGK6_L 0.1903505780647854e+00
#define C_WGK7_L 0.2044329400752989e+00
#define C_WGK8_L 0.2094821410847278e+00

#define C_WG1_L 0.1294849661688697e+00
#define C_WG2_L 0.2797053914892767e+00
#define C_WG3_L 0.3818300505051189e+00
#define C_WG4_L 0.4179591836734694e+00
#else
#define C_02EP3_L 0.2d+03
#define C_1P5_L   1.5d+00
#define C_05EP2_L 0.5d+02

#define C_XGK1_L 0.9914553711208126d+00
#define C_XGK2_L 0.9491079123427585d+00
#define C_XGK3_L 0.8648644233597691d+00
#define C_XGK4_L 0.7415311855993944d+00
#define C_XGK5_L 0.5860872354676911d+00
#define C_XGK6_L 0.4058451513773972d+00
#define C_XGK7_L 0.2077849550078985d+00

#define C_WGK1_L 0.2293532201052922d-01
#define C_WGK2_L 0.6309209262997855d-01
#define C_WGK3_L 0.1047900103222502d+00
#define C_WGK4_L 0.1406532597155259d+00
#define C_WGK5_L 0.1690047266392679d+00
#define C_WGK6_L 0.1903505780647854d+00
#define C_WGK7_L 0.2044329400752989d+00
#define C_WGK8_L 0.2094821410847278d+00

#define C_WG1_L 0.1294849661688697d+00
#define C_WG2_L 0.2797053914892767d+00
#define C_WG3_L 0.3818300505051189d+00
#define C_WG4_L 0.4179591836734694d+00
#endif

!***begin prologue  qk15
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a2
!***keywords  15-point gauss-kronrod rules
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!   de doncker,elise,appl. math. & progr. div - k.u.leuven
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
!               result is computed by applying the 15-point
!               kronrod rule (resk) obtained by optimal addition
!               of abscissae to the7-point gauss rule(resg).
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
!***routines called  M_MACH
!***end prologue  qk15
!
      TREAL a,absc,abserr,b,centr,dhlgth,epmach,f,fc,fsum,fval1,fval2,&
        fv1,fv2,hlgth,resabs,resasc,resg,resk,reskh,result,M_MACH,&
        uflow,&
        wg,wgk,xgk,p
      TINTEGER j,jtw,jtwm1
      TINTEGER i1,i4
      external f
!
      dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)
!
!   the abscissae and weights are given for the interval (-1,1).
!   because of symmetry only the positive abscissae and their
!   corresponding weights are given.
!
!   xgk    - abscissae of the 15-point kronrod rule
!            xgk(2), xgk(4), ...  abscissae of the 7-point
!            gauss rule
!            xgk(1), xgk(3), ...  abscissae which are optimally
!            added to the 7-point gauss rule
!
!   wgk    - weights of the 15-point kronrod rule
!
!   wg     - weights of the 7-point gauss rule
!
      data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/&
           C_XGK1_L,   C_XGK2_L,&
           C_XGK3_L,   C_XGK4_L,&
           C_XGK5_L,   C_XGK6_L,&
           C_XGK7_L,   C_0_R  /
      data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/&
           C_WGK1_L,   C_WGK2_L,&
           C_WGK3_L,   C_WGK4_L,&
           C_WGK5_L,   C_WGK6_L,&
           C_WGK7_L,   C_WGK8_L/
      data wg(1),wg(2),wg(3),wg(4)/&
           C_WG1_L,   C_WG2_L,&
           C_WG3_L,   C_WG4_L/
!
!
!   list of major variables
!   -----------------------
!
!   centr  - mid point of the interval
!   hlgth  - half-length of the interval
!   absc   - abscissa
!   fval*  - function value
!   resg   - result of the 7-point gauss formula
!   resk   - result of the 15-point kronrod formula
!   reskh  - approximation to the mean value of f over (a,b),
!            i.e. to i/(b-a)
!
!   machine dependent constants
!   ---------------------------
!
!   epmach is the largest relative spacing.
!   uflow is the smallest positive magnitude.
!
!***first executable statement  qk15
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
      centr = C_05_R*(a+b)
      hlgth = C_05_R*(b-a)
      dhlgth = abs(hlgth)
!
!   compute the 15-point kronrod approximation to
!   the integral, and estimate the absolute error.
!
      fc = f(centr,p)
      resg = fc*wg(4)
      resk = fc*wgk(8)
      resabs = abs(resk)
      do 10 j=1,3
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
      do 15 j = 1,4
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
      resasc = wgk(8)*abs(fc-reskh)
      do 20 j=1,7
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
   20 continue
      result = resk*hlgth
      resabs = resabs*dhlgth
      resasc = resasc*dhlgth
      abserr = abs((resk-resg)*hlgth)
      if(resasc.ne. C_0_R .and.abserr.ne. C_0_R )&
        abserr = resasc*MIN(C_1_R,&
        (C_02EP3_L*abserr/resasc)**C_1P5_L)
      if(resabs.gt.uflow/(C_05EP2_L*epmach)) abserr = MAX&
        ((epmach*C_05EP2_L)*resabs,abserr)
      return
      end

