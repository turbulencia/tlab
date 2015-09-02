      subroutine qage(f,p,a,b,epsabs,epsrel,key,limit,result,abserr,&
         neval,ier,alist,blist,rlist,elist,iord,last)

      IMPLICIT NONE

#include "types.h"

#ifdef SINGLE_PREC
#define C_05EP2_L  0.5e+02
#define C_05EN14_L 0.5e-14
#define C_01EN4_L  0.1e-04
#define C_099_L    0.99e+00
#define C_01EP4_L  0.1e+04
#define C_100_L    100.0e+0
#else
#define C_05EP2_L  0.5d+02
#define C_05EN14_L 0.5d-14
#define C_01EN4_L  0.1d-04
#define C_099_L    0.99d+00
#define C_01EP4_L  0.1d+04
#define C_100_L    100.0d+0
#endif

!***begin prologue  qage
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
!***keywords  automatic integrator, general-purpose,
!     integrand examinator, globally adaptive,
!     gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
!   de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!    definite integral   i = integral of f over (a,b),
!    hopefully satisfying following claim for accuracy
!    abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
!***description
!
!computation of a definite integral
!standard fortran subroutine
!real version
!
!parameters
! on entry
!    f      - real
!             function subprogram defining the integrand
!             function f(x). the actual name for f needs to be
!             declared e x t e r n a l in the driver program.
!
!    a      - real
!             lower limit of integration
!
!    b      - real
!             upper limit of integration
!
!    epsabs - real
!             absolute accuracy requested
!    epsrel - real
!             relative accuracy requested
!             if  epsabs.le.0
!             and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!             the routine will end with ier = 6.
!
!    key    - integer
!             key for choice of local integration rule
!             a gauss-kronrod pair is used with
!                  7 - 15 points if key.lt.2,
!                 10 - 21 points if key = 2,
!                 15 - 31 points if key = 3,
!                 20 - 41 points if key = 4,
!                 25 - 51 points if key = 5,
!                 30 - 61 points if key.gt.5.
!
!    limit  - integer
!             gives an upperbound on the number of subintervals
!             in the partition of (a,b), limit.ge.1.
!
! on return
!    result - real
!             approximation to the integral
!
!    abserr - real
!             estimate of the modulus of the absolute error,
!             which should equal or exceed abs(i-result)
!
!    neval  - integer
!             number of integrand evaluations
!
!    ier    - integer
!             ier = 0 normal and reliable termination of the
!                     routine. it is assumed that the requested
!                     accuracy has been achieved.
!             ier.gt.0 abnormal termination of the routine
!                     the estimates for result and error are
!                     less reliable. it is assumed that the
!                     requested accuracy has not been achieved.
!    error messages
!             ier = 1 maximum number of subdivisions allowed
!                     has been achieved. one can allow more
!                     subdivisions by increasing the value
!                     of limit.
!                     however, if this yields no improvement it
!                     is rather advised to analyze the integrand
!                     in order to determine the integration
!                     difficulties. if the position of a local
!                     difficulty can be determined(e.g.
!                     singularity, discontinuity within the
!                     interval) one will probably gain from
!                     splitting up the interval at this point
!                     and calling the integrator on the
!                     subranges. if possible, an appropriate
!                     special-purpose integrator should be used
!                     which is designed for handling the type of
!                     difficulty involved.
!                 = 2 the occurrence of roundoff error is
!                     detected, which prevents the requested
!                     tolerance from being achieved.
!                 = 3 extremely bad integrand behaviour occurs
!                     at some points of the integration
!                     interval.
!                 = 6 the input is invalid, because
!                     (epsabs.le.0 and
!                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                     result, abserr, neval, last, rlist(1) ,
!                     elist(1) and iord(1) are set to zero.
!                     alist(1) and blist(1) are set to a and b
!                     respectively.
!
!    alist   - real
!              vector of dimension at least limit, the first
!               last  elements of which are the left
!              end points of the subintervals in the partition
!              of the given integration range (a,b)
!
!    blist   - real
!              vector of dimension at least limit, the first
!               last  elements of which are the right
!              end points of the subintervals in the partition
!              of the given integration range (a,b)
!
!    rlist   - real
!              vector of dimension at least limit, the first
!               last  elements of which are the
!              integral approximations on the subintervals
!
!    elist   - real
!              vector of dimension at least limit, the first
!               last  elements of which are the moduli of the
!              absolute error estimates on the subintervals
!
!    iord    - integer
!              vector of dimension at least limit, the first k
!              elements of which are pointers to the
!              error estimates over the subintervals,
!              such that elist(iord(1)), ...,
!              elist(iord(k)) form a decreasing sequence,
!              with k = last if last.le.(limit/2+2), and
!              k = limit+1-last otherwise
!
!    last    - integer
!              number of subintervals actually produced in the
!              subdivision process
!
!***references  (none)
!***routines called  qk15,qk21,qk31,qk41,qk51,qk61,qpsrt,d1mach
!***end prologue  qage
!
      TREAL a,abserr,alist,area,area1,area12,area2,a1,a2,b,blist,&
        b1,b2,defabs,defab1,defab2,M_MACH,elist,epmach,p,&
        epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f,&
        resabs,result,rlist,uflow
      TINTEGER ier,iord,iroff1,iroff2,k,key,keyf,last,&
        limit,maxerr,neval,nrmax
      TINTEGER i1,i4
!
      dimension alist(limit),blist(limit),elist(limit),iord(limit),&
        rlist(limit)
!
      external f
!
!    list of major variables
!    -----------------------
!
!   alist     - list of left end points of all subintervals
!               considered up to now
!   blist     - list of right end points of all subintervals
!               considered up to now
!   rlist(i)  - approximation to the integral over
!              (alist(i),blist(i))
!   elist(i)  - error estimate applying to rlist(i)
!   maxerr    - pointer to the interval with largest
!               error estimate
!   errmax    - elist(maxerr)
!   area      - sum of the integrals over the subintervals
!   errsum    - sum of the errors over the subintervals
!   errbnd    - requested accuracy max(epsabs,epsrel*
!               abs(result))
!   *****1    - variable for the left subinterval
!   *****2    - variable for the right subinterval
!   last      - index for subdivision
!
!
!   machine dependent constants
!   ---------------------------
!
!   epmach  is the largest relative spacing.
!   uflow  is the smallest positive magnitude.
!
!***first executable statement  qage
      i1 = 1
      i4 = 4
      epmach = M_MACH(i4)
      uflow = M_MACH(i1)
!
!   test on validity of parameters
!   ------------------------------
!
      ier = 0
      neval = 0
      last = 0
      result = C_0_R
      abserr = C_0_R
      alist(1) = a
      blist(1) = b
      rlist(1) = C_0_R
      elist(1) = C_0_R
      iord(1) = 0

      if(epsabs.le.C_0_R.and.&
        epsrel.lt.MAX(C_05EP2_L*epmach,C_05EN14_L)) ier = 6
      if(ier.eq.6) go to 999
!
!   first approximation to the integral
!   -----------------------------------
!
      keyf = key
      if(key.le.0) keyf = 1
      if(key.ge.7) keyf = 6
      neval = 0
      if(keyf.eq.1) call qk15(f,p,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.2) call qk21(f,p,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.3) call qk31(f,p,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.4) call qk41(f,p,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.5) call qk51(f,p,a,b,result,abserr,defabs,resabs)
      if(keyf.eq.6) call qk61(f,p,a,b,result,abserr,defabs,resabs)
      last = 1
      rlist(1) = result
      elist(1) = abserr
      iord(1) = 1
!
!   test on accuracy.
!
      errbnd = MAX(epsabs,epsrel*abs(result))
      if(abserr.le.C_05EP2_L*epmach*defabs.and.abserr.gt.&
        errbnd) ier = 2
      if(limit.eq.1) ier = 1
      if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs)&
        .or.abserr.eq.C_0_R) go to 60
!
!   initialization
!   --------------
!
!
      errmax = abserr
      maxerr = 1
      area = result
      errsum = abserr
      nrmax = 1
      iroff1 = 0
      iroff2 = 0
!
!   main do-loop
!   ------------
!
      do 30 last = 2,limit
!
!   bisect the subinterval with the largest error estimate.
!
        a1 = alist(maxerr)
        b1 = C_05_R*(alist(maxerr)+blist(maxerr))
        a2 = b1
        b2 = blist(maxerr)
        if(keyf.eq.1) call qk15(f,p,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.2) call qk21(f,p,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.3) call qk31(f,p,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.4) call qk41(f,p,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.5) call qk51(f,p,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.6) call qk61(f,p,a1,b1,area1,error1,resabs,defab1)
        if(keyf.eq.1) call qk15(f,p,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.2) call qk21(f,p,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.3) call qk31(f,p,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.4) call qk41(f,p,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.5) call qk51(f,p,a2,b2,area2,error2,resabs,defab2)
        if(keyf.eq.6) call qk61(f,p,a2,b2,area2,error2,resabs,defab2)
!
!   improve previous approximations to integral
!   and error and test for accuracy.
!
        neval = neval+1
        area12 = area1+area2
        erro12 = error1+error2
        errsum = errsum+erro12-errmax
        area = area+area12-rlist(maxerr)
        if(defab1.eq.error1.or.defab2.eq.error2) go to 5
        if(abs(rlist(maxerr)-area12).le.C_01EN4_L*abs(area12)&
        .and.erro12.ge.C_099_L*errmax) iroff1 = iroff1+1
        if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
    5   rlist(maxerr) = area1
        rlist(last) = area2
        errbnd = MAX(epsabs,epsrel*abs(area))
        if(errsum.le.errbnd) go to 8
!
!   test for roundoff error and eventually
!   set error flag.
!
        if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
!
!   set error flag in the case that the number of
!   subintervals equals limit.
!
        if(last.eq.limit) ier = 1
!
!   set error flag in the case of bad integrand behaviour
!   at a point of the integration range.
!
        if(MAX(abs(a1),abs(b2)).le.(C_1_R+C_100_L*&
             epmach)*(abs(a2)+C_01EP4_L*uflow)) ier = 3
!
!   append the newly-created intervals to the list.
!
    8   if(error2.gt.error1) go to 10
        alist(last) = a2
        blist(maxerr) = b1
        blist(last) = b2
        elist(maxerr) = error1
        elist(last) = error2
        go to 20
   10   alist(maxerr) = a2
        alist(last) = a1
        blist(last) = b1
        rlist(maxerr) = area2
        rlist(last) = area1
        elist(maxerr) = error2
        elist(last) = error1
!
!   call subroutine qpsrt to maintain the descending ordering
!   in the list of error estimates and select the
!   subinterval with the largest error estimate (to be
!   bisected next).
!
   20   call qpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
! ***jump out of do-loop
        if(ier.ne.0.or.errsum.le.errbnd) go to 40
   30 continue
!
!   compute final result.
!   ---------------------
!
   40 result = C_0_R
      do 50 k=1,last
        result = result+rlist(k)
   50 continue
      abserr = errsum
   60 if(keyf.ne.1) neval = (10*keyf+1)*(2*neval+1)
      if(keyf.eq.1) neval = 30*neval+15
  999 return
      end
