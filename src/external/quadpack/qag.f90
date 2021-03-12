      subroutine qag(f,p,a,b,epsabs,epsrel,key,result,abserr,neval,ier,&
          limit,lenw,last,iwork,work)

      IMPLICIT NONE

#include "types.h"

!***begin prologue  qag
!***date written   800101   (yymmdd)
!***revision date  830518   (yymmdd)
!***category no.  h2a1a1
!***keywords  automatic integrator, general-purpose,
!     integrand examinator, globally adaptive,
!     gauss-kronrod
!***author  piessens,robert,appl. math. & progr. div - k.u.leuven
!   de doncker,elise,appl. math. & progr. div. - k.u.leuven
!***purpose  the routine calculates an approximation result to a given
!    definite integral i = integral of f over (a,b),
!    hopefully satisfying following claim for accuracy
!    abs(i-result)le.max(epsabs,epsrel*abs(i)).
!***description
!
!computation of a definite integral
!standard fortran subroutine
!real version
!
!    f      - real
!             function subprogam defining the integrand
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
!               7 - 15 points if key.lt.2,
!              10 - 21 points if key = 2,
!              15 - 31 points if key = 3,
!              20 - 41 points if key = 4,
!              25 - 51 points if key = 5,
!              30 - 61 points if key.gt.5.
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
!              error messages
!             ier = 1 maximum number of subdivisions allowed
!                     has been achieved. one can allow more
!                     subdivisions by increasing the value of
!                     limit (and taking the according dimension
!                     adjustments into account). however, if
!                     this yield no improvement it is advised
!                     to analyze the integrand in order to
!                     determine the integration difficulaties.
!                     if the position of a local difficulty can
!                     be determined (i.e.singularity,
!                     discontinuity within the interval) one
!                     will probably gain from splitting up the
!                     interval at this point and calling the
!                     integrator on the subranges. if possible,
!                     an appropriate special-purpose integrator
!                     should be used which is designed for
!                     handling the type of difficulty involved.
!                 = 2 the occurrence of roundoff error is
!                     detected, which prevents the requested
!                     tolerance from being achieved.
!                 = 3 extremely bad integrand behaviour occurs
!                     at some points of the integration
!                     interval.
!                 = 6 the input is invalid, because
!                     (epsabs.le.0 and
!                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                     or limit.lt.1 or lenw.lt.limit*4.
!                     result, abserr, neval, last are set
!                     to zero.
!                     except when lenw is invalid, iwork(1),
!                     work(limit*2+1) and work(limit*3+1) are
!                     set to zero, work(1) is set to a and
!                     work(limit+1) to b.
!
! dimensioning parameters
!    limit - integer
!            dimensioning parameter for iwork
!            limit determines the maximum number of subintervals
!            in the partition of the given integration interval
!            (a,b), limit.ge.1.
!            if limit.lt.1, the routine will end with ier = 6.
!
!    lenw  - integer
!            dimensioning parameter for work
!            lenw must be at least limit*4.
!            if lenw.lt.limit*4, the routine will end with
!            ier = 6.
!
!    last  - integer
!            on return, last equals the number of subintervals
!            produced in the subdivision process, which
!            determines the number of significant elements
!            actually in the work arrays.
!
! work arrays
!    iwork - integer
!            vector of dimension at least limit, the first k
!            elements of which contain pointers to the error
!            estimates over the subintervals, such that
!            work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
!            form a decreasing sequence with k = last if
!            last.le.(limit/2+2), and k = limit+1-last otherwise
!
!    work  - real
!            vector of dimension at least lenw
!            on return
!            work(1), ..., work(last) contain the left end
!            points of the subintervals in the partition of
!             (a,b),
!            work(limit+1), ..., work(limit+last) contain the
!             right end points,
!            work(limit*2+1), ..., work(limit*2+last) contain
!             the integral approximations over the subintervals,
!            work(limit*3+1), ..., work(limit*3+last) contain
!             the error estimates.
!
!***references  (none)
!***routines called  qage,xerror
!***end prologue  qag
!
      TREAL a,abserr,b,epsabs,epsrel,f,result,work,p
      TINTEGER ier,iwork,key,lenw,limit,lvl,l1,l2,l3,neval, last
!
      dimension iwork(limit),work(lenw)
!
      external f
!
! check validity of lenw.
!
!***first executable statement  qag
      ier = 6
      neval = 0
      last = 0
      result = C_0_R
      abserr = C_0_R
      if(limit.lt.1.or.lenw.lt.limit*4) go to 10
!
! prepare call for qage.
!
      l1 = limit+1
      l2 = limit+l1
      l3 = limit+l2
!
      call qage(f,p,a,b,epsabs,epsrel,key,limit,result,abserr,neval,&
        ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
! call error handler if necessary.
!
      lvl = 0
10    if(ier.eq.6) lvl = 1
      if(ier.ne.0) call xerror(26habnormal return from  qag ,&
        26,ier,lvl)
      return
      end
