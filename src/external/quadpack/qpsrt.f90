      subroutine qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
      IMPLICIT NONE

#include "types.h"
!***begin prologue  qpsrt
!***refer to  qage,qagie,qagpe,qagse,qawce,qawse,qawoe
!***routines called  (none)
!***keywords  sequential sorting
!***description
!
! 1.        qpsrt
!   ordering routine
!      standard fortran subroutine
!      real version
!
! 2.        purpose
!      this routine maintains the descending ordering
!      in the list of the local error estimates resulting from
!      the interval subdivision process. at each call two error
!      estimates are inserted using the sequential search
!      method, top-down for the largest error estimate
!      and bottom-up for the smallest error estimate.
!
! 3.        calling sequence
!      call qpsrt(limit,last,maxerr,ermax,elist,iord,nrmax)
!
!   parameters (meaning at output)
!      limit  - integer
!               maximum number of error estimates the list
!               can contain
!
!      last   - integer
!               number of error estimates currently
!               in the list
!
!      maxerr - integer
!               maxerr points to the nrmax-th largest error
!               estimate currently in the list
!
!      ermax  - real
!               nrmax-th largest error estimate
!               ermax = elist(maxerr)
!
!      elist  - real
!               vector of dimension last containing
!               the error estimates
!
!      iord   - integer
!               vector of dimension last, the first k
!               elements of which contain pointers
!               to the error estimates, such that
!               elist(iord(1)),... , elist(iord(k))
!               form a decreasing sequence, with
!               k = last if last.le.(limit/2+2), and
!               k = limit+1-last otherwise
!
!      nrmax  - integer
!               maxerr = iord(nrmax)
!
! 4.        no subroutines or functions needed
!***end prologue  qpsrt
!
      TREAL elist,ermax,errmax,errmin
      TINTEGER i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last,limit,maxerr,&
        nrmax
      dimension elist(last),iord(last)
!
!   check whether the list contains more than
!   two error estimates.
!
!***first executable statement  qpsrt
      if(last.gt.2) go to 10
      iord(1) = 1
      iord(2) = 2
      go to 90
!
!   this part of the routine is only executed
!   if, due to a difficult integrand, subdivision
!   increased the error estimate. in the normal case
!   the insert procedure should start after the
!   nrmax-th largest error estimate.
!
   10 errmax = elist(maxerr)
      if(nrmax.eq.1) go to 30
      ido = nrmax-1
      do 20 i = 1,ido
        isucc = iord(nrmax-1)
! ***jump out of do-loop
        if(errmax.le.elist(isucc)) go to 30
        iord(nrmax) = isucc
        nrmax = nrmax-1
   20    continue
!
!   compute the number of elements in the list to
!   be maintained in descending order. this number
!   depends on the number of subdivisions still
!   allowed.
!
   30 jupbn = last
      if(last.gt.(limit/2+2)) jupbn = limit+3-last
      errmin = elist(last)
!
!   insert errmax by traversing the list top-down,
!   starting comparison from the element elist(iord(nrmax+1)).
!
      jbnd = jupbn-1
      ibeg = nrmax+1
      if(ibeg.gt.jbnd) go to 50
      do 40 i=ibeg,jbnd
        isucc = iord(i)
! ***jump out of do-loop
        if(errmax.ge.elist(isucc)) go to 60
        iord(i-1) = isucc
   40 continue
   50 iord(jbnd) = maxerr
      iord(jupbn) = last
      go to 90
!
!   insert errmin by traversing the list bottom-up.
!
   60 iord(i-1) = maxerr
      k = jbnd
      do 70 j=i,jbnd
        isucc = iord(k)
! ***jump out of do-loop
        if(errmin.lt.elist(isucc)) go to 80
        iord(k+1) = isucc
        k = k-1
   70 continue
      iord(i) = last
      go to 90
   80 iord(k+1) = last
!
!   set maxerr and ermax.
!
   90 maxerr = iord(nrmax)
      ermax = elist(maxerr)
      return
      end
