      subroutine fpcyt1(a,n,nn)
! (l u)-decomposition of a cyclic tridiagonal matrix with the non-zero
! elements stored as follows
!
!| a(1,2) a(1,3)                                    a(1,1)  |
!| a(2,1) a(2,2) a(2,3)                                     |
!|        a(3,1) a(3,2) a(3,3)                              |
!|               ...............                            |
!|                               a(n-1,1) a(n-1,2) a(n-1,3) |
!| a(n,3)                                  a(n,1)   a(n,2)  |
!
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TINTEGER n,nn
!  ..array arguments..
      TREAL a(nn,6)
!  ..local scalars..
      TREAL aa,beta,gamma,sum,teta,v,one
      TINTEGER i,n1,n2
!  ..
!  set constant
      one = C_1_R
      n2 = n-2
      beta = one/a(1,2)
      gamma = a(n,3)
      teta = a(1,1)*beta
      a(1,4) = beta
      a(1,5) = gamma
      a(1,6) = teta
      sum = gamma*teta
      do 10 i=2,n2
         v = a(i-1,3)*beta
         aa = a(i,1)
         beta = one/(a(i,2)-aa*v)
         gamma = -gamma*v
         teta = -teta*aa*beta
         a(i,4) = beta
         a(i,5) = gamma
         a(i,6) = teta
         sum = sum+gamma*teta
  10  continue
      n1 = n-1
      v = a(n2,3)*beta
      aa = a(n1,1)
      beta = one/(a(n1,2)-aa*v)
      gamma = a(n,1)-gamma*v
      teta = (a(n1,3)-teta*aa)*beta
      a(n1,4) = beta
      a(n1,5) = gamma
      a(n1,6) = teta
      a(n,4) = one/(a(n,2)-(sum+gamma*teta))
      return
      end
