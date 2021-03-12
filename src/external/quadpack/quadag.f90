      SUBROUTINE QUADAG(f,p,a,b,epsabs,epsrel,key,result,abserr,neval,&
           ier, limit,lenw,last,iwork,work)

      IMPLICIT NONE

#include "types.h"

      TREAL a,abserr,b,epsabs,epsrel,f,result
      TREAL work, p
      TINTEGER key, limit, lenw, last, iwork, ier, neval 
      external f

      CALL qag(f,p,a,b,epsabs,epsrel,key,result,abserr,neval,ier,&
           limit,lenw,last,iwork,work)

      RETURN
      END
