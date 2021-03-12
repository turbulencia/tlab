      SUBROUTINE QUADNG(f,p,a,b,epsabs,epsrel,result,abserr,neval,&
           ier)

      IMPLICIT NONE

#include "types.h"

      TREAL a,abserr,b,epsabs,epsrel,f,result
      TREAL ier,neval,p
      external f

      CALL qng(f,p,a,b,epsabs,epsrel,result,abserr,neval,ier)

      RETURN
      END
