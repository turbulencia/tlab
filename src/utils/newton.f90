      SUBROUTINE NEWTON(FR, p, x0, imax, eps, niter, root)

      IMPLICIT NONE

#include "types.h"
      
      TREAL FR, p(*), x0, eps, root
      TINTEGER imax, niter
      TREAL x, delta
      TINTEGER i

      EXTERNAL FR

      x = x0
      DO i=1, imax
#ifdef _DEBUG
!     PRINT *, i, x
#endif
         delta = FR(x,p)
         x = x + delta
         IF ( ABS(delta) .LT. eps ) GOTO 11
      ENDDO

 11   root = x
      niter = i

      RETURN
      END
