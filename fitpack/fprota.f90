      subroutine fprota(cos,sin,a,b)
!  subroutine fprota applies a givens rotation to a and b.
!  ..

      IMPLICIT NONE

#include "types.h"

!  ..scalar arguments..
      TREAL cos,sin,a,b
! ..local scalars..
      TREAL stor1,stor2
!  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end
