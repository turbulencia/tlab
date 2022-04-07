#include "types.h"
#include "dns_const.h"

SUBROUTINE FI_RANDOM(random, nx,ny,nz, h, tmp)

  USE TLAB_TYPES,  ONLY : term_dt

  IMPLICIT NONE
  
  TYPE(term_dt),                INTENT(IN)    :: random
  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx,ny,nz),   INTENT(INOUT) :: h
  TREAL, DIMENSION(nx,ny,nz)                  :: tmp

  ! -----------------------------------------------------------------------
  TREAL dummy

  dummy=random%parameters(1)

  CALL RANDOM_NUMBER(tmp) 

  tmp = dummy*(tmp*2.0 - 1.0)
  h = h * (1+tmp)

  RETURN
  
END SUBROUTINE FI_RANDOM
