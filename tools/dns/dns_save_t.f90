SUBROUTINE DNS_SAVE_T( icount, t_dat )
  
#include "types.h"

  USE DNS_GLOBAL, ONLY : nspa_rest, itime, rtime
  
  IMPLICIT NONE
  
  TINTEGER icount
  TREAL t_dat(2,nspa_rest)
  
  t_dat(1,icount) = M_REAL(itime)
  t_dat(2,icount) = rtime
  
  RETURN
END SUBROUTINE DNS_SAVE_T
      
