SUBROUTINE DNS_STOP(ierc)
  
#include "types.h"

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE
  
  INTEGER ierc

  CHARACTER*32 line
  
  WRITE(line,*) ierc
  line = 'Error Code = '//TRIM(ADJUSTL(line))
  CALL IO_WRITE_ASCII(efile,line)
  
  CALL DNS_END(1)
  
  STOP
END SUBROUTINE DNS_STOP
