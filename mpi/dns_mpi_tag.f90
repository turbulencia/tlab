#include "types.h"
  
SUBROUTINE DNS_MPI_TAGUPDT
  
  USE DNS_MPI, ONLY : ims_tag

  IMPLICIT NONE
  
  ims_tag = ims_tag+1
  
  IF ( ims_tag .GT. 32000 ) THEN
     CALL DNS_MPI_TAGRESET
  ENDIF
  
  RETURN
END SUBROUTINE DNS_MPI_TAGUPDT

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_TAGRESET
  
  USE DNS_MPI, ONLY : ims_tag

  IMPLICIT NONE
  
  ims_tag = 0
  
  RETURN
END SUBROUTINE DNS_MPI_TAGRESET
    
