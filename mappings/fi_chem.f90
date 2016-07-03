#include "types.h"
#include "dns_const.h"

SUBROUTINE FI_CHEM(chemistry, nx,ny,nz, is, s, source)

  USE DNS_TYPES, ONLY : term_structure

  IMPLICIT NONE

  TYPE(term_structure),         INTENT(IN)  :: chemistry
  TINTEGER,                     INTENT(IN)  :: nx,ny,nz, is
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: source

! -----------------------------------------------------------------------
!########################################################################
  IF ( chemistry%type .EQ. EQNS_CHEM_QUADRATIC ) THEN

     IF ( is .EQ. 2 ) source = chemistry%parameters(is) * s(:,2) *s(:,3)
     IF ( is .EQ. 3 ) source = chemistry%parameters(is) * s(:,2) *s(:,3)
     IF ( is .EQ. 4 ) source = chemistry%parameters(is) * s(:,2) *s(:,3)
        
  ENDIF  
  
  RETURN
END SUBROUTINE FI_CHEM
