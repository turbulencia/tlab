#include "types.h"
#include "dns_const.h"

SUBROUTINE FI_SUBSIDENCE(subsidence, nx,ny,nz, s, source, wrk2d,wrk3d)

  USE DNS_TYPES,  ONLY : term_dt
  USE DNS_GLOBAL, ONLY : g

  IMPLICIT NONE

#include "integers.h"

  TYPE(term_dt),                INTENT(IN)    :: subsidence
  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT)   :: source
  TREAL, DIMENSION(*),          INTENT(INOUT) :: wrk2d, wrk3d

! -----------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL W_LOC

  TINTEGER idummy   ! To use old wrappers to calculate derivatives
  TREAL    dummy(1)

!########################################################################
  SELECT CASE( subsidence%type )

  CASE( EQNS_SUB_CONSTANT )
     
     CALL PARTIAL_Y(idummy, nx,ny,nz, idummy, dummy, s, source, i0,i0, dummy,wrk2d,wrk3d)

     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        W_LOC = g(2)%nodes(is) *subsidence%parameters(1)
        
        DO i = 1,nx
           ij = ij +1
           
           source(ij) = source(ij) *W_LOC
           
        ENDDO
     
     ENDDO

  END SELECT
  
  RETURN
END SUBROUTINE FI_SUBSIDENCE
