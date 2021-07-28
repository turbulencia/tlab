#include "types.h"
#include "dns_const.h"

SUBROUTINE FI_SUBSIDENCE(subsidence, nx,ny,nz, s, source, wrk1d,wrk2d,wrk3d)

  USE DNS_TYPES,  ONLY : term_dt
  USE TLAB_VARS, ONLY : g

  IMPLICIT NONE

#include "integers.h"

  TYPE(term_dt),                INTENT(IN)    :: subsidence
  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT)   :: source
  TREAL, DIMENSION(ny,2),       INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(*),          INTENT(INOUT) :: wrk2d, wrk3d

! -----------------------------------------------------------------------
  TINTEGER ij, i, jk, is, bcs(2,2)
  TREAL W_LOC

!########################################################################
  bcs = 0
  
  SELECT CASE( subsidence%type )

  CASE( EQNS_SUB_CONSTANT_LOCAL )
     
     CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), s, source, wrk3d, wrk2d,wrk3d)

     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        W_LOC = g(2)%nodes(is) *subsidence%parameters(1)
        
        DO i = 1,nx
           ij = ij +1
           
           source(ij) = source(ij) *W_LOC
           
        ENDDO
     
     ENDDO

  CASE( EQNS_SUB_CONSTANT_GLOBAL )

     CALL AVG1V2D_V(nx,ny,nz, i1, s, wrk1d(1,1), wrk1d(1,2)) ! Calculate averaged scalar into wrk1d(1,1)

     CALL OPR_PARTIAL_Y(OPR_P1, i1,ny,i1, bcs, g(2), wrk1d(1,1), wrk1d(1,2), wrk3d, wrk2d,wrk3d)

     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        W_LOC = g(2)%nodes(is) *subsidence%parameters(1) *wrk1d(is,2)
        
        DO i = 1,nx
           ij = ij +1
           
           source(ij) = W_LOC
           
        ENDDO
     
     ENDDO
  END SELECT
  
  RETURN
END SUBROUTINE FI_SUBSIDENCE
