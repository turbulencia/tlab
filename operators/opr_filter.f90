#include "types.h"
#include "dns_const.h"

SUBROUTINE OPR_FILTER(nx,ny,nz, f, u, wrk1d,wrk2d,txc)
        
  USE DNS_TYPES,  ONLY : filter_dt
  USE DNS_GLOBAL, ONLY : isize_txc_field
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TYPE(filter_dt), DIMENSION(3),                 INTENT(IN)    :: f
  TREAL,           DIMENSION(nx*ny*nz),          INTENT(INOUT) :: u   ! Inplace operation
  TREAL,           DIMENSION(*),                 INTENT(INOUT) :: wrk1d, wrk2d
  TREAL,           DIMENSION(isize_txc_field,*), INTENT(INOUT) :: txc ! 2 if ADM type, 1 otherwise
  
! ###################################################################
  IF ( f(1)%type .NE. DNS_FILTER_NONE ) THEN
     CALL OPR_FILTER_X(nx,ny,nz, f(1), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF

  IF ( f(2)%type .NE. DNS_FILTER_NONE ) THEN
     CALL OPR_FILTER_Y(nx,ny,nz, f(2), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF
     
  IF ( f(3)%type .NE. DNS_FILTER_NONE ) THEN
     CALL OPR_FILTER_Z(nx,ny,nz, f(3), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF

  RETURN
END SUBROUTINE OPR_FILTER

