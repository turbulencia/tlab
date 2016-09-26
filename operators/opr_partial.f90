#include "types.h"
#include "dns_const.h"

SUBROUTINE OPR_PARTIAL(imode_fdm, nlines, g, u,up, bcs_min,bcs_max, wrk2d)

  USE DNS_TYPES, ONLY : grid_structure
  
  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: imode_fdm
  TINTEGER,                        INTENT(IN)    :: nlines  ! # of lines to be solved
  TINTEGER,                        INTENT(IN)    :: bcs_min ! BC derivative: 0 biased, non-zero
  TINTEGER,                        INTENT(IN)    :: bcs_max !                1 forced to zero
  TYPE(grid_structure),            INTENT(IN)    :: g
  TREAL, DIMENSION(nlines*g%size), INTENT(IN)    :: u
  TREAL, DIMENSION(nlines*g%size), INTENT(OUT)   :: up
  TREAL, DIMENSION(nlines),        INTENT(INOUT) :: wrk2d

! -------------------------------------------------------------------
  TINTEGER ip

! ###################################################################
  IF ( g%periodic ) THEN
     SELECT CASE( imode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4P_RHS(g%size,nlines, u, up)
        
     CASE( FDM_COM6_JACOBIAN, FDM_COM6_DIRECT ) ! Direct = Jacobian because uniform grid
        CALL FDM_C1N6P_RHS(g%size,nlines, u, up)
        
     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8P_RHS(g%size,nlines, u, up)
        
     END SELECT
     
     ! ip  = inb_grid_1 - 1
     ! CALL TRIDPSS(g%size,nlines, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3),dx(1,ip+4),dx(1,ip+5), up,wrk2d)
     CALL TRIDPSS(g%size,nlines, g%lu1(1,1),g%lu1(1,2),g%lu1(1,3),g%lu1(1,4),g%lu1(1,5), up,wrk2d)

! -------------------------------------------------------------------
  ELSE
     SELECT CASE( imode_fdm )
        
     CASE( FDM_COM4_JACOBIAN )
        CALL FDM_C1N4_RHS(g%size,nlines, bcs_min,bcs_max, u, up)

     CASE( FDM_COM6_JACOBIAN )
        CALL FDM_C1N6_RHS(g%size,nlines, bcs_min,bcs_max, u, up)

     CASE( FDM_COM8_JACOBIAN )
        CALL FDM_C1N8_RHS(g%size,nlines, bcs_min,bcs_max, u, up)

     CASE( FDM_COM6_DIRECT   ) ! Not yet implemented
        CALL FDM_C1N6_RHS(g%size,nlines, bcs_min,bcs_max, u, up)

     END SELECT
     
     ! ip = inb_grid_1 + (bcs_min + bcs_max*2)*3 - 1
     ! CALL TRIDSS(g%size,nlines, dx(1,ip+1),dx(1,ip+2),dx(1,ip+3), up)
     ip = (bcs_min + bcs_max*2)*3 
     CALL TRIDSS(g%size,nlines, g%lu1(1,ip+1),g%lu1(1,ip+2),g%lu1(1,ip+3), up)

  ENDIF
  
  RETURN
END SUBROUTINE OPR_PARTIAL

! ###################################################################
! ###################################################################
! Add here OPR_DDX, OPR_DDY, and OPR_DDZ
