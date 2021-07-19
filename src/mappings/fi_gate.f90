#include "types.h"
#include "dns_const.h"

SUBROUTINE FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
                   nx,ny,nz, igate_size, gate_threshold, q,s, txc, gate, wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : g

  IMPLICIT NONE

  TINTEGER,                          INTENT(IN)    :: opt_cond, opt_cond_relative, opt_cond_scal
  TINTEGER,                          INTENT(IN)    :: nx,ny,nz, igate_size
  TREAL,                             INTENT(INOUT) :: gate_threshold(igate_size)
  TREAL,      DIMENSION(nx*ny*nz,*), INTENT(IN)    :: q,s
  TREAL,      DIMENSION(nx*ny*nz,5), INTENT(INOUT) :: txc
  INTEGER(1), DIMENSION(nx*ny*nz),   INTENT(OUT)   :: gate
  TREAL,      DIMENSION(*),          INTENT(INOUT) :: wrk2d,wrk3d

! -----------------------------------------------------------------------
  TINTEGER ij,n
  TREAL umin,umax, gate_threshold_loc(igate_size)
  TINTEGER bcs(2,2)

! #######################################################################
! Preparing indicator field

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  IF      ( opt_cond .EQ. 2 ) THEN ! Based on scalar
     txc(:,1) = s(:,opt_cond_scal)

  ELSE IF ( opt_cond .EQ. 3 ) THEN ! Based on vorticity
     CALL FI_VORTICITY(nx,ny,nz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3), wrk2d,wrk3d)

  ELSE IF ( opt_cond .EQ. 4 ) THEN ! Based on scalar gradient
     CALL FI_GRADIENT(nx,ny,nz, s, txc(1,1),txc(1,2),txc(1,3), wrk2d,wrk3d)

  ELSE IF ( opt_cond .EQ. 5 ) THEN ! Based on vertical velocity
     txc(:,1) = q(:,2)

  ELSE IF ( opt_cond .EQ. 6 .OR. opt_cond .EQ. 7 ) THEN ! Based on scalar fluctuation
     txc(:,1) = s(:,opt_cond_scal)
     CALL FI_FLUCTUATION_INPLACE(nx,ny,nz, txc(1,1))

  ELSE IF ( opt_cond .EQ. 8 ) THEN ! Based on potential vorticity
     CALL FI_CURL(nx,ny,nz, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)

     txc(:,4) = s(:,opt_cond_scal)
     CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), txc(1,4),txc(1,5), wrk3d, wrk2d,wrk3d)
     txc(:,1) =           txc(:,1)*txc(:,5)
     CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), txc(1,4),txc(1,5), wrk3d, wrk2d,wrk3d)
     txc(:,1) = txc(:,1) +txc(:,2)*txc(:,5)
     CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), txc(1,4),txc(1,5), wrk3d, wrk2d,wrk3d)
     txc(:,1) = txc(:,1) +txc(:,3)*txc(:,5)

     txc(:,1) = txc(:,1)*txc(:,1)
     txc(:,1) = LOG(txc(:,1)+C_SMALL_R)

  ENDIF

! #######################################################################
! define gate field

  IF ( opt_cond .EQ. 7 ) THEN ! double conditioning; flux
     DO ij = 1,nx*ny*nz
        IF      ( txc(ij,1) .GT. C_0_R .AND. q(ij,2) .GE. C_0_R ) THEN; gate(ij) = 1;
        ELSE IF ( txc(ij,1) .LE. C_0_R .AND. q(ij,2) .GT. C_0_R ) THEN; gate(ij) = 2;
        ELSE IF ( txc(ij,1) .LT. C_0_R .AND. q(ij,2) .LE. C_0_R ) THEN; gate(ij) = 3;
        ELSE;                                                           gate(ij) = 4; ENDIF
     ENDDO

  ELSE                             ! Local file
     IF ( opt_cond_relative .EQ. 1 ) THEN ! case of threshold relative to maximum
        CALL MINMAX(nx,ny,nz, txc(1,1), umin,umax)
        gate_threshold_loc = gate_threshold *(umax-umin)

     ELSE
        gate_threshold_loc = gate_threshold

     ENDIF

! Loop over all thresholds (# thresholds is 1 less than # gate levels)
! Thresholds are in ascending order
     DO ij = 1,nx*ny*nz
        DO n = 1,igate_size-1
           IF ( txc(ij,1) .LT. gate_threshold_loc(n) ) EXIT
        ENDDO
        gate(ij) = INT(n,KIND=1) ! note that gate can get -- correctly -- the value igate_size
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE FI_GATE
