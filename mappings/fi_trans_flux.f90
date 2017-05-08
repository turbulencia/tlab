#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2015/06/25 - A de Lozar
!#              Created
!# 2016/06/25 - A de Lozar
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the transport terms due to settling in an airwater mixture.
!# There are two options for the settling model, a simplified calculation and
!# a more exact one (differences of 10^-4)
!#
!########################################################################
SUBROUTINE FI_TRANS_FLUX(transport, flag_grad, nx,ny,nz, is, s,trans, tmp, wrk2d,wrk3d)

  USE DNS_TYPES,  ONLY : term_dt
  USE DNS_GLOBAL, ONLY : g
  
  IMPLICIT NONE

  TYPE(term_dt),                INTENT(IN)    :: transport
  TINTEGER,                     INTENT(IN)    :: nx,ny,nz, flag_grad
  TINTEGER,                     INTENT(IN)    :: is
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)    :: s
  TREAL, DIMENSION(nx*ny*nz,1), INTENT(OUT)   :: trans ! Transport component. It could have eventually three directions
  TREAL, DIMENSION(nx*ny*nz,1), INTENT(INOUT) :: tmp   ! To avoid re-calculations when repetedly calling this routine
  TREAL, DIMENSION(*),          INTENT(INOUT) :: wrk2d,wrk3d

! -----------------------------------------------------------------------
  TREAL dummy, exponent
  TINTEGER is_ref, bcs
  
!########################################################################
  bcs = 0
  
  exponent = transport%auxiliar(1)
  is_ref   = transport%scalar(1)

  IF     ( transport%type .EQ. EQNS_TRANS_AIRWATERSIMPLIFIED ) THEN
     IF ( flag_grad .EQ. 1 ) THEN
        CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), s(1,is_ref), tmp, wrk3d, wrk2d,wrk3d)
        IF ( exponent .GT. C_0_R ) tmp(:,1) = tmp(:,1) *(s(:,is_ref)**exponent)
     ENDIF

     dummy = transport%parameters(is) *( C_1_R + exponent )
     trans(:,1) = dummy*tmp(:,1)

  ELSEIF ( transport%type .EQ. EQNS_TRANS_AIRWATER ) THEN
     dummy = C_1_R + exponent

     IF ( exponent .GT. C_0_R ) THEN
        tmp(:,1) = (transport%parameters(is) - transport%parameters(5)*s(:,is)) *(s(:,is_ref)**dummy)
     ELSE
        tmp(:,1) = (transport%parameters(is) - transport%parameters(5)*s(:,is)) * s(:,is_ref)
     ENDIF
     CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), tmp(1,1), trans(1,1), wrk3d, wrk2d,wrk3d)

  ENDIF
  
  RETURN
END SUBROUTINE FI_TRANS_FLUX
