#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2012/07/16 Cedrick Ansorge - created 
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT(dte,etime, kex,kim, kco, &
     q,hq,s,hs, txc, vaux, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : isize_field, isize_txc_field
  USE DNS_GLOBAL,    ONLY : imode_eqns

  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : VA_BCS_HT, VA_BCS_HB
  USE DNS_LOCAL,  ONLY : rkm_mode

  IMPLICIT NONE

  TREAL,                              INTENT(IN)    :: dte,etime 
  TREAL,                              INTENT(IN)    :: kex,kim, kco
  TREAL, DIMENSION(isize_field, *),   INTENT(INOUT) :: q,hq, s,hs
  TREAL, DIMENSION(isize_txc_field,9),INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk1d,wrk2d,wrk3d, vaux
! -----------------------------------------------------------------------

! #######################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT')
#endif

! #######################################################################
! Evaluate standard RHS of incompressible equations
! #######################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE ) THEN
     IF      ( rkm_mode .EQ. RKM_IMP3_DIFFUSION ) THEN
        CALL RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2(&
             dte, kex,kim,kco,  &
             q, hq, q(:,1),q(:,2),q(:,3), hq(1,1),hq(1,2),hq(1,3), s,hs, &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), &
             vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),vaux, wrk1d,wrk2d,wrk3d)
! pressure-correction algorithm; to be checked
        ! CALL RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_3(&
        !      dte, kex,kim,kco,  &
        !      q, hq, q(:,1),q(:,2),q(:,3), hq(1,1),hq(1,2),hq(1,3), s,hs, &
        !      txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), txc(1,8), &
        !      vaux(vindex(VA_BCS_HB)),vaux(vindex(VA_BCS_HT)),&
        !      vaux,wrk1d,wrk2d,wrk3d)
     ELSE 
        CALL IO_WRITE_ASCII(efile,'TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT. Undeveloped formulation.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)

     ENDIF

! #######################################################################
! Evaluate standard RHS of anelastic equations
! #######################################################################
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL IO_WRITE_ASCII(efile,'TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT. Undeveloped anelastic formulation.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)

  ENDIF
  
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT')
#endif

  RETURN

END SUBROUTINE TIME_SUBSTEP_INCOMPRESSIBLE_IMPLICIT
