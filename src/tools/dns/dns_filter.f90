#include "types.h"
#include "dns_const.h"

!########################################################################
!########################################################################
SUBROUTINE DNS_FILTER()

  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, inb_flow,inb_scal, isize_field
  USE TLAB_VARS,    ONLY : imode_eqns, imode_sim
  USE TLAB_VARS,    ONLY : itime, rtime
  USE TLAB_VARS,    ONLY : FilterDomain
  USE TLAB_VARS,    ONLY : g, area
  USE TLAB_ARRAYS
  USE TLAB_PROCS
  USE OPR_FILTERS
  USE DNS_LOCAL,     ONLY : ilimit_scal, s_bound_min, s_bound_max
  USE DNS_LOCAL,     ONLY : nitera_stats_spa, nitera_first,nitera_stats
  USE STATISTICS
  USE AVGS, ONLY: AVG_IK_V

  IMPLICIT NONE

#include "integers.h"

  ! -----------------------------------------------------------------------
  TINTEGER iq,is,ij
  CHARACTER*250 line
  CHARACTER*64 fname, varnames(1), groupnames(1)

  ! #######################################################################
  WRITE(line,*) itime; line = 'Filtering fields at It'//TRIM(ADJUSTL(line))//'.'
  CALL TLAB_WRITE_ASCII(lfile,line)

  ! -------------------------------------------------------------------
  ! Statistics
#define Tke0(j)   mean(j,1)
#define Eps0(j)   mean(j,2)
#define Tke1(j)   mean(j,3)
#define Eps1(j)   mean(j,4)

  varnames(1) = 'TkeBefore EpsBefore TkeAfter EpsAfter'
  groupnames(1) = ''

  IF ( MOD(itime-nitera_first,nitera_stats_spa) == 0 ) THEN   ! Accumulate statistics in spatially evolving cases
    CALL AVG_TKE_ZT_REDUCE(q(1,5), q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), txc(1,3), txc(1,4), &
        txc(1,5), txc(1,6), txc(1,7), mean_flow, wrk2d)
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. MOD(itime-nitera_first,nitera_stats) == 0 ) THEN
    CALL FI_RTKE(imax,jmax,kmax,q, wrk1d,wrk3d)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Tke0(1), wrk1d, area)
    CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, txc, g(1)%jac,g(3)%jac, Eps0(1), wrk1d, area)
  ENDIF

  ! -------------------------------------------------------------------
  ! filtering

  ! Might be better to filter the pressure instead of the energy in compressible flows
  ! IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
  !    iq_loc = (/ 5,1,2,3,6 /) ! Filtered variables: rho, u,v,w, p
  ! ELSE
  !    iq_loc = (/ 1,2,3 /)
  ! ENDIF

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! contruct fields per unit volume
    DO iq = 1,inb_flow-1
      q(:,iq) = q(:,iq) *q(:,inb_flow)
    ENDDO
    DO is = 1,inb_scal
      s(:,is) = s(:,is) *q(:,inb_flow)
    ENDDO
  ENDIF

  DO iq = 1,inb_flow
    CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, q(1,iq), wrk1d,wrk2d,txc)
  ENDDO
  DO is = 1,inb_scal
    CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, s(1,is), wrk1d,wrk2d,txc)
  ENDDO

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! re-contruct fields per unit mass
    DO iq = 1,inb_flow-1
      q(:,iq) = q(:,iq) /q(:,inb_flow)
    ENDDO
    DO is = 1,inb_scal
      s(:,is) = s(:,is) /q(:,inb_flow)
    ENDDO
  ENDIF

  ! -------------------------------------------------------------------
  IF ( inb_scal .GT. 0 .AND. ilimit_scal .EQ. 1 ) THEN
     DO is =1,inb_scal
        DO ij = 1,isize_field
           s(ij,is) = MIN(MAX(s(ij,is),s_bound_min(is)), s_bound_max(is))
        ENDDO
     ENDDO
  ENDIF

  ! -------------------------------------------------------------------
  ! statistics
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. MOD(itime-nitera_first,nitera_stats) == 0 ) THEN
    CALL FI_RTKE(imax,jmax,kmax, q, wrk1d,wrk3d)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac,g(3)%jac, Tke1(1), wrk1d, area)
    CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
    CALL AVG_IK_V(imax,jmax,kmax, jmax, txc, g(1)%jac,g(3)%jac, Eps1(1), wrk1d, area)

    WRITE(fname,*) itime; fname='kin'//TRIM(ADJUSTL(fname))
    CALL IO_WRITE_AVERAGES( fname, itime,rtime, jmax,4,1, g(2)%nodes, varnames,groupnames, mean )

  ENDIF

  ! -------------------------------------------------------------------
  ! recalculation of diagnostic variables
  CALL FI_DIAGNOSTIC( imax,jmax,kmax, q,s, wrk3d )

  RETURN
END SUBROUTINE DNS_FILTER
