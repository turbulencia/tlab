#include "types.h"
#include "dns_const.h"

!########################################################################
!########################################################################
SUBROUTINE DNS_FILTER()

  USE DNS_CONSTANTS, ONLY : lfile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow,inb_scal, isize_field
  USE DNS_GLOBAL,    ONLY : imode_eqns, imode_sim
  USE DNS_GLOBAL,    ONLY : itime, rtime
  USE DNS_GLOBAL,    ONLY : FilterDomain
  USE TLAB_ARRAYS
  USE DNS_LOCAL,     ONLY : ilimit_scal, s_bound_min, s_bound_max
  USE DNS_LOCAL,     ONLY : nitera_stats_spa, nitera_first,nitera_stats
  USE STATISTICS

  IMPLICIT NONE

#include "integers.h"

  ! -----------------------------------------------------------------------
  LOGICAL flag_save
  TINTEGER iq,is,ij
  CHARACTER*250 line

  ! #######################################################################
  WRITE(line,*) itime; line = 'Filtering fields at It'//TRIM(ADJUSTL(line))//'.'
  CALL IO_WRITE_ASCII(lfile,line)

  IF ( MOD(itime-nitera_first,nitera_stats)  == 0 ) THEN; flag_save = .TRUE.
  ELSE;                                                   flag_save = .FALSE.
  ENDIF
  ! Might be better to filter the pressure instead of the energy in compressible flows
  ! IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
  !    iq_loc = (/ 5,1,2,3,6 /) ! Filtered variables: rho, u,v,w, p
  ! ELSE
  !    iq_loc = (/ 1,2,3 /)
  ! ENDIF

  ! Statistics
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. nitera_stats_spa .GT. 0 ) THEN
    CALL AVG_TKE_ZT_REDUCE(q(1,5), q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), txc(1,3), txc(1,4), &
        txc(1,5), txc(1,6), txc(1,7), mean_flow, wrk2d)
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. stats_filter .AND. flag_save ) THEN
    CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
    CALL AVG_ENERGY_XZ(i1, itime, rtime, imax,jmax,kmax,q, txc, mean, wrk3d)
  ENDIF

  ! filtering
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

  IF ( ilimit_scal .EQ. 1 ) THEN
    DO ij = 1,isize_field
      s(ij,is) = MIN(MAX(s(ij,is),s_bound_min(is)), s_bound_max(is))
    ENDDO
  ENDIF

  ! statistics
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. stats_filter .AND. flag_save ) THEN
    CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
    CALL AVG_ENERGY_XZ(i2, itime, rtime, imax,jmax,kmax, q, txc, mean, wrk3d)
  ENDIF

  ! recalculation of diagnostic variables
  CALL FI_DIAGNOSTIC( imax,jmax,kmax, q,s, wrk3d )

  RETURN
END SUBROUTINE DNS_FILTER

! #######################################################################
! #######################################################################
#define rR(j)     mean2d(j,1)
#define fU(j)     mean2d(j,2)
#define fV(j)     mean2d(j,3)
#define fW(j)     mean2d(j,4)
#define Kin_0(j)  mean2d(j,5)
#define Kin_1(j)  mean2d(j,6)
#define Eps_0(j)  mean2d(j,7)
#define Eps_1(j)  mean2d(j,8)

SUBROUTINE AVG_ENERGY_XZ(iflag, itime, rtime, nx,ny,nz, q, eps, mean2d, wrk3d)

  USE DNS_GLOBAL, ONLY : imode_eqns
  USE DNS_GLOBAL, ONLY : g, area

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iflag, itime, nx,ny,nz
  TREAL rtime
  TREAL, DIMENSION(nx,ny,nz,*), INTENT(IN)    :: q
  TREAL, DIMENSION(nx,ny,nz),   INTENT(IN)    :: eps
  TREAL, DIMENSION(nx,ny,nz),   INTENT(INOUT) :: wrk3d
  TREAL, DIMENSION(ny,8),       INTENT(INOUT) :: mean2d

  ! -----------------------------------------------------------------------
  TINTEGER i, j, k
  TREAL AVG_IK
  CHARACTER*64 fname, varnames(1), groupnames(1)
#ifdef USE_MPI
  INTEGER ims_err, ims_pro
#endif

  ! #######################################################################
  DO j = 1,ny

    IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
      rR(j) = AVG_IK(nx, ny, nz, j, q(1,1,1,5), g(1)%jac, g(3)%jac, area)
      DO k = 1,nz; DO i = 1,nx
        wrk3d(i,1,k) = q(i,j,k,5) *q(i,j,k,1)
        wrk3d(i,2,k) = q(i,j,k,5) *q(i,j,k,2)
        wrk3d(i,3,k) = q(i,j,k,5) *q(i,j,k,3)
      ENDDO; ENDDO

    ELSE
      rR(j) = C_1_R
      DO k=1, nz; DO i=1, nx
        wrk3d(i,1,k) = q(i,j,k,1)
        wrk3d(i,2,k) = q(i,j,k,2)
        wrk3d(i,3,k) = q(i,j,k,3)
      ENDDO; ENDDO

    ENDIF

    fU(j) = AVG_IK(nx,ny,nz, i1, wrk3d, g(1)%jac,g(3)%jac, area) /rR(j)
    fV(j) = AVG_IK(nx,ny,nz, i2, wrk3d, g(1)%jac,g(3)%jac, area) /rR(j)
    fW(j) = AVG_IK(nx,ny,nz, i3, wrk3d, g(1)%jac,g(3)%jac, area) /rR(j)

    DO k = 1,nz; DO i = 1,nx
      wrk3d(i,1,k) = rR(j) *((q(i,j,k,1)-fU(j))**2 + (q(i,j,k,2)-fV(j))**2 + (q(i,j,k,3)-fW(j))**2) *C_05_R
    ENDDO; ENDDO

    IF ( iflag .EQ. 1 ) THEN
      Kin_0(j) = AVG_IK(nx,ny,nz, i1, wrk3d, g(1)%jac,g(3)%jac, area)
      Eps_0(j) = AVG_IK(nx,ny,nz, j,  eps,   g(1)%jac,g(3)%jac, area)
    ELSE
      Kin_1(j) = AVG_IK(nx,ny,nz, i1, wrk3d, g(1)%jac,g(3)%jac, area)
      Eps_1(j) = AVG_IK(nx,ny,nz, j,  eps,   g(1)%jac,g(3)%jac, area)
    ENDIF

  ENDDO

  IF ( iflag .EQ. 2 ) THEN
    varnames(1) = 'TkeBefore TkeAfter EpsBefore EpsAfter'
    groupnames(1) = ''
    WRITE(fname,*) itime; fname='kin'//TRIM(ADJUSTL(fname))
    CALL IO_WRITE_AVERAGES( fname, itime,rtime, ny,4,1, g(2)%nodes, varnames,groupnames, mean2d(1,5) )
  ENDIF

  RETURN
END SUBROUTINE AVG_ENERGY_XZ
