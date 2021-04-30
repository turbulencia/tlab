#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#include "dns_const_mpi.h"
#include "avgij_map.h"

!########################################################################
!#
!# Performing the time integration over a given number of steps
!#
!########################################################################
SUBROUTINE TIME_INTEGRATION(q,hq, s,hs, q_inf,s_inf, txc, wrk1d,wrk2d,wrk3d, &
    l_q, l_hq, l_txc, l_comm)

  USE DNS_CONSTANTS, ONLY : tag_flow, tag_scal, tag_part, tag_traj, lfile
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field
  USE DNS_GLOBAL, ONLY : isize_particle
  USE DNS_GLOBAL, ONLY : imode_sim
  USE DNS_GLOBAL, ONLY : icalc_flow, icalc_scal, icalc_part
  USE DNS_GLOBAL, ONLY : visc
  USE DNS_GLOBAL, ONLY : itime, rtime
  USE DNS_LOCAL
  USE DNS_TOWER
  USE LAGRANGE_GLOBAL, ONLY : itrajectory, l_g
  USE BOUNDARY_INFLOW, ONLY : BOUNDARY_INFLOW_FILTER, FilterInflowStep
  USE BOUNDARY_BCS, ONLY : BcsFlowImin, BcsScalImin
  USE STATISTICS
  USE PARTICLE_TRAJECTORIES
  USE AVG_SCAL_ZT
  USE PLANES
#ifdef LES
  USE LES_GLOBAL, ONLY : iles
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(isize_field,*) :: q,hq, s,hs
  TREAL, DIMENSION(*)             :: txc
  TREAL, DIMENSION(*)             :: q_inf, s_inf
  TREAL, DIMENSION(*)             :: wrk1d, wrk2d, wrk3d

  TREAL, DIMENSION(isize_particle,*) :: l_q, l_hq, l_txc
  TREAL, DIMENSION(*)                :: l_comm

  ! -------------------------------------------------------------------
  CHARACTER*32 fname
  CHARACTER*250 line1
  LOGICAL flag_save

  ! ###################################################################
  ! Loop on iterations: itime counter
  ! ###################################################################
  itime = nitera_first

  WRITE(line1,*) itime; line1 = 'Starting time integration at It'//TRIM(ADJUSTL(line1))//'.'
  CALL IO_WRITE_ASCII(lfile,line1)

  DO WHILE ( itime < nitera_last )

    CALL TIME_RUNGEKUTTA(q,hq, s,hs, q_inf,s_inf, txc, wrk1d,wrk2d,wrk3d, l_q, l_hq, l_txc, l_comm)

    itime = itime + 1
    rtime = rtime + dtime

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,FilterDomainStep) == 0 ) THEN
      IF ( MOD(itime-nitera_first,nitera_stats)  == 0 ) THEN; flag_save = .TRUE.
      ELSE;                                                     flag_save = .FALSE.
      ENDIF
      CALL DNS_FILTER(flag_save, q,s, txc, wrk1d,wrk2d,wrk3d)
    ENDIF

    ! This should be integrated into the inflow buffer, as the filter contribution
    IF ( MOD(itime-nitera_first,FilterInflowStep) == 0 ) THEN ! Inflow filter in spatial mode
      CALL BOUNDARY_INFLOW_FILTER(BcsFlowImin%ref, BcsScalImin%ref, q,s, txc, wrk1d,wrk2d,wrk3d)
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( iviscchg == 1 ) THEN ! Change viscosity if necessary
      visc = visc - dtime*visctime
      IF ( ( (visc < viscstop) .AND. (viscstart > viscstop) )   .OR. &
          ( (visc > viscstop) .AND. (viscstart < viscstop) ) ) THEN
        iviscchg = 0; visc = viscstop
      ENDIF
    ENDIF

    ! -----------------------------------------------------------------------
    CALL TIME_COURANT(q,s, wrk3d)

    ! ###################################################################
    ! The rest: Logging, postprocessing and saving
    ! ###################################################################
    IF ( MOD(itime-nitera_first,nitera_log) == 0 .OR. INT(logs_data(1)) /= 0 ) THEN ! Log files
      CALL DNS_LOGS(i2)
#ifdef LES
      IF ( iles == 1 ) CALL LES_LOGS(i2)
#endif
    ENDIF

    ! -----------------------------------------------------------------------
    ! Accumulate statistics in spatially evolving cases
    IF ( MOD(itime-nitera_first,nitera_stats_spa) == 0 ) THEN
      CALL AVG_FLOW_ZT_REDUCE(q, hq,txc, mean_flow, wrk2d,wrk3d)
      IF ( icalc_scal == 1 ) THEN
        CALL AVG_SCAL_ZT_REDUCE(q,s, hq,txc, mean_scal, wrk2d,wrk3d)
      ENDIF
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( tower_mode == 1 ) THEN
      CALL DNS_TOWER_ACCUMULATE(q,1,wrk1d)
      CALL DNS_TOWER_ACCUMULATE(s,2,wrk1d)
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( itrajectory /= LAG_TRAJECTORY_NONE ) THEN
      CALL PARTICLE_TRAJECTORIES_ACCUMULATE(q,s, txc, l_g,l_q,l_hq,l_txc,l_comm, wrk2d,wrk3d)
    END IF

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,nitera_stats) == 0 ) THEN ! Calculate statistics
      IF     ( imode_sim == DNS_MODE_TEMPORAL ) THEN
        CALL STATISTICS_TEMPORAL_LAYER(q,s,hq, txc, wrk1d,wrk2d,wrk3d)
        IF ( icalc_part == 1 ) THEN
          CALL STATS_TEMPORAL_LAGRANGIAN(q,s,hq, l_q,l_txc,l_comm, txc, mean, wrk1d,wrk2d,wrk3d)
        ENDIF
      ELSE IF ( imode_sim == DNS_MODE_SPATIAL ) THEN
        CALL STATISTICS_SPATIAL_LAYER(txc, wrk1d,wrk2d)
      ENDIF
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,nitera_save) == 0 .OR. &      ! Save restart files
        itime == nitera_last .OR. INT(logs_data(1)) /= 0 ) THEN ! Secure that one restart file is saved

      IF ( icalc_flow == 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, inb_flow, isize_field, q, wrk3d)
      ENDIF
      IF ( icalc_scal == 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, isize_field, s, wrk3d)
      ENDIF

      IF ( tower_mode == 1 ) THEN
        CALL DNS_TOWER_WRITE(wrk3d)
      ENDIF

      IF ( icalc_part == 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
        CALL IO_WRITE_PARTICLE(fname, l_g, l_q)
        IF ( itrajectory /= LAG_TRAJECTORY_NONE ) THEN
          WRITE(fname,*) itime; fname =TRIM(ADJUSTL(tag_traj))//TRIM(ADJUSTL(fname))
          CALL PARTICLE_TRAJECTORIES_WRITE(fname)
        END IF
      END IF

      IF ( imode_sim == DNS_MODE_SPATIAL .AND. nitera_stats_spa > 0 ) THEN ! Spatial; running averages
        WRITE(fname,*) itime; fname='st'//TRIM(ADJUSTL(fname))
        CALL IO_WRITE_AVG_SPATIAL(fname, mean_flow, mean_scal)
      ENDIF

    ENDIF

    ! -----------------------------------------------------------------------
    IF ( MOD(itime-nitera_first,nitera_pln) == 0 ) THEN ! Save planes
      CALL PLANES_SAVE( q,s, hq,txc, wrk1d,wrk2d,wrk3d )
    ENDIF

    ! -----------------------------------------------------------------------
    IF ( INT(logs_data(1)) /= 0 ) CALL DNS_STOP(INT(logs_data(1)))

  ENDDO

  RETURN
END SUBROUTINE TIME_INTEGRATION
