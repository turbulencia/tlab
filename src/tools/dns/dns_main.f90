#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "DNS"

PROGRAM DNS

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE TLAB_ARRAYS
  USE LAGRANGE_GLOBAL
  USE LAGRANGE_ARRAYS
  USE DNS_LOCAL
  USE DNS_ARRAYS
  USE TIME
  USE DNS_TOWER
  USE PLANES
  USE BOUNDARY_INFLOW
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS
  USE STATISTICS
  USE PARTICLE_TRAJECTORIES

  IMPLICIT NONE
  SAVE

#include "integers.h"

  ! -------------------------------------------------------------------
  CHARACTER*32 fname
  CHARACTER*128 str
  TINTEGER ig

  ! ###################################################################
  CALL DNS_START()

  CALL DNS_READ_GLOBAL(ifile)
  IF ( icalc_part == 1 ) THEN
    CALL PARTICLE_READ_GLOBAL(ifile)
  END IF
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(ifile)
#endif
  CALL DNS_READ_LOCAL(ifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
  IF ( imode_rhs == EQNS_RHS_NONBLOCKING ) CALL DNS_NB3DFFT_INITIALIZE
#endif

  ! #######################################################################
  ! Allocating memory space
  ! #######################################################################
  CALL TLAB_ALLOCATE(C_FILE_LOC)

  CALL PARTICLE_ALLOCATE(C_FILE_LOC)

  CALL DNS_ALLOCATE()

  CALL STATISTICS_INITIALIZE()

  CALL PLANES_INITIALIZE()

  IF ( tower_mode == 1 ) THEN
    CALL DNS_TOWER_INITIALIZE(tower_stride)
  END IF

  ! ###################################################################
  ! Read the grid
  ! ###################################################################
#include "dns_read_grid.h"

  ! ###################################################################
  ! Initialize operators and reference data
  ! ###################################################################
  DO ig = 1,3
    CALL OPR_FILTER_INITIALIZE( g(ig), FilterDomain(ig), wrk1d )
  END DO

  IF ( ifourier == 1 ) THEN
    CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  END IF

  CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)

  CALL FI_PROFILES_INITIALIZE(wrk1d)

  ! ###################################################################
  ! Read fields
  ! ###################################################################
  itime = nitera_first

  visc_stop  = visc ! Value read in ifile

  IF ( icalc_scal == 1 ) THEN
    WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
    CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
  END IF

  WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
  CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, inb_flow, i0, isize_wrk3d, q, wrk3d)

  IF ( icalc_part == 1 ) THEN
    WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
    CALL IO_READ_PARTICLE(fname, l_g, l_q)
    IF ( itrajectory /= LAG_TRAJECTORY_NONE ) THEN
      CALL PARTICLE_TRAJECTORIES_INITIALIZE(nitera_save, nitera_last)
    END IF
  END IF

  IF ( imode_sim == DNS_MODE_SPATIAL .AND. nitera_stats_spa > 0 ) THEN
    WRITE(fname,*) nitera_first; fname = 'st'//TRIM(ADJUSTL(fname))
    CALL IO_READ_AVG_SPATIAL(fname, mean_flow, mean_scal)
  END IF

  ! ###################################################################
  ! Initialize diagnostic thermodynamic quantities
  ! ###################################################################
  CALL FI_DIAGNOSTIC( imax,jmax,kmax, q,s, wrk3d )

  ! ###################################################################
  ! Initialize change in viscosity
  ! ###################################################################
  flag_viscosity = .FALSE.
  IF ( visc /= visc_stop ) THEN
    WRITE(str,*) visc
    CALL IO_WRITE_ASCII(lfile,'Changing original viscosity '//TRIM(ADJUSTL(str))//' to new value.')
    IF ( visc_time > C_0_R ) THEN
      visc_rate = ( visc_stop -visc ) /visc_time
      visc_time = rtime +visc_time                 ! Stop when this time is reached
      flag_viscosity = .TRUE.
    ELSE
      visc = visc_stop
    END IF
  END IF

  ! ###################################################################
  ! Check
  ! ###################################################################
  logs_data(1) = 0 ! Status
  CALL DNS_CONTROL(i0, q,s, txc, wrk2d,wrk3d)

  ! ###################################################################
  ! Initialize particle simumulation
  ! ###################################################################
  IF ( icalc_part == 1 ) THEN
    CALL PARTICLE_INITIALIZE()
  END IF

  ! ###################################################################
  ! Initialize data for boundary conditions
  ! ###################################################################
  CALL BOUNDARY_BUFFER_INITIALIZE(q,s, txc, wrk3d)

  CALL BOUNDARY_BCS_INITIALIZE(wrk3d)

  IF ( imode_sim == DNS_MODE_SPATIAL ) THEN
    CALL BOUNDARY_INFLOW_INITIALIZE(rtime, txc, wrk1d,wrk2d,wrk3d)
  END IF

  ! ###################################################################
  ! Initialize time marching scheme
  ! ###################################################################
  CALL TIME_INITIALIZE()
  CALL TIME_COURANT(q, wrk3d)

  ! ###################################################################
  ! Initialize logfiles
  ! ###################################################################
  CALL DNS_LOGS(i1) ! headers
  CALL DNS_LOGS(i2) ! first line

  ! ###################################################################
  ! Do simulation: Integrate equations
  ! ###################################################################
  itime = nitera_first

  CALL TIME_INTEGRATION()

  ! ###################################################################
  CALL DNS_STOP(INT(logs_data(1)))

CONTAINS
  ! #######################################################################
  ! #######################################################################
  SUBROUTINE DNS_ALLOCATE()
  IMPLICIT NONE

  TINTEGER ierr
  CHARACTER*128 line

  ! ###################################################################
  WRITE(str,*) inb_flow; line = 'Allocating array rhs flow of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(hq(isize_field,inb_flow),    stat=ierr)
  IF ( ierr /= 0 ) THEN
    CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_q.')
    CALL DNS_STOP(DNS_ERROR_ALLOC)
  END IF

  WRITE(str,*) inb_scal; line = 'Allocating array rhs scal of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(hs(isize_field,inb_scal),    stat=ierr)
  IF ( ierr /= 0 ) THEN
    CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_s.')
    CALL DNS_STOP(DNS_ERROR_ALLOC)
  END IF

  ! -------------------------------------------------------------------
  IF ( icalc_part == 1 ) THEN
    WRITE(str,*) isize_l_comm; line = 'Allocating array l_comm of size '//TRIM(ADJUSTL(str))
    CALL IO_WRITE_ASCII(lfile,line)
    ALLOCATE(l_comm(isize_l_comm), stat=ierr)
    IF ( ierr /= 0 ) THEN
      CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_comm.')
      CALL DNS_STOP(DNS_ERROR_ALLOC)
    END IF

    WRITE(str,*) isize_particle; line = 'Allocating array l_hq of size '//TRIM(ADJUSTL(str))//'x'
    WRITE(str,*) inb_part; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    CALL IO_WRITE_ASCII(lfile,line)
    ALLOCATE(l_hq(isize_particle,inb_part),stat=ierr)
    IF ( ierr /= 0 ) THEN
      CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_hq.')
      CALL DNS_STOP(DNS_ERROR_ALLOC)
    END IF

  END IF

  RETURN
  END SUBROUTINE DNS_ALLOCATE
END PROGRAM DNS
