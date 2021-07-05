#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "avgij_map.h"

#define C_FILE_LOC "DNS"

PROGRAM DNS

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE LAGRANGE_GLOBAL
  USE DNS_LOCAL
  USE TIME
  USE DNS_TOWER
  USE PLANES
  USE BOUNDARY_INFLOW
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS
  USE STATISTICS
  USE PARTICLE_TRAJECTORIES
#ifdef LES
  USE LES_GLOBAL
#endif

  IMPLICIT NONE
  SAVE

#include "integers.h"

  ! -------------------------------------------------------------------
  TREAL, DIMENSION(:,:), ALLOCATABLE :: x,y,z     ! Grid and associated arrays

  TREAL, DIMENSION(:,:), ALLOCATABLE :: q,s       ! Eulerian fields
  TREAL, DIMENSION(:,:), ALLOCATABLE :: txc       ! Temporary space for Eulerian fields

  TREAL, DIMENSION(:,:), ALLOCATABLE :: l_q       ! Lagrangian fields
  TREAL, DIMENSION(:,:), ALLOCATABLE :: l_txc     ! Temporary space for Lagrangian fields

  TREAL, DIMENSION(:),   ALLOCATABLE :: wrk1d,wrk2d,wrk3d   ! Work arrays (scratch space)

  TREAL, DIMENSION(:,:), ALLOCATABLE :: h_q,h_s   ! Right-hand sides Eulerian fields
  TREAL, DIMENSION(:,:), ALLOCATABLE :: l_hq      ! Right-hand sides Lagrangian fields
  TREAL, DIMENSION(:),   ALLOCATABLE :: l_comm    ! Communication space for Lagrangian fields

  TARGET q, x,y,z

  ! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: e, rho, p, T

  CHARACTER*32 fname
  CHARACTER*128 str, line
  TINTEGER idummy, ig, ierr

  ! ###################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(ifile)
  IF ( icalc_part == 1 ) THEN
    CALL PARTICLE_READ_GLOBAL(ifile)
  ENDIF
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(ifile)
#endif
#ifdef LES
  CALL LES_READ_INI(ifile)
#endif
  CALL DNS_READ_LOCAL(ifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
  IF ( imode_rhs == EQNS_RHS_NONBLOCKING ) CALL DNS_NB3DFFT_INITIALIZE
#endif

  ! #######################################################################
  ! Memory management
  ! #######################################################################
  ! txc
  inb_txc = 9
  IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
    inb_txc = 6
    IF ( rkm_mode == RKM_IMP3_DIFFUSION ) inb_txc = inb_txc+1
  ELSE IF ( imode_eqns == DNS_EQNS_INTERNAL       .AND. &
      iadvection == EQNS_SKEWSYMMETRIC      .AND. &
      iviscous   == EQNS_EXPLICIT                 ) THEN
    inb_txc = 6
  ENDIF
  IF ( imixture == MIXT_TYPE_AIRWATER .AND. damkohler(3) > C_0_R ) inb_txc = inb_txc + 1

  IF ( imode_sim == DNS_MODE_SPATIAL ) THEN ! because of the statistics
    inb_txc = MAX(inb_txc,7)

    IF ( stats_averages ) THEN
      idummy = MAX(MA_MOMENTUM_SIZE,MS_SCALAR_SIZE)
      IF ( MOD( nstatavg*jmax*idummy , isize_txc_field ) > 0 ) THEN
        idummy = nstatavg*jmax*idummy / isize_txc_field + 1
      ELSE
        idummy = nstatavg*jmax*idummy / isize_txc_field
      ENDIF
      inb_txc = MAX(inb_txc,idummy)
    ENDIF

  ENDIF

#ifdef USE_PSFFT
  IF ( imode_rhs == EQNS_RHS_NONBLOCKING ) inb_txc = MAX(inb_txc,15)
#endif

  ! wkr3d
  isize_wrk3d = MAX(imax,g_inf(1)%size)*MAX(jmax,g_inf(2)%size)*kmax
  isize_wrk3d = MAX(isize_wrk3d,isize_txc_field)
  IF ( icalc_part == 1) THEN
    isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
    isize_wrk3d = MAX(isize_wrk3d,(jmax*(kmax+1)*inb_particle_interp*2))
    isize_wrk3d = MAX(isize_wrk3d,(jmax*(imax+1)*inb_particle_interp*2))
  END IF
  IF ( tower_mode == 1 ) THEN
    isize_wrk3d = MAX(isize_wrk3d,nitera_save*(g(2)%size+2))
  ENDIF

#ifdef LES
#ifdef USE_MPI
  IF ( ims_npro > 1 ) THEN ! wrk3d for OZ filter in PARALLEL mode may require more space in LES
    isize_wrk3d = MAX(isize_wrk3d,imax*jmax*(isgs_f0size + isgs_f1size + kmax))
  ENDIF
#endif
#endif

  ! -------------------------------------------------------------------
  ! Allocating basic memory space
  ! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

#include "dns_alloc_arrays.h"
  IF ( icalc_part == 1 ) THEN
#include "dns_alloc_larrays.h"
  ENDIF

  ! -------------------------------------------------------------------
  ! Allocating other stuff
  ! -------------------------------------------------------------------
  CALL DNS_ALLOCATE()

  CALL STATISTICS_INITIALIZE()
  CALL PLANES_INITIALIZE()
  IF ( tower_mode == 1 ) THEN
    CALL DNS_TOWER_INITIALIZE(tower_stride)
  ENDIF

  ! ###################################################################
  ! Read the grid
  ! ###################################################################
#include "dns_read_grid.h"

  ! ###################################################################
  ! Initialize operators and reference data
  ! ###################################################################
  ! Filters
  DO ig = 1,3
    CALL OPR_FILTER_INITIALIZE( g(ig), FilterDomain(ig), wrk1d )
  END DO

  ! Spectral Poisson solver
  IF ( ifourier == 1 ) THEN
    CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  ! Check operators
  CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)

  ! Thermodynamic quantities
  CALL FI_PROFILES_INITIALIZE(wrk1d)

  ! ###################################################################
  ! Read fields
  ! ###################################################################
  itime = nitera_first

  visc_stop  = visc ! Value read in ifile

  IF ( icalc_scal == 1 ) THEN
    WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
    CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
  ENDIF

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
  ENDIF

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
    ENDIF
  END IF

  ! ###################################################################
  ! Initialize Lagrangian stuff
  ! ###################################################################
  IF ( icalc_part == 1 ) THEN
    ! set boundarys for residence time pdf
    IF ( ilagrange == LAG_TYPE_BIL_CLOUD_4 ) THEN
      l_y_lambda =  (g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(1)%ymean - C_2_R
      l_y_base =   ((g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(1)%ymean -(g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(3)%ymean )/C_2_R &
          +  (g(2)%nodes(jmax)-g(2)%nodes(1)) *sbg(3)%ymean
      IF (residence_reset == 1) THEN
        l_q(:,6:7) = C_0_R
      ENDIF
    ENDIF
  END IF

  ! ###################################################################
  ! Initialize thermodynamic quantities
  ! ###################################################################
  IF ( imode_eqns == DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns == DNS_EQNS_ANELASTIC ) THEN
    IF      ( imixture == MIXT_TYPE_AIRWATER .AND. damkohler(3) <= C_0_R ) THEN ! Calculate q_l
      CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2), s(1,1), epbackground,pbackground)

    ELSE IF ( imixture == MIXT_TYPE_AIRWATER_LINEAR                        ) THEN
      CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))

    ENDIF

  ELSE
    e   => q(:,4)
    rho => q(:,5)
    p   => q(:,6)
    T   => q(:,7)

    CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
    CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
    IF ( itransport == EQNS_TRANS_SUTHERLAND .OR. itransport == EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, q(:,8))

#ifdef CHEMISTRY
    IF ( ireactive /= CHEM_NONE ) THEN ! Calculate TGFM if reactive case
      CALL THERMO_GAMMA(imax,jmax,kmax, s, T, wrk3d)
      CALL CHEM_BURKESCHUMANN(imax,jmax,kmax, s, T, wrk3d, h_q(1,1))
    ENDIF
#endif

  ENDIF

  ! ###################################################################
  ! Check
  ! ###################################################################
  logs_data(1) = 0 ! Status
  CALL DNS_CONTROL(i0, q,s, txc, wrk2d,wrk3d)

  ! ###################################################################
  ! Initialize data for boundary conditions
  ! ###################################################################
  CALL BOUNDARY_BUFFER_INITIALIZE(q,s, txc, wrk3d)

  CALL BOUNDARY_BCS_INITIALIZE(wrk3d)

  IF ( imode_sim == DNS_MODE_SPATIAL ) THEN
    CALL BOUNDARY_INFLOW_INITIALIZE(rtime, txc, wrk1d,wrk2d,wrk3d)
  ENDIF

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

  CALL TIME_INTEGRATION(q,h_q, s,h_s, txc, wrk1d,wrk2d,wrk3d, l_q, l_hq, l_txc, l_comm)

  ! ###################################################################
#ifdef USE_FFTW
  IF ( ifourier == 1 ) THEN
    CALL dfftw_destroy_plan(fft_plan_fx)
    CALL dfftw_destroy_plan(fft_plan_bx)
    IF ( g(3)%size > 1 ) THEN
      CALL dfftw_destroy_plan(fft_plan_fz)
      CALL dfftw_destroy_plan(fft_plan_bz)
    ENDIF
  ENDIF
#endif

  CALL DNS_STOP(INT(logs_data(1)))

CONTAINS
  SUBROUTINE DNS_ALLOCATE()
  IMPLICIT NONE

  ! ###################################################################
  WRITE(str,*) inb_flow; line = 'Allocating array rhs flow of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(h_q(isize_field,inb_flow),    stat=ierr)
  IF ( ierr /= 0 ) THEN
    CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_q.')
    CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  WRITE(str,*) inb_scal; line = 'Allocating array rhs scal of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(h_s(isize_field,inb_scal),    stat=ierr)
  IF ( ierr /= 0 ) THEN
    CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for h_s.')
    CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  ! -------------------------------------------------------------------
  IF ( icalc_part == 1 ) THEN
    WRITE(str,*) isize_l_comm; line = 'Allocating array l_comm of size '//TRIM(ADJUSTL(str))
    CALL IO_WRITE_ASCII(lfile,line)
    ALLOCATE(l_comm(isize_l_comm), stat=ierr)
    IF ( ierr /= 0 ) THEN
      CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_comm.')
      CALL DNS_STOP(DNS_ERROR_ALLOC)
    ENDIF

    WRITE(str,*) isize_particle; line = 'Allocating array l_hq of size '//TRIM(ADJUSTL(str))//'x'
    WRITE(str,*) inb_part; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    CALL IO_WRITE_ASCII(lfile,line)
    ALLOCATE(l_hq(isize_particle,inb_part),stat=ierr)
    IF ( ierr /= 0 ) THEN
      CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_hq.')
      CALL DNS_STOP(DNS_ERROR_ALLOC)
    ENDIF

  ENDIF

  RETURN
  END SUBROUTINE DNS_ALLOCATE
END PROGRAM DNS
