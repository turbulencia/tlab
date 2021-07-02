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
  USE DNS_TOWER
  USE DNS_IBM
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

#include "integers.h"

  ! -------------------------------------------------------------------
  ! Grid and associated arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z

  ! Flow/Scalar variables and RHS space
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: q,s, h_q,h_s, txc

  ! Particle data
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q, l_hq, l_txc
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: l_comm

  ! Work arrays
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk1d,wrk2d,wrk3d

  ! Inflow arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x_inf, y_inf, z_inf
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE         :: q_inf, s_inf

  TARGET q

  ! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: e, rho, p, T

  CHARACTER*32 fname, inifile
  CHARACTER*128 str, line
  TINTEGER idummy, ig
  TINTEGER ierr, isize_wrk3d, isize_loc
  LOGICAL ibm_allocated

  ! ###################################################################
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_part == 1 ) THEN
    CALL PARTICLE_READ_GLOBAL(inifile)
  ENDIF
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(inifile)
#endif
#ifdef LES
  CALL LES_READ_INI(inifile)
#endif
  CALL DNS_READ_LOCAL(inifile)

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
  IF ( imode_rhs == EQNS_RHS_NONBLOCKING ) CALL DNS_NB3DFFT_INITIALIZE
#endif

  ! #######################################################################
  ! Memory management
  ! #######################################################################
  isize_loc = MAX(g_inf(1)%size,MAX(g_inf(2)%size,g_inf(3)%size))
  isize_wrk1d = MAX(isize_wrk1d,isize_loc)

  isize_loc = MAX(g_inf(1)%size*g_inf(2)%size,MAX(g_inf(1)%size*g_inf(3)%size,g_inf(2)%size*g_inf(3)%size))
  isize_wrk2d = MAX(isize_wrk2d, isize_loc)

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

#ifdef LES
  IF ( iles == 1 ) THEN ! this number needs to be revised
    isize_loc = 13

    IF ( iles_type_regu == LES_REGU_SMGDYN .OR. iles_type_regu == LES_REGU_SMGDYNRMS ) THEN
      IF ( iles_type_tran == LES_TRAN_NONE ) isize_loc = 10
      isize_loc = isize_loc + 4 ! space for aux_sg in dynamic smagorinsky
    ENDIF
    IF ( iles_type_chem == LES_CHEM_QUASIBS ) isize_loc = isize_loc + 4 ! space for chi

    idummy =  isize_wrk1d*isize_wrk1d*3 / isize_txc_field
    IF ( MOD( isize_wrk1d*isize_wrk1d*3 , isize_txc_field ) > 0 ) idummy = idummy+1
    isize_loc = MAX(isize_loc,idummy)

    inb_txc = MAX(isize_loc, inb_txc)

  ENDIF
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

  ALLOCATE(x_inf(g_inf(1)%size,g_inf(1)%inb_grid)) ! Inflow fields for spatial simulations
  ALLOCATE(y_inf(g_inf(2)%size,g_inf(2)%inb_grid))
  ALLOCATE(z_inf(g_inf(3)%size,g_inf(3)%inb_grid))
  ALLOCATE(q_inf(g_inf(1)%size*g_inf(2)%size*kmax,inb_flow_array))
  ALLOCATE(s_inf(g_inf(1)%size*g_inf(2)%size*kmax,inb_scal_array))

  ! Rhs
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

! IBM
  IF ( imode_ibm == 1 ) THEN
    ibm_allocated = .FALSE.
    CALL ALLOCATE_IBM(ibm_allocated)
  ENDIF

  ! Lagrangian part
  IF ( icalc_part == 1 ) THEN
#include "dns_alloc_larrays.h"

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

  ! -------------------------------------------------------------------
  ! Allocating other stuff
  ! -------------------------------------------------------------------
  CALL STATISTICS_INITIALIZE()
  CALL PLANES_INITIALIZE()
  IF ( tower_mode == 1 ) THEN
    CALL DNS_TOWER_INITIALIZE(tower_stride)
  ENDIF

  ! ###################################################################
  ! Read the grid
  ! ###################################################################
#include "dns_read_grid.h"

  IF ( g_inf(1)%size > 1 ) THEN ! Inflow fields for spatial simulations
    CALL IO_READ_GRID('grid.inf', g_inf(1)%size, g_inf(2)%size, g_inf(3)%size,  &
        g_inf(1)%scale,g_inf(2)%scale,g_inf(3)%scale, x_inf,y_inf,z_inf)
    CALL FDM_INITIALIZE(x_inf, g_inf(1), wrk1d)
    g_inf(2)%nodes => y_inf(:,1)
    g_inf(3)%nodes => z_inf(:,1)
  ENDIF

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

  visc_stop  = visc ! Value read in inifile

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
    CALL BOUNDARY_INFLOW_INITIALIZE(rtime, q_inf,s_inf, txc, wrk2d,wrk3d)
  ENDIF

  ! ###################################################################
  ! Initialize IBM
  ! ###################################################################
  IF ( imode_ibm == 1 ) THEN
    CALL INITIALIZE_GEOMETRY(txc, wrk3d)
    ! CALL INITIALIZE_IBM()
  ENDIF  

  ! ###################################################################
  ! Initialize LES
  ! ###################################################################
#ifdef LES
  IF ( iles == 1 ) CALL LES_INI(q,s,h_q,h_s, txc, vaux, wrk1d,wrk2d,wrk3d)
#endif

  ! ###################################################################
  ! Initialize time step dt
  ! ###################################################################
  CALL TIME_COURANT(q,s, wrk3d)

  ! ###################################################################
  ! Initialize logfiles
  ! ###################################################################
  CALL DNS_LOGS(i1) ! headers
  CALL DNS_LOGS(i2) ! first line

  ! ###################################################################
  ! Do simulation: Integrate equations
  ! ###################################################################
  itime = nitera_first

  CALL TIME_INTEGRATION(q,h_q, s,h_s, q_inf,s_inf, txc, wrk1d,wrk2d,wrk3d, &
      l_q, l_hq, l_txc, l_comm)

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
END PROGRAM DNS
