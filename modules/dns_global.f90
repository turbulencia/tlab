#include "types.h"

MODULE DNS_CONSTANTS
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_VARS = 50
  TINTEGER, PARAMETER :: MAX_PROF = 10
  TINTEGER, PARAMETER :: MAX_JETS =  5

  TINTEGER, PARAMETER :: MAX_NSP = 10 ! Species in the mixture

!  TINTEGER, PARAMETER :: MAX_AVG_TEMPORAL  = 115
  TINTEGER, PARAMETER :: MAX_AVG_TEMPORAL  = 230
  TINTEGER, PARAMETER :: MAX_AVG_SPATIAL   = 228
  TINTEGER, PARAMETER :: MAX_STATS_SPATIAL = 100 ! Running statistics

  CHARACTER*32, PARAMETER :: gfile = 'grid'
  CHARACTER*32, PARAMETER :: ofile = 'dns.out'
  CHARACTER*32, PARAMETER :: lfile = 'dns.log'
  CHARACTER*32, PARAMETER :: efile = 'dns.err'
  CHARACTER*32, PARAMETER :: wfile = 'dns.war'
  CHARACTER*32, PARAMETER :: tfile = 'dns.trc'

  CHARACTER*32, PARAMETER :: tag_flow ='flow.'
  CHARACTER*32, PARAMETER :: tag_scal ='scal.' 
  
END MODULE DNS_CONSTANTS

MODULE DNS_GLOBAL
  USE DNS_TYPES,     ONLY : grid_structure, subarray_structure, term_structure
  USE DNS_CONSTANTS, ONLY : MAX_VARS, MAX_PROF, MAX_JETS, MAX_NSP
  USE DNS_CONSTANTS, ONLY : MAX_STATS_SPATIAL
  IMPLICIT NONE
  SAVE

! ###################################################################
! OpenMP
! ###################################################################
  TINTEGER :: dns_omp_numThreads 
  TINTEGER :: dns_omp_error

! ###################################################################
! General options
! ###################################################################
  TINTEGER :: icalc_flow, icalc_scal, icalc_particle
  TINTEGER :: imode_sim                ! type of simulation (spatial, temporal)
  TINTEGER :: imode_flow               ! type of geometry
  TINTEGER :: imode_files              ! files format
  TINTEGER :: imode_verbosity          ! level of verbosity used in log files
  TINTEGER :: imode_eqns               ! set of equations to be solved 
  TINTEGER :: iadvection, iviscous, idiffusion, icoriolis, ibodyforce ! formulation
  TINTEGER :: ifourier
  TINTEGER :: itransport, ireactive

  TINTEGER :: imode_fdm                ! finite-difference method for spatial operators

! ###################################################################
! Iteration
! ###################################################################
  TINTEGER :: itime                    ! iteration number
  TREAL    :: rtime                    ! physical time

! ###################################################################
! Arrays size
! ###################################################################
! grid
  TINTEGER :: imax_total,jmax_total,kmax_total
  TINTEGER :: inb_grid, inb_grid_1, inb_grid_2, inb_grid_3
  TINTEGER :: imax,jmax,kmax              ! locally per processor

! fields
  TINTEGER :: isize_field
  TINTEGER :: inb_flow, inb_flow_array    ! transported & array space
  TINTEGER :: inb_scal, inb_scal_array    ! transported & array space
  TINTEGER :: inb_vars                    ! simply inb_flow + inb_scal

! auxiliary arrays
  TINTEGER :: isize_wrk1d,     inb_wrk1d  ! 1D arrays
  TINTEGER :: isize_wrk2d,     inb_wrk2d  ! 2D arrays
  TINTEGER :: isize_txc_field, inb_txc    ! 3D arrays for intermediate calculations
  TINTEGER :: isize_txc                       ! total space, >= isize_txc_field*inb_txc
  TINTEGER :: isize_txc_dimx, isize_txc_dimz  ! partition for MPI data transposition

! Particle arrays
  TINTEGER :: isize_particle           ! particle number per processor
  TINTEGER :: inb_particle             ! number of particle properties
  TINTEGER :: inb_particle_txc         ! particle arrays for intermediate calculations

! temporal statistics
  TINTEGER :: inb_mean_temporal, inb_mean_spatial

! spatial statistics
  TINTEGER :: nspa_rest, nspa_step

! subarray information (offset)
  TYPE(subarray_structure), DIMENSION(10) :: io_aux

! ###################################################################
! Grid
! ###################################################################
  TYPE(grid_structure), DIMENSION(3) :: grid
  TINTEGER :: iunifx,iunify,iunifz               ! uniform
  TINTEGER :: i1bc,j1bc,k1bc                     ! biased
  TREAL    :: scalex,scaley,scalez

  TREAL    :: area,volume
  
! ###################################################################
! Profiles
! ###################################################################
! Flow: velocities
  TINTEGER :: iprof_u
  TREAL    :: mean_u, delta_u, thick_u,   & ! Velocity profile parameters
              ycoor_u,                    & ! Relative reference position 
              prof_u(MAX_PROF)        ,   & ! Further parameters
              diam_u, jet_u(MAX_JETS)       ! Jet-related geometry

  TREAL    :: mean_v, mean_w      ! in these two directions, only mean values

! Flow: thermodynamics
  TREAL    :: p_init              ! reference pressure

  TINTEGER :: iprof_rho
  TREAL    :: mean_rho, delta_rho, thick_rho, ycoor_rho, &
              prof_rho(MAX_PROF), diam_rho, jet_rho(MAX_JETS)

  TINTEGER :: iprof_tem
  TREAL    :: mean_tem, delta_tem, thick_tem, ycoor_tem, &
              prof_tem(MAX_PROF), diam_tem, jet_tem(MAX_JETS)

! Scalars
  TINTEGER :: iprof_i(MAX_NSP)
  TREAL    :: mean_i(MAX_NSP), delta_i(MAX_NSP), thick_i(MAX_NSP), ycoor_i(MAX_NSP), &
              prof_i(MAX_PROF,MAX_NSP), diam_i(MAX_NSP), jet_i(MAX_JETS,MAX_NSP)
  
! ###################################################################
! Body force vector and buoyancy function parameters
! ###################################################################
  TREAL    :: body_vector(3)                       ! vector
  TREAL    :: body_param(MAX_PROF)                 ! buoyancy function parameters
  TINTEGER :: ibodyforce_x,ibodyforce_y,ibodyforce_z
  TYPE(term_structure) :: body

! ###################################################################
! Rotation parameters
! ###################################################################
  TREAL    :: rotn_vector(3)
  TREAL    :: rotn_param(MAX_PROF)
  TINTEGER :: icoriolis_x,icoriolis_y,icoriolis_z
  TYPE(term_structure) :: coriolis

! ###################################################################
  TYPE(term_structure) :: radiation ! Radiation parameters
  TYPE(term_structure) :: transport ! Transport parameters
  TYPE(term_structure) :: chemistry ! Chemistry parameters

! ###################################################################
! Nondimensional numbers
! ###################################################################
  TREAL    :: reynolds, prandtl, schmidt(MAX_NSP) ! molecular transport
  TREAL    :: mach                                ! compressibility
  TREAL    :: damkohler(MAX_NSP)                  ! reaction
  TREAL    :: froude                              ! body force
  TREAL    :: rossby                              ! Coriolis force
  TREAL    :: stokes                              ! Stokes number of liquid particles
  TREAL    :: settling                            ! sedimentation parameter for liquid particle

  TREAL    :: visc                                ! simply 1/reynolds

! ###################################################################
! FFTW
! ###################################################################
  INTEGER(8) :: fft_plan_fx, fft_plan_bx, fft_plan_fx_bcs
  INTEGER(8) :: fft_plan_fy, fft_plan_by
  INTEGER(8) :: fft_plan_fz, fft_plan_bz

  TINTEGER :: fft_reordering 

! ###################################################################
! Jet Statistics
! ###################################################################
  TINTEGER :: nstatavg, statavg(MAX_STATS_SPATIAL), nstatavg_points, istattimeorg, &
              nstatlin, statlin_i(MAX_STATS_SPATIAL), statlin_j(MAX_STATS_SPATIAL), &
              nstatpln, statpln(MAX_STATS_SPATIAL), nstatplnextra, nstatplnvars, &
              iupdate_stat, istat_min_ver, istat_maj_ver
  TREAL    :: rstattimeorg

END MODULE DNS_GLOBAL
