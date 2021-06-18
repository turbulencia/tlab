#include "types.h"
#include "dns_const.h"

MODULE DNS_CONSTANTS
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MajorVersion = 7
  TINTEGER, PARAMETER :: MinorVersion = 0

  TINTEGER, PARAMETER :: MAX_VARS = 20
  TINTEGER, PARAMETER :: MAX_PROF = 10
  TINTEGER, PARAMETER :: MAX_JETS =  5

  TINTEGER, PARAMETER :: MAX_NSP = 10 ! Species in the mixture

  TINTEGER, PARAMETER :: MAX_AVG_TEMPORAL  = 230
  TINTEGER, PARAMETER :: MAX_STATS_SPATIAL = 100 ! Running statistics

  CHARACTER*32, PARAMETER :: gfile = 'grid'
  CHARACTER*32, PARAMETER :: ofile = 'dns.out'
  CHARACTER*32, PARAMETER :: lfile = 'dns.log'
  CHARACTER*32, PARAMETER :: efile = 'dns.err'
  CHARACTER*32, PARAMETER :: wfile = 'dns.war'
  CHARACTER*32, PARAMETER :: tfile = 'dns.trc'

  CHARACTER*32, PARAMETER :: tag_flow ='flow.'
  CHARACTER*32, PARAMETER :: tag_scal ='scal.'
  CHARACTER*32, PARAMETER :: tag_part ='part.'
  CHARACTER*32, PARAMETER :: tag_traj ='traj.'

END MODULE DNS_CONSTANTS

MODULE DNS_GLOBAL
  USE DNS_TYPES,     ONLY : grid_dt, filter_dt, subarray_dt, term_dt, background_dt
  USE DNS_CONSTANTS, ONLY : MAX_VARS, MAX_NSP
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
  TINTEGER :: icalc_flow, icalc_scal, icalc_part
  TINTEGER :: imode_sim                ! type of simulation (spatial, temporal)
  TINTEGER :: imode_files              ! files format
  TINTEGER :: imode_verbosity          ! level of verbosity used in log files
  TINTEGER :: imode_eqns               ! set of equations to be solved
  TINTEGER :: iadvection, iviscous, idiffusion ! formulation
  TINTEGER :: ifourier
  TINTEGER :: itransport

  TINTEGER :: imode_fdm                ! finite-difference method for spatial operators

  TINTEGER :: imode_ibm                ! IBM mode

! ###################################################################
! Iteration
! ###################################################################
  TINTEGER :: itime                    ! iteration number
  TREAL    :: rtime                    ! physical time

! ###################################################################
! Arrays size
! ###################################################################
! fields
  TINTEGER :: imax,jmax,kmax, isize_field ! locally per processor
  TINTEGER :: inb_flow, inb_flow_array    ! transported & array space
  TINTEGER :: inb_scal, inb_scal_array    ! transported & array space

! auxiliary arrays
  TINTEGER :: isize_wrk1d,     inb_wrk1d  ! 1D arrays
  TINTEGER :: isize_wrk2d,     inb_wrk2d  ! 2D arrays
  TINTEGER :: isize_txc_field, inb_txc    ! 3D arrays for intermediate calculations
  TINTEGER :: isize_txc_dimx, isize_txc_dimz  ! partition for MPI data transposition

! Particle arrays
  TINTEGER :: isize_particle              ! max number of particles per processor
  TINTEGER :: inb_part, inb_part_array
  TINTEGER :: inb_part_txc

! subarray information (offset)
  TYPE(subarray_dt), DIMENSION(IO_SUBARRAY_SIZE) :: io_aux

! ###################################################################
  TYPE(grid_dt), DIMENSION(3) :: g     ! Grid information along 3 directions
  TREAL :: area,volume                 ! Horizontal area and volume

! ###################################################################
  TYPE(background_dt) :: qbg(3)        ! Velocity background
  TYPE(background_dt) :: sbg(MAX_NSP)  ! Scalars backgrounds
  TYPE(background_dt) :: pbg, rbg, tbg ! Pressure, density, temperature backgrounds

  TREAL, DIMENSION(:), ALLOCATABLE :: pbackground, tbackground, rbackground, ribackground
  TREAL, DIMENSION(:), ALLOCATABLE :: bbackground, epbackground

! ###################################################################
  TYPE(term_dt) :: buoyancy   ! Buoyancy parameters
  TYPE(term_dt) :: coriolis   ! Coriolis parameters
  TYPE(term_dt) :: radiation  ! Radiation parameters
  TYPE(term_dt) :: transport  ! Transport parameters
  TYPE(term_dt) :: chemistry  ! Chemistry parameters
  TYPE(term_dt) :: subsidence ! Large-scale parameters

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

  TREAL    :: visc                                ! 1/reynolds

! ###########################################################
! Filters
! ###########################################################
  TYPE(filter_dt), DIMENSION(3)        :: FilterDomain
  LOGICAL,         DIMENSION(MAX_VARS) :: FilterDomainActive
  TINTEGER,        DIMENSION(MAX_VARS) :: FilterDomainBcsFlow, FilterDomainBcsScal

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
  TINTEGER :: nstatavg, statavg(MAX_STATS_SPATIAL), & ! Ox planes at which to accumulate statistics
              nstatavg_points, &                      ! number of accumulated points
              istattimeorg, &                         ! time at which accumulation of statistics started
              istat_min_ver, istat_maj_ver
  TREAL    :: rstattimeorg

END MODULE DNS_GLOBAL
