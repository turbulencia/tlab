#include "types.h"
#include "dns_const.h"

MODULE TLAB_VARS
  USE TLAB_TYPES,     ONLY : grid_dt, filter_dt, subarray_dt, term_dt, background_dt
  USE TLAB_CONSTANTS, ONLY : MAX_VARS, MAX_NSP
  USE TLAB_CONSTANTS, ONLY : MAX_STATS_SPATIAL
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
  TINTEGER :: imode_precision_files
  TINTEGER :: imode_verbosity          ! level of verbosity used in log files
  TINTEGER :: imode_eqns               ! set of equations to be solved
  TINTEGER :: iadvection, iviscous, idiffusion ! formulation
  TINTEGER :: ifourier
  TINTEGER :: itransport
  TINTEGER :: istagger, ivfilter       ! horizontal staggering of pressure 
                                       ! vertical   filtering  of pressure
  TREAL    :: vfilter_param            ! vertical filter parameter

  TINTEGER :: imode_fdm                ! finite-difference method for spatial operators

! ###################################################################
! Iteration
! ###################################################################
  TINTEGER :: itime                    ! iteration number
  TREAL    :: rtime                    ! physical time

! ###################################################################
! Arrays size
! ###################################################################
! fields
  TINTEGER :: imax,jmax,kmax, isize_field     ! locally per processor
  TINTEGER :: inb_flow, inb_flow_array        ! transported & array space
  TINTEGER :: inb_scal, inb_scal_array        ! transported & array space

! auxiliary arrays
  TINTEGER :: isize_wrk1d,     inb_wrk1d      ! 1D arrays
  TINTEGER :: isize_wrk2d,     inb_wrk2d      ! 2D arrays
  TINTEGER :: isize_wrk3d                     ! 2D arrays
  TINTEGER :: isize_txc_field, inb_txc        ! 3D arrays for intermediate calculations
  TINTEGER :: isize_txc_dimx, isize_txc_dimz  ! partition for MPI data transposition

! Particle arrays
  TINTEGER :: isize_particle                  ! max number of particles per processor
  TINTEGER :: inb_part, inb_part_array
  TINTEGER :: inb_part_txc

! subarray information (offset)
  TYPE(subarray_dt), DIMENSION(IO_SUBARRAY_SIZE) :: io_aux

! ###################################################################
  TYPE(grid_dt), DIMENSION(3) :: g      ! Grid information along 3 directions
  TREAL :: area                         ! Horizontal area and volume

! ###################################################################
  TYPE(background_dt) :: qbg(3)         ! Velocity background
  TYPE(background_dt) :: sbg(MAX_NSP)   ! Scalars backgrounds
  TYPE(background_dt) :: pbg, rbg, tbg  ! Pressure, density, temperature backgrounds

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

! ###################################################################
! Compact parameters (1st derivative of 6th-order pentadiagonal)
! ###################################################################
  TREAL    :: C1N6M_ALPHA,  C1N6M_BETA
  TREAL    :: C1N6M_ALPHA2, C1N6M_BETA2
  TREAL    :: C1N6M_A,   C1N6M_B,   C1N6M_C
  TREAL    :: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
  TREAL    ::            C1N6M_BD2, C1N6M_CD3

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
              istattimeorg                            ! time at which accumulation of statistics started
  TREAL    :: rstattimeorg

END MODULE TLAB_VARS
