#include "dns_const.h"

module TLAB_VARS
    use TLAB_TYPES, only: grid_dt, filter_dt, subarray_dt, term_dt, profiles_dt
    use TLAB_CONSTANTS, only: MAX_VARS, MAX_NSP, wp, wi, sp
    use TLAB_CONSTANTS, only: MAX_STATS_SPATIAL
    implicit none
    save

! ###################################################################
! OpenMP
! ###################################################################
    integer :: dns_omp_numThreads
    integer :: dns_omp_error

! ###################################################################
! General options
! ###################################################################
    integer :: icalc_flow, icalc_scal
    integer :: imode_sim                ! type of simulation (spatial, temporal)
    integer :: imode_files              ! files format
    integer :: imode_precision_files    ! whether restart files in single or double precision
    integer :: imode_verbosity = 1      ! level of verbosity used in log files
    integer :: imode_eqns               ! set of equations to be solved
    integer :: iadvection, iviscous, idiffusion,  itransport ! formulation
    integer :: ifourier
    integer :: istagger, ivfilter       ! horizontal staggering of pressure

    real(wp) :: vfilter_param           ! vertical filter parameter of pressure

    integer :: imode_fdm                ! finite-difference method for spatial operators

    integer :: imode_ibm                ! IBM mode

! ###################################################################
! Iteration
! ###################################################################
    integer(wi) :: itime                    ! iteration number
    real(wp) :: rtime                       ! physical time

! ###################################################################
! Arrays size
! ###################################################################
! fields
    integer(wi) :: imax, jmax, kmax, isize_field    ! locally per processor
    integer(wi) :: inb_flow, inb_flow_array         ! transported & array space
    integer(wi) :: inb_scal, inb_scal_array         ! transported & array space

! auxiliary arrays
    integer(wi) :: isize_wrk1d, inb_wrk1d           ! 1D arrays
    integer(wi) :: isize_wrk2d, inb_wrk2d           ! 2D arrays
    integer(wi) :: isize_wrk3d                      ! 2D arrays
    integer(wi) :: isize_txc_field, inb_txc         ! 3D arrays for intermediate calculations
    integer(wi) :: isize_txc_dimx, isize_txc_dimz   ! partition for MPI data transposition

! subarray information (offset)
    type(subarray_dt), dimension(IO_SUBARRAY_SIZE) :: io_aux

! ###################################################################
    type(grid_dt), dimension(3) :: g            ! Grid information along 3 directions
    real(wp) :: area                            ! Horizontal area and volume

! ###################################################################
    type(profiles_dt) :: qbg(3)             ! Velocity background
    type(profiles_dt) :: sbg(MAX_NSP)       ! Scalars backgrounds
    type(profiles_dt) :: pbg, rbg, tbg, hbg ! Pressure, density, temperature, enthalpy backgrounds

    real(wp), dimension(:), allocatable :: pbackground, tbackground, rbackground, ribackground
    real(wp), dimension(:), allocatable :: bbackground, epbackground

! ###################################################################
    type(term_dt) :: buoyancy   ! Buoyancy parameters
    type(term_dt) :: coriolis   ! Coriolis parameters
    type(term_dt) :: radiation  ! Radiation parameters
    type(term_dt) :: transport  ! Transport parameters
    type(term_dt) :: chemistry  ! Chemistry parameters
    type(term_dt) :: subsidence ! Large-scale parameters
    type(term_dt) :: random     ! Random Forcing parameters

! ###################################################################
! Nondimensional numbers
! ###################################################################
    real(wp) :: reynolds, prandtl, schmidt(MAX_NSP) ! molecular transport
    real(wp) :: mach                                ! compressibility
    real(wp) :: damkohler(MAX_NSP)                  ! reaction
    real(wp) :: froude                              ! body force
    real(wp) :: rossby                              ! Coriolis force
    real(wp) :: stokes                              ! inertial effects
    real(wp) :: settling                            ! sedimentation effects

    real(wp) :: visc                                ! 1/reynolds, to save computation time

! ###################################################################
! Compact parameters (1st derivative of 6th-order pentadiagonal)
! ###################################################################
    real(wp) :: C1N6M_ALPHA, C1N6M_BETA
    real(wp) :: C1N6M_ALPHA2, C1N6M_BETA2
    real(wp) :: C1N6M_A, C1N6M_B, C1N6M_C
    real(wp) :: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
    real(wp) :: C1N6M_BD2, C1N6M_CD3

! ###########################################################
! Filters
! ###########################################################
    type(filter_dt) :: FilterDomain(3)
    logical :: FilterDomainActive(MAX_VARS)
    integer :: FilterDomainBcsFlow(MAX_VARS), FilterDomainBcsScal(MAX_VARS)

    type(filter_dt) :: Dealiasing(3)
    type(filter_dt) :: vprefil(3)

! ###################################################################
! Jet Statistics
! ###################################################################
    integer :: nstatavg, statavg(MAX_STATS_SPATIAL) ! Ox planes at which to accumulate statistics
    integer :: nstatavg_points                      ! number of accumulated points
    integer :: istattimeorg                         ! time at which accumulation of statistics started
    real(wp) :: rstattimeorg

end module TLAB_VARS
