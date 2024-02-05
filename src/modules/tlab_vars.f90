module TLAB_VARS
    use TLAB_TYPES, only: grid_dt, filter_dt, term_dt, profiles_dt
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
    integer :: imode_sim                ! type of simulation (spatial, temporal)
    integer :: imode_files              ! files format
    integer :: imode_precision_files    ! whether restart files in single or double precision
    integer :: imode_verbosity = 1      ! level of verbosity used in log files
    integer :: imode_eqns               ! set of equations to be solved
    integer :: iadvection, iviscous, idiffusion, itransport ! formulation

    integer :: imode_elliptic           ! finite-difference method for pressure-Poisson and Helmholtz equations

    logical :: flow_on = .true.         ! calculate flow parts of the code
    logical :: scal_on = .true.         ! calculate scal parts of the code
    logical :: fourier_on = .false.     ! using FFT libraries
    logical :: stagger_on = .false.     ! horizontal staggering of pressure

    integer :: imode_ibm                ! IBM mode

! ###################################################################
! Iteration
! ###################################################################
    integer(wi) :: itime                ! iteration number
    real(wp) :: rtime                   ! physical time

! ###################################################################
! Arrays sizes
! ###################################################################
! fields
    integer(wi) :: imax, jmax, kmax                 ! number of grid nodes per direction locally per processor
    integer(wi) :: isize_field                      ! =imax*jmax*kmax, 3D fields sizes locally per processor
    integer(wi) :: inb_flow                         ! # of prognostic 3d flow fields (flow evolution equations)
    integer(wi) :: inb_flow_array                   ! >= inb_flow, # of prognostic and diagnostic 3d flow arrays
    integer(wi) :: inb_scal                         ! # of prognostic 3d scal fields (scal evolution equations)
    integer(wi) :: inb_scal_array                   ! >= inb_scal, # of prognostic and diagnostic 3d scal arrays

! auxiliary arrays
    integer(wi) :: isize_wrk1d, inb_wrk1d           ! 1D scratch arrays
    integer(wi) :: isize_wrk2d, inb_wrk2d           ! 2D scratch arrays
    integer(wi) :: isize_wrk3d                      ! 3D scratch arrays
    integer(wi) :: isize_txc_field, inb_txc         ! 3D arrays for intermediate calculations
    integer(wi) :: isize_txc_dimx, isize_txc_dimz   ! partition for MPI data transposition

! ###################################################################
    type(grid_dt), dimension(3) :: g        ! Grid information along 3 directions
    real(wp) :: area                        ! Horizontal area and volume

! ###################################################################
    type(profiles_dt) :: qbg(3)                 ! Velocity background information
    type(profiles_dt) :: sbg(MAX_NSP)           ! Scalars background information
    type(profiles_dt) :: pbg, rbg, tbg, hbg     ! Pressure, density, temperature, enthalpy background information

    real(wp), allocatable :: sbackground(:, :)  ! Scalar reference profiles
    real(wp), allocatable :: bbackground(:)     ! Buoyancy

! ###################################################################
    type(term_dt) :: buoyancy               ! Buoyancy parameters
    type(term_dt) :: coriolis               ! Coriolis parameters
    type(term_dt) :: radiation              ! Radiation parameters
    type(term_dt) :: transport              ! Transport parameters
    type(term_dt) :: chemistry              ! Chemistry parameters
    type(term_dt) :: subsidence             ! Large-scale parameters
    type(term_dt) :: random                 ! Random Forcing parameters

! ###################################################################
! Nondimensional numbers
! ###################################################################
    real(wp) :: visc, prandtl, schmidt(MAX_NSP)     ! molecular transport
    real(wp) :: mach                                ! compressibility
    real(wp) :: damkohler(MAX_NSP)                  ! reaction
    real(wp) :: froude                              ! body force
    real(wp) :: rossby                              ! Coriolis force
    real(wp) :: stokes                              ! inertial effects
    real(wp) :: settling                            ! sedimentation effects

! ###########################################################
! Filters
! ###########################################################
    type(filter_dt) :: FilterDomain(3)
    logical :: FilterDomainActive(MAX_VARS)
    integer :: FilterDomainBcsFlow(MAX_VARS), FilterDomainBcsScal(MAX_VARS)

    type(filter_dt) :: Dealiasing(3)
    type(filter_dt) :: PressureFilter(3)

! ###################################################################
! Jet Statistic
! ###################################################################
    integer :: nstatavg, statavg(MAX_STATS_SPATIAL) ! Ox planes at which to accumulate statistics
    integer :: nstatavg_points                      ! number of accumulated points
    integer :: istattimeorg                         ! time at which accumulation of statistics started
    real(wp) :: rstattimeorg

end module TLAB_VARS
