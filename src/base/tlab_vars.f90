module TLAB_VARS
    use TLab_Types, only: filter_dt, term_dt, profiles_dt
    use TLab_Constants, only: MAX_VARS, wp, wi, sp
    implicit none
    save

! ###################################################################
! General options
! ###################################################################
    integer :: imode_sim                ! type of simulation (spatial, temporal)
    integer :: imode_files              ! files format
    integer :: imode_precision_files    ! whether restart files in single or double precision
    integer :: imode_verbosity = 1      ! level of verbosity used in log files

    logical :: flow_on = .true.         ! calculate flow parts of the code
    logical :: scal_on = .true.         ! calculate scal parts of the code
    logical :: fourier_on = .false.     ! using FFT libraries
    logical :: stagger_on = .false.     ! horizontal staggering of pressure

    integer :: imode_eqns                       ! set of equations to be solved: internal energy, total energy, anelastic, Boussinesq
    integer :: iadvection, iviscous, idiffusion ! formulation of the Burgers operator

! ###################################################################
! Iteration
! ###################################################################
    integer(wi) :: itime                ! iteration number
    real(wp) :: rtime                   ! physical time

! ###################################################################
! Arrays sizes
! ###################################################################
! fields
    integer(wi) :: imax, jmax, kmax     ! number of grid nodes per direction locally per processor
    integer(wi) :: isize_field          ! =imax*jmax*kmax, 3D fields sizes locally per processor
    integer(wi) :: inb_flow             ! # of prognostic 3d flow fields (flow evolution equations)
    integer(wi) :: inb_flow_array       ! >= inb_flow, # of prognostic and diagnostic 3d flow arrays
    integer(wi) :: inb_scal             ! # of prognostic 3d scal fields (scal evolution equations)
    integer(wi) :: inb_scal_array       ! >= inb_scal, # of prognostic and diagnostic 3d scal arrays

! auxiliary arrays
    integer(wi) :: isize_wrk1d, inb_wrk1d           ! 1D scratch arrays
    integer(wi) :: isize_wrk2d, inb_wrk2d           ! 2D scratch arrays
    integer(wi) :: isize_wrk3d                      ! 3D scratch array (only 1)
    integer(wi) :: isize_txc_field, inb_txc         ! 3D arrays for intermediate calculations
    integer(wi) :: isize_txc_dimx, isize_txc_dimz   ! partition for MPI data transposition

! ###################################################################
! information to set up bcs, ics, and reference background profiles
! ###################################################################
    type(profiles_dt) :: qbg(3)                     ! Velocity
    type(profiles_dt) :: sbg(MAX_VARS)              ! Scalars
    type(profiles_dt) :: pbg, rbg, tbg, hbg         ! Pressure, density, temperature, enthalpy

! ###################################################################
! phenomena in addition to the navier-stokes equations
! ###################################################################
    type(term_dt) :: buoyancy                       ! Buoyancy parameters
    type(term_dt) :: coriolis                       ! Coriolis parameters
    type(term_dt) :: subsidence                     ! Large-scale parameters

! ###################################################################
! Nondimensional numbers
! ###################################################################
    real(wp) :: visc, prandtl, schmidt(MAX_VARS)    ! molecular transport
    real(wp) :: mach                                ! compressibility
    real(wp) :: damkohler(MAX_VARS)                 ! reaction
    real(wp) :: froude                              ! gravity force
    real(wp) :: rossby                              ! Coriolis force
    real(wp) :: stokes                              ! particle inertial effects
    real(wp) :: settling                            ! sedimentation effects

! ###########################################################
! Filters
! ###########################################################
    type(filter_dt) :: FilterDomain(3)
    logical :: FilterDomainActive(MAX_VARS)
    integer :: FilterDomainBcsFlow(MAX_VARS), FilterDomainBcsScal(MAX_VARS)

    type(filter_dt) :: Dealiasing(3)
    type(filter_dt) :: PressureFilter(3)

end module TLAB_VARS

! ###################################################################
! Jet Statistic
! ###################################################################
module TLab_Spatial
    use TLab_Types, only: wp
    implicit none
    save

    integer, parameter :: MAX_STATS_SPATIAL = 100 ! Running statistics

    integer :: nstatavg, statavg(MAX_STATS_SPATIAL) ! Ox planes at which to accumulate statistics
    integer :: nstatavg_points                      ! number of accumulated points
    integer :: istattimeorg                         ! time at which accumulation of statistics started
    real(wp) :: rstattimeorg

end module TLAB_Spatial
