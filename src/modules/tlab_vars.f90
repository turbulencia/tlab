#include "types.h"
#include "dns_const.h"

module TLAB_VARS
    use TLAB_TYPES, only: grid_dt, filter_dt, subarray_dt, term_dt, profiles_dt
    use TLAB_CONSTANTS, only: MAX_VARS, MAX_NSP
    use TLAB_CONSTANTS, only: MAX_STATS_SPATIAL
    implicit none
    save

! ###################################################################
! OpenMP
! ###################################################################
    TINTEGER :: dns_omp_numThreads
    TINTEGER :: dns_omp_error

! ###################################################################
! General options
! ###################################################################
    TINTEGER :: icalc_flow, icalc_scal
    TINTEGER :: imode_sim                ! type of simulation (spatial, temporal)
    TINTEGER :: imode_files              ! files format
    TINTEGER :: imode_precision_files
    TINTEGER :: imode_verbosity = 1      ! level of verbosity used in log files
    TINTEGER :: imode_eqns               ! set of equations to be solved
    TINTEGER :: iadvection, iviscous, idiffusion ! formulation
    TINTEGER :: ifourier
    TINTEGER :: itransport
    TINTEGER :: istagger, ivfilter       ! horizontal staggering of pressure
    ! vertical   filtering  of pressure
    TREAL :: vfilter_param            ! vertical filter parameter

    TINTEGER :: imode_fdm                ! finite-difference method for spatial operators

    TINTEGER :: imode_ibm                ! IBM mode

! ###################################################################
! Iteration
! ###################################################################
    TINTEGER :: itime                    ! iteration number
    TREAL :: rtime                    ! physical time

! ###################################################################
! Arrays size
! ###################################################################
! fields
    TINTEGER :: imax, jmax, kmax, isize_field     ! locally per processor
    TINTEGER :: inb_flow, inb_flow_array        ! transported & array space
    TINTEGER :: inb_scal, inb_scal_array        ! transported & array space

! auxiliary arrays
    TINTEGER :: isize_wrk1d, inb_wrk1d      ! 1D arrays
    TINTEGER :: isize_wrk2d, inb_wrk2d      ! 2D arrays
    TINTEGER :: isize_wrk3d                     ! 2D arrays
    TINTEGER :: isize_txc_field, inb_txc        ! 3D arrays for intermediate calculations
    TINTEGER :: isize_txc_dimx, isize_txc_dimz  ! partition for MPI data transposition

! subarray information (offset)
    type(subarray_dt), dimension(IO_SUBARRAY_SIZE) :: io_aux

! ###################################################################
    type(grid_dt), dimension(3) :: g      ! Grid information along 3 directions
    TREAL :: area                         ! Horizontal area and volume

! ###################################################################
    type(profiles_dt) :: qbg(3)         ! Velocity background
    type(profiles_dt) :: sbg(MAX_NSP)   ! Scalars backgrounds
    type(profiles_dt) :: pbg, rbg, tbg  ! Pressure, density, temperature backgrounds

    TREAL, dimension(:), allocatable :: pbackground, tbackground, rbackground, ribackground
    TREAL, dimension(:), allocatable :: bbackground, epbackground

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
    TREAL :: reynolds, prandtl, schmidt(MAX_NSP) ! molecular transport
    TREAL :: mach                                ! compressibility
    TREAL :: damkohler(MAX_NSP)                  ! reaction
    TREAL :: froude                              ! body force
    TREAL :: rossby                              ! Coriolis force
    TREAL :: stokes                              ! Stokes number of liquid particles
    TREAL :: settling                            ! sedimentation parameter for liquid particle

    TREAL :: visc                                ! 1/reynolds

! ###################################################################
! Compact parameters (1st derivative of 6th-order pentadiagonal)
! ###################################################################
    TREAL :: C1N6M_ALPHA, C1N6M_BETA
    TREAL :: C1N6M_ALPHA2, C1N6M_BETA2
    TREAL :: C1N6M_A, C1N6M_B, C1N6M_C
    TREAL :: C1N6M_AD2, C1N6M_BD4, C1N6M_CD6
    TREAL :: C1N6M_BD2, C1N6M_CD3

! ###########################################################
! Filters
! ###########################################################
    type(filter_dt), dimension(3) :: FilterDomain
    logical, dimension(MAX_VARS) :: FilterDomainActive
    TINTEGER, dimension(MAX_VARS) :: FilterDomainBcsFlow, FilterDomainBcsScal

    type(filter_dt), dimension(3) :: Dealiasing

! ###################################################################
! FFTW
! ###################################################################
    integer(8) :: fft_plan_fx, fft_plan_bx, fft_plan_fx_bcs
    integer(8) :: fft_plan_fy, fft_plan_by, fft_plan_fy1d, fft_plan_by1d
    integer(8) :: fft_plan_fz, fft_plan_bz

    TINTEGER :: fft_reordering

! ###################################################################
! Jet Statistics
! ###################################################################
    TINTEGER :: nstatavg, statavg(MAX_STATS_SPATIAL), & ! Ox planes at which to accumulate statistics
        nstatavg_points, &                      ! number of accumulated points
        istattimeorg                            ! time at which accumulation of statistics started
    TREAL :: rstattimeorg

end module TLAB_VARS
