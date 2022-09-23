module PARTICLE_VARS
    use TLAB_CONSTANTS, only: wp, wi, longi, MAX_PARS, MAX_NSP
    use PARTICLE_TYPES
    implicit none
    save

    ! Possible values of imode_part
    integer, parameter :: PART_TYPE_NONE = 0
    integer, parameter :: PART_TYPE_TRACER = 1
    integer, parameter :: PART_TYPE_SIMPLE_SETT = 2
    integer, parameter :: PART_TYPE_BIL_CLOUD_3 = 3
    integer, parameter :: PART_TYPE_BIL_CLOUD_4 = 4

    ! Posible values of imode_traj
    integer, parameter :: TRAJ_TYPE_NONE = 0
    integer, parameter :: TRAJ_TYPE_FIRST = 1
    integer, parameter :: TRAJ_TYPE_LARGEST = 2
    integer, parameter :: TRAJ_TYPE_VORTICITY = 3

    integer(wi)       :: imode_part                   ! type if particle formulation, e.g., tracer, inertia...
    integer(longi)    :: isize_part_total             ! total # of particles

    integer(wi)       :: isize_part                   ! maximum # of particles per processor (to allocate memory space)
    integer(wi)       :: inb_part_array               ! # of particle properties in arrays (prognostic & diagnostic)
    integer(wi)       :: inb_part                     ! # of particle properties in Runge-Kutta (prognostic)
    integer(wi)       :: inb_part_txc                 ! # of particle auxiliary properties for intermediate calculations
    integer(wi)       :: inb_part_interp              ! # of interpolated fields into lagrangian framework
    integer(wi)       :: isize_l_comm                 ! memory space for the halo regions

#ifdef USE_MPI
    integer(wi)   :: isize_pbuffer                    ! space for communication of halo regions
#endif

    ! Initialization
    integer(wi)   :: part_ini_mode       ! Type of initialization
    real(wp)      :: part_ini_ymean      ! Mean position where particles positions will be initialize
    real(wp)      :: part_ini_thick      ! Width of initial particle distribution

    ! Trajectory
    integer(wi)   :: imode_traj          ! Type of trajectories
    integer(wi)   :: isize_traj          ! # of saved trajectories
    integer(wi)   :: inb_traj            ! # of properties saved along trajectories

    ! Calculation of residence times
    integer(wi)   :: residence_reset     !if reseidence l_q should be reset
    real(wp)      :: l_y_lambda          !y coordinate where approx radiation begins for residence times (set in dns_main)
    real(wp)      :: l_y_base            !set to be 1/3 of cloud domain between two bouyancy stratification for residence times

    ! Calculation of pdfs
    integer(wi)   :: icalc_part_pdf      ! if calculation of pdf for particles
    real(wp)      :: particle_pdf_subdomain(6)
    real(wp)      :: particle_pdf_max
    real(wp)      :: particle_pdf_interval

    ! Auxiliary   data
    real(wp)      :: particle_param(MAX_PARS) ! lagrange function parameters
    character*32  :: particle_spname(MAX_NSP) !Name of different lagrange species

end module PARTICLE_VARS
