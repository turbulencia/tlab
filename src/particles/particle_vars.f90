module PARTICLE_VARS
    use TLAB_CONSTANTS, only: wp, wi, longi, MAX_PARS, MAX_NSP
    use TLAB_TYPES, only: profiles_dt, term_dt
    use PARTICLE_TYPES
    implicit none
    save

    ! Possible values of part%type
    integer, parameter :: PART_TYPE_NONE = 0
    integer, parameter :: PART_TYPE_TRACER = 1
    integer, parameter :: PART_TYPE_SIMPLE_SETT = 2
    integer, parameter :: PART_TYPE_BIL_CLOUD_3 = 3
    integer, parameter :: PART_TYPE_BIL_CLOUD_4 = 4

    ! Posible values of imode_traj
    integer, parameter :: TRAJ_TYPE_NONE = 0
    integer, parameter :: TRAJ_TYPE_BASIC = 1
    integer, parameter :: TRAJ_TYPE_VORTICITY = 2

    type(term_dt)     :: part                         ! particle formulation, e.g., tracer, inertia... Maybe new derived type

    character(len=32) :: part_spname(MAX_NSP)

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
    type(profiles_dt) :: IniP                         ! Information about the initialization 
    integer, parameter :: PART_INITYPE_SCALAR = 101   ! Special type of particle initialization not included in default profile data

    ! Trajectory
    integer(wi)   :: imode_traj             ! Type of trajectory information that is saved
    integer(wi)   :: isize_traj             ! # of saved trajectories
    integer(wi)   :: inb_traj               ! # of properties saved along trajectories
    character(len=32) :: traj_filename      ! file with the particle tags to be tracked; if void, then the first isize_traj particles are used

    ! Calculation of residence times
    integer(wi)   :: residence_reset     !if reseidence l_q should be reset
    real(wp)      :: l_y_lambda          !y coordinate where approx radiation begins for residence times (set in dns_main)
    real(wp)      :: l_y_base            !set to be 1/3 of cloud domain between two bouyancy stratification for residence times

    ! Calculation of pdfs
    logical       :: particle_pdf_calc      ! if calculation of pdf for particles
    real(wp)      :: particle_pdf_subdomain(6)
    real(wp)      :: particle_pdf_max
    real(wp)      :: particle_pdf_interval

end module PARTICLE_VARS
