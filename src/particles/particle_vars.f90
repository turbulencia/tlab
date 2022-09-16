module PARTICLE_VARS
    use TLAB_TYPES, only: wp, wi, longi
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

    integer(wi) :: imode_part                   ! type if particle formulation, e.g., tracer, inertia...
    integer(longi) :: isize_part_total          ! total # of particles

    type(particle_dt) :: l_g                    ! particle tags and Oy-node information in local processor

    integer(wi) :: isize_part                   ! maximum # of particles per processor (to allocate memory space)
    integer(wi) :: inb_part_array               ! # of particle properties in arrays (prognostic & diagnostic)
    integer(wi) :: inb_part                     ! # of particle properties in Runge-Kutta (prognostic)
    integer(wi) :: inb_part_txc                 ! # of particle auxiliary properties for intermediate calculations
    integer(wi) :: inb_part_interp              ! # of interpolated fields into lagrangian framework

#ifdef USE_MPI
    integer(wi), dimension(:), allocatable :: ims_size_p ! vector with all # of particles per processor
#endif

    integer(wi) :: isize_l_comm, isize_pbuffer

    ! Initialization
    integer(wi) :: part_ini_mode    ! Type of initialization
    real(wp) :: part_ini_ymean      ! Mean position where particles positions will be initialize
    real(wp) :: part_ini_thick      ! Width of initial particle distribution

    ! Trajectory
    integer(wi) :: imode_traj       ! Type of trajectories
    integer(wi) :: isize_traj       ! # of saved trajectories
    integer(wi) :: inb_traj         ! # of properties saved along trajectories

    ! Calculation of residence times
    integer(wi) :: residence_reset  !if reseidence l_q should be reset
    real(wp) :: l_y_lambda !y coordinate where approx radiation begins for residence times (set in dns_main)
    real(wp) :: l_y_base   !set to be 1/3 of cloud domain between two bouyancy stratification for residence times

    ! Calculation of pdfs
    integer(wi) :: icalc_part_pdf    ! if calculation of pdf for particles
    real(wp) :: particle_pdf_subdomain(6)
    real(wp) :: particle_pdf_max
    real(wp) :: particle_pdf_interval

    ! Auxiliary data
    integer(wi), parameter :: MAX_LAGPARAM = 10 ! Maximum size of Lagrange Parameters

    real(wp) :: particle_param(MAX_LAGPARAM)                 ! lagrange function parameters
    character*32, dimension(15) :: particle_spname            !Name of different lagrange species

end module PARTICLE_VARS
