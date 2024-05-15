module TLAB_CONSTANTS
    implicit none
    save

    integer, parameter :: MajorVersion = 7
    integer, parameter :: MinorVersion = 0

    integer, parameter :: MAX_PARS = 10
    integer, parameter :: MAX_VARS = 10
    integer, parameter :: MAX_MODES = 20
    integer, parameter :: MAX_PROF = 10
    integer, parameter :: MAX_JETS = 5
    integer, parameter :: MAX_AVG_TEMPORAL = 235
    integer, parameter :: MAX_STATS_SPATIAL = 100 ! Running statistics
    integer, parameter :: MAX_PATH_LENGTH = 128

    character(len=*), parameter :: gfile = 'grid'
    character(len=*), parameter :: ifile = 'tlab.ini'
    character(len=*), parameter :: lfile = 'tlab.log'
    character(len=*), parameter :: efile = 'tlab.err'
    character(len=*), parameter :: wfile = 'tlab.war'
    character(len=*), parameter :: tfile = 'dns.trc'

    character(len=*), parameter :: tag_flow = 'flow.'
    character(len=*), parameter :: tag_scal = 'scal.'
    character(len=*), parameter :: tag_part = 'part.'
    character(len=*), parameter :: tag_traj = 'traj.'

    character(len=*), parameter :: fmt_r = 'e13.5e3'

    ! from https://fortran-lang.org/en/learn/best_practices/floating_point/
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
! !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
! integer, parameter :: sp = selected_real_kind(6, 37)
! !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
! integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: wp = dp             ! working precision

! !> Char length for integers, range -2⁷ to 2⁷-1; 8 bits
! integer, parameter :: i1 = selected_int_kind(2)
! !> Short length for integers, range -2¹⁵ to 2¹⁵-1; 16 bits
! integer, parameter :: i2 = selected_int_kind(4)
!> Length of default integers, range -2³¹ to 2³¹-1; 32 bits
    integer, parameter :: i4_ = selected_int_kind(9)            ! i4 was already used...
! !> Long length for integers, range -2⁶³ to 2⁶³-1; 64 bits
! integer, parameter :: i8 = selected_int_kind(18)
    integer, parameter :: i8_ = selected_int_kind(18)           ! i8 was already used...
    integer, parameter :: wi = i4_                ! working integer type
    integer, parameter :: longi = i8_             ! long integer type; different variable name to avoid errors

    integer, parameter :: sizeofreal = sizeof(1.0_wp)
    integer, parameter :: sizeofint = sizeof(1_wi)
    integer, parameter :: sizeoflongint = sizeof(1_longi)

    real(wp), parameter :: pi_wp = 3.14159265358979323846_wp
    real(wp), parameter :: small_wp = 1.0e-20_wp
    real(wp), parameter :: big_wp = 1.0e20_wp

    integer, parameter :: BCS_PERIODIC = -1
    integer, parameter :: BCS_DD = 0     ! Dirichlet/Dirichlet
    integer, parameter :: BCS_ND = 1     ! Neumann/Dirichlet
    integer, parameter :: BCS_DN = 2     ! Dirichlet/Neumann
    integer, parameter :: BCS_NN = 3     ! Neumann/Neumann

    integer, parameter :: BCS_NONE = 0   ! No special treatment of boundaries
    integer, parameter :: BCS_MIN = 1    ! Special treatment at the lower interval limit
    integer, parameter :: BCS_MAX = 2    ! Special treatment at the upper interval limit
    integer, parameter :: BCS_BOTH = 3

end module TLAB_CONSTANTS
