module TLAB_VARS
    use TLab_Types, only: term_dt
    use TLab_Constants, only: wp, wi, sp
    implicit none
    save

! ###################################################################
! General options
! ###################################################################
    integer :: imode_sim                ! type of simulation (spatial, temporal)

    logical :: flow_on = .true.         ! calculate flow parts of the code
    logical :: scal_on = .true.         ! calculate scal parts of the code
    logical :: fourier_on = .false.     ! using FFT libraries
    logical :: stagger_on = .false.     ! horizontal staggering of pressure

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
! phenomena in addition to the navier-stokes equations
! ###################################################################
    type(term_dt) :: coriolis                       ! Coriolis parameters

    ! nondimensional parameters; to be moved to navierstokes
    real(wp), public :: mach                                ! compressibility

end module TLAB_VARS
