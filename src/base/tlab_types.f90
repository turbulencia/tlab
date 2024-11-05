module TLab_Types
    use TLab_Constants
    implicit none
    save

    type pointers_dt
        sequence
        character(len=32) :: tag
        real(wp), pointer :: field(:)
    end type pointers_dt

    type pointers3d_dt
        sequence
        character(len=32) :: tag
        real(wp), pointer :: field(:, :, :)
    end type pointers3d_dt

    type term_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type term_dt

    type :: profiles_dt                             ! I wonder if this should be in module profiles, which needs to change dependecies...
        sequence
        integer type
        integer :: padding = 0_i4_
        logical :: relative = .true.                ! use reference spatial position relative to the extent of the domain
        real(wp) :: mean = 0.0_wp                   ! mean value of f
        real(wp) :: delta = 1.0_wp                  ! increment of f
        real(wp) :: ymean = 0.0_wp                  ! reference spatial position at which f changes      
        real(wp) :: ymean_rel = 0.5_wp              ! same but relative to the extent of the domain
        real(wp) :: thick = 1.0_wp                  ! spatial interval over which f changes
        real(wp) :: lslope = 0.0_wp                 ! slope of f below the ymean
        real(wp) :: uslope = 0.0_wp                 ! slope of f above ymean
        real(wp) :: diam = 0.0_wp                   ! diameter
        real(wp) :: parameters(MAX_PARS) = 0.0_wp   ! additional parameters
    end type profiles_dt

    type filter_dt
        sequence
        integer type, ipadding
        integer(wi) size, inb_filter
        logical uniform, periodic, lpadding(2)
        real(wp) parameters(MAX_PARS)
        integer BcsMin, BcsMax                  ! boundary conditions
        integer repeat
        integer mpitype
        real(wp), allocatable :: coeffs(:, :)    ! filted coefficients
    end type filter_dt

    type discrete_dt
        sequence
        integer type, size
        integer, dimension(MAX_MODES) :: modex, modez
        real(wp), dimension(MAX_MODES) :: amplitude, phasex, phasez
        real(wp), dimension(MAX_PARS) :: parameters
    end type discrete_dt

end module TLab_Types
