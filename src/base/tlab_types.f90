module TLab_Types
    use TLab_Constants, only: wp, wi, i4_, longi, MAX_VARS, MAX_PARS, MAX_MODES
    implicit none
    save

    type term_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type term_dt

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
