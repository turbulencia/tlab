module TLAB_TYPES
    use TLAB_CONSTANTS
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

    type grid_dt
        sequence
        character*8 name
        integer(wi) size, inb_grid
        integer mode_fdm                    ! finite-difference method for spatial operators
        logical uniform, periodic, anelastic
        real(wp) scale
        real(wp), pointer :: nodes(:)
        real(wp), pointer :: jac(:, :)      ! pointer to Jacobians
        real(wp), pointer :: lu0i(:, :)     ! pointer to LU decomposition for interpolation
        real(wp), pointer :: lu1(:, :)      ! pointer to LU decomposition for 1. derivative
        real(wp), pointer :: lu1i(:, :)     ! pointer to LU decomposition for 1. derivative inc. interp.
        real(wp), pointer :: lu2(:, :)      ! pointer to LU decomposition for 2. derivative
        real(wp), pointer :: lu2d(:, :)     ! pointer to LU decomposition for 2. derivative inc. diffusion
        real(wp), pointer :: mwn(:, :)      ! pointer to modified wavenumbers
        real(wp), pointer :: rhoinv(:)      ! pointer to density correction in anelastic
    end type grid_dt

    type profiles_dt
        sequence
        integer type, padding
        logical relative
        real(wp) mean, delta, ymean, ymean_rel, thick, lslope, uslope, diam
        real(wp) parameters(MAX_PARS)
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
        real(wp), allocatable :: coeffs(:,:)    ! filted coefficients
    end type filter_dt

    type discrete_dt
        sequence
        integer type, size
        integer, dimension(MAX_MODES) :: modex, modez
        real(wp), dimension(MAX_MODES) :: amplitude, phasex, phasez
        real(wp), dimension(MAX_PARS) :: parameters
    end type discrete_dt

end module TLAB_TYPES
