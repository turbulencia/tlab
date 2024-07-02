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
        integer mode_fdm1                   ! finite-difference method for 1. order derivative
        integer mode_fdm2                   ! finite-difference method for 2. order derivative
        logical uniform, periodic, anelastic
        logical :: need_1der = .false.      ! In Jacobian formulation, I need 1. order derivative for the 2. order if non-uniform
        integer nb_diag_1(2)                ! # of left and right diagonals 1. order derivative (max 5/7)
        integer nb_diag_2(2)                ! # of left and right diagonals 2. order derivative (max 5/7)
        real(wp) scale
        real(wp) :: rhs1_b(4, 7), rhs1_t(4, 7)          ! RHS data for Neumann boundary conditions, 1. order derivative max. # of diagonals is 7, # rows is 7/2+1
        real(wp) :: rhsr_b(5, 0:7), rhsr_t(0:4, 8)      ! RHS data for reduced boundary conditions; max. # of diagonals is 7, # rows is 7/2+1
        real(wp) :: rhsi_b(5*2, 0:7), rhsi_t(0:9, 8)    ! RHS data for integration, 2x bcs
        real(wp), pointer :: nodes(:)
        real(wp), pointer :: jac(:, :)      ! pointer to Jacobians
        !
        real(wp), pointer :: lhs1(:, :)     ! pointer to LHS for 1. derivative
        real(wp), pointer :: rhs1(:, :)     ! pointer to RHS for 1. derivative
        real(wp), pointer :: lu1(:, :)      ! pointer to LU decomposition for 1. derivative
        real(wp), pointer :: lu0i(:, :)     ! pointer to LU decomposition for interpolation
        real(wp), pointer :: lu1i(:, :)     ! pointer to LU decomposition for 1. derivative inc. interp.
        real(wp), pointer :: lhsi(:, :)     ! pointer to LHS for 1. order integration
        real(wp), pointer :: rhsi(:, :)     ! pointer to RHS for 1. order integration
        real(wp), pointer :: mwn1(:)        ! pointer to modified wavenumbers
        !
        real(wp), pointer :: lhs2(:, :)     ! pointer to LHS for 2. derivative
        real(wp), pointer :: rhs2(:, :)     ! pointer to RHS for 2. derivative
        real(wp), pointer :: lu2(:, :)      ! pointer to LU decomposition for 2. derivative
        real(wp), pointer :: lu2d(:, :)     ! pointer to LU decomposition for 2. derivative inc. diffusion
        real(wp), pointer :: mwn2(:)        ! pointer to modified wavenumbers
        !
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
        real(wp), allocatable :: coeffs(:, :)    ! filted coefficients
    end type filter_dt

    type discrete_dt
        sequence
        integer type, size
        integer, dimension(MAX_MODES) :: modex, modez
        real(wp), dimension(MAX_MODES) :: amplitude, phasex, phasez
        real(wp), dimension(MAX_PARS) :: parameters
    end type discrete_dt

    type decomposition_dt
        sequence
        character(32) :: name
    end type decomposition_dt

    type phaseavg_dt
        sequence
        integer phaseavg
        integer stride
    end type phaseavg_dt

end module TLAB_TYPES
