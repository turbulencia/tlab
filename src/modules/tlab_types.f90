#include "dns_const.h"

module TLAB_TYPES
    use TLAB_CONSTANTS
#ifdef USE_MPI
    use MPI
#endif
    implicit none
    save

    type discrete_dt
        sequence
        integer type, size
        integer, dimension(MAX_MODES) :: modex, modez
        real(wp), dimension(MAX_MODES) :: amplitude, phasex, phasez
        real(wp), dimension(MAX_PARS) :: parameters
    end type discrete_dt

    type profiles_dt
        sequence
        integer type, padding
        logical relative
        real(wp) mean, delta, ymean, ymean_rel, thick, lslope, uslope, diam
        real(wp), dimension(MAX_PARS) :: parameters
    end type profiles_dt

    type term_dt
        sequence
        integer type
        integer, dimension(MAX_PARS) :: scalar     ! fields defining this term
        logical, dimension(MAX_PARS) :: active, lpadding(3)     ! fields affected by this term
        real(wp), dimension(MAX_PARS) :: parameters
        real(wp), dimension(MAX_PARS) :: auxiliar
        real(wp), dimension(3) :: vector
    end type term_dt

    type grid_dt
        sequence
        character*8 name
        integer(wi) size, inb_grid
        integer mode_fdm
        logical uniform, periodic, anelastic
        real(wp) scale
        real(wp), dimension(:), pointer :: nodes
        real(wp), dimension(:, :), pointer :: jac   ! pointer to Jacobians
        real(wp), dimension(:, :), pointer :: lu0i  ! pointer to LU decomposition for interpolation
        real(wp), dimension(:, :), pointer :: lu1   ! pointer to LU decomposition for 1. derivative
        real(wp), dimension(:, :), pointer :: lu1i  ! pointer to LU decomposition for 1. derivative inc. interp.
        real(wp), dimension(:, :), pointer :: lu2   ! pointer to LU decomposition for 2. derivative
        real(wp), dimension(:, :), pointer :: lu2d  ! pointer to LU decomposition for 2. derivative inc. diffusion
        real(wp), dimension(:, :), pointer :: mwn   ! pointer to modified wavenumbers
        real(wp), dimension(:), pointer :: rhoinv! pointer to density correction in anelastic
    end type grid_dt

    type filter_dt
        sequence
        integer type, ipadding
        integer(wi) size, inb_filter
        logical uniform, periodic, lpadding(2)
        real(wp), dimension(MAX_PARS) :: parameters
        integer BcsMin, BcsMax                  ! boundary conditions
        integer repeat
        integer mpitype
        real(wp), dimension(:, :), allocatable :: coeffs ! pointer to coefficients
    end type filter_dt

    type pointers_dt
        sequence
        character*32 :: tag
        real(wp), dimension(:), pointer :: field
    end type pointers_dt

    type pointers3d_dt
        sequence
        real(wp), dimension(:, :, :), pointer :: field
    end type pointers3d_dt

end module TLAB_TYPES
