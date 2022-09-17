#include "types.h"

module TLAB_TYPES
    use TLAB_CONSTANTS
#ifdef USE_MPI
    use MPI
#endif
    implicit none
    save

    type discrete_dt
        sequence
        TINTEGER type, size
        TINTEGER, dimension(MAX_MODES) :: modex, modez
        TREAL, dimension(MAX_MODES) :: amplitude, phasex, phasez
        TREAL, dimension(MAX_PARS) :: parameters
    end type discrete_dt

    type profiles_dt
        sequence
        TINTEGER type, padding
        logical relative
        TREAL mean, delta, ymean, ymean_rel, thick, lslope, uslope, diam
        TREAL, dimension(MAX_PARS) :: parameters
    end type profiles_dt

    type term_dt
        sequence
        TINTEGER type
        TINTEGER, dimension(MAX_PARS) :: scalar     ! fields defining this term
        logical, dimension(MAX_PARS) :: active, lpadding(3)     ! fields affected by this term
        TREAL, dimension(MAX_PARS) :: parameters
        TREAL, dimension(MAX_PARS) :: auxiliar
        TREAL, dimension(3) :: vector
    end type term_dt

    type grid_dt
        sequence
        character*8 name
        TINTEGER size, inb_grid, mode_fdm
        logical uniform, periodic, anelastic
        TREAL scale
        TREAL, dimension(:), pointer :: nodes
        TREAL, dimension(:, :), pointer :: jac   ! pointer to Jacobians
        TREAL, dimension(:, :), pointer :: lu0i  ! pointer to LU decomposition for interpolation
        TREAL, dimension(:, :), pointer :: lu1   ! pointer to LU decomposition for 1. derivative
        TREAL, dimension(:, :), pointer :: lu1i  ! pointer to LU decomposition for 1. derivative inc. interp.
        TREAL, dimension(:, :), pointer :: lu2   ! pointer to LU decomposition for 2. derivative
        TREAL, dimension(:, :), pointer :: lu2d  ! pointer to LU decomposition for 2. derivative inc. diffusion
        TREAL, dimension(:, :), pointer :: mwn   ! pointer to modified wavenumbers
        TREAL, dimension(:), pointer :: rhoinv! pointer to density correction in anelastic
    end type grid_dt

    type filter_dt
        sequence
        TINTEGER type, ipadding
        TINTEGER size, inb_filter
        logical uniform, periodic, lpadding(2)
        TREAL, dimension(MAX_PARS) :: parameters
        TINTEGER BcsMin, BcsMax                  ! boundary conditions
        TINTEGER repeat
        TINTEGER mpitype
        TREAL, dimension(:, :), allocatable :: coeffs ! pointer to coefficients
    end type filter_dt

    type pointers_dt
        sequence
        character*32 :: tag
        TREAL, dimension(:), pointer :: field
    end type pointers_dt

    type pointers3d_dt
        sequence
        TREAL, dimension(:, :, :), pointer :: field
    end type pointers3d_dt

    type subarray_dt
        sequence
#ifdef USE_MPI
        logical active, lpadding(3)
        integer communicator
        integer subarray
        integer(KIND=MPI_OFFSET_KIND) offset
#else
        integer offset
#endif
    end type subarray_dt

end module TLAB_TYPES
