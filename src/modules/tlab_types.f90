#include "types.h"

module TLAB_TYPES
#ifdef USE_MPI
    use MPI
#endif
    implicit none
    save

    TINTEGER, parameter :: sp = KIND(1.0)
    TINTEGER, parameter :: dp = KIND(1.0d0)
    TINTEGER, parameter :: cp = dp             ! code precision
    ! !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
    ! integer, parameter :: sp = selected_real_kind(6, 37)
    ! !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
    ! integer, parameter :: dp = selected_real_kind(15, 307)

    TINTEGER, parameter :: MAX_PARS = 10
    TINTEGER, parameter :: MAX_VARS = 20
    TINTEGER, parameter :: MAX_MODES = 20

    type discrete_dt
        sequence
        TINTEGER type, size
        TINTEGER, dimension(MAX_MODES) :: modex, modez
        TREAL, dimension(MAX_MODES) :: amplitude, phasex, phasez
        TREAL, dimension(MAX_PARS) :: parameters
    end type discrete_dt

    type background_dt
        sequence
        TINTEGER type, padding
        TREAL reference, mean, delta, ymean, thick, diam
        TREAL, dimension(MAX_PARS) :: parameters
    end type background_dt

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
