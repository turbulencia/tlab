#include "types.h"

module TLAB_TYPES
#ifdef USE_MPI
    use MPI
#endif
    implicit none
    save

    ! from https://fortran-lang.org/en/learn/best_practices/floating_point/
    integer, parameter :: sp = KIND(1.0)
    integer, parameter :: dp = KIND(1.0d0)
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
