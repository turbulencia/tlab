#include "types.h"

MODULE TLAB_TYPES
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_PARS = 10
  TINTEGER, PARAMETER :: MAX_VARS = 20
  TINTEGER, PARAMETER :: MAX_MODES= 20

  TYPE discrete_dt
     SEQUENCE
     TINTEGER type, size
     TINTEGER, DIMENSION(MAX_MODES) :: modex, modez
     TREAL,    DIMENSION(MAX_MODES) :: amplitude, phasex, phasez
     TREAL,    DIMENSION(MAX_PARS)  :: parameters
  END TYPE discrete_dt

  TYPE background_dt
     SEQUENCE
     TINTEGER type, padding
     TREAL reference, mean, delta, ymean, thick, diam
     TREAL, DIMENSION(MAX_PARS) :: parameters
  END TYPE background_dt

  TYPE term_dt
     SEQUENCE
     TINTEGER type
     TINTEGER, DIMENSION(MAX_PARS) :: scalar     ! fields defining this term
     LOGICAL,  DIMENSION(MAX_PARS) :: active, lpadding(3)     ! fields affected by this term
     TREAL,    DIMENSION(MAX_PARS) :: parameters
     TREAL,    DIMENSION(MAX_PARS) :: auxiliar
     TREAL,    DIMENSION(3)        :: vector
  END TYPE term_dt

  TYPE grid_dt
     SEQUENCE
     CHARACTER*8 name
     TINTEGER size, inb_grid, mode_fdm
     LOGICAL uniform, periodic, anelastic
     TREAL scale, fixed_scale
     TREAL, DIMENSION(:),   POINTER :: nodes
     TREAL, DIMENSION(:,:), POINTER :: jac   ! pointer to Jacobians
     TREAL, DIMENSION(:,:), POINTER :: lu1   ! pointer to LU decomposition for 1. derivative
     TREAL, DIMENSION(:,:), POINTER :: lu2   ! pointer to LU decomposition for 2. derivative
     TREAL, DIMENSION(:,:), POINTER :: lu2d  ! pointer to LU decomposition for 2. derivative inc. diffusion
     TREAL, DIMENSION(:,:), POINTER :: mwn   ! pointer to modified wavenumbers
     TREAL, DIMENSION(:),   POINTER :: rhoinv! pointer to density correction in anelastic
  END TYPE grid_dt

  TYPE filter_dt
     SEQUENCE
     TINTEGER type, ipadding
     TINTEGER size, inb_filter
     LOGICAL uniform, periodic, lpadding(2)
     TREAL,    DIMENSION(MAX_PARS) :: parameters
     TINTEGER BcsMin, BcsMax                  ! boundary conditions
     TINTEGER repeat
     TINTEGER mpitype
     TREAL, DIMENSION(:,:), ALLOCATABLE :: coeffs ! pointer to coefficients
  END TYPE filter_dt

  TYPE ibm_geo_dt
     SEQUENCE
     CHARACTER(32) :: name
     TINTEGER      :: number, height, width, length
     LOGICAL       :: mirrored
  END TYPE ibm_geo_dt

  TYPE pointers_dt
     SEQUENCE
     CHARACTER*32                 :: tag
     TREAL, DIMENSION(:), POINTER :: field
  END TYPE pointers_dt

  TYPE pointers3d_dt
     SEQUENCE
     TREAL, DIMENSION(:,:,:), POINTER :: field
  END TYPE pointers3d_dt

#ifdef USE_MPI
#include "mpif.h"
#endif

  TYPE subarray_dt
     SEQUENCE
#ifdef USE_MPI
     LOGICAL active, lpadding(3)
     INTEGER communicator
     INTEGER subarray
     INTEGER(KIND=MPI_OFFSET_KIND) offset
#else
     INTEGER offset
#endif
  END type subarray_dt

END MODULE TLAB_TYPES
