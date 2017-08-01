#include "types.h"

MODULE DNS_TYPES
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_PARS = 10
  TINTEGER, PARAMETER :: MAX_VARS = 20

  TYPE background_dt
     SEQUENCE
     TINTEGER type
     TREAL reference, mean, delta, ymean, thick, diam
     TREAL, DIMENSION(MAX_PARS) :: parameters
  END TYPE background_dt

  TYPE term_dt
     SEQUENCE
     TINTEGER type
     TINTEGER, DIMENSION(MAX_PARS) :: scalar     ! fields defining this term
     LOGICAL,  DIMENSION(MAX_PARS) :: active     ! fields affected by this term
     TREAL,    DIMENSION(MAX_PARS) :: parameters
     TREAL,    DIMENSION(MAX_PARS) :: auxiliar
     TREAL,    DIMENSION(3)        :: vector
  END TYPE term_dt

  TYPE grid_dt
     SEQUENCE
     CHARACTER*8 name
     TINTEGER size, inb_grid, mode_fdm
     LOGICAL uniform, periodic
     TREAL scale
     TREAL, DIMENSION(:),   POINTER :: nodes
     TREAL, DIMENSION(:,:), POINTER :: jac   ! pointer to Jacobians
     TREAL, DIMENSION(:,:), POINTER :: lu1   ! pointer to LU decomposition for 1. derivative
     TREAL, DIMENSION(:,:), POINTER :: lu2   ! pointer to LU decomposition for 2. derivative
     TREAL, DIMENSION(:,:), POINTER :: lu2d  ! pointer to LU decomposition for 2. derivative inc. diffusion
     TREAL, DIMENSION(:,:), POINTER :: mwn   ! pointer to modified wavenumbers
  END TYPE grid_dt

  TYPE pointers_dt
     SEQUENCE
     TREAL, DIMENSION(:), POINTER :: field
  END TYPE pointers_dt

#ifdef USE_MPI
#include "mpif.h"
#endif

  TYPE subarray_dt
     SEQUENCE
     LOGICAL active
     INTEGER communicator
     INTEGER subarray
#ifdef USE_MPI
     INTEGER(KIND=MPI_OFFSET_KIND) offset
#else
     INTEGER offset
#endif
  END type subarray_dt

  TYPE filter_dt
     SEQUENCE
     TINTEGER type, size, inb_filter, delta
     TREAL alpha
     TINTEGER bcs_min, bcs_max
     LOGICAL uniform, periodic
     TINTEGER mpitype
     TREAL, DIMENSION(:,:), POINTER :: coeffs ! pointer to coefficients
  END TYPE filter_dt
  
  TYPE buffer_dt
     SEQUENCE
     TINTEGER type                              ! relaxation, filter...
     TINTEGER size                              ! # points in buffer layer
     LOGICAL active(MAX_VARS), hard
     TREAL strength(MAX_VARS), sigma(MAX_VARS), hardvalues(MAX_VARS)  ! parameters in relaxation term
     TREAL, ALLOCATABLE, DIMENSION(:,:)     :: tau  ! pointer to relaxation timescale
     TREAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: ref  ! pointer to reference fields
  END TYPE buffer_dt
  
END MODULE DNS_TYPES
