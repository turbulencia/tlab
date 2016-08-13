#include "types.h"

MODULE DNS_TYPES
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: MAX_PARS = 10

  TYPE grid_structure
     SEQUENCE
     CHARACTER*8 name
     TINTEGER size
     LOGICAL periodic
     TREAL scale
     TREAL, DIMENSION(:),   POINTER :: nodes
     TREAL, DIMENSION(:,:), POINTER :: aux
  END TYPE grid_structure

  TYPE pointers_structure
     SEQUENCE
     TREAL, DIMENSION(:), POINTER :: field
  END TYPE pointers_structure

#ifdef USE_MPI
#include "mpif.h"
#endif

  TYPE subarray_structure
     SEQUENCE
     LOGICAL active
     INTEGER communicator
     INTEGER subarray
#ifdef USE_MPI
     INTEGER(KIND=MPI_OFFSET_KIND) offset
#else
     INTEGER offset
#endif
  END type subarray_structure
  
  TYPE term_structure
     SEQUENCE
     TINTEGER type
     TINTEGER, DIMENSION(MAX_PARS) :: scalar     ! fields defining this term
     LOGICAL,  DIMENSION(MAX_PARS) :: active     ! fields affected by this term
     TREAL,    DIMENSION(MAX_PARS) :: parameters
     TREAL,    DIMENSION(MAX_PARS) :: auxiliar
     TREAL,    DIMENSION(3)        :: vector
  END TYPE term_structure

END MODULE DNS_TYPES
