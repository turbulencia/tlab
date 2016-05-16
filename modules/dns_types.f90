#include "types.h"

MODULE DNS_TYPES
  IMPLICIT NONE
  SAVE

  TYPE grid_structure
     SEQUENCE
     CHARACTER*8 name
     TINTEGER size
     LOGICAL periodic
     TREAL scale
     TREAL, DIMENSION(:),   ALLOCATABLE :: nodes
     TREAL, DIMENSION(:,:), ALLOCATABLE :: aux
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
  
END MODULE DNS_TYPES
