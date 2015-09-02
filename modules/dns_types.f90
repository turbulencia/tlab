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

END MODULE DNS_TYPES
