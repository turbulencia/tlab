#include "types.h"

MODULE TLAB_ARRAYS
  IMPLICIT NONE
  SAVE
  PRIVATE

  TREAL, ALLOCATABLE, PUBLIC :: x(:,:),y(:,:),z(:,:)  ! Grid and associated arrays
  TREAL, ALLOCATABLE, PUBLIC :: q(:,:)                ! Eulerian fields, flow vartiables
  TREAL, ALLOCATABLE, PUBLIC :: s(:,:)                ! Eulerian fields, scalar variables
  TREAL, ALLOCATABLE, PUBLIC :: txc(:,:)              ! Temporary space for Eulerian fields
  TREAL, ALLOCATABLE, PUBLIC :: wrk1d(:,:)            ! Work arrays (scratch space)
  TREAL, ALLOCATABLE, PUBLIC :: wrk2d(:,:)            ! Work arrays (scratch space)
  TREAL, ALLOCATABLE, PUBLIC :: wrk3d(:)              ! Work arrays (scratch space)

  TARGET x, y, z
  TARGET q, s, txc

END MODULE TLAB_ARRAYS
