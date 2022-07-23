#include "types.h"

module TLAB_ARRAYS
    implicit none
    save
    private

    TREAL, allocatable, public :: x(:, :), y(:, :), z(:, :)     ! Grid and associated arrays
    TREAL, allocatable, public :: q(:, :)                       ! Eulerian fields, flow vartiables
    TREAL, allocatable, public :: s(:, :)                       ! Eulerian fields, scalar variables
    TREAL, allocatable, public :: txc(:, :)                     ! Temporary space for Eulerian fields
    TREAL, allocatable, public :: wrk1d(:, :)                   ! Work arrays (scratch space)
    TREAL, allocatable, public :: wrk2d(:, :)                   ! Work arrays (scratch space)
    TREAL, allocatable, public :: wrk3d(:)                      ! Work arrays (scratch space)

    target x, y, z
    target q, s, txc

end module TLAB_ARRAYS
