module TLAB_ARRAYS
    use TLAB_CONSTANTS, only: wp
    implicit none
    save
    private

    real(wp), allocatable, public :: x(:, :), y(:, :), z(:, :)     ! Grid and associated arrays
    real(wp), allocatable, public :: q(:, :)                       ! Eulerian fields, flow vartiables
    real(wp), allocatable, public :: s(:, :)                       ! Eulerian fields, scalar variables
    real(wp), allocatable, public :: txc(:, :)                     ! Temporary space for Eulerian fields
    real(wp), allocatable, public :: wrk1d(:, :)                   ! Work arrays (scratch space)
    real(wp), allocatable, public :: wrk2d(:, :)                   ! Work arrays (scratch space)
    real(wp), allocatable, public :: wrk3d(:)                      ! Work arrays (scratch space)

    target x, y, z
    target q, s, txc, wrk3d

end module TLAB_ARRAYS
