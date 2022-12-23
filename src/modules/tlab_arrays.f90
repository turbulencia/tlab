module TLAB_ARRAYS
    use TLAB_CONSTANTS, only: wp
    implicit none
    save

    real(wp), allocatable :: x(:, :), y(:, :), z(:, :)     ! Grid and associated arrays
    real(wp), allocatable :: q(:, :)                       ! Eulerian fields, flow vartiables
    real(wp), allocatable :: s(:, :)                       ! Eulerian fields, scalar variables
    real(wp), allocatable :: txc(:, :)                     ! Temporary space for Eulerian fields
    real(wp), allocatable :: wrk1d(:, :)                   ! Work arrays (scratch space)
    real(wp), allocatable :: wrk2d(:, :)                   ! Work arrays (scratch space)
    real(wp), allocatable :: wrk3d(:)                      ! Work arrays (scratch space)

    target x, y, z
    target q, s, txc, wrk3d

end module TLAB_ARRAYS
