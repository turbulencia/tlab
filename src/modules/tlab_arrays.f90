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

module TLAB_POINTERS
    use TLAB_CONSTANTS, only: wp
    use TLAB_ARRAYS
    implicit none

    real(wp), pointer :: u(:) => null()
    real(wp), pointer :: v(:) => null()
    real(wp), pointer :: w(:) => null()
    real(wp), pointer :: tmp1(:) => null()
    real(wp), pointer :: tmp2(:) => null()
    real(wp), pointer :: tmp3(:) => null()
    real(wp), pointer :: tmp4(:) => null()
    real(wp), pointer :: tmp5(:) => null()
    real(wp), pointer :: tmp6(:) => null()
end module TLAB_POINTERS