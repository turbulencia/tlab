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
    real(wp), allocatable :: wrkdea(:,:)                   ! Work arrays for dealiasing (scratch space)

    target x, y, z
    target q, s, txc, wrk1d, wrk2d, wrk3d, wrkdea

end module TLAB_ARRAYS

module TLAB_POINTERS
    use TLAB_CONSTANTS, only: wp
    implicit none

    real(wp), pointer :: u(:) => null()
    real(wp), pointer :: v(:) => null()
    real(wp), pointer :: w(:) => null()
    real(wp), pointer :: e(:) => null()
    real(wp), pointer :: rho(:) => null()
    real(wp), pointer :: p(:) => null()
    real(wp), pointer :: T(:) => null()
    real(wp), pointer :: vis(:) => null()

    real(wp), pointer :: tmp1(:) => null()
    real(wp), pointer :: tmp2(:) => null()
    real(wp), pointer :: tmp3(:) => null()
    real(wp), pointer :: tmp4(:) => null()
    real(wp), pointer :: tmp5(:) => null()
    real(wp), pointer :: tmp6(:) => null()
    real(wp), pointer :: tmp7(:) => null()
    real(wp), pointer :: tmp8(:) => null()
    real(wp), pointer :: tmp9(:) => null()

end module TLAB_POINTERS

module TLAB_POINTERS_3D
    use TLAB_CONSTANTS, only: wp
    implicit none

    real(wp), pointer :: u(:, :, :) => null()
    real(wp), pointer :: v(:, :, :) => null()
    real(wp), pointer :: w(:, :, :) => null()
    real(wp), pointer :: e(:, :, :) => null()
    real(wp), pointer :: rho(:, :, :) => null()
    real(wp), pointer :: p(:, :, :) => null()
    real(wp), pointer :: T(:, :, :) => null()
    real(wp), pointer :: vis(:, :, :) => null()

    real(wp), pointer :: p_q(:, :, :, :) => null()
    real(wp), pointer :: p_s(:, :, :, :) => null()
    real(wp), pointer :: p_wrk1d(:, :) => null()
    real(wp), pointer :: p_wrk2d(:, :, :) => null()
    real(wp), pointer :: p_wrk3d(:, :, :) => null()

    real(wp), pointer :: tmp1(:, :, :) => null()
    real(wp), pointer :: tmp2(:, :, :) => null()
    real(wp), pointer :: tmp3(:, :, :) => null()
    real(wp), pointer :: tmp4(:, :, :) => null()
    real(wp), pointer :: tmp5(:, :, :) => null()
    real(wp), pointer :: tmp6(:, :, :) => null()
    real(wp), pointer :: tmp7(:, :, :) => null()
    real(wp), pointer :: tmp8(:, :, :) => null()
    real(wp), pointer :: tmp9(:, :, :) => null()

end module TLAB_POINTERS_3D

module TLAB_POINTERS_C
    use TLAB_CONSTANTS, only: wp
    implicit none

    complex(wp), pointer :: c_wrk1d(:, :) => null()
    complex(wp), pointer :: c_wrk3d(:, :) => null()

end module TLAB_POINTERS_C
