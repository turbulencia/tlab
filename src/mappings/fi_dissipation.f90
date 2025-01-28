#include "dns_const.h"

!########################################################################
!# Calculates turbulent dissipation per unit volume \rho \epsilon = \tau'_{ij}u'_{ij}
!# It assumes constant viscosity
!# It is not multiplied by the viscosity here, to reduce cross-dependencies.
!########################################################################
subroutine FI_DISSIPATION(nx, ny, nz, u, v, w, eps, tmp1, tmp2, tmp3, tmp4)
    use TLab_Constants, only: wp, wi
    use FDM, only: g
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: p_wrk3d
    use Averages, only: AVG_IK_V
    use OPR_PARTIAL
    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
    real(wp), dimension(nx, ny, nz), intent(OUT) :: eps
    real(wp), dimension(nx, ny, nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

! -------------------------------------------------------------------
    integer(wi) bcs(2, 2), j

! ###################################################################
    bcs = 0

! Diagonal terms
    call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp3)
    tmp4 = (tmp1 + tmp2 + tmp3)*2.0_wp/3.0_wp

! 11
    p_wrk3d = 2.0_wp*tmp1 - tmp4 ! )*vis
    ! if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, wrk1d(1, 1), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
        end do
    ! end if
    eps = p_wrk3d*tmp1

! 22
    p_wrk3d = 2.0_wp*tmp2 - tmp4 ! )*vis
    ! if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, wrk1d(1, 1), wrk1d(1, 2))
        call AVG_IK_V(nx, ny, nz, ny, v, wrk1d(1, 3), wrk1d(1, 2))
        call OPR_PARTIAL_Y(OPR_P1, 1, ny, 1, bcs, g(2), wrk1d(1, 3), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
            tmp2(:, j, :) = tmp2(:, j, :) - wrk1d(j, 2)
        end do
    ! end if
    eps = eps + p_wrk3d*tmp2

! 33
    p_wrk3d = 2.0_wp*tmp3 - tmp4 ! )*vis
    ! if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, wrk1d(1, 1), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
        end do
    ! end if
    eps = eps + p_wrk3d*tmp3

! Off-diagonal terms
! 12
    call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp1)
    call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp2)

    p_wrk3d = tmp1 + tmp2 ! )*vis
    ! if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, wrk1d(1, 1), wrk1d(1, 2))
        call AVG_IK_V(nx, ny, nz, ny, u, wrk1d(1, 3), wrk1d(1, 2))
        call OPR_PARTIAL_Y(OPR_P1, 1, ny, 1, bcs, g(2), wrk1d(1, 3), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
            tmp1(:, j, :) = tmp1(:, j, :) - wrk1d(j, 2)
        end do
    ! end if
    eps = eps + p_wrk3d*(tmp1 + tmp2)

! 13 term
    call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp1)
    call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp2)

    p_wrk3d = tmp1 + tmp2 ! )*vis
    ! if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, wrk1d(1, 1), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
        end do
    ! end if
    eps = eps + p_wrk3d*(tmp1 + tmp2)

! 23 term
    call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
    call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)

    p_wrk3d = tmp1 + tmp2 ! )*vis
    ! if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, wrk1d(1, 1), wrk1d(1, 2))
        call AVG_IK_V(nx, ny, nz, ny, w, wrk1d(1, 3), wrk1d(1, 2))
        call OPR_PARTIAL_Y(OPR_P1, 1, ny, 1, bcs, g(2), wrk1d(1, 3), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
            tmp1(:, j, :) = tmp1(:, j, :) - wrk1d(j, 2)
        end do
    ! end if
    eps = eps + p_wrk3d*(tmp1 + tmp2)

    return
end subroutine FI_DISSIPATION

!########################################################################
! Reynolds fluctuations of array a
!########################################################################
subroutine FI_FLUCTUATION_INPLACE(nx, ny, nz, a)
    use TLab_Constants, only: wp, wi
    use Averages, only: AVG_IK

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(INOUT) :: a(nx, ny, nz)

    ! -------------------------------------------------------------------
    real(wp) dummy
    integer(wi) j

    ! ###################################################################
    do j = 1, ny
        dummy = AVG_IK(nx, ny, nz, j, a)
        a(:, j, :) = a(:, j, :) - dummy
    end do

    return
end subroutine FI_FLUCTUATION_INPLACE
