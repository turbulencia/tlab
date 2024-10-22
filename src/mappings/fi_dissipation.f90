#include "dns_const.h"

!########################################################################
!# Calculates turbulent dissipation per unit volume \rho \epsilon = \tau'_{ij}u'_{ij}
!# It assumes constant visocsity
!########################################################################
subroutine FI_DISSIPATION(flag, nx, ny, nz, u, v, w, eps, tmp1, tmp2, tmp3, tmp4)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: g
    use TLAB_VARS, only: area, visc
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: p_wrk3d
    use AVGS, only: AVG_IK_V
    use OPR_PARTIAL
    implicit none

    integer(wi), intent(IN) :: flag ! 0 for tau_ji  * u_i,j
    ! 1 for tau'_ij * u'_i,j
    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
    real(wp), dimension(nx, ny, nz), intent(OUT) :: eps
    real(wp), dimension(nx, ny, nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

! -------------------------------------------------------------------
    integer(wi) bcs(2, 2), j
    integer, parameter :: i1 = 1

! ###################################################################
    bcs = 0

! Diagonal terms
    call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp3)
    tmp4 = (tmp1 + tmp2 + tmp3)*2.0_wp/3.0_wp

! 11
    p_wrk3d = 2.0_wp*tmp1 - tmp4 ! )*vis
    if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 1), wrk1d(1, 2), area)
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
        end do
    end if
    eps = p_wrk3d*tmp1

! 22
    p_wrk3d = 2.0_wp*tmp2 - tmp4 ! )*vis
    if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 1), wrk1d(1, 2), area)
        call AVG_IK_V(nx, ny, nz, ny, v, g(1)%jac, g(3)%jac, wrk1d(1, 3), wrk1d(1, 2), area)
        call OPR_PARTIAL_Y(OPR_P1, i1, ny, i1, bcs, g(2), wrk1d(1, 3), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
            tmp2(:, j, :) = tmp2(:, j, :) - wrk1d(j, 2)
        end do
    end if
    eps = eps + p_wrk3d*tmp2

! 33
    p_wrk3d = 2.0_wp*tmp3 - tmp4 ! )*vis
    if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 1), wrk1d(1, 2), area)
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
        end do
    end if
    eps = eps + p_wrk3d*tmp3

! Off-diagonal terms
! 12
    call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp1)
    call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp2)

    p_wrk3d = tmp1 + tmp2 ! )*vis
    if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 1), wrk1d(1, 2), area)
        call AVG_IK_V(nx, ny, nz, ny, u, g(1)%jac, g(3)%jac, wrk1d(1, 3), wrk1d(1, 2), area)
        call OPR_PARTIAL_Y(OPR_P1, i1, ny, i1, bcs, g(2), wrk1d(1, 3), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
            tmp1(:, j, :) = tmp1(:, j, :) - wrk1d(j, 2)
        end do
    end if
    eps = eps + p_wrk3d*(tmp1 + tmp2)

! 13 term
    call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp1)
    call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp2)

    p_wrk3d = tmp1 + tmp2 ! )*vis
    if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 1), wrk1d(1, 2), area)
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
        end do
    end if
    eps = eps + p_wrk3d*(tmp1 + tmp2)

! 23 term
    call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
    call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)

    p_wrk3d = tmp1 + tmp2 ! )*vis
    if (flag == 1) then
        call AVG_IK_V(nx, ny, nz, ny, p_wrk3d, g(1)%jac, g(3)%jac, wrk1d(1, 1), wrk1d(1, 2), area)
        call AVG_IK_V(nx, ny, nz, ny, w, g(1)%jac, g(3)%jac, wrk1d(1, 3), wrk1d(1, 2), area)
        call OPR_PARTIAL_Y(OPR_P1, i1, ny, i1, bcs, g(2), wrk1d(1, 3), wrk1d(1, 2))
        do j = 1, ny
            p_wrk3d(:, j, :) = p_wrk3d(:, j, :) - wrk1d(j, 1)
            tmp1(:, j, :) = tmp1(:, j, :) - wrk1d(j, 2)
        end do
    end if
    eps = eps + p_wrk3d*(tmp1 + tmp2)

! Final calculation
    eps = eps*visc

    return
end subroutine FI_DISSIPATION

! #######################################################################
! Calculate kinetic energy of fluctuating field per unit volume
! #######################################################################
#define rR(j)     wrk1d(j,1)
#define fU(j)     wrk1d(j,2)
#define fV(j)     wrk1d(j,3)
#define fW(j)     wrk1d(j,4)
#define aux(j)    wrk1d(j,5)

subroutine FI_RTKE(nx, ny, nz, q, ke)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: imode_eqns, inb_flow
    use TLAB_VARS, only: g, area
    use TLab_Arrays, only: wrk1d
    use THERMO_ANELASTIC, only : rbackground
    use AVGS, only: AVG_IK_V

    implicit none

    integer(wi) nx, ny, nz
    real(wp), intent(in) :: q(nx, ny, nz, inb_flow)
    real(wp), intent(out) :: ke(nx, ny, nz)

    ! -----------------------------------------------------------------------
    integer(wi) j

    ! #######################################################################
    select case (imode_eqns)
    case (DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL)
        call AVG_IK_V(nx, ny, nz, ny, q(1, 1, 1, 5), g(1)%jac, g(3)%jac, rR(1), aux(1), area)

        ke = q(:, :, :, 5)*q(:, :, :, 1)
        call AVG_IK_V(nx, ny, nz, ny, ke, g(1)%jac, g(3)%jac, fU(1), aux(1), area)
        fU(:) = fU(:)/rR(:)

        ke = q(:, :, :, 5)*q(:, :, :, 2)
        call AVG_IK_V(nx, ny, nz, ny, ke, g(1)%jac, g(3)%jac, fV(1), aux(1), area)
        fV(:) = fV(:)/rR(:)

        ke = q(:, :, :, 5)*q(:, :, :, 3)
        call AVG_IK_V(nx, ny, nz, ny, ke, g(1)%jac, g(3)%jac, fW(1), aux(1), area)
        fW(:) = fW(:)/rR(:)

    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        rR(:) = rbackground(:)
        call AVG_IK_V(nx, ny, nz, ny, q(1, 1, 1, 1), g(1)%jac, g(3)%jac, fU(1), aux(1), area)
        call AVG_IK_V(nx, ny, nz, ny, q(1, 1, 1, 2), g(1)%jac, g(3)%jac, fV(1), aux(1), area)
        call AVG_IK_V(nx, ny, nz, ny, q(1, 1, 1, 3), g(1)%jac, g(3)%jac, fW(1), aux(1), area)

    end select

    do j = 1, ny
        ke(:, j, :) = 0.5_wp*rR(j)*((q(:, j, :, 1) - fU(j))**2 + (q(:, j, :, 2) - fV(j))**2 + (q(:, j, :, 3) - fW(j))**2)
    end do

    return
end subroutine FI_RTKE

!########################################################################
! Reynolds fluctuations of array a
!########################################################################
subroutine FI_FLUCTUATION_INPLACE(nx, ny, nz, a)
    use TLab_Constants, only: wp, wi
    use TLAB_VARS, only: g, area
    use AVGS, only: AVG_IK

    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(INOUT) :: a(nx, ny, nz)

    ! -------------------------------------------------------------------
    real(wp) dummy
    integer(wi) j

    ! ###################################################################
    do j = 1, ny
        dummy = AVG_IK(nx, ny, nz, j, a, g(1)%jac, g(3)%jac, area)
        a(:, j, :) = a(:, j, :) - dummy
    end do

    return
end subroutine FI_FLUCTUATION_INPLACE
