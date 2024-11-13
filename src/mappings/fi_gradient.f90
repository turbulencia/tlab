#include "dns_const.h"

!########################################################################
!# Calculate the magnitude of the scalar gradient as given
!# by G_i G_i, where G_i is ds/dx_i
!# and terms in its evolution equation.
!########################################################################
module FI_GRADIENT_EQN
    use TLab_Constants, only: wp, wi
    use FDM, only: g
    use IBM_VARS, only: imode_ibm, ibm_partial
    use OPR_PARTIAL
    implicit none
    private

    integer(wi), parameter :: bcs(2, 2) = 0

    public :: FI_GRADIENT
    public :: FI_GRADIENT_PRODUCTION
    public :: FI_GRADIENT_DIFFUSION

contains
!########################################################################
! Calculate the magnitude of the scalar gradient
!########################################################################
    subroutine FI_GRADIENT(nx, ny, nz, s, result, tmp1)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)

! ###################################################################
        ! IBM   (if .true., OPR_PARTIAL_X/Y/Z uses modified fields for derivatives)
        if (imode_ibm == 1) ibm_partial = .true.

        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), s, result)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), s, tmp1)
        result = result*result + tmp1*tmp1
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), s, tmp1)
        result = result + tmp1*tmp1

        if (imode_ibm == 1) then
            ibm_partial = .false.
            call IBM_BCS_FIELD(result)
        end if

        return
    end subroutine FI_GRADIENT

!########################################################################
! Calculate the scalar gradient production term as given by -(G_i G_j s_ij)
!########################################################################
    subroutine FI_GRADIENT_PRODUCTION(nx, ny, nz, s, u, v, w, result, grad_x, grad_y, grad_z, tmp1, tmp2)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(in) :: s, u, v, w
        real(wp), dimension(nx*ny*nz), intent(out) :: result
        real(wp), dimension(nx*ny*nz), intent(inout) :: grad_x, grad_y, grad_z
        real(wp), dimension(nx*ny*nz), intent(inout) :: tmp1, tmp2

! ###################################################################
! Vorticity vector
! ###################################################################
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), s, grad_x)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), s, grad_y)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), s, grad_z)

! ###################################################################
! Production term
! ###################################################################
! Ux, Vy, Wz
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
        result = tmp1*grad_x*grad_x + tmp2*grad_y*grad_y

        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp2)
        result = result + tmp2*grad_z*grad_z

! Uy, Vx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        result = result + (tmp1 + tmp2)*grad_x*grad_y

! Uz, Wx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp2)
        result = result + (tmp1 + tmp2)*grad_x*grad_z

! Vz, Wy
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        result = -(result + (tmp1 + tmp2)*grad_y*grad_z)

        return
    end subroutine FI_GRADIENT_PRODUCTION

!########################################################################
! Calculate the gradient diffusion term as given by G_i lap G_i
! The diffusivity D is not multiplied here.
!########################################################################
    subroutine FI_GRADIENT_DIFFUSION(nx, ny, nz, s, result, grad, tmp1, tmp2, tmp3, tmp4)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(in) :: s
        real(wp), dimension(nx*ny*nz), intent(out) :: result
        real(wp), dimension(nx*ny*nz), intent(inout) :: grad
        real(wp), dimension(nx*ny*nz), intent(inout) :: tmp1, tmp2, tmp3, tmp4

! ###################################################################
! G_x
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), s, grad)

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), grad, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), grad, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), grad, tmp1, tmp4)
        result = grad*(tmp1 + tmp2 + tmp3)

! G_y
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), s, grad)

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), grad, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), grad, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), grad, tmp1, tmp4)
        result = result + grad*(tmp1 + tmp2 + tmp3)

! G_z
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), s, grad)

        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), grad, tmp3, tmp4)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), grad, tmp2, tmp4)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), grad, tmp1, tmp4)
        result = result + grad*(tmp1 + tmp2 + tmp3)

        return
    end subroutine FI_GRADIENT_DIFFUSION

end module FI_GRADIENT_EQN
