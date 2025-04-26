#include "dns_const.h"

!########################################################################
!# Calculate the magnitude of the vorticity as given by w_i w_i:
!# (v,x-u,y)^2 + (u,z-w,x)^2 + (w,y-v,z)^2
!# and terms in its evolution equation.
!########################################################################
module FI_VORTICITY_EQN
    use TLab_Constants, only: wp, wi
    use FDM, only: g
    use IBM_VARS, only: imode_ibm, ibm_partial
    use OPR_Partial
    implicit none
    private

    integer(wi), parameter :: bcs(2, 2) = 0

    public :: FI_VORTICITY
    public :: FI_VORTICITY_PRODUCTION
    public :: FI_VORTICITY_DIFFUSION
    public :: FI_VORTICITY_BAROCLINIC

contains

!########################################################################
! Calculate the magnitude of the vorticity
!########################################################################
    subroutine FI_VORTICITY(nx, ny, nz, u, v, w, result, tmp1, tmp2)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2

! ###################################################################
        ! IBM   (if .true., OPR_Partial_X/Y/Z uses modified fields for derivatives)
        if (imode_ibm == 1) ibm_partial = .true.

! v,x-u,y
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        result = (tmp1 - tmp2)*(tmp1 - tmp2)

! u,z-w,x
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp1)
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp2)
        result = result + (tmp1 - tmp2)*(tmp1 - tmp2)

! w,y-v,z
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        result = result + (tmp1 - tmp2)*(tmp1 - tmp2)

        if (imode_ibm == 1) then
            ibm_partial = .false.
            call IBM_BCS_FIELD(result)
        end if

        return
    end subroutine FI_VORTICITY

!########################################################################
! Calculate the vorticity production term as given by w_i w_j s_ij
!########################################################################
    subroutine FI_VORTICITY_PRODUCTION(nx, ny, nz, u, v, w, result, vort_x, vort_y, vort_z, tmp1, tmp2)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: vort_x, vort_y, vort_z
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2

! ###################################################################
! Vorticity vector
! ###################################################################
! v,x-u,y
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        vort_z = tmp1 - tmp2

! u,z-w,x
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp1)
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp2)
        vort_y = tmp1 - tmp2

! w,y-v,z
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        vort_x = tmp1 - tmp2

! ###################################################################
! Production term
! ###################################################################
! Ux, Vy, Wz
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1)
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2)
        result = tmp1*vort_x*vort_x + tmp2*vort_y*vort_y

        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp2)
        result = result + tmp2*vort_z*vort_z

! Uy, Vx
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        result = result + (tmp1 + tmp2)*vort_x*vort_y

! Uz, Wx
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1)
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp2)
        result = result + (tmp1 + tmp2)*vort_x*vort_z

! Vz, Wy
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        result = result + (tmp1 + tmp2)*vort_y*vort_z

        return
    end subroutine FI_VORTICITY_PRODUCTION

!########################################################################
! Calculate the vorticity diffusion term as given by w_i lap w_i
! The kinematic viscosity \nu is not multiplied here.
!########################################################################
    subroutine FI_VORTICITY_DIFFUSION(nx, ny, nz, u, v, w, result, vort, tmp1, tmp2, tmp3, tmp4)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: vort
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

! ###################################################################
! -------------------------------------------------------------------
! W_z = V,x - U,y
! -------------------------------------------------------------------
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1)
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2)
        vort = tmp1 - tmp2

        call OPR_Partial_Z(OPR_P2, nx, ny, nz, bcs, g(3), vort, tmp3, tmp4)
        call OPR_Partial_Y(OPR_P2, nx, ny, nz, bcs, g(2), vort, tmp2, tmp4)
        call OPR_Partial_X(OPR_P2, nx, ny, nz, bcs, g(1), vort, tmp1, tmp4)
        result = vort*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! W_y = U,z - W,x
! -------------------------------------------------------------------
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp1)
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp2)
        vort = tmp1 - tmp2

        call OPR_Partial_Z(OPR_P2, nx, ny, nz, bcs, g(3), vort, tmp3, tmp4)
        call OPR_Partial_Y(OPR_P2, nx, ny, nz, bcs, g(2), vort, tmp2, tmp4)
        call OPR_Partial_X(OPR_P2, nx, ny, nz, bcs, g(1), vort, tmp1, tmp4)
        result = result + vort*(tmp1 + tmp2 + tmp3)

! -------------------------------------------------------------------
! W_z = W,y - V,z
! -------------------------------------------------------------------
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1)
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2)
        vort = tmp1 - tmp2

        call OPR_Partial_Z(OPR_P2, nx, ny, nz, bcs, g(3), vort, tmp3, tmp4)
        call OPR_Partial_Y(OPR_P2, nx, ny, nz, bcs, g(2), vort, tmp2, tmp4)
        call OPR_Partial_X(OPR_P2, nx, ny, nz, bcs, g(1), vort, tmp1, tmp4)
        result = result + vort*(tmp1 + tmp2 + tmp3)

        return
    end subroutine FI_VORTICITY_DIFFUSION

!########################################################################
! Calculate the baroclinic production of vorticity
! (grad\rho x grad p)/\rho^2
!########################################################################
    subroutine FI_VORTICITY_BAROCLINIC(nx, ny, nz, r, p, result, tmp1, tmp2)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: r, p
        real(wp), dimension(nx*ny*nz, 3), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2

! ###################################################################
! r,yp,z-r,zp,y
! ###################################################################
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), r, result(1, 1))
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), p, tmp1)

        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), p, result(1, 2))
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), r, tmp2)

        result(:, 1) = (result(:, 1)*tmp1(:) - result(:, 2)*tmp2(:))/(r(:)*r(:))

! ###################################################################
! r,zp,x-r,xp,z
! p,z and r,z are in tmp1 and tmp2 respectively
! ###################################################################
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), p, result(1, 3))
        result(:, 2) = result(:, 3)*tmp2(:)

        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), r, tmp2)
        result(:, 2) = (result(:, 2) - tmp1(:)*tmp2(:))/(r(:)*r(:))

! ###################################################################
! r,xp,y-r,yp,x
! p,x and r,x are in result3 and tmp2 respectively
! ###################################################################
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), r, tmp1)
        result(:, 3) = result(:, 3)*tmp1(:)

        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), p, tmp1)
        result(:, 3) = (tmp1(:)*tmp2(:) - result(:, 3))/(r(:)*r(:))

        return
    end subroutine FI_VORTICITY_BAROCLINIC

end module FI_VORTICITY_EQN
