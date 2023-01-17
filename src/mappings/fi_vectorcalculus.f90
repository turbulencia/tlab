#include "dns_const.h"

module FI_VECTORCALCULUS
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imode_ibm
    use TLAB_ARRAYS, only: wrk2d, wrk3d
    use IBM_VARS, only: ibm_partial
    use OPR_PARTIAL
    implicit none
    private

    public :: FI_CURL
    public :: FI_SOLENOIDAL
    public :: FI_INVARIANT_P, FI_INVARIANT_P_STAG, FI_INVARIANT_Q, FI_INVARIANT_R
    public :: FI_ISOSURFACE_ANGLE, FI_ISOSURFACE_CURVATURE

contains
!########################################################################
! Calculate the curl of the vector (u,v,w) in Cartesian coordinates
!########################################################################
    subroutine FI_CURL(nx, ny, nz, u, v, w, wx, wy, wz, tmp)
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz), v(nx*ny*nz), w(nx*ny*nz)
        real(wp), intent(out) :: wx(nx*ny*nz), wy(nx*ny*nz), wz(nx*ny*nz), tmp(nx*ny*nz)

! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)

! ###################################################################
        bcs = 0

! IBM   (if .true., OPR_PARTIAL_X/Y/Z uses modified fields for derivatives)
        if (imode_ibm == 1) ibm_partial = .true.

! v,x-u,y
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, wz, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp, wrk3d, wrk2d, wrk3d)
        wz = wz - tmp

! u,z-w,x
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, wy, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp, wrk3d, wrk2d, wrk3d)
        wy = wy - tmp

! w,y-v,z
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, wx, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp, wrk3d, wrk2d, wrk3d)
        wx = wx - tmp

        tmp = wx*wx + wy*wy + wz*wz

        if (imode_ibm == 1) then
            ibm_partial = .false.
            call IBM_BCS_FIELD(wx)
            call IBM_BCS_FIELD(wy)
            call IBM_BCS_FIELD(wz)
            call IBM_BCS_FIELD(tmp)
        end if

        return
    end subroutine FI_CURL

!########################################################################
!# remove divergence part of a vector field a=(u,v,w)
!#
!# Calculate scalar field phi s.t. lap phi = -div a, with BCs phi = 0
!# at top and bottom (i.e. zero tangential component of vector grad phi)
!# Then, add grad phi to vector a, where n grad phi = 0 at top and bottom.
!#
!# The BCs are such that a and a + grad phi are the same at top and bottom
!#
!########################################################################
    subroutine FI_SOLENOIDAL(nx, ny, nz, u, v, w, tmp1, tmp2, tmp3)
        use TLAB_POINTERS_3D, only: p_wrk2d
        use OPR_ELLIPTIC

        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), intent(INOUT) :: u(nx, ny, nz), v(nx, ny, nz), w(nx, ny, nz)
        real(wp), intent(INOUT) :: tmp1(nx, ny, nz), tmp2(nx, ny, nz), tmp3(nx, ny, nz)

! -------------------------------------------------------------------
        integer(wi) ibc, bcs(2, 2)

! ###################################################################
        bcs = 0

!  IF ( iwall .EQ. 1) THEN; ibc = 4
!  ELSE;                    ibc = 3; ENDIF
        ibc = 3

! -------------------------------------------------------------------
! Solve lap(phi) = - div(u)
! -------------------------------------------------------------------
        call FI_INVARIANT_P(nx, ny, nz, u, v, w, tmp1, tmp2)

        if (g(1)%periodic .and. g(3)%periodic) then ! Doubly periodic in xOz
            p_wrk2d(:, :, 1:2) = 0.0_wp  ! bcs
            call OPR_POISSON_FXZ(nx, ny, nz, g, ibc, tmp1, tmp2, tmp3, p_wrk2d(:, :, 1), p_wrk2d(:, :, 2))

        else                                          ! General treatment
#ifdef USE_CGLOC
! Need to define global variable with ipos,jpos,kpos,ci,cj,ck,
            tmp2 = -tmp1            ! change of forcing term sign
            call CGPOISSON(i1, nx, ny, nz, g(3)%size, tmp1, tmp2, tmp3, tmp4, ipos, jpos, kpos, ci, cj, ck, wrk2d)
#endif
        end if

! -------------------------------------------------------------------
! Eliminate solenoidal part of u by adding grad(phi)
! -------------------------------------------------------------------
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
        u = u + tmp2
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
        v = v + tmp2
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
        w = w + tmp2

        return
    end subroutine FI_SOLENOIDAL

!########################################################################
! First invariant of the velocity gradient tensor, div u
!########################################################################
    subroutine FI_INVARIANT_P(nx, ny, nz, u, v, w, result, tmp1)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)

        ! ###################################################################
        bcs = 0

        ! -------------------------------------------------------------------
        ! IBM   (if .true., OPR_PARTIAL_X/Y/Z uses modified fields for derivatives)
        if (imode_ibm == 1) ibm_partial = .true.

        ! -------------------------------------------------------------------

        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, result, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp1, wrk3d, wrk2d, wrk3d)

        result = result + tmp1
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp1, wrk3d, wrk2d, wrk3d)

        result = -(result + tmp1)

        ! -------------------------------------------------------------------
        if (imode_ibm == 1) ibm_partial = .false.

        return
    end subroutine FI_INVARIANT_P

    !########################################################################
    ! First invariant on horizontal pressure nodes
    ! (caution: div(u)=0 condition only holds on pressure nodes,
    !           in IBM-mode without splines, because of combined schemes)
    !########################################################################
    subroutine FI_INVARIANT_P_STAG(nx, ny, nz, u, v, w, result, tmp1, tmp2)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)

        ! ###################################################################
        bcs = 0

        ! dudx
        call OPR_PARTIAL_X(OPR_P1_INT_VP, nx, ny, nz, bcs, g(1), u, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, nx, ny, nz, bcs, g(3), tmp1, result, wrk3d, wrk2d, wrk3d)
        ! dvdy
        call OPR_PARTIAL_X(OPR_P0_INT_VP, nx, ny, nz, bcs, g(1), v, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, nx, ny, nz, bcs, g(3), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
        result = result + tmp1
        ! dwdz
        call OPR_PARTIAL_X(OPR_P0_INT_VP, nx, ny, nz, bcs, g(1), w, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1_INT_VP, nx, ny, nz, bcs, g(3), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
        result = -(result + tmp1)

        return
    end subroutine FI_INVARIANT_P_STAG

    !########################################################################
    ! Second invariant of the velocity gradient tensor, Q, as defined
    ! by Chong et al. (1990). Regions with high positive values of Q are
    ! vorticity dominated, regions with high negative values of Q are
    ! strain dominated.
    ! If incompressible, P=0 and then Q=(w^2-2s^2)/4, where P is the first
    ! invariant, P=-div(v).
    !########################################################################
    subroutine FI_INVARIANT_Q(nx, ny, nz, u, v, w, result, tmp1, tmp2, tmp3)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2, tmp3

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)

        ! ###################################################################
        bcs = 0

        ! -------------------------------------------------------------------
        ! diagonal terms
        ! -------------------------------------------------------------------
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp3, wrk3d, wrk2d, wrk3d)
        result = tmp1*tmp2 + tmp2*tmp3 + tmp3*tmp1

        ! -------------------------------------------------------------------
        ! off-diagonal terms
        ! -------------------------------------------------------------------
        ! UyVx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp2, wrk3d, wrk2d, wrk3d)
        result = result - tmp1*tmp2

        ! UzWx
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp2, wrk3d, wrk2d, wrk3d)
        result = result - tmp1*tmp2

        ! WyVz
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp2, wrk3d, wrk2d, wrk3d)
        result = result - tmp1*tmp2

        return
    end subroutine FI_INVARIANT_Q

    !########################################################################
    ! Calculate third invariant of the velocity gradient tensor (-determinant)
    ! The derivatives v,y and v,z are repeated to avoid more arrays
    !########################################################################
    subroutine FI_INVARIANT_R(nx, ny, nz, u, v, w, result, tmp1, tmp2, tmp3, tmp4, tmp5)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: u, v, w
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4, tmp5

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)

        ! ###################################################################
        bcs = 0

        ! ###################################################################
        ! Term u,x (v,y w,z-w,y v,z)
        ! ###################################################################
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), w, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), w, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp3, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp4, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), u, tmp5, wrk3d, wrk2d, wrk3d)
        result = tmp5*(tmp1*tmp3 - tmp2*tmp4)

        ! ###################################################################
        ! Term v,x (u,z w,y-w,z u,y)
        ! ###################################################################
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), u, tmp3, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), u, tmp4, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), v, tmp5, wrk3d, wrk2d, wrk3d)
        result = result + tmp5*(tmp2*tmp3 - tmp1*tmp4)

        ! ###################################################################
        ! Term v,x (u,z w,y-w,z u,y)
        ! ###################################################################
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), v, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), v, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), w, tmp5, wrk3d, wrk2d, wrk3d)
        result = result + tmp5*(tmp1*tmp4 - tmp2*tmp3)

        result = -result ! set the right sign

        return
    end subroutine FI_INVARIANT_R

    !########################################################################
    !# Calculate the angle between isosurfaces of a and b at each point;
    !# the gradient vectors are calculated first and then the angle between them.
    !########################################################################
    subroutine FI_ISOSURFACE_ANGLE(nx, ny, nz, a, b, result, tmp1, tmp2, tmp3, tmp4)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: a, b
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2), ij

        ! ###################################################################
        bcs = 0

        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), a, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), b, tmp2, wrk3d, wrk2d, wrk3d)
        result = tmp1*tmp2
        tmp3 = tmp1*tmp1
        tmp4 = tmp2*tmp2

        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), a, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), b, tmp2, wrk3d, wrk2d, wrk3d)
        result = result + tmp1*tmp2
        tmp3 = tmp3 + tmp1*tmp1
        tmp4 = tmp4 + tmp2*tmp2

        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), a, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), b, tmp2, wrk3d, wrk2d, wrk3d)
        result = result + tmp1*tmp2
        tmp3 = tmp3 + tmp1*tmp1
        tmp4 = tmp4 + tmp2*tmp2

        do ij = 1, nx*ny*nz
            if (tmp3(ij) > 0.0_wp .and. tmp4(ij) > 0.0_wp) then
                result(ij) = result(ij)/sqrt(tmp3(ij)*tmp4(ij))
            end if
        end do

        return
    end subroutine FI_ISOSURFACE_ANGLE

    !########################################################################
    !# Calculate the curvature \kappa of a scalar field s and -div(n), where n
    !# is the normal vector n = grad(s)/|grad(s)|.
    !# it is expanded and calculated as
    !# \kappa=-( lap(s) - n*grad( n*grad(s) ) )/|grad(s)| to use the second-derivative FD operator
    !########################################################################
    subroutine FI_ISOSURFACE_CURVATURE(nx, ny, nz, s, result, tmp1, tmp2, tmp3, tmp4)
        integer(wi), intent(IN) :: nx, ny, nz
        real(wp), dimension(nx*ny*nz), intent(IN) :: s
        real(wp), dimension(nx*ny*nz), intent(OUT) :: result
        real(wp), dimension(nx*ny*nz), intent(INOUT) :: tmp1, tmp2, tmp3, tmp4 ! tmp1 contains G\kappa

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)

        ! ###################################################################
        bcs = 0

        ! -------------------------------------------------------------------
        ! |grad(s)|
        ! -------------------------------------------------------------------
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), s, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), s, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), s, tmp3, wrk3d, wrk2d, wrk3d)
        tmp4 = sqrt(tmp1*tmp1 + tmp2*tmp2 + tmp3*tmp3)

        ! -------------------------------------------------------------------
        ! first derivative terms
        ! -------------------------------------------------------------------
        call OPR_PARTIAL_X(OPR_P1, nx, ny, nz, bcs, g(1), tmp4, result, wrk3d, wrk2d, wrk3d)
        result = result*tmp1
        ! tmp1 is now free

        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp4, tmp1, wrk3d, wrk2d, wrk3d)
        result = result + tmp1*tmp2
        ! tmp2 is now free

        call OPR_PARTIAL_Z(OPR_P1, nx, ny, nz, bcs, g(3), tmp4, tmp2, wrk3d, wrk2d, wrk3d)
        result = (result + tmp2*tmp3)/tmp4
        ! tmp3 is now free

        ! -------------------------------------------------------------------
        ! second derivative terms and final calculations
        ! -------------------------------------------------------------------
        call OPR_PARTIAL_Z(OPR_P2, nx, ny, nz, bcs, g(3), s, tmp1, tmp3, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P2, nx, ny, nz, bcs, g(2), s, tmp2, tmp3, wrk2d, wrk3d)
        result = result - (tmp1 + tmp2)
        call OPR_PARTIAL_X(OPR_P2, nx, ny, nz, bcs, g(1), s, tmp1, tmp3, wrk2d, wrk3d)
        tmp1 = result - tmp1
        result = tmp1/tmp4

        return
    end subroutine FI_ISOSURFACE_CURVATURE

end module FI_VECTORCALCULUS
