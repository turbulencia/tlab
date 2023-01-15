#include "dns_const.h"

!########################################################################
!#
!# remove divergence part of a vector field a=(u,v,w)
!#
!# Calculate scalar field phi s.t. lap phi = -div a, with BCs phi = 0
!# at top and bottom (i.e. zero tangential component of vector grad phi)
!# Then, add grad phi to vector a, where n grad phi = 0 at top and bottom.
!#
!# The BCs are such that a and a + grad phi are the same at top and bottom
!#
!########################################################################
subroutine FI_SOLENOIDAL(iwall, nx, ny, nz, u, v, w, tmp1, tmp2, tmp3)
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: g
    use TLAB_ARRAYS, only: wrk2d, wrk3d
    use TLAB_POINTERS_3D, only: p_wrk2d
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    implicit none

    integer(wi), intent(IN) :: iwall, nx, ny, nz
    real(wp), dimension(nx, ny, nz), intent(INOUT) :: u, v, w
    real(wp), dimension(nx, ny, nz), intent(INOUT) :: tmp1, tmp2, tmp3

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
    call FI_INVARIANT_P(nx, ny, nz, u, v, w, tmp1, tmp2, wrk2d, wrk3d)

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
