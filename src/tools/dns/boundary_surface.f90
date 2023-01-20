#include "types.h"
#include "dns_const.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE BOUNDARY_SURFACE_J
! Calculates and updates interactive surface boundary condition
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BOUNDARY_SURFACE_J(is, bcs, s, hs, tmp1, tmp2, aux)
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_CONSTANTS, only: lfile
    use TLAB_VARS, only: imax, jmax, kmax, g
    use TLAB_VARS, only: isize_field
    use TLAB_VARS, only: visc, schmidt
    use TLAB_ARRAYS, only: wrk2d, wrk3d
    use BOUNDARY_BCS, only: BcsScalJmin, BcsScalJmax
    use AVGS, only: AVG1V2D
    use OPR_PARTIAL

    implicit none

#include "integers.h"

    TINTEGER is
    TINTEGER, dimension(2, 2), intent(IN) :: bcs          ! Boundary conditions from derivative operator
    TREAL, dimension(isize_field, *) :: s, hs
    TREAL, dimension(isize_field) :: tmp1, tmp2
    TREAL, dimension(imax, kmax, 6), target :: aux

    TINTEGER nxy, ip, k
    TREAL, dimension(:, :), pointer :: hfx, hfx_anom
    TREAL :: diff, hfx_avg

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING SUBROUTINE BOUNDARY_SURFACE_J')
#endif
    diff = visc/schmidt(is)
    nxy = imax*jmax

    ! vertical derivative of scalar for flux at the boundaries
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:, is), tmp1, wrk3d, wrk2d, wrk3d)

    ! ------------------------------------------------------------
    ! Bottom Boundary
    ! ------------------------------------------------------------
    if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) then
        hfx => aux(:, :, 1)
        hfx_anom => aux(:, :, 2)
        ip = 1
        do k = 1, kmax    ! Calculate the surface flux
            hfx(:, k) = diff*tmp1(ip:ip + imax - 1); ip = ip + nxy
        end do
        hfx_avg = diff*AVG1V2D(imax, jmax, kmax, 1, 1, tmp1)
        hfx_anom = hfx - hfx_avg
        BcsScalJmin%ref(:, :, is) = BcsScalJmin%ref(:, :, is) + BcsScalJmin%cpl(is)*hfx_anom
    end if

    ! ------------------------------------------------------------
    ! Top Boundary
    ! ------------------------------------------------------------
    if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) then
        hfx => aux(:, :, 3)
        hfx_anom => aux(:, :, 4)
        ip = imax*(jmax - 1) + 1
        do k = 1, kmax; ! Calculate the surface flux
            hfx(:, k) = -diff*tmp1(ip:ip + imax - 1); ip = ip + nxy; 
        end do
        hfx_avg = diff*AVG1V2D(imax, jmax, kmax, 1, 1, tmp1)
        hfx_anom = hfx - hfx_avg
        BcsScalJmax%ref(:, :, is) = BcsScalJmax%ref(:, :, is) + BcsScalJmax%cpl(is)*hfx_anom
    end if

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(TFILE, 'LEAVING SUBROUTINE BOUNDAR_SURFACE_J')
#endif

    return

end subroutine BOUNDARY_SURFACE_J
