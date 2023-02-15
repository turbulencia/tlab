#include "dns_const.h"

!########################################################################
!#
!# Momentum equations, nonlinear term in divergence form and the
!# viscous term explicit. 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# Scalar needed for the buoyancy term
!#
!########################################################################
subroutine RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3()
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: g
    use TLAB_VARS, only: visc
    use TLAB_ARRAYS, only: q, wrk2d, wrk3d
    use TLAB_POINTERS, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    use DNS_ARRAYS, only: hq
    use TIME, only: dte
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS, only: BcsFlowJmin, BcsFlowJmax
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use BOUNDARY_BCS

    implicit none

! -----------------------------------------------------------------------
    integer(wi) iq, iq_max
    integer(wi) ibc, bcs(2, 2)
    real(wp) alpha

    real(wp), dimension(:, :, :), pointer :: p_bcs

! #######################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

! #######################################################################
! Diffusion and convection terms in momentum equations
! #######################################################################
    if (g(3)%size > 1) then
        iq_max = 3
    else
        iq_max = 2
    end if

    do iq = 1, iq_max
        call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), q(:, iq), tmp6, tmp3, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), q(:, iq), tmp5, tmp2, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), q(:, iq), tmp4, tmp1, wrk2d, wrk3d)
        hq(:, iq) = hq(:, iq) + visc*(tmp6 + tmp5 + tmp4)

        tmp6 = q(:, iq)*w
        tmp5 = q(:, iq)*v
        tmp4 = q(:, iq)*u
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp6, tmp3, wrk3d, wrk2d, wrk3d)
!  For iq = 2, i.e., in the wall-normal direction
!  BCs s.t. to make sure that product vd/dy(v) at the boundary is zero because v is zero.
!  bcs_loc = 1
!  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs_loc, g(2), tmp5, tmp2, wrk3d, wrk2d,wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp5, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1, wrk3d, wrk2d, wrk3d)
        hq(:, iq) = hq(:, iq) - (tmp3 + tmp2 + tmp1)

    end do

! -----------------------------------------------------------------------
! Dilatation term for Bcs
! -----------------------------------------------------------------------
!  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
!  DO k = 1,kmax
!     DO i = 1,imax
!        ij = i                 + imax*jmax*(k-1) ! bottom
!        h1(ij) = h1(ij) + u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
!        h1(ij) = h1(ij) + u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     ENDDO
!  ENDDO

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
    if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
        call BOUNDARY_BUFFER_RELAX_FLOW()
    end if

! #######################################################################
! Pressure term
! #######################################################################
! -----------------------------------------------------------------------
! Poisson equation
! -----------------------------------------------------------------------
    alpha = 1.0_wp/dte
    tmp1 = hq(:, 1) + alpha*u
    tmp2 = hq(:, 2) + alpha*v
    tmp3 = hq(:, 3) + alpha*w
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp4, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp5, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp6, wrk3d, wrk2d, wrk3d)

! forcing term
    tmp1 = tmp6 + tmp5 + tmp4

! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, 2)
    BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)
    BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)

! pressure in tmp1, Oy derivative in tmp3
    ibc = 3
    call OPR_POISSON_FXZ(imax, jmax, kmax, g, ibc, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)

! horizontal derivatives
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4, wrk3d, wrk2d, wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
    hq(:, 1) = hq(:, 1) - tmp2
    hq(:, 2) = hq(:, 2) - tmp3
    hq(:, 3) = hq(:, 3) - tmp4

! #######################################################################
! Boundary conditions
! #######################################################################
    BcsFlowJmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsFlowJmax%ref = 0.0_wp ! Scalar BCs initialized at start of routine

    do iq = 1, iq_max
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hq(1, iq), &
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), tmp1)
        end if

        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, iq)
        p_bcs(:, 1, :) = BcsFlowJmin%ref(:, :, iq)
        p_bcs(:, jmax, :) = BcsFlowJmax%ref(:, :, iq)

    end do

    return
end subroutine RHS_FLOW_GLOBAL_INCOMPRESSIBLE_3
