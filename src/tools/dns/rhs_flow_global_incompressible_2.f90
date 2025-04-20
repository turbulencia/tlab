#include "dns_const.h"

!########################################################################
!#
!# Momentum equations, nonlinear term in skew-symmetric form and the
!# viscous term explicit. 9 2nd order + 9+9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# Scalar needed for the buoyancy term
!#
!########################################################################
subroutine RHS_FLOW_GLOBAL_INCOMPRESSIBLE_2()
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax
    use FDM, only: g
    use NavierStokes, only: visc
    use TLab_WorkFlow, only: stagger_on
    use TLab_Arrays, only: q
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    use DNS_ARRAYS, only: hq
    use TIME, only: dte
    use IBM_VARS, only: imode_ibm, ibm_partial
    use OPR_PARTIAL
    use OPR_Elliptic
    use BOUNDARY_BUFFER
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
! Preliminaries for IBM use
! (OPR_PARTIAL_X/Y/Z uses modified fields for derivatives)
! #######################################################################
    if (imode_ibm == 1) ibm_partial = .true.

! #######################################################################
! Diffusion and convection terms in momentum equations
! #######################################################################
    if (g(3)%size > 1) then
        iq_max = 3
    else
        iq_max = 2
    end if

    do iq = 1, iq_max
        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), q(:, iq), tmp6, tmp3)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), q(:, iq), tmp5, tmp2)
        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), q(:, iq), tmp4, tmp1)
        hq(:, iq) = hq(:, iq) + visc*(tmp6 + tmp5 + tmp4) - 0.5_wp*(w*tmp3 + v*tmp2 + u*tmp1)

        ! Adding the conservative-form terms to complete skew-symmetric formulation
        tmp1 = u*q(:, iq)
        tmp2 = v*q(:, iq)
        tmp3 = w*q(:, iq)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp5)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp6)
        hq(:, iq) = hq(:, iq) - 0.5_wp*(tmp6 + tmp5 + tmp4)

    end do

! -----------------------------------------------------------------------
! In case there is dilatation
! -----------------------------------------------------------------------
!  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1)
!  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2)
!  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3)
!  DO k = 1,kmax
!     DO i = 1,imax
!        ij = i                 + imax*jmax*(k-1) ! bottom
!        h1(ij) = h1(ij) + 0.5_wp*u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + 0.5_wp*w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        ij = i + imax*(jmax-1) + imax*jmax*(k-1) ! top
!        h1(ij) = h1(ij) + 0.5_wp*u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!        h3(ij) = h3(ij) + 0.5_wp*w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     ENDDO
!  ENDDO
!  DO ij = 1,isize_field
!     h1(ij) = h1(ij) - 0.5_wp*u(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     h2(ij) = h2(ij) - 0.5_wp*v(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!     h3(ij) = h3(ij) - 0.5_wp*w(ij)*( tmp3(ij) + tmp2(ij) + tmp1(ij) )
!  ENDDO
! no big impact of this part, but in the BCs I need it because I already have
! add 1/2 of it.

! #######################################################################
! IBM
! #######################################################################
    if (imode_ibm == 1) ibm_partial = .false. ! until here, IBM is used for flow fields

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
    if (imode_ibm == 1) then
        call IBM_BCS_FIELD(tmp1)
        call IBM_BCS_FIELD(tmp2)
        call IBM_BCS_FIELD(tmp3)
    end if
    if (stagger_on) then ! staggering on horizontal pressure nodes
        !  Ox derivative
        call OPR_PARTIAL_X(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(1), tmp1, tmp5)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)
        !  Oy derivative
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp2, tmp6)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp6, tmp2)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp5)
        !  Oz derivative
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp3, tmp1)
        call OPR_PARTIAL_Z(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(3), tmp1, tmp6)
    else
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp5)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp6)
    end if

! forcing term
    tmp1 = tmp6 + tmp5 + tmp4

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
! Stagger also Bcs
    if (imode_ibm == 1) call IBM_BCS_FIELD(hq(:, 2))
    if (stagger_on) then ! todo: only need to stagger upper/lower boundary plane, not full h2-array
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), hq(:, 2), tmp5)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)
        if (imode_ibm == 1) call IBM_BCS_FIELD_STAGGER(tmp4)
        p_bcs(1:imax, 1:jmax, 1:kmax) => tmp4(1:imax*jmax*kmax)
    else
        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, 2)
    end if

    BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)
    BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)

! pressure in tmp1, Oy derivative in tmp3
    call OPR_Poisson(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)

    if (stagger_on) then
        !  vertical pressure derivative   dpdy - back on horizontal velocity nodes
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp3, tmp5)
        call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp3)
        !  horizontal pressure derivative dpdz - back on horizontal velocity nodes
        call OPR_PARTIAL_Z(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5)
        call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp4)
        !  horizontal pressure derivative dpdx - back on horizontal velocity nodes
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5)
        call OPR_PARTIAL_X(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp2)
    else
        !  horizontal pressure derivatives
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4)
    end if

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
        if (imode_ibm == 1) call IBM_BCS_FIELD(hq(1, iq)) ! set tendency in solid to zero

        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, iq)
        p_bcs(:, 1, :) = BcsFlowJmin%ref(:, :, iq)
        p_bcs(:, jmax, :) = BcsFlowJmax%ref(:, :, iq)

    end do

    return
end subroutine RHS_FLOW_GLOBAL_INCOMPRESSIBLE_2
