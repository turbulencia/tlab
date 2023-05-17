#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2012/07/10 - C. Ansorge
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Derived from RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1
!#
!# Momentum equations, nonlinear term in convective form and the
!# viscous term semi-implicit: 9 1st order derivatives.
!# The explicit laplacian is eliminated by augmenting the array on which
!# the Helmholtz equation is solved.
!#
!# Pressure term requires 3 1st order derivatives
!# BCs need 2 1st order derivatives in Oy
!# Scalar needed for the buoyancy term
!#
!########################################################################
subroutine RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_3(kex, kim, kco, &
                                                q, hq, u, v, w, h1, h2, h3, s, hs, &
                                                tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
    use TLAB_CONSTANTS
#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: isize_field, isize_txc_field, inb_scal, inb_flow
    use TLAB_VARS, only: scal_on
    use TLAB_VARS, only: visc, schmidt, rossby, ipressure
    use TLAB_VARS, only: buoyancy, coriolis
    use TLAB_VARS, only: bbackground
    use TLAB_ARRAYS, only: wrk1d, wrk2d, wrk3d
    use TIME, only: dte
    use DNS_LOCAL, only: remove_divergence
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use FI_SOURCES, only: FI_BUOYANCY

    implicit none

    real(wp) kex, kim, kco
    real(wp), dimension(isize_field, *) :: q, hq
    real(wp), dimension(isize_field), intent(INOUT) :: u, v, w, h1, h2, h3
    real(wp), dimension(isize_field, inb_scal), intent(INOUT) :: s, hs
    real(wp), dimension(isize_txc_field), intent(OUT) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

    target tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, h2, u, v, w

! -----------------------------------------------------------------------
    integer(wi) iq, is, ij, k, nxy, ip, ip_b, ip_t
    integer(wi) ibc, bcs(2, 2)
    real(wp) dummy, visc_exp, visc_imp, visc_tot, diff, alpha, beta, kef, aug

    real(wp), dimension(:), pointer :: p_bcs
    real(wp), dimension(imax, kmax, 4:6) :: bcs_hb, bcs_ht

    integer, parameter :: i0 = 0

    kef = kex/kim
    aug = 1.0_wp + kef

! #######################################################################
    nxy = imax*jmax

    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    visc_exp = kex*visc
    visc_imp = kim*visc
    visc_tot = visc

! SAVE old tendencies until end of implicit substep of integration

    tmp4(1:isize_field) = hq(:, 1)
    tmp5(1:isize_field) = hq(:, 2)
    tmp6(1:isize_field) = hq(:, 3)

! ################################################################################
! SAVE VALUES AT BOUNDARIES
! ################################################################################

    ip = 1
    do k = 1, kmax
        ip_b = (k - 1)*imax*jmax + 1; 
        ip_t = ip_b + (jmax - 1)*imax
        ! bcs_hb(:,k,1) = u(ip_b:ip_b+imax-1); bcs_hb(:,k,2) = w(ip_b:ip_b+imax-1)
        ! bcs_ht(:,k,1) = u(ip_t:ip_t+imax-1); bcs_ht(:,k,2) = w(ip_t:ip_t+imax-1)
        BcsFlowJmin%ref(1:imax, k, 1) = u(ip_b:ip_b + imax - 1); BcsFlowJmin%ref(1:imax, k, 3) = w(ip_b:ip_b + imax - 1)
        BcsFlowJmax%ref(1:imax, k, 1) = u(ip_t:ip_t + imax - 1); BcsFlowJmax%ref(1:imax, k, 3) = w(ip_t:ip_t + imax - 1)
        do is = 1, inb_scal
!        bcs_sb(1:imax,k,is) = s(ip_b:ip_b+imax-1,is); bcs_st(1:imax,k,is) = s(ip_t:ip_t+imax-1,is)
            BcsScalJmin%ref(1:imax, k, is) = s(ip_b:ip_b + imax - 1, is); BcsScalJmax%ref(1:imax, k, is) = s(ip_t:ip_t + imax - 1, is)
        end do
    end do

! ################################################################################
! EXPLICIT PART OF SUBSTEP - change viscosity accordingly
! ################################################################################
    visc = visc_exp

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
    ! h3   contains explicit nonlinear w-tendency (nonlinear operator N(u_n))
    if (g(3)%size > 1) then
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), w, h3)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), w, tmp3, tmp2) ! tmp2 used for BCs below

        do k = 1, kmax
            ip = (k - 1)*imax*jmax + 1; 
!        bcs_hb(1:imax,k,6)=bcs_hb(1:imax,k,2) + visc_tot*dte*tmp2(ip:ip+imax-1) - dte*1.0_wp;
            bcs_hb(1:imax, k, 6) = BcsFlowJmin%ref(1:imax, k, 3) + visc_tot*dte*tmp2(ip:ip + imax - 1) - dte*1.0_wp; 
            ip = ip + (jmax - 1)*imax; 
!        bcs_ht(1:imax,k,6)=bcs_ht(1:imax,k,2) + visc_tot*dte*tmp2(ip:ip+imax-1) - dte*1.0_wp;
            bcs_ht(1:imax, k, 6) = BcsFlowJmax%ref(1:imax, k, 3) + visc_tot*dte*tmp2(ip:ip + imax - 1) - dte*1.0_wp; 
        end do

        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)
        do ij = 1, isize_field
            h3(ij) = -u(ij)*h3(ij) - v(ij)*tmp2(ij) - w(ij)*tmp3(ij)
        end do

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
        if (buoyancy%active(3)) then
            wrk1d(:, 1) = 0.0_wp
            call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, wrk1d)
            dummy = buoyancy%vector(3)
            do ij = 1, isize_field
                h3(ij) = h3(ij) + dummy*wrk3d(ij)
            end do
        end if

! -----------------------------------------------------------------------
! Coriolis (so far, rotation only in the Oy direction)
! -----------------------------------------------------------------------
        if (coriolis%type == EQNS_COR_NORMALIZED) then
            dummy = 1.0_wp/rossby
            do ij = 1, isize_field
                h3(ij) = h3(ij) + dummy*(u(ij) - 1.0_wp)
            end do
        end if

    end if

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! #######################################################################
! h1 contains explicit nonlinear u-tendency
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, h1)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), u, tmp3, tmp2) ! tmp2 used for BCs below

    do k = 1, kmax
        ip = (k - 1)*imax*jmax + 1; 
!     bcs_hb(1:imax,k,4)=bcs_hb(1:imax,k,1) + visc_tot*dte*tmp2(ip:ip+imax-1)
        bcs_hb(1:imax, k, 4) = BcsFlowJmin%ref(1:imax, k, 1) + visc_tot*dte*tmp2(ip:ip + imax - 1)
        ip = ip + (jmax - 1)*imax; 
!     bcs_ht(1:imax,k,4)=bcs_ht(1:imax,k,1) + visc_tot*dte*tmp2(ip:ip+imax-1)
        bcs_ht(1:imax, k, 4) = BcsFlowJmax%ref(1:imax, k, 1) + visc_tot*dte*tmp2(ip:ip + imax - 1)
    end do

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), u, tmp3)
    do ij = 1, isize_field
        h1(ij) = -u(ij)*h1(ij) - v(ij)*tmp2(ij) - w(ij)*tmp3(ij)
    end do
! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
    if (buoyancy%active(1)) then
        wrk1d(:, 1) = 0.0_wp
        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, wrk1d)
        dummy = buoyancy%vector(1)
        do ij = 1, isize_field
            h1(ij) = h1(ij) + dummy*wrk3d(ij)
        end do
    end if

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction.
! -----------------------------------------------------------------------

    if (coriolis%type == EQNS_COR_NORMALIZED) then
        dummy = 1.0_wp/rossby
        do ij = 1, isize_field
            h1(ij) = h1(ij) - dummy*w(ij)
        end do
    end if

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
!   h2 contains explicit nonlinear v-tendency
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), v, h2)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), v, tmp3, tmp2) ! tmp2 used for BCs below

    do k = 1, kmax
        ip = (k - 1)*imax*jmax + 1; 
        bcs_hb(1:imax, k, 5) = 0.0_wp + visc_tot*dte*tmp2(ip:ip + imax - 1)
        ip = ip + (jmax - 1)*imax; 
        bcs_ht(1:imax, k, 5) = 0.0_wp + visc_tot*dte*tmp2(ip:ip + imax - 1)
    end do

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), v, tmp3)
    do ij = 1, isize_field
        h2(ij) = -u(ij)*h2(ij) - v(ij)*tmp2(ij) - w(ij)*tmp3(ij)
    end do
! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
    if (buoyancy%active(2)) then
        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, bbackground)
        dummy = buoyancy%vector(2)
        do ij = 1, isize_field
            h2(ij) = h2(ij) - w(ij)*tmp3(ij) + dummy*wrk3d(ij)
        end do
    else
        do ij = 1, isize_field
            h2(ij) = h2(ij) - w(ij)*tmp3(ij)
        end do
    end if

! #######################################################################
! Impose buffer zone as relaxation terms (Flow)
! #######################################################################
    if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
        call BOUNDARY_BUFFER_RELAX_FLOW()
    end if

! old pressure tendencies
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp8, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp8, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp8, tmp3)
    do ij = 1, isize_field
        h1(ij) = h1(ij) - tmp1(ij)
        h2(ij) = h2(ij) - tmp2(ij)
        h3(ij) = h3(ij) - tmp3(ij)
    end do

! Explicit part of time integration
! tmp4-6 contain tendencies, h1-3 contain old tendencies
    do ij = 1, isize_field
        tmp1(ij) = u(ij)*aug + dte*(h1(ij) + kco*tmp4(ij)); 
        tmp2(ij) = v(ij)*aug + dte*(h2(ij) + kco*tmp5(ij)); 
        tmp3(ij) = w(ij)*aug + dte*(h3(ij) + kco*tmp6(ij)); 
    end do
! tmp1-tmp3:  u,v,w updated with explicit scheme
! h1-h3:      explicit tendencies from this substep
! u-w:        old velocities

! ################################################################################
! ADVECTION-DIFFUSION FOR SCALAR
! ################################################################################
! hs  -> explicit tendency without diffusion
    if (scal_on) then
        do is = 1, inb_scal
            diff = visc_exp/schmidt(is)
            tmp4 = hs(:, is)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(1, is), hs(1, is))
            call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), s(1, is), tmp6, tmp5)
!        CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s,       tmp6,     tmp5)
            do k = 1, kmax
                ! ip=(k-1)*imax*jmax+1; bcs_hb(1:imax,k,3)= 0.0_wp + visc_tot*dte*tmp6(ip:ip+imax-1)
                ! ip=ip+(jmax-1)*imax;  bcs_ht(1:imax,k,3)= 0.0_wp + visc_tot*dte*tmp6(ip:ip+imax-1)
                ip = (k - 1)*imax*jmax + 1; BcsFlowJmin%ref(1:imax, k, 2) = 0.0_wp + visc_tot*dte*tmp6(ip:ip + imax - 1)
                ip = ip + (jmax - 1)*imax; BcsFlowJmax%ref(1:imax, k, 2) = 0.0_wp + visc_tot*dte*tmp6(ip:ip + imax - 1)
            end do
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(1, is), tmp6)
            do ij = 1, isize_field; 
                hs(ij, is) = -u(ij)*hs(ij, is) - v(ij)*tmp5(ij) - w(ij)*tmp6(ij)
            end do

! #######################################################################
! Impose buffer zone as relaxation terms (Scalar #is)
! #######################################################################
            if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
                call BOUNDARY_BUFFER_RELAX_SCAL_I(is, s(1, is), hs(1, is))
            end if

! --------------------------------------------------
! Explicit Time Stepping for scalar #is
! --------------------------------------------------
            do ij = 1, isize_field
                tmp7(ij) = s(ij, is)*aug + dte*(hs(ij, is) + kco*tmp6(ij)); 
            end do
! --------------------------------------------------
! semi-Implicit Diffusion for Scalar #is
! --------------------------------------------------
            diff = visc_imp/schmidt(is)
            alpha = dte*diff
            beta = -1.0_wp/alpha

            ip_b = 1; ip_t = 1 + imax*kmax
            do k = 1, kmax
!!!!!! I THINK THIS IS WRONG AND SHOULD USE bcs_hb(:,:,3) and bcs_ht(:,:,3) instead of bcs_sb and bcs_st
                ! ip=(k-1)*imax*jmax+1; wrk2d(ip_b:ip_b+imax-1)=-alpha*aug*bcs_sb(1:imax,k,is); ip_b=ip_b+imax
                ! ip=ip+(jmax-1)*imax;  wrk2d(ip_t:ip_t+imax-1)=-alpha*aug*bcs_st(1:imax,k,is); ip_t=ip_t+imax
!!!!!! I THINK THIS IS WRONG AND SHOULD USE BcsFlowJmin%ref(1:imax,k,2)
                ! ip=(k-1)*imax*jmax+1; wrk2d(ip_b:ip_b+imax-1)=-alpha*aug*BcsScalJmin%ref(1:imax,k,is); ip_b=ip_b+imax
                ! ip=ip+(jmax-1)*imax;  wrk2d(ip_t:ip_t+imax-1)=-alpha*aug*BcsScalJmax%ref(1:imax,k,is); ip_t=ip_t+imax
                ip = (k - 1)*imax*jmax + 1; wrk2d(ip_b:ip_b + imax - 1, 1) = -alpha*aug*BcsFlowJmin%ref(1:imax, k, 2); ip_b = ip_b + imax
                ip = ip + (jmax - 1)*imax; wrk2d(ip_t:ip_t + imax - 1, 1) = -alpha*aug*BcsFlowJmax%ref(1:imax, k, 2); ip_t = ip_t + imax
            end do
            ip_b = 1; ip_t = 1 + imax*kmax
            call OPR_HELMHOLTZ_FXZ(imax, jmax, kmax, g, 0, beta, &
                                   tmp7, tmp5, tmp6, wrk2d(ip_b, 1), wrk2d(ip_t, 1))

            do ij = 1, isize_field
                s(ij, is) = tmp7(ij)*beta - kef*s(ij, is)
            end do
        end do
    end if

! ################################################################################
! IMPLICIT PART OF SUBSTEP - change viscosity accordingly
! ################################################################################
! TO BE CLEANED: in this section wrk2d is used for the BCS of the implicit solver

    visc = visc_imp
    alpha = dte*visc_imp
    beta = -1.0_wp/alpha

    bcs_hb(1:imax, 1:kmax, 4:6) = -alpha*aug*bcs_hb(1:imax, 1:kmax, 4:6)
    bcs_ht(1:imax, 1:kmax, 4:6) = -alpha*aug*bcs_ht(1:imax, 1:kmax, 4:6)
    call OPR_HELMHOLTZ_FXZ(imax, jmax, kmax, g, 0, beta, &
                           tmp1, tmp5, tmp6, bcs_hb(1, 1, 4), bcs_ht(1, 1, 4))
    call OPR_HELMHOLTZ_FXZ(imax, jmax, kmax, g, 0, beta, &
                           tmp2, tmp5, tmp6, bcs_hb(1, 1, 5), bcs_ht(1, 1, 5))
    call OPR_HELMHOLTZ_FXZ(imax, jmax, kmax, g, 0, beta, &
                           tmp3, tmp5, tmp6, bcs_hb(1, 1, 6), bcs_ht(1, 1, 6))

    do ij = 1, isize_field
        u(ij) = tmp1(ij)*beta - kef*u(ij)
        v(ij) = tmp2(ij)*beta - kef*v(ij)
        w(ij) = tmp3(ij)*beta - kef*w(ij)
    end do

! ################################################################################
! END OF IMPLICIT PART OF SUBSTEP - set viscosity back to actual value
! ################################################################################
    visc = visc_tot

! #######################################################################
! Pressure term
! #######################################################################
    if (remove_divergence) then ! remove residual divergence
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)
    end if

! -----------------------------------------------------------------------
    do ij = 1, isize_field
        tmp1(ij) = (tmp1(ij) + tmp2(ij) + tmp3(ij))/dte ! forcing term in tmp1
    end do

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    ip_b = 1
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        p_bcs => v(ip_b:); BcsFlowJmin%ref(1:imax, k, 2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
        p_bcs => v(ip_t:); BcsFlowJmax%ref(1:imax, k, 2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
    end do

! pressure in tmp8, Oy derivative in tmp3
    select case (ipressure)
    case (FDM_COM6_JACOBIAN)
        call OPR_POISSON_FXZ(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
        call OPR_POISSON_FXZ_D(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    end select

! update pressure with correction from pressure solver
    do ij = 1, isize_field
        tmp8(ij) = tmp8(ij) + tmp1(ij)
    end do

! horizontal derivatives
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
    do ij = 1, isize_field
        u(ij) = u(ij) - dte*tmp2(ij)
        v(ij) = v(ij) - dte*tmp3(ij)
        w(ij) = w(ij) - dte*tmp4(ij)
    end do

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
! Default is Dirichlet -> values at boundary are kept const;
! BcsFlowJmin/Jmax for u and w initialized to old BC values above
    BcsFlowJmin%ref(:, :, 2) = 0.0_wp
    BcsFlowJmax%ref(:, :, 2) = 0.0_wp
    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), q(1, iq), &
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), tmp1)
        end if
    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsScalJmax%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), s(1, is), &
                                        BcsScalJmin%ref(1, 1, is), BcsScalJmax%ref(1, 1, is), tmp1)
        end if
    end do

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin
! -----------------------------------------------------------------------
    ip_b = 1
    do k = 1, kmax
        u(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 1)
        v(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 2)
        w(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 3)
        do is = 1, inb_scal
            s(ip_b:ip_b + imax - 1, is) = BcsScalJmin%ref(1:imax, k, is)
        end do
        ip_b = ip_b + nxy
    end do

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        u(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 1)
        v(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 2)
        w(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 3)
        do is = 1, inb_scal
            s(ip_t:ip_t + imax - 1, is) = BcsScalJmax%ref(1:imax, k, is)
        end do
        ip_t = ip_t + nxy
    end do
! ! -----------------------------------------------------------------------
! ! Preliminaries
! ! -----------------------------------------------------------------------

!   ibc = 0
! ! Default is Dirichlet -> values at boundary are kept const;
! ! bcs_hb and bcs_ht initialized to old BC values above
!   IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
!   IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
!   IF ( ibc .GT. 0 ) THEN
!      ! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
!      CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), u, &
!           bcs_hb(1,1,1),bcs_ht(1,1,1), tmp1)
!      CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), w, &
!           bcs_hb(1,1,2),bcs_ht(1,1,2), tmp1)
!   ENDIF

! ! -----------------------------------------------------------------------
! ! Impose bottom BCs at Jmin
! ! -----------------------------------------------------------------------
!   ip_b =                 1
!   DO k = 1,kmax
!      u(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
!      v(ip_b:ip_b+imax-1) = 0.0_wp               ! no penetration
!      w(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2);
!      IF ( scal_on ) THEN
!         DO is=1,inb_scal
!            s(ip_b:ip_b+imax-1,is) = bcs_sb(1:imax,k,is)
!         ENDDO
!      ENDIF
!      ip_b = ip_b + nxy
!   ENDDO

! ! -----------------------------------------------------------------------
! ! Impose top BCs at Jmax
! ! -----------------------------------------------------------------------
!   ip_t = imax*(jmax-1) + 1
!   DO k = 1,kmax
!      u(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
!      v(ip_t:ip_t+imax-1) = 0.0_wp               ! no penetration
!      w(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2);
!      IF ( scal_on ) THEN
!         DO is=1,inb_scal
!            s(ip_t:ip_t+imax-1,is) = bcs_st(1:imax,k,is)
!         ENDDO
!      ENDIF
!      ip_t = ip_t + nxy
!   ENDDO

    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_3
