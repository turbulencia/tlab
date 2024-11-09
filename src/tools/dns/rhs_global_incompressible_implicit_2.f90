#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2012/07/10 - C. Ansorge
!#              Created
!#
!########################################################################
!#
!# Derived from RHS_GLOBAL_INCOMPRESSIBLE_1
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
subroutine RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2(kex, kim, kco)
#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use TLab_Constants
    use FDM, only: g
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: inb_flow, inb_scal
    use TLAB_VARS, only: scal_on
    use TLAB_VARS, only: visc, schmidt
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Arrays
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
    use DNS_ARRAYS
    use TIME, only: dte
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use FI_SOURCES
    implicit none

    real(wp), intent(in) :: kex, kim, kco

    ! -----------------------------------------------------------------------
    integer(wi) iq, is, iq_max
    integer(wi) ibc, bcs(2, 2)
    real(wp) dummy, visc_exp, visc_imp, visc_tot, alpha, beta, kef, aug, vdummy!, u_geo, w_geo

    real(wp), dimension(imax, kmax, 4) :: bcs_locb, bcs_loct

    integer, parameter :: i0 = 0

    real(wp), dimension(:, :, :), pointer :: p_bcs

! #######################################################################
    do is = 1, inb_scal
        if (BcsScalJmin%SfcType(is) == DNS_SFC_STATIC .and. &
            BcsScalJmax%SfcType(is) == DNS_SFC_STATIC) then
            continue
        else
            call TLab_Write_ASCII(efile, 'Only static surface implemented in implicit mode')
        end if
    end do

    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    visc_exp = kex*visc
    visc_imp = kim*visc
    visc_tot = visc

    kef = kex/kim
    aug = 1.0_wp + kef
    alpha = dte*visc_imp

!  vdummy= -alpha*visc_tot*dte
!  vdummy= -alpha*visc_exp*dte
    vdummy = 0.0_wp

! ################################################################################
! SAVE VALUES AT BOUNDARIES. Recovered at the end of the routine
! ################################################################################
    iq = 1
    p_bcs(1:imax, 1:jmax, 1:kmax) => q(1:imax*jmax*kmax, iq)
    BcsFlowJmin%ref(:, :, iq) = p_bcs(:, 1, :)
    BcsFlowJmax%ref(:, :, iq) = p_bcs(:, jmax, :)

    iq = 3
    p_bcs(1:imax, 1:jmax, 1:kmax) => q(1:imax*jmax*kmax, iq)
    BcsFlowJmin%ref(:, :, iq) = p_bcs(:, 1, :)
    BcsFlowJmax%ref(:, :, iq) = p_bcs(:, jmax, :)

    do is = 1, inb_scal
        p_bcs(1:imax, 1:jmax, 1:kmax) => s(1:imax*jmax*kmax, is)
        BcsScalJmin%ref(:, :, is) = p_bcs(:, 1, :)
        BcsScalJmax%ref(:, :, is) = p_bcs(:, jmax, :)
    end do

! #######################################################################
! Diffusion and convection terms
! #######################################################################
    tmp6 = hq(:, 3) ! SAVE old tendencies until end of implicit substep
    tmp5 = hq(:, 2) ! SAVE old tendencies until end of implicit substep
    tmp4 = hq(:, 1) ! SAVE old tendencies until end of implicit substep

    if (g(3)%size > 1) then
        iq_max = 3
    else
        iq_max = 2
    end if

    do iq = 1, iq_max
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(:, iq), hq(:, iq))
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(:, iq), tmp3)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), q(:, iq), tmp1, tmp2) ! tmp1 used for BCs below
        hq(:, iq) = -u*hq(:, iq) - v*tmp2 - w*tmp3      ! hq(:,iq) contains explicit nonlinear tendency

        ! Diffusion part of the BCs for intermediate velocity
        p_bcs(1:imax, 1:jmax, 1:kmax) => tmp1(1:imax*jmax*kmax)
        bcs_locb(:, :, iq) = vdummy*p_bcs(:, 1, :)      ! bottom
        bcs_loct(:, :, iq) = vdummy*p_bcs(:, jmax, :)   ! top

    end do

! #######################################################################
    call FI_SOURCES_FLOW(q, s, hq, tmp1)

! #######################################################################
! Impose buffer zone as relaxation terms (Flow)
! #######################################################################
    if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
        call BOUNDARY_BUFFER_RELAX_FLOW()
    end if

! #######################################################################
! Explicit part of time integration for velocities
! #######################################################################
! tmp4-6 contain tendencies, h1-3 contain old tendencies
    tmp1 = u*aug + dte*(hq(:, 1) + kco*tmp4)
    tmp2 = v*aug + dte*(hq(:, 2) + kco*tmp5)
    tmp3 = w*aug + dte*(hq(:, 3) + kco*tmp6)

! tmp1-tmp3:  u,v,w updated with explicit scheme
! h1-h3:      explicit tendencies from this substep
! u-w:        old velocities

! -----------------------------------------------------------------------
! Remaining part of the BCs for intermediate velocity
! -----------------------------------------------------------------------
    ! This block implements the BCs retaining a dt correction according to the
    ! notes but it leads to first-order convergence rate instead of second-order
    ! and even instability if the diffusion terms are retained.
    !
    ! ip_b =                 1
    ! ip_t = imax*(jmax-1) + 1
    ! DO k = 1,kmax
    !    bcs_locb(:,k,1) = bcs_locb(:,k,1) - alpha*tmp1(ip_b:) ! bottom
    !    bcs_loct(:,k,1) = bcs_loct(:,k,1) - alpha*tmp1(ip_t:) ! top

    !    bcs_locb(:,k,2) = bcs_locb(:,k,2) - alpha*tmp2(ip_b:) ! bottom
    !    bcs_loct(:,k,2) = bcs_loct(:,k,2) - alpha*tmp2(ip_t:) ! top

    !    bcs_locb(:,k,3) = bcs_locb(:,k,3) - alpha*tmp3(ip_b:) ! bottom
    !    bcs_loct(:,k,3) = bcs_loct(:,k,3) - alpha*tmp3(ip_t:) ! top

    !    ip_b = ip_b + nxy
    !    ip_t = ip_t + nxy
    ! ENDDO

    ! This block retains only the leading order terms and provides second-order
    ! convergence rate.
    ! bcs_locb(:,:,1) = -alpha*aug*bcs_hb(:,:,1) ! bottom
    ! bcs_loct(:,:,1) = -alpha*aug*bcs_ht(:,:,1) ! top
    bcs_locb(:, :, 1) = -alpha*aug*BcsFlowJmin%ref(:, :, 1) ! bottom
    bcs_loct(:, :, 1) = -alpha*aug*BcsFlowJmax%ref(:, :, 1) ! top
    bcs_locb(:, :, 2) = 0.0_wp                    ! bottom
    bcs_loct(:, :, 2) = 0.0_wp                    ! top
    bcs_locb(:, :, 3) = -alpha*aug*BcsFlowJmin%ref(:, :, 3) ! bottom
    bcs_loct(:, :, 3) = -alpha*aug*BcsFlowJmax%ref(:, :, 3) ! top

! ################################################################################
! ADVECTION-DIFFUSION FOR SCALAR
! We need this code segment here because we need the old velocities
! ################################################################################
! hs  -> explicit tendency without diffusion
    if (scal_on) then

        do is = 1, inb_scal
            tmp7 = hs(:, is) ! save old values

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(1, is), hs(1, is))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(1, is), tmp6)
            call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), s(1, is), tmp4, tmp5)

            hs(:, is) = -u*hs(:, is) - v*tmp5 - w*tmp6

! -----------------------------------------------------------------------
! Diffusion part of the BCs for intermediate scalar
! -----------------------------------------------------------------------
            p_bcs(1:imax, 1:jmax, 1:kmax) => tmp4(1:imax*jmax*kmax)
            bcs_locb(:, :, 4) = (vdummy/schmidt(is))*p_bcs(:, 1, :) ! bottom
            bcs_loct(:, :, 4) = (vdummy/schmidt(is))*p_bcs(:, jmax, :) ! top

! #######################################################################
! Impose buffer zone as relaxation terms (Scalar #is)
! #######################################################################
            if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
                call BOUNDARY_BUFFER_RELAX_SCAL_I(is, s(1, is), hs(1, is))
            end if

! #######################################################################
! Explicit Time Stepping for scalar #is
! #######################################################################
            tmp4 = s(:, is)*aug + dte*(hs(:, is) + kco*tmp7)

! -----------------------------------------------------------------------
! Remaining part of the BCs for intermediate scalar
! -----------------------------------------------------------------------
            ! See velocity part
            ! ip_b =                 1
            ! ip_t = imax*(jmax-1) + 1
            ! DO k = 1,kmax
            !    bcs_locb(:,k,4) = bcs_locb(:,k,4) - (alpha/schmidt(is))*tmp4(ip_b:); ip_b = ip_b + nxy ! bottom
            !    bcs_loct(:,k,4) = bcs_loct(:,k,4) - (alpha/schmidt(is))*tmp4(ip_t:); ip_t = ip_t + nxy ! top
            ! ENDDO

            ! bcs_locb(:,:,4) = -(alpha/schmidt(is))*aug*bcs_hb(:,:,inb_flow+is)
            ! bcs_loct(:,:,4) = -(alpha/schmidt(is))*aug*bcs_ht(:,:,inb_flow+is)
            bcs_locb(:, :, 4) = -(alpha/schmidt(is))*aug*BcsScalJmin%ref(:, :, is)
            bcs_loct(:, :, 4) = -(alpha/schmidt(is))*aug*BcsScalJmax%ref(:, :, is)

! #######################################################################
! semi-Implicit Diffusion for Scalar #is
! #######################################################################
            beta = -1.0_wp/(alpha/schmidt(is))

            ! call OPR_Helmholtz_FourierXZ_Direct(imax, jmax, kmax, g, 0, beta, &
            !                          tmp4, tmp6, tmp7, bcs_locb(1, 1, 4), bcs_loct(1, 1, 4))
            call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, tmp4, tmp6, tmp7, bcs_locb(1, 1, 4), bcs_loct(1, 1, 4))

            s(:, is) = beta*tmp4 - kef*s(:, is)

        end do
    end if

! ################################################################################
! Implicit part of time integration for velocities
! ################################################################################
    beta = -1.0_wp/alpha

    ! call OPR_Helmholtz_FourierXZ_Direct(imax, jmax, kmax, g, 0, beta, &
    !                                     tmp1, tmp5, tmp6, bcs_locb(1, 1, 1), bcs_loct(1, 1, 1))
    ! call OPR_Helmholtz_FourierXZ_Direct(imax, jmax, kmax, g, 0, beta, &
    !                                     tmp2, tmp5, tmp6, bcs_locb(1, 1, 2), bcs_loct(1, 1, 2))
    ! call OPR_Helmholtz_FourierXZ_Direct(imax, jmax, kmax, g, 0, beta, &
    !                                     tmp3, tmp5, tmp6, bcs_locb(1, 1, 3), bcs_loct(1, 1, 3))
    call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, tmp1, tmp5, tmp6, bcs_locb(1, 1, 1), bcs_loct(1, 1, 1))
    call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, tmp2, tmp5, tmp6, bcs_locb(1, 1, 2), bcs_loct(1, 1, 2))
    call OPR_Helmholtz(imax, jmax, kmax, g, 0, beta, tmp3, tmp5, tmp6, bcs_locb(1, 1, 3), bcs_loct(1, 1, 3))

    u = beta*tmp1 - kef*u
    v = beta*tmp2 - kef*v
    w = beta*tmp3 - kef*w

! #######################################################################
! Pressure term (actually solving for dte*p)
! #######################################################################
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), u, tmp1)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), v, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), w, tmp3)

! -----------------------------------------------------------------------
    tmp1 = tmp1 + tmp2 + tmp3 ! forcing term in tmp1

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    iq = 2
    p_bcs(1:imax, 1:jmax, 1:kmax) => q(1:imax*jmax*kmax, iq)
    BcsFlowJmin%ref(:, :, iq) = p_bcs(:, 1, :)
    BcsFlowJmax%ref(:, :, iq) = p_bcs(:, jmax, :)

! pressure in tmp1, Oy derivative in tmp3
    ! select case (imode_elliptic)
    ! case (FDM_COM6_JACOBIAN)
    !     call OPR_Poisson_FourierXZ_Factorize(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    ! case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
    !     call OPR_Poisson_FourierXZ_Direct(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    ! end select
    call OPR_Poisson(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)

! horizontal derivatives
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
    dummy = 1.0_wp/dte
    u = u - tmp2; hq(:, 1) = hq(:, 1) - dummy*tmp2
    v = v - tmp3; hq(:, 2) = hq(:, 2) - dummy*tmp3
    w = w - tmp4; hq(:, 3) = hq(:, 3) - dummy*tmp4

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

        p_bcs(1:imax, 1:jmax, 1:kmax) => q(1:imax*jmax*kmax, iq)
        p_bcs(:, 1, :) = BcsFlowJmin%ref(:, :, iq)
        p_bcs(:, jmax, :) = BcsFlowJmax%ref(:, :, iq)

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

        p_bcs(1:imax, 1:jmax, 1:kmax) => s(1:imax*jmax*kmax, is)
        p_bcs(:, 1, :) = BcsScalJmin%ref(:, :, is)
        p_bcs(:, jmax, :) = BcsScalJmax%ref(:, :, is)

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
!      IF (  .NE. 0scal_on ) THEN
!         DO is=1,inb_scal
!            s(ip_b:ip_b+imax-1,is) = bcs_hb(1:imax,k,inb_flow+is)
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
!            s(ip_t:ip_t+imax-1,is) = bcs_ht(1:imax,k,inb_flow+is)
!         ENDDO
!      ENDIF
!      ip_t = ip_t + nxy
!   ENDDO

    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2
