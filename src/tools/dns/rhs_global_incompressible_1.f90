#include "dns_const.h"

!########################################################################
!#
!# Evolution equations, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!# It is written such that u and w transposes are calculated first for the
!# Ox and Oz momentum equations, stored in tmp5 and tmp6 and then used as needed.
!# This saves 2 MPI transpositions.
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
subroutine RHS_GLOBAL_INCOMPRESSIBLE_1()
#ifdef USE_OPENMP
    use OMP_LIB
#endif
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
#endif
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: imode_ibm
    use TLAB_VARS, only: imode_eqns, istagger
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g
    use TLAB_VARS, only: rbackground, ribackground
    use TLAB_ARRAYS
    use TLAB_POINTERS, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    use DNS_ARRAYS
    use DNS_LOCAL, only: remove_divergence
    use DNS_LOCAL, only: use_tower
    use TIME, only: rkm_substep, rkm_endstep, dte
    use DNS_TOWER
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use IBM_VARS, only: imode_ibm_scal, ibm_burgers
    use OPR_PARTIAL
    use OPR_BURGERS
    use OPR_ELLIPTIC

    implicit none

! -----------------------------------------------------------------------
    integer(wi) iq, is, ij
    integer ibc, bcs(2, 2)
    real(wp) dummy
    integer, parameter :: i3 = 3

    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning

    real(wp), dimension(:, :, :), pointer :: p_bcs

#ifdef USE_ESSL
    integer ilen
#endif

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
#endif

! #######################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

#ifdef USE_ESSL
    ilen = isize_field
#endif

! #######################################################################
! Preliminaries for Scalar BC
! (flow BCs initialized below as they are used for pressure in between)
! #######################################################################
! Default is zero
    BcsScalJmin%ref(:, :, :) = 0.0_wp
    BcsScalJmax%ref(:, :, :) = 0.0_wp

! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
    if (any(BcsScalJmin%SfcType(1:inb_scal) == DNS_SFC_LINEAR) .or. any(BcsScalJmax%SfcType(1:inb_scal) == DNS_SFC_LINEAR)) then
        do is = 1, inb_scal
            p_bcs(1:imax, 1:jmax, 1:kmax) => hs(1:imax*jmax*kmax, is)
            if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) BcsScalJmin%ref(:, :, is) = p_bcs(:, 1, :)
            if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) BcsScalJmax%ref(:, :, is) = p_bcs(:, jmax, :)
        end do
    end if

! #######################################################################
! Preliminaries for IBM use
! (OPR_BURGERS_X/Y/Z uses modified fields for derivatives)
! #######################################################################
    if (imode_ibm == 1) ibm_burgers = .true.

! #######################################################################
! Ox diffusion and convection terms in Ox momentum eqn
! Initializing tmp5 for the rest of terms
! #######################################################################
    call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), u, u, u, tmp1, tmp5) ! store u transposed in tmp5
    hq(:, 1) = hq(:, 1) + tmp1(:)

! #######################################################################
! Oy diffusion and convection terms in Oy momentum eqn
! Initializing tmp4 for the rest of terms
! #######################################################################
    call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), v, v, v, tmp2, tmp4) ! store v transposed in tmp4
    hq(:, 2) = hq(:, 2) + tmp2(:)

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
    if (g(3)%size > 1) then

        call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), w, u, tmp5, tmp1, tmp6) ! tmp5 contains u transposed
        call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), w, v, tmp4, tmp2, tmp6) ! tmp4 contains v transposed
        call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), w, w, w, tmp3, tmp6) ! store w transposed in tmp6

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        call DNS_OMP_PARTITION(isize_field, srt, end, siz)
        do ij = srt, end
            hq(ij, 3) = hq(ij, 3) + tmp1(ij) + tmp2(ij) + tmp3(ij)
        end do
!$omp end parallel

    end if

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
    call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), v, u, tmp5, tmp1, tmp2) ! tmp5 contains u transposed
    call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), v, w, tmp6, tmp3, tmp2) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call DNS_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 2) = hq(ij, 2) + tmp1(ij) + tmp3(ij)
    end do
!$omp end parallel

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! The term u u''-u u' has been already added in the beginning
! #######################################################################
    call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), u, v, tmp4, tmp2, tmp1) ! tmp4 contains v transposed
    call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), u, w, tmp6, tmp3, tmp1) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call DNS_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 1) = hq(ij, 1) + tmp2(ij) + tmp3(ij)
    end do
!$omp end parallel

! #######################################################################
! IBM
! #######################################################################
    if (imode_ibm == 1) then
        ibm_burgers = .false. ! until here, IBM is used for flow fields
        if (imode_ibm_scal == 1) then ! IBM usage for scalar field
            ! (requirenments: only possible with objects on bottom boundary
            !  with homogeneous temperature in solid regions)
            ibm_burgers = .true.
        end if
    end if

! #######################################################################
! Diffusion and convection terms in scalar eqns
! #######################################################################
    do is = 1, inb_scal

        call OPR_BURGERS_Y(OPR_B_U_IN, is, imax, jmax, kmax, bcs, g(2), s(1, is), v, tmp4, tmp2, tmp1) ! tmp4 contains v transposed
        hs(:, is) = hs(:, is) + tmp2

        call OPR_BURGERS_X(OPR_B_U_IN, is, imax, jmax, kmax, bcs, g(1), s(1, is), u, tmp5, tmp1, tmp2) ! tmp5 contains u transposed
        call OPR_BURGERS_Z(OPR_B_U_IN, is, imax, jmax, kmax, bcs, g(3), s(1, is), w, tmp6, tmp3, tmp2) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        call DNS_OMP_PARTITION(isize_field, srt, end, siz)
        do ij = srt, end
            hs(ij, is) = hs(ij, is) + tmp1(ij) + tmp3(ij)
        end do
!$omp end parallel

    end do

! #######################################################################
! IBM
! #######################################################################
! IBM usage for scalar field, done
    if (imode_ibm_scal == 1) ibm_burgers = .false.

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
    if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
        call BOUNDARY_BUFFER_RELAX_FLOW()
    end if

! #######################################################################
! Pressure term
! #######################################################################
    if (remove_divergence) then ! remove residual divergence

#ifdef USE_ESSL
!$omp parallel default( shared )&
!$omp private( ilen, dummy, srt,end,siz )
#else
!$omp parallel default( shared )&
!$omp private( ij,   dummy, srt,end,siz )
#endif

        call DNS_OMP_PARTITION(isize_field, srt, end, siz)
        dummy = 1.0_wp/dte

#ifdef USE_ESSL
        ilen = siz
        call DZAXPY(ilen, dummy, v(srt), 1, hq(srt, 2), 1, tmp2(srt), 1)
        call DZAXPY(ilen, dummy, u(srt), 1, hq(srt, 1), 1, tmp3(srt), 1)
        call DZAXPY(ilen, dummy, w(srt), 1, hq(srt, 3), 1, tmp4(srt), 1)

#else
        do ij = srt, end
            tmp2(ij) = hq(ij, 2) + v(ij)*dummy
            tmp3(ij) = hq(ij, 1) + u(ij)*dummy
            tmp4(ij) = hq(ij, 3) + w(ij)*dummy
        end do

#endif
!$omp end parallel

        if (imode_ibm == 1) then
            call IBM_BCS_FIELD(tmp2)
            call IBM_BCS_FIELD(tmp3)
            call IBM_BCS_FIELD(tmp4)
        end if
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp2)
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp3)
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp4)
        end if
        if (istagger == 1) then ! staggering on horizontal pressure nodes
            !  Oy derivative
            call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp2, tmp5, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp5, tmp2, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
            !  Ox derivative
            call OPR_PARTIAL_X(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(1), tmp3, tmp5, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp2, wrk3d, wrk2d, wrk3d)
            !  Oz derivative
            call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp4, tmp5, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp3, wrk3d, wrk2d, wrk3d)
        else
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp2, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3, wrk3d, wrk2d, wrk3d)
        end if

    else
        if (imode_ibm == 1) then
            call IBM_BCS_FIELD(hq(:, 2))
            call IBM_BCS_FIELD(hq(:, 1))
            call IBM_BCS_FIELD(hq(:, 3))
        end if
        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 2), tmp2)
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 1), tmp3)
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 3), tmp4)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp2, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3, wrk3d, wrk2d, wrk3d)
        else
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), hq(:, 2), tmp1, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), hq(:, 1), tmp2, wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), hq(:, 3), tmp3, wrk3d, wrk2d, wrk3d)
        end if

    end if

! -----------------------------------------------------------------------
!$omp parallel default( shared ) private( ij,srt,end,siz )
    call DNS_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
    end do
!$omp end parallel

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
! Stagger also Bcs
    if (imode_ibm == 1) call IBM_BCS_FIELD(hq(:, 2))
    if (istagger == 1) then ! todo: only need to stagger upper/lower boundary plane, not full h2-array
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), hq(:, 2), tmp5, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp4, wrk3d, wrk2d, wrk3d)
        if (imode_ibm == 1) call IBM_BCS_FIELD_STAGGER(tmp4)
        p_bcs(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    else
        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, 2)
    end if

    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)*rbackground(1)
        BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)*rbackground(g(2)%size)
    else
        BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)
        BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)
    end if

! pressure in tmp1, Oy derivative in tmp3
    call OPR_POISSON_FXZ(imax, jmax, kmax, g, i3, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)

! Saving pressure for towers to tmp array
    if (use_tower .and. rkm_substep == rkm_endstep) then
        call DNS_TOWER_ACCUMULATE(tmp1, 4, wrk1d)
    end if

    if (istagger == 1) then
        !  vertical pressure derivative   dpdy - back on horizontal velocity nodes
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp3, tmp5, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp3, wrk3d, wrk2d, wrk3d)
        !  horizontal pressure derivative dpdz - back on horizontal velocity nodes
        call OPR_PARTIAL_Z(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp4, wrk3d, wrk2d, wrk3d)
        !  horizontal pressure derivative dpdx - back on horizontal velocity nodes
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp2, wrk3d, wrk2d, wrk3d)
    else
        !  horizontal pressure derivatives
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4, wrk3d, wrk2d, wrk3d)
    end if

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax, jmax, kmax, ribackground, tmp2, hq(:, 1))
        call THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax, jmax, kmax, ribackground, tmp3, hq(:, 2))
        call THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax, jmax, kmax, ribackground, tmp4, hq(:, 3))

    else
#ifdef USE_ESSL
!$omp parallel default( shared ) &
!$omp private( ilen, srt,end,siz,dummy )
#else
!$omp parallel default( shared ) &
!$omp private( ij,   srt,end,siz,dummy )
#endif
        call DNS_OMP_PARTITION(isize_field, srt, end, siz)

#ifdef USE_ESSL
        ilen = siz
        dummy = -1.0_wp
        call DAXPY(ilen, dummy, tmp2(srt), 1, hq(srt, 1), 1)
        call DAXPY(ilen, dummy, tmp3(srt), 1, hq(srt, 2), 1)
        call DAXPY(ilen, dummy, tmp4(srt), 1, hq(srt, 3), 1)
#else
        do ij = srt, end
            hq(ij, 1) = hq(ij, 1) - tmp2(ij)
            hq(ij, 2) = hq(ij, 2) - tmp3(ij)
            hq(ij, 3) = hq(ij, 3) - tmp4(ij)
        end do
#endif
!$omp end parallel
    end if

! #######################################################################
! Boundary conditions
! #######################################################################
    BcsFlowJmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsFlowJmax%ref = 0.0_wp ! Scalar BCs initialized at start of routine

    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hq(1, iq), &
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), wrk1d, tmp1, wrk3d)
        end if
        if (imode_ibm == 1) call IBM_BCS_FIELD(hq(1, iq)) ! set tendency in solid to zero

        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, iq)
        p_bcs(:, 1, :) = BcsFlowJmin%ref(:, :, iq)
        p_bcs(:, jmax, :) = BcsFlowJmax%ref(:, :, iq)

    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsScalJmax%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hs(1, is), &
                                        BcsScalJmin%ref(1, 1, is), BcsScalJmax%ref(1, 1, is), wrk1d, tmp1, wrk3d)
        end if

        if (BcsScalJmin%type(is) /= DNS_SFC_STATIC .or. &
            BcsScalJmax%type(is) /= DNS_SFC_STATIC) then
            call BOUNDARY_SURFACE_J(is, bcs, s, hs, tmp1, tmp2, tmp3)
        end if
        if (imode_ibm == 1) call IBM_BCS_FIELD(hs(1, is)) ! set tendency in solid to zero

        p_bcs(1:imax, 1:jmax, 1:kmax) => hs(1:imax*jmax*kmax, is)
        p_bcs(:, 1, :) = BcsScalJmin%ref(:, :, is)
        p_bcs(:, jmax, :) = BcsScalJmax%ref(:, :, is)

    end do

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
#endif

    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_1
