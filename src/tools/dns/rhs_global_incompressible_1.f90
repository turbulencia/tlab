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
    use TLab_Constants, only: tfile
#endif
    use TLab_Constants, only: wp, wi, BCS_NN
    use NavierStokes, only: nse_eqns
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use FDM, only: g
    use TLAB_VARS, only: stagger_on
    use TLAB_VARS, only: itime
    use TLab_Arrays
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
    use THERMO_ANELASTIC
    use TLab_OpenMP
    use DNS_ARRAYS
    use DNS_LOCAL, only: remove_divergence
    use DNS_LOCAL, only: use_tower
    use DNS_LOCAL, only: nitera_first, nitera_save
    use TIME, only: rkm_substep, rkm_endstep, dte
    use DNS_TOWER
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use IBM_VARS, only: imode_ibm, imode_ibm_scal, ibm_burgers
    use OPR_PARTIAL
    use OPR_Burgers
    use OPR_ELLIPTIC
    use OPR_FILTERS
    use AVG_PHASE

    implicit none

    ! -----------------------------------------------------------------------
    integer(wi) iq, is, ij
    integer ibc, bcs(2, 2)
    real(wp) dummy

    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning

    real(wp), dimension(:, :, :), pointer :: p_bcs

#ifdef USE_ESSL
    integer ilen
#endif

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'ENTERING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
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
    ! Diffusion and advection terms
    ! #######################################################################
    ! Preliminaries for IBM use
    ! (OPR_Burgers_X/Y/Z uses modified fields for derivatives)
    if (imode_ibm == 1) ibm_burgers = .true.

    ! Diagonal terms and transposed velocity arrays
    call OPR_Burgers_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, u, u, tmp1, tmp4) ! store u transposed in tmp4
    call OPR_Burgers_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, v, v, tmp2, tmp5) ! store v transposed in tmp5
    call OPR_Burgers_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, w, w, tmp3, tmp6) ! store w transposed in tmp6

    ! Ox momentum equation
    call OPR_Burgers_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, u, v, tmp7, tmp9, tmp5) ! tmp5 contains v transposed
    call OPR_Burgers_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, u, w, tmp8, tmp9, tmp6) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 1) = hq(ij, 1) + tmp1(ij) + tmp7(ij) + tmp8(ij)
    end do
!$omp end parallel

    ! Oy momentum equation
    call OPR_Burgers_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, v, u, tmp7, tmp9, tmp4) ! tmp4 contains u transposed
    call OPR_Burgers_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, v, w, tmp8, tmp9, tmp6) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 2) = hq(ij, 2) + tmp2(ij) + tmp7(ij) + tmp8(ij)
    end do
!$omp end parallel

    ! Oz momentum equation
    call OPR_Burgers_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, w, u, tmp7, tmp9, tmp4) ! tmp4 contains u transposed
    call OPR_Burgers_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, w, v, tmp8, tmp9, tmp5) ! tmp5 contains v transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        hq(ij, 3) = hq(ij, 3) + tmp3(ij) + tmp7(ij) + tmp8(ij)
    end do
!$omp end parallel

    ! IBM
    if (imode_ibm == 1) then
        ibm_burgers = .false. ! until here, IBM is used for flow fields
        if (imode_ibm_scal == 1) then ! IBM usage for scalar field
            ! (requirenments: only possible with objects on bottom boundary
            !  with homogeneous temperature in solid regions)
            ibm_burgers = .true.
        end if
    end if

    ! Scalar equations
    do is = 1, inb_scal
        call OPR_Burgers_X(OPR_B_U_IN, is, imax, jmax, kmax, bcs, s(1, is), u, tmp1, tmp9, tmp4) ! tmp4 contains u transposed
        call OPR_Burgers_Y(OPR_B_U_IN, is, imax, jmax, kmax, bcs, s(1, is), v, tmp2, tmp9, tmp5) ! tmp5 contains v transposed
        call OPR_Burgers_Z(OPR_B_U_IN, is, imax, jmax, kmax, bcs, s(1, is), w, tmp3, tmp9, tmp6) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        call TLab_OMP_PARTITION(isize_field, srt, end, siz)
        do ij = srt, end
            hs(ij, is) = hs(ij, is) + tmp1(ij) + tmp2(ij) + tmp3(ij)
        end do
!$omp end parallel

    end do

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

        call TLab_OMP_PARTITION(isize_field, srt, end, siz)
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
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp2)
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp3)
            call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp4)
        end if
        if (stagger_on) then ! staggering on horizontal pressure nodes
            !  Oy derivative
            call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp2, tmp5)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp5, tmp2)
            call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp1)
            !  Ox derivative
            call OPR_PARTIAL_X(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(1), tmp3, tmp5)
            call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp2)
            !  Oz derivative
            call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp4, tmp5)
            call OPR_PARTIAL_Z(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp3)
        else
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp1)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp2)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
        end if

    else
        if (imode_ibm == 1) then
            call IBM_BCS_FIELD(hq(:, 2))
            call IBM_BCS_FIELD(hq(:, 1))
            call IBM_BCS_FIELD(hq(:, 3))
        end if
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 2), tmp2)
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 1), tmp3)
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, hq(:, 3), tmp4)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp1)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp2)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
        else
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), hq(:, 2), tmp1)
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), hq(:, 1), tmp2)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), hq(:, 3), tmp3)
        end if

    end if

    ! -----------------------------------------------------------------------
!$omp parallel default( shared ) private( ij,srt,end,siz )
    call TLab_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
    end do
!$omp end parallel

    ! -----------------------------------------------------------------------
    ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    ! Stagger also Bcs
    if (imode_ibm == 1) call IBM_BCS_FIELD(hq(:, 2))
    if (stagger_on) then ! todo: only need to stagger upper/lower boundary plane, not full h2-array
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), hq(:, 2), tmp5)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)
        if (imode_ibm == 1) call IBM_BCS_FIELD_STAGGER(tmp4)
        p_bcs(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    else
        p_bcs(1:imax, 1:jmax, 1:kmax) => hq(1:imax*jmax*kmax, 2)
    end if

    if (nse_eqns == DNS_EQNS_ANELASTIC) then
        BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)*rbackground(1)
        BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)*rbackground(g(2)%size)
    else
        BcsFlowJmin%ref(:, :, 2) = p_bcs(:, 1, :)
        BcsFlowJmax%ref(:, :, 2) = p_bcs(:, jmax, :)
    end if

    ! pressure in tmp1, Oy derivative in tmp3
    ! select case (imode_elliptic)
    ! case (FDM_COM6_JACOBIAN)
    !     call OPR_Poisson_FourierXZ_Factorize(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    ! case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
    !     call OPR_Poisson_FourierXZ_Direct(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)
    ! end select
    call OPR_Poisson(imax, jmax, kmax, g, BCS_NN, tmp1, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), tmp3)

    ! filter pressure p and its vertical gradient dpdy
    if (any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        call OPR_FILTER(imax, jmax, kmax, PressureFilter, tmp1, txc(1:isize_field,4:6))
        call OPR_FILTER(imax, jmax, kmax, PressureFilter, tmp3, txc(1:isize_field,4:6))
    end if

    ! Saving pressure for towers to tmp array
    if (rkm_substep == rkm_endstep) then
        if (stagger_on .and. ( use_tower .or. PhAvg%active )) then ! Stagger pressure field back on velocity grid (only for towers)
            call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), tmp1, tmp5)
            call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp5, tmp4)
        endif
        if ( use_tower ) &
            call DNS_TOWER_ACCUMULATE(tmp4, 4, wrk1d)
        if ( PhAvg%active) then   
            if (mod((itime+1),PhAvg%stride) == 0)  then
                call AvgPhaseSpace(wrk2d, 1, (itime+1)/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, tmp4)
            end if
        end if
    end if

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
    if (nse_eqns == DNS_EQNS_ANELASTIC) then
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
        call TLab_OMP_PARTITION(isize_field, srt, end, siz)

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
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), tmp1)
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
                                        BcsScalJmin%ref(1, 1, is), BcsScalJmax%ref(1, 1, is), tmp1)
        end if

        if (BcsScalJmin%type(is) /= DNS_SFC_STATIC .or. &
            BcsScalJmax%type(is) /= DNS_SFC_STATIC) then
            call BOUNDARY_BCS_SURFACE_Y(is, bcs, s, hs, tmp1, tmp2)
        end if
        if (imode_ibm == 1) call IBM_BCS_FIELD(hs(1, is)) ! set tendency in solid to zero

        p_bcs(1:imax, 1:jmax, 1:kmax) => hs(1:imax*jmax*kmax, is)
        p_bcs(:, 1, :) = BcsScalJmin%ref(:, :, is)
        p_bcs(:, jmax, :) = BcsScalJmax%ref(:, :, is)

    end do

#ifdef TRACE_ON
    call TLab_Write_ASCII(tfile, 'LEAVING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
#endif

    return
end subroutine RHS_GLOBAL_INCOMPRESSIBLE_1
