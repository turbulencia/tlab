#include "dns_const.h"

!########################################################################
!#
!# Calculate the pressure field from a divergence free velocity field and a body force.
!#
!########################################################################
subroutine FI_PRESSURE_BOUSSINESQ(q, s, p, tmp1, tmp2, tmp)
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: imode_eqns, imode_ibm, istagger
    use TLAB_VARS, only: rbackground, PressureFilter
    use TLAB_ARRAYS, only: wrk2d, wrk3d
    use TLAB_POINTERS_3D, only: p_wrk2d
    use IBM_VARS, only: ibm_burgers
    use OPR_PARTIAL
    use OPR_BURGERS
    use OPR_ELLIPTIC
    use FI_SOURCES
    use OPR_FILTERS

    implicit none

    real(wp), intent(in) :: q(isize_field, 3)
    real(wp), intent(in) :: s(isize_field, *)
    real(wp), intent(out) :: p(isize_field)
    real(wp), intent(inout) :: tmp1(isize_field), tmp2(isize_field)
    real(wp), intent(inout) :: tmp(isize_field, 3)

    target q, tmp
! -----------------------------------------------------------------------
    integer(wi) bcs(2, 2)
    integer, parameter :: i3 = 3

! Pointers to existing allocated space
    real(wp), dimension(:), pointer :: u, v, w
    real(wp), dimension(:), pointer :: tmp3, tmp4, tmp5
    real(wp), dimension(:, :, :), pointer :: p_bcs

! #######################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    p = 0.0_wp
    tmp = 0.0_wp

! Define pointers
    u => q(:, 1)
    v => q(:, 2)
    w => q(:, 3)

! #######################################################################
! Sources
    call FI_SOURCES_FLOW(q, s, tmp, tmp1)

    tmp3 => tmp(:, 1)
    tmp4 => tmp(:, 2)
    tmp5 => tmp(:, 3)

! If IBM, then use modified fields for derivatives
    if (imode_ibm == 1) ibm_burgers = .true.

! Advection and diffusion terms
    call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), u, u, u, p, tmp1) ! store u transposed in tmp1
    tmp3 = tmp3 + p
    call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), v, u, tmp1, p, tmp2) ! tmp1 contains u transposed
    tmp4 = tmp4 + p
    call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), w, u, tmp1, p, tmp2) ! tmp1 contains u transposed
    tmp5 = tmp5 + p

    call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), v, v, v, p, tmp1) ! store v transposed in tmp1
    tmp4 = tmp4 + p
    call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), u, v, tmp1, p, tmp2) ! tmp1 contains v transposed
    tmp3 = tmp3 + p
    call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), w, v, tmp1, p, tmp2) ! tmp1 contains v transposed
    tmp5 = tmp5 + p

    call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), w, w, w, p, tmp1) ! store w transposed in tmp1
    tmp5 = tmp5 + p
    call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), v, w, tmp1, p, tmp2) ! tmp1 contains w transposed
    tmp4 = tmp4 + p
    call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), u, w, tmp1, p, tmp2) ! tmp1 contains w transposed
    tmp3 = tmp3 + p

! If IBM, set flag back to false
    if (imode_ibm == 1) ibm_burgers = .false.

! Set p-field back to zero
    p = 0.0_wp

! Apply IBM BCs
    if (imode_ibm == 1) then
        call IBM_BCS_FIELD(tmp3)
        call IBM_BCS_FIELD(tmp4)
        call IBM_BCS_FIELD(tmp5)
    end if

! Calculate forcing term Ox
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp3)
    end if
    if (istagger == 1) then
        call OPR_PARTIAL_X(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(1), tmp3, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
    else
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp1, wrk3d, wrk2d, wrk3d)
    end if
    p = p + tmp1

! Calculate forcing term Oy
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp4)
    end if
    if (istagger == 1) then
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp4, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp3, tmp1, wrk3d, wrk2d, wrk3d)
    else
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp1, wrk3d, wrk2d, wrk3d)
    end if
    p = p + tmp1

! Calculate forcing term Oz
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp5)
    end if
    if (istagger == 1) then
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp5, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
    else
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp5, tmp1, wrk3d, wrk2d, wrk3d)
    end if
    p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    if (istagger == 1) then ! todo: only need to stagger upper/lower boundary plane, not full h2-array
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp4, tmp5, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp4, wrk3d, wrk2d, wrk3d)
        if (imode_ibm == 1) call IBM_BCS_FIELD_STAGGER(tmp4)
    end if
    p_bcs(1:imax, 1:jmax, 1:kmax) => tmp4(1:imax*jmax*kmax)
    p_wrk2d(:, :, 1) = p_bcs(:, 1, :)
    p_wrk2d(:, :, 2) = p_bcs(:, jmax, :)

! Pressure field in p
    call OPR_POISSON_FXZ(imax, jmax, kmax, g, i3, p, tmp1, tmp2, p_wrk2d(:, :, 1), p_wrk2d(:, :, 2))

! filter pressure and dpdy
    if      (PressureFilter(1)%type /= DNS_FILTER_NONE) then
        call OPR_FILTER_X(imax, jmax, kmax, PressureFilter(1), tmp1); call OPR_FILTER_X(imax, jmax, kmax, PressureFilter(1), tmp3)
    else if (PressureFilter(2)%type /= DNS_FILTER_NONE) then
        call OPR_FILTER_Y(imax, jmax, kmax, PressureFilter(2), tmp1); call OPR_FILTER_Y(imax, jmax, kmax, PressureFilter(2), tmp3)
    else if (PressureFilter(3)%type /= DNS_FILTER_NONE) then
        call OPR_FILTER_Z(imax, jmax, kmax, PressureFilter(3), tmp1); call OPR_FILTER_Z(imax, jmax, kmax, PressureFilter(3), tmp3)
    end if

! Stagger pressure field p back on velocity grid
    if (istagger == 1) then
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), p, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp1, p, wrk3d, wrk2d, wrk3d)
    end if

    nullify (u, v, w, tmp3, tmp4, tmp5, p_bcs)

    return
end subroutine FI_PRESSURE_BOUSSINESQ
