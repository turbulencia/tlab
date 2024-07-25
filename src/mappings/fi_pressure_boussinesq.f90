#include "dns_const.h"

!########################################################################
!#
!# Calculate the pressure field from a divergence free velocity field and a body force.
!#
!########################################################################
subroutine FI_PRESSURE_BOUSSINESQ(q, s, p, tmp1, tmp2, tmp)
    use TLAB_CONSTANTS, only: wp, wi, BCS_NN
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: imode_eqns, imode_ibm, imode_elliptic
    use TLAB_VARS, only: PressureFilter, stagger_on
    use TLAB_VARS, only: buoyancy, coriolis
    use TLAB_VARS, only: inb_txc, pdecomposition
    use TLAB_ARRAYS, only: wrk1d
    use TLAB_POINTERS_3D, only: p_wrk2d
    use THERMO_ANELASTIC
    use IBM_VARS, only: ibm_burgers
    use OPR_PARTIAL
    use OPR_BURGERS
    use OPR_ELLIPTIC
    use FI_SOURCES
    use OPR_FILTERS
    use IO_FIELDS

    implicit none

    real(wp), intent(in) :: q(isize_field, 3)
    real(wp), intent(in) :: s(isize_field, *)
    real(wp), intent(out) :: p(isize_field)
    real(wp), intent(inout) :: tmp1(isize_field), tmp2(isize_field)
    real(wp), intent(inout) :: tmp(isize_field, inb_txc - 3)

    target q, tmp, s
! -----------------------------------------------------------------------
    integer(wi) bcs(2, 2)
    integer(wi) iq
    integer(wi) siz, srt, end
    integer(wi) :: i

! -----------------------------------------------------------------------
    real(wp) dummy

! -----------------------------------------------------------------------
#ifdef USE_BLAS
    integer ILEN
#endif
! -----------------------------------------------------------------------

! Pointers to existing allocated space
    real(wp), dimension(:), pointer :: u, v, w
    real(wp), dimension(:), pointer :: tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
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
    if (pdecomposition%name == 'total') then 
        call FI_SOURCES_FLOW(q, s, tmp, tmp1)
    end if

    tmp3 => tmp(:, 1)
    tmp4 => tmp(:, 2)
    tmp5 => tmp(:, 3)

    if (pdecomposition%name /= 'total') then
        tmp6 => tmp(:, 4)
        tmp7 => tmp(:, 5)
        tmp8 => tmp(:, 6)
        tmp9 => tmp(:, 7)
        tmp3 = 0.0_wp
        tmp4 = 0.0_wp
        tmp5 = 0.0_wp
        tmp6 = 0.0_wp
        tmp7 = 0.0_wp
        tmp8 = 0.0_wp
        tmp9 = 0.0_wp
    end if

! If IBM, then use modified fields for derivatives
    if (imode_ibm == 1) ibm_burgers = .true.

    if (pdecomposition%name == 'advdiffu' .OR. pdecomposition%name == 'total' .OR. pdecomposition%name == 'advction') then
        !  Advection and diffusion terms
        call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), u, u, p, tmp1) ! store u transposed in tmp1
        tmp3 = tmp3 + p
        call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), v, u, p, tmp2, tmp1) ! tmp1 contains u transposed
        tmp4 = tmp4 + p
        call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), w, u, p, tmp2, tmp1) ! tmp1 contains u transposed
        tmp5 = tmp5 + p

        call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), v, v, p, tmp1) ! store v transposed in tmp1
        tmp4 = tmp4 + p
        call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), u, v, p, tmp2, tmp1) ! tmp1 contains v transposed
        tmp3 = tmp3 + p
        call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), w, v, p, tmp2, tmp1) ! tmp1 contains v transposed
        tmp5 = tmp5 + p

        call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), w, w, p, tmp1) ! store w transposed in tmp1
        tmp5 = tmp5 + p
        call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), v, w, p, tmp2, tmp1) ! tmp1 contains w transposed
        tmp4 = tmp4 + p
        call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), u, w, p, tmp2, tmp1) ! tmp1 contains w transposed
        tmp3 = tmp3 + p

    end if

    if (pdecomposition%name == 'advction' .OR. pdecomposition%name == 'difusion') then
        tmp9 = 0.0_wp
        ! Sepereating Diffusion
        ! NSE X-Comp
        call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), u, tmp9, p, tmp1)
        tmp6 = tmp6 + p   ! Diffusion d2u/dx2
        call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), v, tmp9, p, tmp2, tmp1)
        tmp7 = tmp7 + p
        call OPR_BURGERS_X(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(1), w, tmp9, p, tmp2, tmp1)
        tmp8 = tmp8 + p

        ! NSE Y-Comp
        call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), v, tmp9, p, tmp1)
        tmp7 = tmp7 + p ! Diffusion d2v/dx2 + d2v/dy2 
        call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), u, tmp9, p, tmp2, tmp1)
        tmp6 = tmp6 + p
        call OPR_BURGERS_Y(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(2), w, tmp9, p, tmp2, tmp1)
        tmp8 = tmp8 + p

        ! NSE Z-Comp
        call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), w, tmp9, p, tmp1)
        tmp8 = tmp8 + p  ! Diffusion d2w/dx2 + d2w/dy2 + d2w/dz2
        call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), v, tmp9, p, tmp2, tmp1)
        tmp7 = tmp7 + p
        call OPR_BURGERS_Z(OPR_B_U_IN, 0, imax, jmax, kmax, bcs, g(3), u, tmp9, p, tmp2, tmp1)
        tmp6 = tmp6 + p

    end if

    if (pdecomposition%name == 'advction') then
        tmp3 = tmp3 - tmp6
        tmp4 = tmp4 - tmp7
        tmp5 = tmp5 - tmp8

    else if (pdecomposition%name == 'difusion') then
        tmp3 = tmp6
        tmp4 = tmp7
        tmp5 = tmp8

    end if

    ! Coriolis Forcing term
    if (pdecomposition%name == 'coriolis') then
        call FI_CORIOLIS(coriolis,imax, jmax, kmax, q, tmp)
        tmp3        => tmp(:, 1)
        tmp4        => tmp(:, 2)
        tmp5        => tmp(:, 3)
    end if

    ! Buoyancy Forcing term
    if (pdecomposition%name == 'buoyancy') then
        do iq = 1, 3
            if (buoyancy%type == EQNS_EXPLICIT) then
                call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, tmp1)
            else
                if (buoyancy%active(iq)) then
                    if (iq == 2) then
                        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, tmp1, bbackground)
                    else
                        wrk1d(:, 1) = 0.0_wp
                        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, tmp1, wrk1d)
                    end if

                    call DNS_OMP_PARTITION(isize_field, srt, end, siz)
                    dummy = buoyancy%vector(iq)
#ifdef USE_BLAS
                    ILEN = isize_field
                    call DAXPY(ILEN, dummy, tmp1(srt), 1, tmp(srt, iq), 1)
#else
                    do i = srt, end
                            tmp(i, iq) = tmp(i, iq) + dummy*tmp1(i)
                    end do
#endif
                end if
            end if
        end do

        tmp3        => tmp(:, 1)
        tmp4        => tmp(:, 2)
        tmp5        => tmp(:, 3)

    end if

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
    if (stagger_on) then
        call OPR_PARTIAL_X(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(1), tmp3, tmp2)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp1)
    else
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp1)
    end if
    p = p + tmp1

! Calculate forcing term Oy
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp4)
    end if
    if (stagger_on) then
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp4, tmp2)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp3, tmp1)
    else
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp1)
    end if
    p = p + tmp1

! Calculate forcing term Oz
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call THERMO_ANELASTIC_WEIGHT_INPLACE(imax, jmax, kmax, rbackground, tmp5)
    end if
    if (stagger_on) then
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp5, tmp2)
        call OPR_PARTIAL_Z(OPR_P1_INT_VP, imax, jmax, kmax, bcs, g(3), tmp2, tmp1)
    else
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp5, tmp1)
    end if
    p = p + tmp1

! #######################################################################
! Solve Poisson equation
! #######################################################################
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    if (stagger_on) then ! todo: only need to stagger upper/lower boundary plane, not full h2-array
        call OPR_PARTIAL_X(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(1), tmp4, tmp5)
        call OPR_PARTIAL_Z(OPR_P0_INT_VP, imax, jmax, kmax, bcs, g(3), tmp5, tmp4)
        if (imode_ibm == 1) call IBM_BCS_FIELD_STAGGER(tmp4)
    end if
    p_bcs(1:imax, 1:jmax, 1:kmax) => tmp4(1:imax*jmax*kmax)
    p_wrk2d(:, :, 1) = p_bcs(:, 1, :)
    p_wrk2d(:, :, 2) = p_bcs(:, jmax, :)

! Pressure field in p
    select case (imode_elliptic)
    case (FDM_COM6_JACOBIAN)
        call OPR_POISSON_FXZ(imax, jmax, kmax, g, BCS_NN, p, tmp1, tmp2, p_wrk2d(:, :, 1), p_wrk2d(:, :, 2))

    case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
        call OPR_POISSON_FXZ_D(imax, jmax, kmax, g, BCS_NN, p, tmp1, tmp2, p_wrk2d(:, :, 1), p_wrk2d(:, :, 2))
    end select

    ! filter pressure p
    if (any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        call OPR_FILTER(imax, jmax, kmax, PressureFilter, p, tmp)
    end if

! Stagger pressure field p back on velocity grid
    if (stagger_on) then
        call OPR_PARTIAL_Z(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(3), p, tmp1)
        call OPR_PARTIAL_X(OPR_P0_INT_PV, imax, jmax, kmax, bcs, g(1), tmp1, p)
    end if

    nullify (u, v, w, tmp3, tmp4, tmp5, p_bcs)

    return
end subroutine FI_PRESSURE_BOUSSINESQ
