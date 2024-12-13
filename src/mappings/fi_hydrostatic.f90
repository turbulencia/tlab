#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
! Compute hydrostatic equilibrium from profiles s=(h,q_t).
! Evaluate the integral \int_pbg%ymean^y dx/H(x), where H(x) is the scale height in the system
!########################################################################
subroutine FI_HYDROSTATIC_H(g, s, ep, T, p, wrk1d)
    use TLab_Constants, only: wp, wi, BCS_MIN, BCS_BOTH
    use FDM, only: grid_dt
    use TLAB_VARS, only: inb_scal_array
    use TLAB_VARS, only: imode_eqns
    use TLAB_VARS, only: pbg, damkohler, buoyancy
    use Thermodynamics, only: imixture, GRATIO, scaleheightinv
    use THERMO_ANELASTIC
    use THERMO_AIRWATER
    use THERMO_THERMAL
    use OPR_ODES

    implicit none

    type(grid_dt), intent(in) :: g
    real(wp), dimension(g%size, inb_scal_array), intent(inout) :: s      ! We calculate equilibrium composition
    real(wp), dimension(g%size), intent(out) :: ep, T, p
    real(wp), dimension(g%size, 2), intent(inout) :: wrk1d

    ! -------------------------------------------------------------------
    integer(wi) iter, niter, j, jcenter
    real(wp) dummy

    ! ###################################################################
    ! Get the center
    do j = 1, g%size
        if (g%nodes(j) <= pbg%ymean .and. &
            g%nodes(j + 1) > pbg%ymean) then
            jcenter = j
            exit
        end if
    end do

#define p_aux(i)        wrk1d(i,1)
#define r_aux(i)        wrk1d(i,2)

    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        ep(:) = (g%nodes - pbg%ymean)*GRATIO*scaleheightinv
        epbackground(:) = ep(:)
    end if

! Setting the pressure entry to 1 to get 1/RT
    p_aux(:) = 1.0_wp

    niter = 10

    p(:) = pbg%mean             ! initialize iteration
    if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then       ! Get ql, if necessary
        s(:, 3) = 0.0_wp
    end if
    do iter = 1, niter           ! iterate
        if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
            pbackground(:) = p_aux(:)
            call THERMO_ANELASTIC_DENSITY(1, g%size, 1, s, r_aux(:))            ! Get r_aux=1/RT
            r_aux(:) = -scaleheightinv*r_aux(:)
        else
            call THERMO_AIRWATER_PH_RE(g%size, s(1, 2), p, s(1, 1), T)
            call THERMO_THERMAL_DENSITY(g%size, s(:, 2), p_aux(:), T, r_aux(:)) ! Get r_aux=1/RT
            r_aux(:) = buoyancy%vector(2)*r_aux(:)
        end if

        p(1) = 0.0_wp
        call OPR_Integral1(1, g, r_aux(:), p, BCS_MIN)

        ! Calculate pressure and normalize s.t. p=pbg%mean at y=pbg%ymean
        p(:) = exp(p(:))
        if (abs(pbg%ymean - g%nodes(jcenter)) == 0.0_wp) then
            dummy = p(jcenter)
        else
            dummy = p(jcenter) + (p(jcenter + 1) - p(jcenter)) &
                    /(g%nodes(jcenter + 1) - g%nodes(jcenter))*(pbg%ymean - g%nodes(jcenter))
        end if
        dummy = pbg%mean/dummy
        p(:) = dummy*p(:)

        if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
            pbackground(:) = p(:)
            if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
                call THERMO_ANELASTIC_PH(1, g%size, 1, s(1, 2), s(1, 1))
            else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                call THERMO_AIRWATER_LINEAR(g%size, s, s(:, inb_scal_array))
            end if
            call THERMO_ANELASTIC_TEMPERATURE(1, g%size, 1, s, T)
        else
            if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
                call THERMO_AIRWATER_PH_RE(g%size, s(1, 2), p, s(1, 1), T)
            end if
        end if

    end do

#undef p_aux
#undef r_aux

    return
end subroutine FI_HYDROSTATIC_H
