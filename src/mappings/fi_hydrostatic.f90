#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
! Compute hydrostatic equilibrium from profiles s=(h,q_t).
! Evaluate the integral \int_pbg%ymean^y dx/H(x), where H(x) is the scale height in the system
!########################################################################
subroutine FI_HYDROSTATIC_H(g, s, e, T, p, wrk1d)
    use TLab_Constants, only: wp, wi, BCS_MIN, BCS_BOTH
    use TLab_Types, only: grid_dt
    use TLAB_VARS, only: imode_eqns
    use TLAB_VARS, only: pbg, damkohler, buoyancy
    use Thermodynamics, only: imixture, scaleheight
    use THERMO_ANELASTIC
    use THERMO_AIRWATER
    use THERMO_THERMAL
    use OPR_ODES

    implicit none

    type(grid_dt), intent(IN) :: g
    real(wp), dimension(g%size), intent(IN) :: e
    real(wp), dimension(g%size), intent(OUT) :: T, p
    real(wp), dimension(g%size, *), intent(INOUT) :: s      ! We calculate equilibrium composition
    real(wp), dimension(g%size, *), intent(INOUT) :: wrk1d

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

    ! Setting the pressure entry to 1 to get 1/RT
    p_aux(:) = 1.0_wp

    niter = 10

    p(:) = pbg%mean             ! initialize iteration
    if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then       ! Get ql, if necessary
        s(:, 3) = 0.0_wp
    end if
    do iter = 1, niter           ! iterate
        if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
            epbackground(:) = e(:)
            pbackground(:) = p_aux(:)
            call THERMO_ANELASTIC_DENSITY(1, g%size, 1, s, r_aux(:))   ! Get r_aus=1/RT
            dummy = -1.0_wp/sign(scaleheight, buoyancy%vector(2))
        else
            call THERMO_AIRWATER_PH_RE(g%size, s(1, 2), p, s(1, 1), T)
            call THERMO_THERMAL_DENSITY(g%size, s(:, 2), p_aux(:), T, r_aux(:)) ! Get r_aux=1/RT
            dummy = buoyancy%vector(2)
        end if
        r_aux(:) = dummy*r_aux(:)

        p(1) = 0.0_wp
        call OPR_Integral1(1, g, r_aux(:), p, BCS_MIN)

        ! Calculate pressure and normalize s.t. p=pbg%mean at y=pbg%ymean_rel
        p(:) = exp(p(:))
        if (abs(pbg%ymean - g%nodes(jcenter)) == 0.0_wp) then
            dummy = p(jcenter)
        else
            dummy = p(jcenter) + (p(jcenter + 1) - p(jcenter)) &
                    /(g%nodes(jcenter + 1) - g%nodes(jcenter))*(pbg%ymean - g%nodes(jcenter))
        end if
        dummy = pbg%mean/dummy
        p(:) = dummy*p(:)

        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then ! Get ql, if necessary
            if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
                epbackground(:) = e(:)
                pbackground(:) = p(:)
                call THERMO_ANELASTIC_PH(1, g%size, 1, s(1, 2), s(1, 1))
            else
                call THERMO_AIRWATER_PH_RE(g%size, s(1, 2), p, s(1, 1), T)
            end if
        end if

    end do

#undef p_aux
#undef r_aux

    ! compute equilibrium values of T
    if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
        epbackground(:) = e(:)
        call THERMO_ANELASTIC_TEMPERATURE(1, g%size, 1, s, T)
    end if

    return
end subroutine FI_HYDROSTATIC_H
