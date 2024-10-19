#include "dns_const.h"

subroutine SCAL_MEAN(is, s)

    use TLAB_TYPES, only: wp, wi
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imax, jmax, kmax, inb_wrk2d
    use TLAB_VARS, only: imode_sim
    use TLAB_VARS, only: pbg, rbg, tbg, sbg, qbg
    use TLAB_ARRAYS, only: wrk2d
    use TLAB_POINTERS_3D, only: p_wrk1d
    use THERMO_THERMAL
    use Profiles
    implicit none

    integer(wi) is
    real(wp), intent(OUT) :: s(imax, jmax, kmax)
    ! -------------------------------------------------------------------
    integer(wi) j, k
    real(wp) dummy

    real(wp), pointer :: p_wrk2d(:, :, :)    ! this one is planes OXY instead of the global OXZ

    !########################################################################
    if (imode_sim == DNS_MODE_TEMPORAL) then
        do j = 1, jmax
            dummy = Profiles_Calculate(sbg(is), g(2)%nodes(j))
            s(:, j, :) = dummy + s(:, j, :)
        end do

    else if (imode_sim == DNS_MODE_SPATIAL .and. rbg%type == PROFILE_NONE) then ! temperature/mixture profile are given
        p_wrk2d(1:imax, 1:jmax, 1:inb_wrk2d) => wrk2d(1:imax*jmax*inb_wrk2d, 1)

#define rho_vi(j) p_wrk1d(j,1)
#define u_vi(j)   p_wrk1d(j,2)
#define z_vi(j)   p_wrk1d(j,3)
#define aux1(j)   p_wrk1d(j,4)
#define aux2(j)   p_wrk1d(j,5)
#define aux3(j)   p_wrk1d(j,6)
#define aux4(j)   p_wrk1d(j,7)
#define rho_loc(i,j) p_wrk2d(i,j,1)
#define p_loc(i,j)   p_wrk2d(i,j,2)
#define u_loc(i,j)   p_wrk2d(i,j,3)
#define v_loc(i,j)   p_wrk2d(i,j,4)
#define t_loc(i,j)   p_wrk2d(i,j,5)
        ! Inflow profile of scalar
        do j = 1, jmax
            z_vi(j) = Profiles_Calculate(sbg(is), g(2)%nodes(j))
        end do

        ! Initialize density field
        rho_vi(1:jmax) = 0.0_wp
        do j = 1, jmax
            dummy = Profiles_Calculate(tbg, g(2)%nodes(j))
            ! pilot to be added: ijet_pilot, rjet_pilot_thickness, XIST
            t_loc(:, j) = dummy
        end do
        ! the species array here is wrong for multispecies case !!!
        p_loc(:, :) = pbg%mean
        call THERMO_THERMAL_DENSITY(imax*jmax, s, p_loc(:, :), t_loc(:, :), rho_loc(:, :))

        ! Inflow profile of density
        rho_vi(:) = rho_loc(1, :)

        ! inflow profile of velocity
        u_vi(1:jmax) = 0.0_wp
        do j = 1, jmax
            u_vi(j) = Profiles_Calculate(qbg(1), g(2)%nodes(j))
            ! pilot to be added: ijet_pilot, rjet_pilot_thickness, rjet_pilot_velocity
        end do

        ! 2D distributions of density and velocity
        if (rbg%delta /= 0.0_wp) then
            call FLOW_SPATIAL_DENSITY(imax, jmax, tbg, qbg(1), &
                                      g(1)%nodes, g(2)%nodes, s, p_loc(1, 1), rho_vi(1), u_vi(1), aux1(1), rho_loc(1, 1), &
                                      aux2(1), aux3(1), aux4(1))
        end if
        call FLOW_SPATIAL_VELOCITY(imax, jmax, qbg(1), qbg(1)%diam, &
                                   qbg(1)%parameters(2), qbg(1)%parameters(3), qbg(1)%parameters(4), &
                     g(1)%nodes, g(2)%nodes, rho_vi(1), u_vi(1), rho_loc(1, 1), u_loc(1, 1), v_loc(1, 1), aux1(1), p_wrk2d(1, 1, 6))
        ! 2D distribution of scalar
        call FLOW_SPATIAL_SCALAR(imax, jmax, sbg(is), sbg(is)%diam, sbg(is)%diam, &
                                 sbg(is)%parameters(2), sbg(is)%parameters(3), sbg(is)%parameters(4), &
                                 g(1)%nodes, g(2)%nodes, rho_vi(1), u_vi(1), z_vi(1), rho_loc(1, 1), u_loc(1, 1), s, p_wrk2d(1, 1, 6))
        if (g(3)%size > 1) then
            do k = 2, kmax
                s(:, :, k) = s(:, :, 1)
            end do
        end if

    end if

    return
end subroutine SCAL_MEAN
