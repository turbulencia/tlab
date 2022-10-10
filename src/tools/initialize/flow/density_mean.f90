#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Setting the density field.
!# The default case is when temperature field is set similar to a passive
!# scalar, and the density is then derived taking pressure constant.
!# The inputs are however density maximum/minimum:
!#
!# T = T_0 + s(T_1-T_0), s.t. 0<=s<=1, =>
!#           => rho = W / [ W_0/rho_0 + s(W_1/rho_1-W_0/rho_0) ]
!#
!# This allow to compare the evolution of the passive scalar with that
!# of the temperature.
!#
!# If volumetric forces are present, then rho is computed from
!# dp/dy = rho g, to ensure hydrostatic equilibrium including numerical
!# errors.
!#
!# Spatial case with mulstispecies is not jet done
!#
!########################################################################
subroutine DENSITY_MEAN(rho, p, T, s, txc, wrk1d, wrk2d, wrk3d)
    use TLAB_CONSTANTS, only: wp, wi, efile
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imode_sim, inb_scal, imax, jmax, kmax
    use TLAB_VARS, only: rbg, tbg, sbg, qbg
    use TLAB_VARS, only: buoyancy
    use THERMO_VARS, only: imixture
    use PROFILES
    implicit none

    real(wp), dimension(imax, jmax, kmax), intent(IN) :: p, T
    real(wp), dimension(imax, jmax, kmax), intent(OUT) :: rho
    real(wp), dimension(imax, jmax, kmax, *), intent(OUT) :: s
    real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: txc, wrk3d
    real(wp), dimension(jmax, *), intent(INOUT) :: wrk1d, wrk2d

    ! -------------------------------------------------------------------
    real(wp) dummy
    integer(wi) j, k, is, bcs(2, 2)

    bcs = 0

    ! -------------------------------------------------------------------
    ! Temporal shear layer case without volumetric force:
    ! Calculate density from equation of state
    ! -------------------------------------------------------------------
    if (imode_sim == DNS_MODE_TEMPORAL) then
        if (buoyancy%type == EQNS_NONE) then

#define TEM_MEAN_LOC(i,j,k) wrk3d(i,j,k)
#define RHO_MEAN_LOC(i,j,k) txc(i,j,k)

            ! temperature/mixture profile are given
            if (rbg%type == PROFILE_NONE) then
                do j = 1, jmax
                    dummy = PROFILES_CALCULATE(tbg, g(2)%nodes(j))
                    TEM_MEAN_LOC(:, j, :) = dummy
                end do

                do is = 1, inb_scal
                    do j = 1, jmax
                        dummy = PROFILES_CALCULATE(sbg(is), g(2)%nodes(j))
                        s(:, j, :, is) = dummy
                    end do
                end do

                ! define liquid content in AirWater case: (p,T) given
                if (imixture == MIXT_TYPE_AIRWATER) then
                    call THERMO_AIRWATER_PT(imax, jmax, kmax, s, p, TEM_MEAN_LOC(1, 1, 1))
                end if

                call THERMO_THERMAL_DENSITY(imax, jmax, kmax, s, p, TEM_MEAN_LOC(1, 1, 1), RHO_MEAN_LOC(1, 1, 1))
                rho(:, :, :) = rho(:, :, :) + RHO_MEAN_LOC(:, :, :)

                ! density profile itself is given
            else
                do j = 1, jmax
                    dummy = PROFILES_CALCULATE(rbg, g(2)%nodes(j))
                    rho(:, j, :) = rho(:, j, :) + dummy
                end do

            end if

#undef TEM_MEAN_LOC
#undef RHO_MEAN_LOC

            ! -------------------------------------------------------------------
            ! Temporal shear layer case with volumetric force:
            ! Calculate density from hydrostatic equilibrium
            ! assuming a volumetric force along OY
            ! -------------------------------------------------------------------
        else
            ! AIRWATER case. Routine OPR_PARTIAL_Y introduces small errors in equilibrium
            if (imixture == MIXT_TYPE_AIRWATER) then
                call THERMO_THERMAL_DENSITY(imax, jmax, kmax, s, p, T, rho)

                ! General case
            else
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), p, txc, wrk3d, wrk2d, wrk3d)
                dummy = 1.0_wp/buoyancy%vector(2)
                rho(:, :, :) = rho(:, :, :) + txc(:, :, :)*dummy
            end if

        end if

        ! -------------------------------------------------------------------
        ! Spatial jet
        ! Only if there is a density variation. Constant density is already
        ! initialized in previous routine segment.
        ! -------------------------------------------------------------------
    else if (imode_sim == DNS_MODE_SPATIAL) then
        if (rbg%type == PROFILE_NONE) then ! temperature/mixture profile are given
#define rho_vi(j) wrk1d(j,1)
#define u_vi(j)   wrk1d(j,2)
#define aux1(j)   wrk1d(j,3)
#define aux2(j)   wrk1d(j,4)
#define aux3(j)   wrk1d(j,5)
#define aux4(j)   wrk1d(j,6)
            ! Inflow profile of density
            ! rho_vi(:) = rho(1,:,1) ! I need to update this because rho(1,:,1) is now undefined

            ! Inflow profile of axial velocity
            do j = 1, jmax
                u_vi(j) = PROFILES_CALCULATE(qbg(1), g(2)%nodes(j))
            end do

            ! 2D distribution of density
            call FLOW_SPATIAL_DENSITY(imax, jmax, tbg, qbg(1), &
                              g(1)%nodes, g(2)%nodes, s, p, rho_vi(1), u_vi(1), aux1(1), rho, aux2(1), aux3(1), aux4(1))

            do k = 2, kmax
                rho(:, :, k) = rho(:, :, 1)
            end do

        else ! density profile itself is given
            do j = 1, jmax
                dummy = PROFILES_CALCULATE(rbg, g(2)%nodes(j))
                rho(:, j, :) = rho(:, j, :) + dummy
            end do

        end if
    end if

    return
end subroutine DENSITY_MEAN
