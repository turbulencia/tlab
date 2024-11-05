#include "dns_const.h"
#include "dns_error.h"

module FLOW_MEAN
    use TLab_Constants, only: wp, wi, efile
    use FDM, only: g
    use TLAB_VARS, only: imode_sim, imax, jmax, kmax, inb_scal
    use TLAB_VARS, only: qbg, pbg, rbg, tbg, hbg, sbg
    use TLAB_VARS, only: coriolis, buoyancy
    use TLab_Arrays, only: wrk1d
    use TLab_Pointers_3D, only: p_wrk1d, p_wrk3d
    use TLab_WorkFlow
    use Thermodynamics, only: imixture
    use THERMO_THERMAL
    use THERMO_AIRWATER
    use Profiles
    use OPR_PARTIAL
    implicit none
    private

    public :: VELOCITY_MEAN
    public :: PRESSURE_MEAN
    public :: DENSITY_MEAN

contains

    ! ###################################################################
    subroutine VELOCITY_MEAN(u, v, w)
        real(wp), dimension(imax, jmax, kmax), intent(OUT) :: u, v, w

        ! -------------------------------------------------------------------
        integer(wi) iq, j
        real(wp) calpha, salpha

        !########################################################################
        if (imode_sim == DNS_MODE_TEMPORAL) then

            ! Construct reference profiles into array wrk1d
            do iq = 1, 3
                do j = 1, jmax
                    wrk1d(j, iq) = Profiles_Calculate(qbg(iq), g(2)%nodes(j))
                end do
            end do

            ! Construct velocity field
            if (coriolis%type == EQNS_COR_NORMALIZED) then
                calpha = cos(coriolis%parameters(1)); salpha = sin(coriolis%parameters(1))
                wrk1d(:, 3) = wrk1d(:, 3)*sign(1.0_wp, coriolis%vector(2)) ! right angular velocity vector (Garratt, 1992, p.42)

                do j = 1, jmax
                    u(:, j, :) = u(:, j, :) + wrk1d(j, 1)*calpha + wrk1d(j, 3)*salpha
                    v(:, j, :) = v(:, j, :) + wrk1d(j, 2)
                    w(:, j, :) = w(:, j, :) - wrk1d(j, 1)*salpha + wrk1d(j, 3)*calpha
                end do

            else
                do j = 1, jmax
                    u(:, j, :) = u(:, j, :) + wrk1d(j, 1)
                    v(:, j, :) = v(:, j, :) + wrk1d(j, 2)
                    w(:, j, :) = w(:, j, :) + wrk1d(j, 3)
                end do

            end if

            ! -------------------------------------------------------------------
        else if (imode_sim == DNS_MODE_SPATIAL) then
! I need to pass rho; need to check
! #define rho_vi(j) wrk1d(j,1)
! #define u_vi(j)   wrk1d(j,2)
! #define aux(j)    wrk1d(j,3)
!     DO j = 1,jmax
!       u_vi(j) = Profiles_Calculate( qbg(1), g(2)%nodes(j) )
!     ENDDO
!     rho_vi(:) = rho(1,:,1)
!
!     CALL FLOW_SPATIAL_VELOCITY(imax,jmax, qbg(1), qbg(1)%diam, &
!        qbg(1)%parameters(2), qbg(1)%parameters(3), qbg(1)%parameters(4), &
!     g(1)%nodes, g(2)%nodes, rho_vi(1), u_vi(1), rho, u, v, aux(1), wrk3d)
!     IF ( g(3)%size .GT. 1 ) THEN
!       DO k = 2,kmax
!         u(:,:,k) = u(:,:,1)
!         v(:,:,k) = v(:,:,1)
!       ENDDO
!       w = w + qbg(3)%mean
!     ENDIF
! #undef rho_vi
! #undef u_vi
! #undef aux

        end if

        ! -------------------------------------------------------------------
        if (g(3)%size == 1) w = 0.0_wp

        return
    end subroutine VELOCITY_MEAN

    ! ###################################################################
    subroutine PRESSURE_MEAN(p, T, s)
        real(wp), dimension(imax, jmax, kmax), intent(OUT) :: p
        real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: T
        real(wp), dimension(imax, jmax, kmax, *), intent(INOUT) :: s

        ! -------------------------------------------------------------------
        integer(wi) j
        real(wp) pmin, pmax

        if (buoyancy%type == EQNS_NONE) then
            p = pbg%mean

        else            ! hydrostatic equilibrium

#define p_loc(i)       p_wrk1d(i,1)
#define r_loc(i)       p_wrk1d(i,2)
#define t_loc(i)       p_wrk1d(i,3)
#define ep_loc(i)      p_wrk1d(i,4)
#define h_loc(i)       p_wrk1d(i,5)
#define z1_loc(i)      p_wrk1d(i,6)
#define z2_loc(i)      p_wrk1d(i,7)
#define z3_loc(i)      p_wrk1d(i,8)
#define wrk1d_loc(i)   p_wrk1d(i,9)

            if (hbg%type /= PROFILE_NONE) then
                select case (imixture)
                case (MIXT_TYPE_AIRWATER)
                    do j = 1, jmax
                        z1_loc(j) = Profiles_Calculate(hbg, g(2)%nodes(j))
                        z2_loc(j) = Profiles_Calculate(sbg(1), g(2)%nodes(j))
                    end do
                    ! ep contains the potential energy but it is not used in the compressible formulation
                    call FI_HYDROSTATIC_H(g(2), z1_loc(1), ep_loc(1), t_loc(1), p_loc(1), wrk1d_loc(1))
                    do j = 1, jmax
                        s(:, j, :, 1) = z2_loc(j)
                        s(:, j, :, 2) = z3_loc(j)
                        T(:, j, :) = t_loc(j)
                        p(:, j, :) = p_loc(j)
                    end do

                case default
                    call TLab_Write_ASCII(efile, 'PRESSURE_MEAN. Mixture case undeveloped.')
                    call TLab_Stop(DNS_ERROR_UNDEVELOP)

                end select

            end if

            if (tbg%type /= PROFILE_NONE) then
                call TLab_Write_ASCII(efile, 'PRESSURE_MEAN. Temperature case undeveloped.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            if (rbg%type /= PROFILE_NONE) then
                call TLab_Write_ASCII(efile, 'PRESSURE_MEAN. Density case undeveloped.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

        end if

        ! Control
        call MINMAX(imax, jmax, kmax, p, pmin, pmax)
        if (pmin < 0.0_wp .or. pmax < 0.0_wp) then
            call TLab_Write_ASCII(efile, 'PRESSURE_MEAN. Negative pressure.')
            call TLab_Stop(DNS_ERROR_NEGPRESS)
        end if

        return
    end subroutine PRESSURE_MEAN

    ! ###################################################################
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
    !# dp/dy = rho g, to ensure hydrostatic equilibrium including numerical errors.
    ! ###################################################################
    subroutine DENSITY_MEAN(rho, p, T, s, txc)
        real(wp), dimension(imax, jmax, kmax), intent(IN) :: p, T
        real(wp), dimension(imax, jmax, kmax), intent(OUT) :: rho
        real(wp), dimension(imax, jmax, kmax, *), intent(OUT) :: s
        real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: txc

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

#define TEM_MEAN_LOC(i,j,k) p_wrk3d(i,j,k)
#define RHO_MEAN_LOC(i,j,k) txc(i,j,k)

                ! temperature/mixture profile are given
                if (rbg%type == PROFILE_NONE) then
                    do j = 1, jmax
                        dummy = Profiles_Calculate(tbg, g(2)%nodes(j))
                        TEM_MEAN_LOC(:, j, :) = dummy
                    end do

                    do is = 1, inb_scal
                        do j = 1, jmax
                            dummy = Profiles_Calculate(sbg(is), g(2)%nodes(j))
                            s(:, j, :, is) = dummy
                        end do
                    end do

                    ! define liquid content in AirWater case: (p,T) given
                    if (imixture == MIXT_TYPE_AIRWATER) then
                        call THERMO_AIRWATER_PT(imax*jmax*kmax, s, p, TEM_MEAN_LOC(:, :, :))
                    end if

                    call THERMO_THERMAL_DENSITY(imax*jmax*kmax, s, p, TEM_MEAN_LOC(:, :, :), RHO_MEAN_LOC(:, :, :))
                    rho(:, :, :) = rho(:, :, :) + RHO_MEAN_LOC(:, :, :)

                    ! density profile itself is given
                else
                    do j = 1, jmax
                        dummy = Profiles_Calculate(rbg, g(2)%nodes(j))
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
                    call THERMO_THERMAL_DENSITY(imax*jmax*kmax, s, p, T, rho)

                    ! General case
                else
                    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), p, txc)
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
#define rho_vi(j) p_wrk1d(j,1)
#define u_vi(j)   p_wrk1d(j,2)
#define aux1(j)   p_wrk1d(j,3)
#define aux2(j)   p_wrk1d(j,4)
#define aux3(j)   p_wrk1d(j,5)
#define aux4(j)   p_wrk1d(j,6)
                ! Inflow profile of density
                ! rho_vi(:) = rho(1,:,1) ! I need to update this because rho(1,:,1) is now undefined

                ! Inflow profile of axial velocity
                do j = 1, jmax
                    u_vi(j) = Profiles_Calculate(qbg(1), g(2)%nodes(j))
                end do

                ! 2D distribution of density
                call FLOW_SPATIAL_DENSITY(imax, jmax, tbg, qbg(1), &
                                          g(1)%nodes, g(2)%nodes, s, p, rho_vi(1), u_vi(1), aux1(1), rho, aux2(1), aux3(1), aux4(1))

                do k = 2, kmax
                    rho(:, :, k) = rho(:, :, 1)
                end do

            else ! density profile itself is given
                do j = 1, jmax
                    dummy = Profiles_Calculate(rbg, g(2)%nodes(j))
                    rho(:, j, :) = rho(:, j, :) + dummy
                end do

            end if
        end if

        return
    end subroutine DENSITY_MEAN

end module FLOW_MEAN
