#include "dns_error.h"
#include "dns_const.h"

subroutine PRESSURE_MEAN(p, T, s, wrk1d)
    use TLAB_TYPES, only: cp, ci
    use TLAB_CONSTANTS, only: efile
    use TLAB_VARS, only: g
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: rbg, pbg, tbg, sbg
    use TLAB_VARS, only: buoyancy
    use TLAB_PROCS
    use THERMO_VARS, only: imixture

    implicit none

    real(cp), dimension(imax, jmax, kmax), intent(OUT) :: p
    real(cp), dimension(imax, jmax, kmax), intent(INOUT) :: T
    real(cp), dimension(imax, jmax, kmax, *), intent(INOUT) :: s
    real(cp), dimension(jmax, *), intent(INOUT) :: wrk1d

! -------------------------------------------------------------------
    integer(ci) j, iprof_loc
    real(cp) pmin, pmax, ycenter
    real(cp) PROFILES

    real(cp), dimension(:), pointer :: y, dy

! ###################################################################
! Define pointers
    y => g(2)%nodes; dy => g(2)%jac(:, 1)

! ###################################################################
! Constant pressure
! ###################################################################
    if (buoyancy%type == EQNS_NONE) then
        p = pbg%mean

! ###################################################################
! Hydrostatic equilibrium
! ###################################################################
    else

#define p_loc(i)       wrk1d(i,1)
#define r_loc(i)       wrk1d(i,2)
#define t_loc(i)       wrk1d(i,3)
#define ep_loc(i)      wrk1d(i,4)
#define h_loc(i)       wrk1d(i,5)
#define z1_loc(i)      wrk1d(i,6)
#define z2_loc(i)      wrk1d(i,7)
#define z3_loc(i)      wrk1d(i,8)
#define wrk1d_loc(i)   wrk1d(i,9)

! -------------------------------------------------------------------
! Temperature profile is given
! -------------------------------------------------------------------
        if (rbg%type == PROFILE_NONE) then

! AIRWATER case: temperature/mixture profile is given
            if (imixture == MIXT_TYPE_AIRWATER .and. tbg%type > 0) then
                do j = 1, jmax
                    ycenter = y(1) + g(2)%scale*tbg%ymean
                    t_loc(j) = PROFILES(tbg%type, tbg%thick, tbg%delta, tbg%mean, ycenter, tbg%parameters, y(j))

                    ycenter = y(1) + g(2)%scale*sbg(1)%ymean
               z1_loc(j) = PROFILES(sbg(1)%type, sbg(1)%thick, sbg(1)%delta, sbg(1)%mean, ycenter, sbg(1)%parameters, g(2)%nodes(j))

                end do
                ! CALL FI_HYDROSTATIC_AIRWATER_T&
                !      (y, dy, z1_loc(1), t_loc(1), p_loc(1), r_loc(1), wrk1d_loc(1), wrk2d, wrk3d)
                call TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Hydrostatic equilibrium 1 undeveloped')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
                do j = 1, jmax
                    s(:, j, :, 1) = z1_loc(j)
                    s(:, j, :, 2) = z2_loc(j)
                    T(:, j, :) = t_loc(j)
                end do

! AIRWATER case: enthalpy/mixture profile is given
            else if (imixture == MIXT_TYPE_AIRWATER .and. tbg%type < 0) then
                do j = 1, jmax
                    ycenter = y(1) + g(2)%scale*tbg%ymean
                    iprof_loc = -tbg%type
                    z1_loc(j) = PROFILES(iprof_loc, tbg%thick, tbg%delta, tbg%mean, ycenter, tbg%parameters, y(j))

                    ycenter = y(1) + g(2)%scale*sbg(1)%ymean
               z2_loc(j) = PROFILES(sbg(1)%type, sbg(1)%thick, sbg(1)%delta, sbg(1)%mean, ycenter, sbg(1)%parameters, g(2)%nodes(j))

                end do
!           CALL FI_HYDROSTATIC_H_OLD(jmax, y, z1_loc(1), ep_loc(1), t_loc(1), p_loc(1), wrk1d_loc(1))
                call FI_HYDROSTATIC_H(g(2), z1_loc(1), ep_loc(1), t_loc(1), p_loc(1), wrk1d_loc(1))
                do j = 1, jmax
                    s(:, j, :, 1) = z2_loc(j)
                    s(:, j, :, 2) = z3_loc(j)
                    T(:, j, :) = t_loc(j)
                end do

! General case: temperature/mixture profile is given
            else
                ycenter = y(1) + tbg%ymean*g(2)%scale
!           CALL FI_HYDROSTATIC(i1, jmax, i1, ycenter, y, p_loc(1))
                call TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Hydrostatic equilibrium 2 undeveloped')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)

                do j = 1, jmax
                    p_loc(j) = pbg%mean*EXP(p_loc(j))
                end do

            end if

! -------------------------------------------------------------------
! Density profile is given
! -------------------------------------------------------------------
        else
            call TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Density case undeveloped')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

! -------------------------------------------------------------------
! 3D array. Simple case of g parallel and opposite to OY
! -------------------------------------------------------------------
        do j = 1, jmax
            p(:, j, :) = p_loc(j)
        end do

    end if

! ###################################################################
! Control
! ###################################################################
    call MINMAX(imax, jmax, kmax, p, pmin, pmax)

    if (pmin < 0.0_cp .or. pmax < 0.0_cp) then
        call TLAB_WRITE_ASCII(efile, 'PRESSURE_MEAN. Negative pressure.')
        call TLAB_STOP(DNS_ERROR_NEGPRESS)
    end if

    return
end subroutine PRESSURE_MEAN
