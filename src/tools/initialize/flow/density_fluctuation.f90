#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Setting up a perturbation of the thermodynamic fields by a
!# displacement of the reference center plane.
!#
!# Array s enters with the scalar total field, including fluctuations.
!#
!########################################################################
subroutine DENSITY_FLUCTUATION(code, s, p, rho, T, h, disp)
    use TLAB_TYPES, only: cp, ci
    use TLAB_CONSTANTS, only: efile
    use TLAB_VARS, only: rbg, tbg
    use THERMO_VARS, only: imixture
    use TLAB_VARS, only: rtime ! rtime is overwritten in io_read_fields
    use TLAB_PROCS
    use FLOW_LOCAL

    implicit none

    integer(ci) code

    real(cp), dimension(imax, jmax, kmax) :: T, h, rho, p
    real(cp), dimension(imax, jmax, kmax, *) :: s
    real(cp), dimension(imax, kmax) :: disp

    ! -------------------------------------------------------------------
    real(cp) dummy, ycenter
    real(cp) AVG1V2D, PROFILES
    real(cp) xcenter, amplify

    real(cp), dimension(:), pointer :: x, y, z

    ! ###################################################################
    ! Define pointers
    x => g(1)%nodes
    y => g(2)%nodes
    z => g(3)%nodes

#ifdef USE_MPI
    idsp = ims_offset_i; kdsp = ims_offset_k
#else
    idsp = 0; kdsp = 0
#endif

    ! ###################################################################
    ! Center plane displacement
    ! ###################################################################
    disp = C_0_R

    ! -------------------------------------------------------------------
    ! Broadband case
    ! -------------------------------------------------------------------
    if (code == 4) then
        dummy = rtime   ! rtime is overwritten in io_read_fields
        call IO_READ_FIELDS('scal.rand', IO_SCAL, imax, 1, kmax, 1, 1, disp, s) ! using array s as aux array
        rtime = dummy
        dummy = AVG1V2D(imax, 1, kmax, 1, 1, disp)     ! remove mean
        disp = disp - dummy

        ! -------------------------------------------------------------------
        ! Discrete case
        ! -------------------------------------------------------------------
    else if (code == 5) then
        wx_1 = C_2_R*C_PI_R/g(1)%scale ! Fundamental wavelengths
        wz_1 = C_2_R*C_PI_R/g(3)%scale

        do im = 1, fp%size
            wx = M_REAL(fp%modex(im))*wx_1
            wz = M_REAL(fp%modez(im))*wz_1

            do k = 1, kmax
                disp(:, k) = disp(:, k) + fp%amplitude(im)*COS(wx*x(idsp + 1:idsp + imax) + fp%phasex(im)) &
                             *COS(wz*z(kdsp + k) + fp%phasez(im))
            end do

        end do

    end if

    ! -------------------------------------------------------------------
    ! Modulation
    ! -------------------------------------------------------------------
    if (fp%type == PROFILE_GAUSSIAN .and. fp%parameters(1) > C_0_R) then
        do k = 1, kmax
            do i = 1, imax
                xcenter = x(i) - g(1)%scale*C_05_R - x(1)
                amplify = EXP(-C_05_R*(xcenter/fp%parameters(1))**2)
                disp(i, k) = disp(i, k)*amplify
            end do
        end do
    end if

    ! ###################################################################
    ! Perturbation in the thermodynamic fields
    ! ###################################################################
    if (rbg%type == PROFILE_NONE) then
        prof_loc = tbg

        if (tbg%type > 0) then ! temperature/mixture profile is given
            do k = 1, kmax
                do i = 1, imax
                    ycenter = y(1) + g(2)%scale*tbg%ymean + disp(i, k)
                    prof_loc%delta = tbg%delta + (tbg%uslope - tbg%lslope)*disp(i, k)*g(2)%scale
                    prof_loc%mean = tbg%mean + C_05_R*(tbg%uslope + tbg%lslope)*disp(i, k)*g(2)%scale
                    do j = 1, jmax
                        T(i, j, k) = PROFILES(prof_loc, ycenter, y(j))
                    end do
                end do
            end do

            if (imixture == MIXT_TYPE_AIRWATER) then
                call THERMO_AIRWATER_PT(imax, jmax, kmax, s, p, T)
            end if

        else if (tbg%type < 0) then ! enthalpy/mixture profile is given
            prof_loc%type = -tbg%type

            do k = 1, kmax
                do i = 1, imax
                    prof_loc%delta = tbg%delta + (tbg%uslope - tbg%lslope)*disp(i, k)*g(2)%scale
                    prof_loc%mean = tbg%mean + C_05_R*(tbg%uslope + tbg%lslope)*disp(i, k)*g(2)%scale
                    ycenter = y(1) + g(2)%scale*tbg%ymean + disp(i, k)
                    do j = 1, jmax
                        h(i, j, k) = PROFILES(prof_loc, ycenter, y(j))
                    end do
                end do
            end do

            if (imixture == MIXT_TYPE_AIRWATER) then
                call THERMO_AIRWATER_PH_RE(imax, jmax, kmax, s, p, h, T)
            end if

        end if

        ! compute perturbation in density
        call THERMO_THERMAL_DENSITY(imax, jmax, kmax, s, p, T, rho)

    else ! Defined in terms of the density, to be developed

    end if

    return
end subroutine DENSITY_FLUCTUATION
