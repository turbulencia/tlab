#include "dns_const.h"
#include "dns_error.h"

module Tlab_Background
    use TLab_Constants, only: wp, wi, efile, lfile, wfile, MAX_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Profiles, only: profiles_dt, Profiles_Calculate
    use THERMO_THERMAL
    implicit none
    private

! ###################################################################
! information to set up bcs, ics, and reference background profiles
! ###################################################################
    type(profiles_dt), public :: qbg(3)                     ! Velocity
    type(profiles_dt), public :: sbg(MAX_VARS)              ! Scalars
    type(profiles_dt), public :: pbg, rbg, tbg, hbg         ! Pressure, density, temperature, enthalpy

    public :: TLab_Initialize_Background
    public :: FLOW_SPATIAL_DENSITY, FLOW_SPATIAL_VELOCITY, FLOW_SPATIAL_SCALAR

    ! background, reference profiles
    real(wp), allocatable :: sbackground(:, :)              ! Scalars

contains
!########################################################################
!# Initialize data of reference profiles
!########################################################################
    subroutine TLab_Initialize_Background(inifile)
        use TLab_Pointers_3D, only: p_wrk1d
        use FDM, only: g
        use TLAB_VARS, only: imode_sim
        use TLAB_VARS, only: inb_scal, inb_scal_array
        use NavierStokes, only: froude, schmidt
        use Thermodynamics
        use THERMO_ANELASTIC
        use Profiles, only: Profiles_ReadBlock
        use Profiles, only: PROFILE_NONE, PROFILE_EKMAN_U, PROFILE_EKMAN_U_P, PROFILE_EKMAN_V
        use Gravity, only: buoyancy, bbackground, Gravity_Buoyancy, Gravity_Hydrostatic_Enthalpy

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        character*64 lstr

        real(wp) ploc(2), rloc(2), Tloc(2), sloc(2, 1:MAX_VARS)

        integer(wi) is, j, idummy

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        ! -----------------------------------------------------------------------
        block = 'Scalar'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        do is = 1, MAX_VARS
            write (lstr, *) is
            call Profiles_ReadBlock(bakfile, inifile, block, 'Scalar'//trim(adjustl(lstr)), sbg(is))
        end do

        if (imode_sim == DNS_MODE_SPATIAL) then
            call TLab_Write_ASCII(bakfile, '#ThickA=<value>')
            call TLab_Write_ASCII(bakfile, '#ThickB=<value>')
            call TLab_Write_ASCII(bakfile, '#Flux=<value>')

            do is = 1, MAX_VARS
                write (lstr, *) is
                call ScanFile_Real(bakfile, inifile, block, 'ThickA'//trim(adjustl(lstr)), '0.14', sbg(is)%parameters(2))
                call ScanFile_Real(bakfile, inifile, block, 'ThickB'//trim(adjustl(lstr)), '2.0', sbg(is)%parameters(3))
                call ScanFile_Real(bakfile, inifile, block, 'Flux'//trim(adjustl(lstr)), '0.94', sbg(is)%parameters(4))
            end do
        end if

        ! -----------------------------------------------------------------------
        block = 'Flow'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityX', qbg(1))
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityY', qbg(2))
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityZ', qbg(3))

        ! backwards compatilibity; originally, all velocity data was contained in block 'Velocity' except for the mean value
        call ScanFile_Char(bakfile, inifile, 'Flow', 'ProfileVelocity', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Velocity', qbg(1))
            call TLab_Write_ASCII(wfile, 'Update tag Flow.Velocity to Flow.VelocityX.')
        end if

        ! Consistency check
        if (any([PROFILE_EKMAN_U, PROFILE_EKMAN_U_P] == qbg(1)%type)) then
            qbg(3) = qbg(1)
            qbg(3)%type = PROFILE_EKMAN_V
        end if

        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Pressure', pbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Density', rbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Temperature', tbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Enthalpy', hbg)

        ! consistency check; two and only two are givem TO BE CHECKED BECAUSE PROFILE_NONE is used as constant profile
        if (imode_thermo == THERMO_TYPE_COMPRESSIBLE) then
            idummy = 0
            if (pbg%type == PROFILE_NONE) idummy = idummy + 1
            if (rbg%type == PROFILE_NONE) idummy = idummy + 1
            if (tbg%type == PROFILE_NONE) idummy = idummy + 1
            if (hbg%type == PROFILE_NONE) idummy = idummy + 1
            if (idummy /= 2) then
                call TLab_Write_ASCII(efile, __FILE__//'. Specify only 2 thermodynamic profiles.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        ! -----------------------------------------------------------------------
        if (imode_sim == DNS_MODE_SPATIAL) then     ! Thickness evolutions delta_i/diam_i=a*(x/diam_i+b)
            call TLab_Write_ASCII(bakfile, '#ThickAVelocity=<value>')
            call TLab_Write_ASCII(bakfile, '#ThickBVelocity=<value>')
            call TLab_Write_ASCII(bakfile, '#FluxVelocity=<value>')
            call TLab_Write_ASCII(bakfile, '#ThickADensity=<value>')
            call TLab_Write_ASCII(bakfile, '#ThickBDensity=<value>')
            call TLab_Write_ASCII(bakfile, '#FluxDensity=<value>')
            call TLab_Write_ASCII(bakfile, '#ThickATemperature=<value>')
            call TLab_Write_ASCII(bakfile, '#ThickBTemperature=<value>')
            call TLab_Write_ASCII(bakfile, '#FluxTemperature=<value>')

            ! Bradbury profile is the default (x0=a*b)
            call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickAVelocity', '0.1235', qbg(1)%parameters(2))
            call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickBVelocity', '-0.873', qbg(1)%parameters(3))
            call ScanFile_Real(bakfile, inifile, 'Flow', 'FluxVelocity', '0.96', qbg(1)%parameters(4))

            ! Ramaprian is the default (x0=a*b)
            call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickADensity', '0.14', rbg%parameters(2))
            call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickBDensity', '2.0', rbg%parameters(3))
            call ScanFile_Real(bakfile, inifile, 'Flow', 'FluxDensity', '0.94', rbg%parameters(4))

            call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickATemperature', '0.14', tbg%parameters(2))
            call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickBTemperature', '2.0', tbg%parameters(3))
            call ScanFile_Real(bakfile, inifile, 'Flow', 'FluxTemperature', '0.94', tbg%parameters(4))

        end if

        ! #######################################################################
        ! mean_rho and delta_rho need to be defined, because of old version.
        if (imode_thermo == THERMO_TYPE_COMPRESSIBLE) then
            if (rbg%type == PROFILE_NONE .and. tbg%type /= PROFILE_NONE) then
                rbg = tbg

                Tloc(1) = tbg%mean + 0.5_wp*tbg%delta
                Tloc(2) = tbg%mean - 0.5_wp*tbg%delta
                ploc(1:2) = pbg%mean
                sloc(1, :) = sbg(:)%mean + 0.5_wp*sbg(:)%delta
                sloc(2, :) = sbg(:)%mean - 0.5_wp*sbg(:)%delta
                call THERMO_THERMAL_DENSITY(2, sloc, ploc, Tloc, rloc)
                rbg%mean = 0.5_wp*(rloc(1) + rloc(2))
                rbg%delta = rloc(1) - rloc(2)

            else if (rbg%type == PROFILE_NONE .and. hbg%type /= PROFILE_NONE) then
                rbg%mean = 1.0_wp ! to be done; I need a nonzero value for dns_control.

            else
                tbg = rbg

                rloc(1) = rbg%mean + 0.5_wp*rbg%delta
                rloc(2) = rbg%mean - 0.5_wp*rbg%delta
                ploc(1:2) = pbg%mean
                sloc(1, :) = sbg(:)%mean + 0.5_wp*sbg(:)%delta
                sloc(2, :) = sbg(:)%mean - 0.5_wp*sbg(:)%delta
                call THERMO_THERMAL_TEMPERATURE(2, sloc, ploc, rloc, Tloc)
                tbg%mean = 0.5_wp*(Tloc(1) + Tloc(2))
                tbg%delta = Tloc(1) - Tloc(2)

            end if
        end if

        do is = 1, size(qbg)
            if (qbg(is)%relative) qbg(is)%ymean = g(2)%nodes(1) + g(2)%scale*qbg(is)%ymean_rel
        end do
        if (pbg%relative) pbg%ymean = g(2)%nodes(1) + g(2)%scale*pbg%ymean_rel
        if (rbg%relative) rbg%ymean = g(2)%nodes(1) + g(2)%scale*rbg%ymean_rel
        if (tbg%relative) tbg%ymean = g(2)%nodes(1) + g(2)%scale*tbg%ymean_rel
        if (hbg%relative) hbg%ymean = g(2)%nodes(1) + g(2)%scale*hbg%ymean_rel
        do is = 1, size(sbg)
            if (sbg(is)%relative) sbg(is)%ymean = g(2)%nodes(1) + g(2)%scale*sbg(is)%ymean_rel
        end do

        ! #######################################################################
        ! -----------------------------------------------------------------------
        ! Construct reference  profiles
        allocate (sbackground(g(2)%size, inb_scal_array))   ! scalar profiles

        do is = 1, inb_scal
            do j = 1, g(2)%size
                sbackground(j, is) = Profiles_Calculate(sbg(is), g(2)%nodes(j))
            end do
        end do

        if (imode_thermo == THERMO_TYPE_ANELASTIC) then     ! thermodynamic profiles
            allocate (epbackground(g(2)%size))
            allocate (tbackground(g(2)%size))
            allocate (pbackground(g(2)%size))
            allocate (rbackground(g(2)%size))
            allocate (ribackground(g(2)%size))

            call Gravity_Hydrostatic_Enthalpy(g(2), sbackground, epbackground, tbackground, pbackground, pbg%ymean, pbg%mean, p_wrk1d(:, 1))

            call THERMO_ANELASTIC_DENSITY(1, g(2)%size, 1, sbackground, rbackground)
            ribackground = 1.0_wp/rbackground

        end if

        if (buoyancy%type /= EQNS_NONE) then
            allocate (bbackground(g(2)%size))                   ! buoyancy profiles

            bbackground(:) = 0.0_wp
            if (buoyancy%active(2)) then
                call Gravity_Buoyancy(buoyancy, 1, g(2)%size, 1, sbackground(:, 1), p_wrk1d, bbackground)
                bbackground(:) = p_wrk1d(:, 1)
            end if

        end if

        ! -----------------------------------------------------------------------
        ! Add diagnostic fields to reference profile data, if any
        do is = inb_scal + 1, inb_scal_array ! Add diagnostic fields, if any
            sbg(is) = sbg(1)
            schmidt(is) = schmidt(1)
        end do

        ! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
        if (buoyancy%type /= EQNS_NONE) then
            sbg(is) = sbg(1)
            sbg(is)%mean = (bbackground(1) + bbackground(g(2)%size))/froude
            sbg(is)%delta = abs(bbackground(1) - bbackground(g(2)%size))/froude
            schmidt(is) = schmidt(1)
        end if

        ! theta_l as next scalar
        if (imode_thermo == THERMO_TYPE_ANELASTIC) then
            if (imixture == MIXT_TYPE_AIRWATER) then
                is = is + 1
                call THERMO_ANELASTIC_THETA_L(1, g(2)%size, 1, sbackground, p_wrk1d)
                sbg(is) = sbg(1)
                sbg(is)%mean = (p_wrk1d(1, 1) + p_wrk1d(g(2)%size, 1))*0.5_wp
                sbg(is)%delta = abs(p_wrk1d(1, 1) - p_wrk1d(g(2)%size, 1))
                schmidt(is) = schmidt(1)
            end if
        end if

        ! end if

        return
    end subroutine TLab_Initialize_Background

!########################################################################
!#
!# Calculate the fields rho(x,y), u(x,y) and v(x,y) s.t. the axial momentum flux
!# is conserved and the continuity equation is satisfied.
!#
!# It assumes passive incompressible mixing of the temperature.
!#
!# The inputs are rho_vi, u_vi and the output of the iteration is rho_vo.
!# Array u_vo and tem_vo are auxiliar.
!#
!########################################################################
    subroutine FLOW_SPATIAL_DENSITY(imax, jmax, tbg, ubg, &
                                    x, y, z1, p, rho_vi, u_vi, tem_vi, rho_vo, u_vo, tem_vo, wrk1d)
        integer(wi) imax, jmax
        type(profiles_dt) tbg, ubg
        real(wp) x(imax)
        real(wp) y(jmax)

        real(wp) z1(jmax), p(jmax)
        real(wp) rho_vi(jmax), u_vi(jmax), tem_vi(jmax)
        real(wp) rho_vo(imax, jmax), u_vo(jmax)
        real(wp) tem_vo(jmax)
        real(wp) wrk1d(jmax, *)

! -------------------------------------------------------------------
        real(wp) tol, err

        integer(wi) i, j, n, nmax, ier
        integer, parameter :: i1 = 1

! ###################################################################

! Convergence parameters
        nmax = 30
        tol = 1.0e-6_wp

        tem_vi(1:jmax) = 0.0_wp
        do j = 1, jmax
            tem_vi(j) = Profiles_Calculate(tbg, y(j))
        end do

#define rho_aux(j) wrk1d(j,1)
#define aux(j)     wrk1d(j,2)
! ###################################################################
! Begining loop over x-planes
! ###################################################################
! The initial condition for rho_vo is rho_vi
        do j = 1, jmax
            rho_aux(j) = rho_vi(j)
        end do
        do i = 1, imax
            ier = 1

! Begining iteration for a fixed x-plane
            do n = 1, nmax
! Velocity profile. Array tem_vo used as auxiliar array
                call FLOW_SPATIAL_VELOCITY(i1, jmax, ubg, ubg%diam, &
                                           ubg%parameters(2), ubg%parameters(3), ubg%parameters(4), x(i), y, &
                                           rho_vi, u_vi, rho_aux(1), u_vo, tem_vo, aux(1), wrk1d(1, 3))
! Normalized temperature and density profiles
                call FLOW_SPATIAL_SCALAR(i1, jmax, tbg, tbg%diam, ubg%diam, &
                                         tbg%parameters(2), tbg%parameters(3), tbg%parameters(4), x(i), y, &
                                         rho_vi, u_vi, tem_vi, rho_aux(1), u_vo, tem_vo, aux(1))
                call THERMO_THERMAL_DENSITY(jmax, z1, p, tem_vo, wrk1d(1, 2))
! Convergence criteria (infinity norm)
                do j = 1, jmax
                    wrk1d(j, 3) = abs(wrk1d(j, 2) - rho_aux(j))
                    rho_aux(j) = wrk1d(j, 2)
                end do
                err = maxval(wrk1d(1:jmax, 3), jmax)/maxval(wrk1d(1:jmax, 2))
                if (err < tol) then
                    ier = 0
                    exit
                end if
            end do

! Final check
            if (ier == 1) then
                call TLab_Write_ASCII(wfile, 'FLOW_SPATIAL: nmax reached.')
            end if
            do j = 1, jmax
                rho_vo(i, j) = rho_aux(j)
            end do
        end do

        return
    end subroutine FLOW_SPATIAL_DENSITY

!########################################################################
!# DESCRIPTION
!#
!# Calculate the fields u(x,y) and v(x,y) s.t. the axial momentum flux
!# is conserved and the continuity equation is satisfied.
!#
!# In this subroutine, the density field is a given input.
!#
!# The calculation of v assumes ymean_rel_u equal to 0.5. This section
!# should also be rewritten in terms of OPR_PARTIAL_ and QUAD routines.
!#
!########################################################################
    subroutine FLOW_SPATIAL_VELOCITY(imax, jmax, prof_loc, diam_u, &
                                     jet_u_a, jet_u_b, jet_u_flux, x, y, rho_vi, u_vi, rho, u, v, wrk1d, wrk2d)
        use Integration, only: Int_Simpson

        integer(wi) imax, jmax
        type(profiles_dt) prof_loc
        real(wp) diam_u
        real(wp) jet_u_a, jet_u_b, jet_u_flux
        real(wp) x(imax)
        real(wp) y(jmax)

        real(wp) rho_vi(jmax), u_vi(jmax)
        real(wp) rho(imax, jmax), u(imax, jmax), v(imax, jmax)
        real(wp) wrk1d(jmax, *), wrk2d(imax, jmax)

! -------------------------------------------------------------------
        real(wp) c1, c2
        real(wp) delta, eta, ExcMom_vi, Q1, Q2, U2, UC
        real(wp) dummy, flux_aux, diam_loc
        real(wp) xi_tr, dxi_tr

        integer(wi) i, j, jsym

! ###################################################################

! Bradbury profile
#ifdef SINGLE_PREC
        c1 = -0.6749e+0
        c2 = 0.027e+0
#else
        c1 = -0.6749d+0
        c2 = 0.027d+0
#endif

        U2 = prof_loc%mean - 0.5_wp*prof_loc%delta

! ###################################################################
! Axial velocity, U_c*f(eta)
! ###################################################################
! -------------------------------------------------------------------
! Transition as a tanh profile around xi_tr between (0,2xi_tr)
! in the slope of the half width.
! It is set s.t. the vorticity thickness associated with this
! tanh is half of the distance between xi_tr and the estimated end of
! transition, 2 xi_tr.
! -------------------------------------------------------------------
        xi_tr = 0.5_wp/jet_u_a - jet_u_b
        if (xi_tr < 0.0_wp) then
            call TLab_Write_ASCII(wfile, 'FLOW_SPATIAL_VELOCITY. xi_tr negative.')
        end if
        dxi_tr = xi_tr/8.0_wp

! -------------------------------------------------------------------
! Normalized velocity, f(eta)
! -------------------------------------------------------------------
        do i = 1, imax
            delta = (dxi_tr*log(exp((x(i)/diam_u - xi_tr)/dxi_tr) + 1.0_wp)*jet_u_a + 0.5_wp)*diam_u
            diam_loc = 2.0_wp*delta

! inflow profile
            prof_loc%parameters(5) = diam_loc
            wrk1d(1:jmax, 1) = 0.0_wp
            do j = 1, jmax
                wrk1d(j, 1) = Profiles_Calculate(prof_loc, y(j))
            end do
            UC = wrk1d(jmax/2, 1) - U2

! U-U2=f(y) reference profile, stored in array u.
            do j = 1, jmax
                eta = (y(j) - prof_loc%ymean)/delta
                u(i, j) = exp(c1*eta**2*(1.0_wp + c2*eta**4))
                dummy = 0.5_wp*(1.0_wp + tanh(0.5_wp*(x(i)/diam_u - xi_tr)/dxi_tr))
                u(i, j) = dummy*u(i, j)*UC + (1.0_wp - dummy)*(wrk1d(j, 1) - U2)
            end do

        end do

! -------------------------------------------------------------------
! Magnitude UC for conserving the axial flux (w/ correction jet_u_flux)
! -------------------------------------------------------------------
! Reference momentum excess at the inflow
        do j = 1, jmax
            wrk1d(j, 1) = rho_vi(j)*u_vi(j)*(u_vi(j) - U2)
        end do
        ExcMom_vi = Int_Simpson(wrk1d(1:jmax, 1), y(1:jmax))

        do i = 1, imax
! Correction factor varying between 1 at the inflow and jet_u_flux
! at the outflow
            dummy = 0.5_wp*(1.0_wp + tanh(0.5_wp*(x(i)/diam_u - xi_tr)/dxi_tr))
            flux_aux = dummy*jet_u_flux + (1.0_wp - dummy)*1.0_wp

! Calculating UC.
! Solve second order equation Q1*UC^2+Q2*UC-J=0, positive root.
            do j = 1, jmax
                wrk1d(j, 1) = rho(i, j)*u(i, j)*u(i, j)
                wrk1d(j, 2) = rho(i, j)*u(i, j)
            end do
            Q1 = Int_Simpson(wrk1d(1:jmax, 1), y(1:jmax))
            Q2 = U2*Int_Simpson(wrk1d(1:jmax, 2), y(1:jmax))
            UC = (-Q2 + sqrt(Q2*Q2 + 4.0_wp*Q1*ExcMom_vi*flux_aux))/2.0_wp/Q1

! Scaled velocity
            do j = 1, jmax
                u(i, j) = U2 + UC*u(i, j)
            end do
        end do

! ###################################################################
! Lateral velocity, d(rho*v)/dy=-d(rho*u)/dx
! ###################################################################
        if (imax > 1) then
! Backwards 1st-order derivative, -d(rho*u)/dx (w used as aux array)
            do i = 2, imax
            do j = 1, jmax
                Q1 = rho(i, j)*u(i, j)
                Q2 = rho(i - 1, j)*u(i - 1, j)
                wrk2d(i, j) = -(Q1 - Q2)/(x(i) - x(i - 1))
            end do
            end do
! Backwards 1st-order derivative for i=1
            i = 1
            do j = 1, jmax
                wrk2d(i, j) = wrk2d(i + 1, j)
            end do

! Midpoint integration, rho*v
            do i = 1, imax
                Q1 = wrk2d(i, jmax/2) + wrk2d(i, jmax/2 + 1)
                Q2 = y(jmax/2 + 1) - y(jmax/2)
                v(i, jmax/2 + 1) = 0.5_wp*(0.5_wp*Q1*Q2)

                do j = jmax/2 + 2, jmax
                    Q1 = wrk2d(i, j) + wrk2d(i, j - 1)
                    Q2 = y(j) - y(j - 1)
                    v(i, j) = v(i, j - 1) + 0.5_wp*Q1*Q2
                end do
            end do

! Division by rho and antisymmetric extension
            do i = 1, imax
            do j = jmax/2 + 1, jmax
                v(i, j) = v(i, j)/rho(i, j)
                jsym = jmax - j + 1
                v(i, jsym) = -v(i, j)
            end do
            end do

! If only one plane, set lateral velocity to zero
        else
            do j = 1, jmax
                v(1, j) = 0.0_wp
            end do

        end if

        return
    end subroutine FLOW_SPATIAL_VELOCITY

!########################################################################
!#
!# Calculate the fields z1(x,y) s.t. the axial scalar flux
!# is conserved.
!#
!# In this subroutine, the density and velocity field are a given input.
!#
!########################################################################
    subroutine FLOW_SPATIAL_SCALAR(imax, jmax, prof_loc, &
                                   diam_z, diam_u, jet_z_a, jet_z_b, jet_z_flux, &
                                   x, y, rho_vi, u_vi, z_vi, rho, u, z1, wrk1d)
        use Integration, only: Int_Simpson

        integer(wi) imax, jmax
        real(wp) diam_u, diam_z
        real(wp) jet_z_a, jet_z_b, jet_z_flux
        real(wp) x(imax)
        real(wp) y(jmax)

        real(wp) rho_vi(jmax), u_vi(jmax), z_vi(jmax)
        real(wp) rho(imax, jmax), u(imax, jmax), z1(imax, jmax)
        real(wp) wrk1d(jmax, *)

! -------------------------------------------------------------------
        real(wp) c1, c2
        real(wp) delta, eta, ExcMom_vi, Q1, Z2, ZC, flux_aux
        real(wp) dummy, diam_loc
        real(wp) xi_tr, dxi_tr
        integer(wi) i, j
        type(profiles_dt) prof_loc

! ###################################################################
!   param = 0.0_wp

! Ramaprian85 for the scalar ?
        c1 = -0.6749_wp
        c2 = 0.027_wp

        Z2 = prof_loc%mean - 0.5_wp*prof_loc%delta

! -------------------------------------------------------------------
! Transition as a tanh profile around xi_tr between (0,2xi_tr)
! in the slope of the half width.
! It is set s.t. the vorticity thickness associated with this
! tanh is half of the distance between xi_tr and the estimated end of
! transition, 2 xi_tr.
! -------------------------------------------------------------------
        xi_tr = 0.5_wp*diam_z/diam_u/jet_z_a - jet_z_b
        if (xi_tr < 0.0_wp) then
            call TLab_Write_ASCII(wfile, 'FLOW_SPATIAL_VELOCITY. xi_tr negative.')
        end if
        dxi_tr = xi_tr/8.0_wp

! -------------------------------------------------------------------
! Normalized scalar, f(eta)
! -------------------------------------------------------------------
        do i = 1, imax
            delta = (dxi_tr*log(exp((x(i)/diam_u - xi_tr)/dxi_tr) + 1.0_wp)*jet_z_a + 0.5_wp)*diam_u
            diam_loc = 2.0_wp*delta

! inflow profile
            prof_loc%parameters(5) = diam_loc
            wrk1d(1:jmax, 1) = 0.0_wp
            do j = 1, jmax
                wrk1d(j, 1) = Profiles_Calculate(prof_loc, y(j))
            end do
            ZC = wrk1d(jmax/2, 1) - Z2

! Z-Z2=f(y) reference profile, stored in array z1.
            do j = 1, jmax
                eta = (y(j) - prof_loc%ymean)/delta
                z1(i, j) = exp(c1*eta**2*(1.0_wp + c2*eta**4))
                dummy = 0.5_wp*(1.0_wp + tanh(0.5_wp*(x(i)/diam_u - xi_tr)/dxi_tr))
                z1(i, j) = dummy*z1(i, j)*ZC + (1.0_wp - dummy)*(wrk1d(j, 1) - Z2)
            end do
        end do

! -------------------------------------------------------------------
! Magnitude ZC for conserving the axial flux (w/ correction jet_u_flux)
! -------------------------------------------------------------------
! Reference momentum excess at the inflow
        do j = 1, jmax
            wrk1d(j, 2) = rho_vi(j)*u_vi(j)*(z_vi(j) - Z2)
        end do
        ExcMom_vi = Int_Simpson(wrk1d(1:jmax, 2), y(1:jmax))

        do i = 1, imax
! Correction factor varying between 1 at the inflow and jet_z_flux
! at the outflow
            dummy = 0.5_wp*(1.0_wp + tanh(0.5_wp*(x(i)/diam_u - xi_tr)/dxi_tr))
            flux_aux = dummy*jet_z_flux + (1.0_wp - dummy)*1.0_wp

! Calculating ZC.
! Solve second order equation Q1*ZC-J=0, positive root.
            do j = 1, jmax
                wrk1d(j, 1) = rho(i, j)*u(i, j)*z1(i, j)
            end do
            Q1 = Int_Simpson(wrk1d(1:jmax, 1), y(1:jmax))
            ZC = flux_aux*ExcMom_vi/Q1
            do j = 1, jmax
                z1(i, j) = Z2 + ZC*z1(i, j)
            end do
        end do

        return
    end subroutine FLOW_SPATIAL_SCALAR

end module Tlab_Background
