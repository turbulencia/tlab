#include "types.h"

!########################################################################
!# DESCRIPTION
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
    use TLAB_TYPES, only: profiles_dt
    use TLAB_CONSTANTS, only: wfile
    use TLAB_PROCS

    implicit none

#include "integers.h"

    TINTEGER imax, jmax
    type(profiles_dt) tbg, ubg
    TREAL x(imax)
    TREAL y(jmax)

    TREAL z1(*), p(*)
    TREAL rho_vi(jmax), u_vi(jmax), tem_vi(jmax)
    TREAL rho_vo(imax, jmax), u_vo(jmax)
    TREAL tem_vo(jmax)
    TREAL wrk1d(jmax, *)

! -------------------------------------------------------------------
    TREAL PROFILES
    TREAL tol, err

    TINTEGER i, j, n, nmax, ier

! ###################################################################

! Convergence parameters
    nmax = 30
    tol = C_1EM6_R

    tem_vi(1:jmax) = C_0_R
    do j = 1, jmax
        tem_vi(j) = PROFILES(tbg, y(j))
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
            call THERMO_THERMAL_DENSITY(i1, jmax, i1, z1, p, tem_vo, wrk1d(1, 2))
! Convergence criteria (infinity norm)
            do j = 1, jmax
                wrk1d(j, 3) = ABS(wrk1d(j, 2) - rho_aux(j))
                rho_aux(j) = wrk1d(j, 2)
            end do
            err = MAXVAL(wrk1d(1:jmax, 3), jmax)/MAXVAL(wrk1d(1:jmax, 2))
            if (err < tol) then
                ier = 0
                exit
            end if
        end do

! Final check
        if (ier == 1) then
            call TLAB_WRITE_ASCII(wfile, 'FLOW_SPATIAL: nmax reached.')
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
!# The calculation of v assumes ycoor_u equal to 0.5. This section
!# should also be rewritten in terms of OPR_PARTIAL_ and QUAD routines.
!#
!########################################################################
subroutine FLOW_SPATIAL_VELOCITY(imax, jmax, prof_loc, diam_u, &
     jet_u_a, jet_u_b, jet_u_flux, x, y, rho_vi, u_vi, rho, u, v, wrk1d, wrk2d)
    use TLAB_TYPES, only: profiles_dt
    use TLAB_CONSTANTS, only: efile, wfile
    use TLAB_VARS, only: g
    use TLAB_PROCS

    implicit none

#include "integers.h"

    TINTEGER imax, jmax
    type(profiles_dt) prof_loc
    TREAL diam_u
    TREAL jet_u_a, jet_u_b, jet_u_flux
    TREAL x(imax)
    TREAL y(jmax)

    TREAL rho_vi(jmax), u_vi(jmax)
    TREAL rho(imax, jmax), u(imax, jmax), v(imax, jmax)
    TREAL wrk1d(jmax, *), wrk2d(imax, jmax)

! -------------------------------------------------------------------
    TREAL c1, c2
    TREAL delta, eta, ExcMom_vi, Q1, Q2, U2, UC
    TREAL dummy, flux_aux, diam_loc
    TREAL SIMPSON_NU, PROFILES, ycenter
    TREAL xi_tr, dxi_tr

    TINTEGER i, j, jsym

! ###################################################################

! Bradbury profile
#ifdef SINGLE_PREC
    c1 = -0.6749e+0
    c2 = 0.027e+0
#else
    c1 = -0.6749d+0
    c2 = 0.027d+0
#endif

    U2 = prof_loc%mean - C_05_R*prof_loc%delta
    ycenter = y(1) + g(2)%scale*prof_loc%ymean_rel

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
    xi_tr = C_05_R/jet_u_a - jet_u_b
    if (xi_tr < C_0_R) then
        call TLAB_WRITE_ASCII(wfile, 'FLOW_SPATIAL_VELOCITY. xi_tr negative.')
    end if
    dxi_tr = xi_tr/C_8_R

! -------------------------------------------------------------------
! Normalized velocity, f(eta)
! -------------------------------------------------------------------
    do i = 1, imax
        delta = (dxi_tr*LOG(EXP((x(i)/diam_u - xi_tr)/dxi_tr) + C_1_R)*jet_u_a + C_05_R)*diam_u
        diam_loc = C_2_R*delta

! inflow profile
        prof_loc%parameters(5) = diam_loc
        wrk1d(1:jmax, 1) = C_0_R
        do j = 1, jmax
            wrk1d(j, 1) = PROFILES(prof_loc, y(j))
        end do
        UC = wrk1d(jmax/2, 1) - U2

! U-U2=f(y) reference profile, stored in array u.
        do j = 1, jmax
            eta = (y(j) - ycenter)/delta
            u(i, j) = EXP(c1*eta**2*(C_1_R + c2*eta**4))
            dummy = C_05_R*(C_1_R + TANH(C_05_R*(x(i)/diam_u - xi_tr)/dxi_tr))
            u(i, j) = dummy*u(i, j)*UC + (C_1_R - dummy)*(wrk1d(j, 1) - U2)
        end do

    end do

! -------------------------------------------------------------------
! Magnitude UC for conserving the axial flux (w/ correction jet_u_flux)
! -------------------------------------------------------------------
! Reference momentum excess at the inflow
    do j = 1, jmax
        wrk1d(j, 1) = rho_vi(j)*u_vi(j)*(u_vi(j) - U2)
    end do
    ExcMom_vi = SIMPSON_NU(jmax, wrk1d, y)

    do i = 1, imax
! Correction factor varying between 1 at the inflow and jet_u_flux
! at the outflow
        dummy = C_05_R*(C_1_R + TANH(C_05_R*(x(i)/diam_u - xi_tr)/dxi_tr))
        flux_aux = dummy*jet_u_flux + (C_1_R - dummy)*C_1_R

! Calculating UC.
! Solve second order equation Q1*UC^2+Q2*UC-J=0, positive root.
        do j = 1, jmax
            wrk1d(j, 1) = rho(i, j)*u(i, j)*u(i, j)
            wrk1d(j, 2) = rho(i, j)*u(i, j)
        end do
        Q1 = SIMPSON_NU(jmax, wrk1d(1, 1), y)
        Q2 = U2*SIMPSON_NU(jmax, wrk1d(1, 2), y)
        UC = (-Q2 + SQRT(Q2*Q2 + C_4_R*Q1*ExcMom_vi*flux_aux))/C_2_R/Q1

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
            v(i, jmax/2 + 1) = C_05_R*(C_05_R*Q1*Q2)

            do j = jmax/2 + 2, jmax
                Q1 = wrk2d(i, j) + wrk2d(i, j - 1)
                Q2 = y(j) - y(j - 1)
                v(i, j) = v(i, j - 1) + C_05_R*Q1*Q2
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
            v(1, j) = C_0_R
        end do

    end if

    return
end subroutine FLOW_SPATIAL_VELOCITY

!########################################################################
!# DESCRIPTION
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
    use TLAB_TYPES, only: profiles_dt
    use TLAB_CONSTANTS, only: wfile
    use TLAB_VARS, only: g
    use TLAB_PROCS

    implicit none

#include "integers.h"

    TINTEGER imax, jmax
    TREAL diam_u, diam_z
    TREAL jet_z_a, jet_z_b, jet_z_flux
    TREAL x(imax)
    TREAL y(jmax)

    TREAL rho_vi(jmax), u_vi(jmax), z_vi(jmax)
    TREAL rho(imax, jmax), u(imax, jmax), z1(imax, jmax)
    TREAL wrk1d(jmax, *)

! -------------------------------------------------------------------
    TREAL c1, c2
    TREAL delta, eta, ExcMom_vi, Q1, Z2, ZC, flux_aux
    TREAL dummy, diam_loc
    TREAL SIMPSON_NU, PROFILES, ycenter
    TREAL xi_tr, dxi_tr
    TINTEGER i, j
    type(profiles_dt) prof_loc

! ###################################################################
!   param = C_0_R

! Ramaprian85 for the scalar ?
#ifdef SINGLE_PREC
    c1 = -0.6749e+0
    c2 = 0.027e+0
#else
    c1 = -0.6749d+0
    c2 = 0.027d+0
#endif

    Z2 = prof_loc%mean - C_05_R*prof_loc%delta
    ycenter = y(1) + g(2)%scale*prof_loc%ymean_rel

! -------------------------------------------------------------------
! Transition as a tanh profile around xi_tr between (0,2xi_tr)
! in the slope of the half width.
! It is set s.t. the vorticity thickness associated with this
! tanh is half of the distance between xi_tr and the estimated end of
! transition, 2 xi_tr.
! -------------------------------------------------------------------
    xi_tr = C_05_R*diam_z/diam_u/jet_z_a - jet_z_b
    if (xi_tr < C_0_R) then
        call TLAB_WRITE_ASCII(wfile, 'FLOW_SPATIAL_VELOCITY. xi_tr negative.')
    end if
    dxi_tr = xi_tr/C_8_R

! -------------------------------------------------------------------
! Normalized scalar, f(eta)
! -------------------------------------------------------------------
    do i = 1, imax
        delta = (dxi_tr*LOG(EXP((x(i)/diam_u - xi_tr)/dxi_tr) + C_1_R)*jet_z_a + C_05_R)*diam_u
        diam_loc = C_2_R*delta

! inflow profile
        prof_loc%parameters(5) = diam_loc
        wrk1d(1:jmax, 1) = C_0_R
        do j = 1, jmax
            wrk1d(j, 1) = PROFILES(prof_loc, y(j))
        end do
        ZC = wrk1d(jmax/2, 1) - Z2

! Z-Z2=f(y) reference profile, stored in array z1.
        do j = 1, jmax
            eta = (y(j) - ycenter)/delta
            z1(i, j) = EXP(c1*eta**2*(C_1_R + c2*eta**4))
            dummy = C_05_R*(C_1_R + TANH(C_05_R*(x(i)/diam_u - xi_tr)/dxi_tr))
            z1(i, j) = dummy*z1(i, j)*ZC + (C_1_R - dummy)*(wrk1d(j, 1) - Z2)
        end do
    end do

! -------------------------------------------------------------------
! Magnitude ZC for conserving the axial flux (w/ correction jet_u_flux)
! -------------------------------------------------------------------
! Reference momentum excess at the inflow
    do j = 1, jmax
        wrk1d(j, 2) = rho_vi(j)*u_vi(j)*(z_vi(j) - Z2)
    end do
    ExcMom_vi = SIMPSON_NU(jmax, wrk1d(1, 2), y)

    do i = 1, imax
! Correction factor varying between 1 at the inflow and jet_z_flux
! at the outflow
        dummy = C_05_R*(C_1_R + TANH(C_05_R*(x(i)/diam_u - xi_tr)/dxi_tr))
        flux_aux = dummy*jet_z_flux + (C_1_R - dummy)*C_1_R

! Calculating ZC.
! Solve second order equation Q1*ZC-J=0, positive root.
        do j = 1, jmax
            wrk1d(j, 1) = rho(i, j)*u(i, j)*z1(i, j)
        end do
        Q1 = SIMPSON_NU(jmax, wrk1d(1, 1), y)
        ZC = flux_aux*ExcMom_vi/Q1
        do j = 1, jmax
            z1(i, j) = Z2 + ZC*z1(i, j)
        end do
    end do

    return
end subroutine FLOW_SPATIAL_SCALAR
