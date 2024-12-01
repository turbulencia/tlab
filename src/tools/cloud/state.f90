#include "dns_const.h"

program STATE
    use TLab_Constants, only: wp, wi
    use TLAB_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use Thermodynamics
    use THERMO_THERMAL
    use THERMO_ANELASTIC
    use THERMO_CALORIC

    implicit none

    real(wp) p(1), ps(1), t(1), qs(1), qv(1), qt(1), ql(1), r(1), e(1), h(1), z1(2), dummy(1), dqldqt(1), ep(1), theta_v(1), theta_l(1), theta_e(1), Td(1)
    real(wp) heat1(1), heat2(1), cp1(1), cp2(1), alpha(1), as(1), bs(1)
    real(wp) r1(1), h1(1), s(3)
    integer(wi) iopt

! ###################################################################
    call TLab_Start()

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call Thermodynamics_Initialize_Parameters()
    ep = 0.0_wp
    dsmooth = 1.0_wp
    scaleheightinv = 1.0_wp

    write (*, *) 'Case p-t (1) or d-e (2) or p-h (3)?'
    read (*, *) iopt

    if (iopt == 1) then
        write (*, *) 'temperature (C) ?'
        read (*, *) t
        t = t + 273.15
        write (*, *) 'pressure (hPa) ?'
        read (*, *) p
        p = p *100.0_wp

    else if (iopt == 2) then
        write (*, *) 'density ?'
        read (*, *) r
        write (*, *) 'energy ?'
        read (*, *) e

    else if (iopt == 3) then
        write (*, *) 'enthalpy (kJ/kg)?'
        read (*, *) h
        write (*, *) 'pressure (hPa) ?'
        read (*, *) p
        p = p *100.0_wp

    end if

    allocate (pbackground(1), epbackground(1), rbackground(1))
    epbackground(1) = ep(1)
    pbackground(1) = p(1)
    rbackground(1) = r(1)

    write (*, *) 'water specific humidity (g/kg) ?'
    read (*, *) qt
    qt = qt*0.001_wp

! ###################################################################
    if (iopt == 1) then
        call Thermo_Psat_Polynomial(1, t, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        if (qt(1) > qs(1)) then
            qv = qs*(1 - qt)
            ql = qt - qv
        else
            qv = qt
            ql = 0.0_wp
        end if
        z1(1) = qt(1)
        z1(2) = ql(1)
        call THERMO_CALORIC_ENTHALPY(1, z1, t, h)
        call THERMO_CALORIC_ENERGY(1, z1, t, e)
        call THERMO_THERMAL_DENSITY(1, z1, p, t, r)

        s(1) = h(1); s(2:3) = z1(1:2)
        call THERMO_ANELASTIC_THETA_V(1, 1, 1, s, theta_v)
        call THERMO_ANELASTIC_THETA_L(1, 1, 1, s, theta_l)
        call THERMO_ANELASTIC_THETA_E(1, 1, 1, s, theta_e)
        call THERMO_ANELASTIC_DEWPOINT(1, 1, 1, s, Td, dummy)

    else if (iopt == 2) then
        z1(1) = qt(1)
        call THERMO_CALORIC_TEMPERATURE(1, z1, e, r, T, dqldqt)
        ql = z1(2)
        qv = qt - ql
        qs = qv ! initial condition for next routine
        call THERMO_THERMAL_PRESSURE(1, z1, r, t, p)
        call Thermo_Psat_Polynomial(1, t, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        call THERMO_CALORIC_ENTHALPY(1, z1, t, h)

        s(1) = h(1); s(2:3) = z1(1:2)
        call THERMO_ANELASTIC_THETA_V(1, 1, 1, s, theta_v)
        call THERMO_ANELASTIC_THETA_L(1, 1, 1, s, theta_l)
        call THERMO_ANELASTIC_THETA_E(1, 1, 1, s, theta_e)
        call THERMO_ANELASTIC_DEWPOINT(1, 1, 1, s, Td, dummy)

    else if (iopt == 3) then
        ! h = h/TREF/1.007
        z1(1) = qt(1)
        call THERMO_ANELASTIC_PH(1, 1, 1, z1, h)
        s(1) = h(1); s(2:3) = z1(1:2)
        call THERMO_ANELASTIC_TEMPERATURE(1, 1, 1, s, T)
        ! CALL THERMO_AIRWATER_PH_RE(1, z1, p, h, T)
        ql(1) = z1(2)
        qv = qt - ql

        call Thermo_Psat_Polynomial(1, T, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        call THERMO_THERMAL_DENSITY(1, z1, p, T, r)
        call THERMO_CALORIC_ENERGY(1, z1, T, e)

        call THERMO_ANELASTIC_THETA_V(1, 1, 1, s, theta_v)
        call THERMO_ANELASTIC_THETA_L(1, 1, 1, s, theta_l)
        call THERMO_ANELASTIC_THETA_E(1, 1, 1, s, theta_e)
        call THERMO_ANELASTIC_DEWPOINT(1, 1, 1, s, Td, dummy)

! check
        call THERMO_ANELASTIC_DENSITY(1, 1, 1, s, r1)
!     r2 = p/(T*(1- qt +qv/rd_ov_rv ) )
        call THERMO_CALORIC_ENTHALPY(1, z1, T, h1)

    end if

    write (*, 1000) 'Saturation specific humidity (g/kg):', qs*1.0e3_wp
    write (*, 1000) 'Vapor specific humidity (g/kg).....:', qv*1.0e3_wp
    write (*, 1000) 'Liquid specific humidity (g/kg)....:', ql*1.0e3_wp
    write (*, 1000) 'Density ...........................:', r
    write (*, 1000) 'Pressure (hPa) ....................:', p/100.0_wp
    write (*, 1000) 'Saturation pressure (hPa) .........:', ps/100.0_wp
    write (*, 1000) 'Temperature (K) ...................:', t !- 273.15
    write (*, 1000) 'Dewpoint temperature (K) ..........:', Td
    write (*, 1000) 'Specific heat capacity (J/kg/T)....:', Cd + qt*Cdv + ql*Cvl
    write (*, 1000) 'Specific energy  (J/kg)............:', e
    write (*, 1000) 'Specific enthalpy (J/kg)...........:', h
    write (*, 1000) 'Reference latent heat (J/kg) ......:', -THERMO_AI(6, 1, 3)
    write (*, 1000) 'Latent heat (J/kg) ................:', (-Cl - t*Lvl)
    write (*, 1000) 'Virtual potential T (K) ...........:', theta_v
    write (*, 1000) 'Liquid-water potential T (K) ......:', theta_l
    write (*, 1000) 'Equivalent potential T (K) ........:', theta_e
    if (iopt == 3) then
        write (*, 1000) 'Density ...........................:', r1
        write (*, 1000) 'Specific enthalpy .................:', h1
    end if

! ###################################################################
    write (*, *) ' '
    write (*, *) 'Calculate reversal linear coefficients (1-yes/0-no) ?'
    read (*, *) iopt

    if (iopt == 1 .and. ql(1) > 1.0_wp) then
        heat1 = -Lvl - Cvl*t
        heat2 = heat1*(1.0_wp + qv/(1.0_wp - qt)) - Cdv*t

        cp1 = (1.0_wp - qt)*Cd + qv*THERMO_AI(1, 1, 1) + ql*Cl
        dummy = (heat1**2)*qv/((t**2)*cp1*CRATIO_INV*Rv)
        cp2 = cp1*(1.0_wp + dummy*(1.0_wp + qv/(1.0_wp - qt)/rd_ov_rv))

        alpha = 1.0_wp + heat1*qv/((1.0_wp - qt)*CRATIO_INV*Rd*t)

        as = -alpha/cp2/t
        bs = heat2*as + 1.0_wp/(1.0_wp - qt)
        write (*, *) 'Enthalpy coefficient ..........:', as
        write (*, *) 'Water fraction coefficient ....:', bs

    else if (iopt == 1 .and. ql(1) == 1.0_wp) then
        cp1 = Cd + qt*Cdv

        as = -1.0_wp/cp1/t
        bs = Cdv/cp1 - Rdv/(Rd + qt*Rdv)
        write (*, *) 'Enthalpy coefficient ..........:', as
        write (*, *) 'Water fraction coefficient ....:', bs
    end if

    stop

1000 format(A, G_FORMAT_R)

end program STATE
