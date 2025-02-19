#include "dns_const.h"

program SMOOTH
    use TLab_Constants, only: wp, wi
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use Thermodynamics
    use THERMO_THERMAL
    use THERMO_CALORIC
    use THERMO_AIRWATER
    use Thermo_Anelastic

    implicit none

    real(wp) qt_min(1), qt_max(1), qt_del(1), qt(1), qs(1), dqldqt(1)
    real(wp) z1(2), e(1), rho(1), p(1), T(1), h(1), ep(1), s(3)
    integer(wi) opt

! ###################################################################
    call TLab_Start

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call Thermodynamics_Initialize_Parameters()
    ep = 0.0_wp

    write (*, *) 'Case d-e (1) or d-p (2) or p-h (3) ?'
    read (*, *) opt

    write (*, *) 'Minimum qt ?'
    read (*, *) qt_min
    write (*, *) 'Maximum qt ?'
    read (*, *) qt_max
    write (*, *) 'Increment qt ?'
    read (*, *) qt_del

    if (opt == 1) then
        write (*, *) 'Density value ?'
        read (*, *) rho
        write (*, *) 'Energy value ?'
        read (*, *) e
    else if (opt == 2) then
        write (*, *) 'Density value ?'
        read (*, *) rho
        write (*, *) 'Pressure value ?'
        read (*, *) p
    else if (opt == 3) then
        write (*, *) 'enthalpy?'
        read (*, *) h
        write (*, *) 'pressure?'
        read (*, *) p
    end if

    write (*, *) 'Smoothing factor ?'
    read (*, *) dsmooth

    allocate (pbackground(1), epbackground(1), rbackground(1))
    epbackground(1) = ep(1)
    pbackground(1) = p(1)
    rbackground(1) = rho(1)

! ###################################################################
    open (21, file='vapor.dat')
    write (21, *) '# qt, ql, qv, qs(T), r, T, p, e, h'

    qt = qt_min
    do while (qt(1) <= qt_max(1))

        z1(1) = qt(1)
        if (opt == 1) then
            call THERMO_CALORIC_TEMPERATURE(1, z1, e, rho, T, dqldqt)
            call Thermo_Psat_Polynomial(1, T, qs)
            qs = qs/(rho*T*Rv)
            call THERMO_THERMAL_PRESSURE(1, z1, rho, T, p)
            call THERMO_CALORIC_ENTHALPY(1, z1, T, h)

        else if (opt == 2) then
            call THERMO_AIRWATER_RP(1, z1, p, rho, T, dqldqt)
            call Thermo_Psat_Polynomial(1, T, qs)
            qs = qs/(rho*T*Rv)
            call THERMO_CALORIC_ENERGY(1, z1, T, e)
            call THERMO_CALORIC_ENTHALPY(1, z1, T, h)

        else if (opt == 3) then
            call Thermo_Anelastic_PH(1, 1, 1, z1, h)
            s(1) = h(1); s(2:3) = z1(1:2)
            call Thermo_Anelastic_TEMPERATURE(1, 1, 1, s, T)
!        CALL THERMO_AIRWATER_PH_RE(1, z1, p, h, T)
            call Thermo_Psat_Polynomial(1, T, qs)
            qs = 1.0_wp/(p/qs - 1.0_wp)*rd_ov_rv
            qs = qs/(1.0_wp + qs)
            call THERMO_THERMAL_DENSITY(1, z1, p, T, rho)
            call THERMO_CALORIC_ENERGY(1, z1, T, e)

        end if
        write (21, 1000) qt, z1(2), qt - z1(2), qs, rho, T, p, e, h

        qt = qt + qt_del
    end do

    close (21)

    stop

1000 format(9(G_FORMAT_R))

end program SMOOTH
