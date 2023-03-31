program SMOOTH

#include "types.h"
#include "dns_const.h"

    use TLAB_VARS
    use TLAB_PROCS
    use THERMO_VARS
    use THERMO_THERMAL
    use THERMO_ANELASTIC
    use THERMO_CALORIC

    implicit none

    TREAL qt_min(1), qt_max(1), qt_del(1), qt(1), qs(1), dqldqt(1)
    TREAL z1(2), e(1), rho(1), p(1), T(1), h(1), ep(1), s(3)
    TINTEGER opt
    integer, parameter :: i1 = 1

! ###################################################################
    call TLAB_START

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call THERMO_INITIALIZE
    ep = C_0_R

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

! ###################################################################
    open (21, file='vapor.dat')
    write (21, *) '# qt, ql, qv, qs(T), r, T, p, e, h'

    qt = qt_min
    do while (qt(1) <= qt_max(1))

        z1(1) = qt(1)
        if (opt == 1) then
            call THERMO_CALORIC_TEMPERATURE(1, z1, e, rho, T, dqldqt)
            call THERMO_POLYNOMIAL_PSAT(i1, i1, i1, T, qs)
            qs = qs/(rho*T*WGHT_INV(1))
            call THERMO_THERMAL_PRESSURE(1, z1, rho, T, p)
            call THERMO_CALORIC_ENTHALPY(1, z1, T, h)

        else if (opt == 2) then
            call THERMO_AIRWATER_RP(i1, i1, i1, z1, p, rho, T, dqldqt)
            call THERMO_POLYNOMIAL_PSAT(i1, i1, i1, T, qs)
            qs = qs/(rho*T*WGHT_INV(1))
            call THERMO_CALORIC_ENERGY(1, z1, T, e)
            call THERMO_CALORIC_ENTHALPY(1, z1, T, h)

        else if (opt == 3) then
            call THERMO_AIRWATER_PH(i1, i1, i1, z1, h, ep, p)
            s(1) = h(1); s(2:3) = z1(1:2)
            call THERMO_ANELASTIC_TEMPERATURE(i1, i1, i1, s, ep, T)
!        CALL THERMO_AIRWATER_PH_RE(i1,i1,i1, z1, p, h, T)
            call THERMO_POLYNOMIAL_PSAT(i1, i1, i1, T, qs)
            qs = C_1_R/(MRATIO*p/qs - C_1_R)*rd_ov_rv
            qs = qs/(C_1_R + qs)
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
