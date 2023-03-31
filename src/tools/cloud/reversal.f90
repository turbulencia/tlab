program REVERSAL

#include "types.h"
#include "dns_const.h"

    use TLAB_VARS
    use TLAB_PROCS
    use THERMO_VARS
    use THERMO_THERMAL
    use THERMO_ANELASTIC
    use THERMO_CALORIC

    implicit none

    TREAL qt_1, qt_2, h_1, h_2, qt, h, x, qsat, dqldqt
    TREAL qvqd, ba_ratio, heat1, heat2, alpha, dummy
    TREAL z1(2), e, rho, p, T, ep, s(3)
    TREAL t_1, t_2, c0, c1, c2
    TREAL r_1, r_2, r_max, r_old, x_c, x_max
    TINTEGER n, nmax, iopt, iup

! ###################################################################
    call TLAB_START()

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call THERMO_INITIALIZE()
    dsmooth = C_0_R

    write (*, *) '1 - Density profile from nondimensional state'
    write (*, *) '2 - Density profile from dimensional state'
    write (*, *) '3 - Coefficient as (enthalpy) table'
    write (*, *) '4 - Coefficient bs (water content) table'
    write (*, *) '5 - Coefficients ratio table'
    write (*, *) '6 - Correction temperature table'
    write (*, *) '7 - Correction temperature table (negative branch)'

    read (*, *) iopt

    if (iopt == 1) then
        write (*, *) 'Initial enthalpy, h_1 ?'
        read (*, *) h_1
        write (*, *) 'Final enthalpy, h_2 ?'
        read (*, *) h_2
        write (*, *) 'Initial water, qt_1 ?'
        read (*, *) qt_1
        write (*, *) 'Final water, qt_2 ?'
        read (*, *) qt_2

        write (*, *) 'pressure (bar) ?'
        read (*, *) p

    else if (iopt == 2) then
        write (*, *) 'Lower temperature (C) ?'
        read (*, *) t_1
        t_1 = (t_1 + 273.15)/TREF
        write (*, *) 'Upper temperature (C) ?'
        read (*, *) t_2
        t_2 = (t_2 + 273.15)/TREF
        write (*, *) 'Lower water content (g/kg) ?'
        read (*, *) qt_1
        qt_1 = qt_1*C_1EM3_R
        write (*, *) 'Upper water content (g/kg) ?'
        read (*, *) qt_2
        qt_2 = qt_2*C_1EM3_R

        write (*, *) 'pressure (bar) ?'
        read (*, *) p

! calculate nondimensional limits
        z1(1) = qt_1
        call THERMO_AIRWATER_PT(1, z1, p, t_1)
        call THERMO_CALORIC_ENTHALPY(1, z1, t_1, h_1)
        call THERMO_THERMAL_DENSITY(1, z1, p, t_1, r_1)

        z1(1) = qt_2
        call THERMO_AIRWATER_PT(1, z1, p, t_2)
        call THERMO_CALORIC_ENTHALPY(1, z1, t_2, h_2)
        call THERMO_THERMAL_DENSITY(1, z1, p, t_2, r_2)

    else if (iopt == 3 .or. iopt == 4 .or. iopt == 5 .or. iopt == 6 .or. iopt == 7) then
        write (*, *) 'pressure (bar) ?'
        read (*, *) p

        write (*, *) 'Minimum temperature (C) ?'
        read (*, *) t_1
        t_1 = (t_1 + 273.15)/TREF
        write (*, *) 'Maximum temperature (C) ?'
        read (*, *) t_2
        t_2 = (t_2 + 273.15)/TREF

        if (iopt == 3) then
            write (*, *) 'Coefficient -as ?'
            read (*, *) ba_ratio
        else if (iopt == 4) then
            write (*, *) 'Coefficient -bs ?'
            read (*, *) ba_ratio
        else if (iopt == 5) then
            write (*, *) 'Coefficients ratio ?'
            read (*, *) ba_ratio
        else
            write (*, *) 'Correction temperature difference (K) ?'
            read (*, *) ba_ratio
            ba_ratio = ba_ratio/TREF
        end if

    end if

    write (*, *) 'number of points ?'
    read (*, *) nmax

    open (21, file='reversal.dat')
! ###################################################################
    if (iopt == 1 .or. iopt == 2) then
        write (21, *) '# x, qt, h, ql, qv, qsat(T), r, T, p, e'
        r_max = C_0_R
        r_old = r_1
        x_c = -C_1_R
        iup = 0
        do n = 1, nmax
            x = M_REAL(n - 1)/M_REAL(nmax - 1)

            qt = qt_1 + x*(qt_2 - qt_1)
            h = h_1 + x*(h_2 - h_1)
            ep = C_0_R

            z1(1) = qt
            call THERMO_AIRWATER_PH(i1, i1, i1, z1, h, ep, p)
            s(1) = h; s(2:3) = z1(1:2)
            call THERMO_ANELASTIC_TEMPERATURE(i1, i1, i1, s, ep, T)
            call THERMO_POLYNOMIAL_PSAT(1, T, qsat)
            qsat = C_1_R/(MRATIO*p/qsat - C_1_R)*rd_ov_rv
            qsat = qsat/(C_1_R + qsat)
            call THERMO_THERMAL_DENSITY(1, z1, p, T, rho)
            call THERMO_CALORIC_ENERGY(1, z1, T, e)

            write (21, 1010) x, qt, h, z1(2), qt - z1(2), qsat, rho, T, p, e

            if ((rho - r_old) > 0 .and. iup == 0) iup = 1
            if (rho < r_1 .and. iup == 1 .and. x_c < C_0_R) x_c = x
            if (rho > r_max) then
                r_max = rho
                x_max = x
            end if

        end do

        write (*, *) 'Maximum density  r_max .................:', r_max
        write (*, *) 'Maximum density diference r_max-r_1 ....:', r_max - r_1
        write (*, *) 'Ratio (r_max-r_1)/(r_1-r_2) ............:', (r_max - r_1)/(r_1 - r_2)
        write (*, *) 'Maximum mixing fraction ................:', x_max
        write (*, *) 'Cross-over mixing fraction .............:', x_c

! ###################################################################
    else if (iopt == 3) then
        write (21, *) '# T (C), T/298K, qt (g/kg)'

        do n = 1, nmax
            t = t_1 + (t_2 - t_1)*M_REAL(n - 1)/M_REAL(nmax - 1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t
            heat2 = heat1*(C_1_R + qvqd) - &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t

            alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

            dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))* &
                    qvqd*(C_1_R + qvqd/rd_ov_rv)
            dummy = dummy + &
                    THERMO_AI(1, 1, 2) + qvqd*(THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3)) - THERMO_AI(1, 1, 3)

            dummy = (alpha/(ba_ratio*t) - THERMO_AI(1, 1, 3))/dummy

            qt = C_1_R - dummy

            if (qt < qsat) exit

            write (21, *) t*TREF - 273.15, t, qt*1000

        end do

! look for crossing point with saturation curve using bracketing
        dummy = (t_2 - t_1)/M_REAL(nmax - 1)
        t_2 = t
        t_1 = t - dummy
        do while ((t_2 - t_1)/t_1 > C_1EM6_R)

            t = C_05_R*(t_2 + t_1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t
            heat2 = heat1*(C_1_R + qvqd) - &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t

            alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

            dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))* &
                    qvqd*(C_1_R + qvqd/rd_ov_rv)
            dummy = dummy + &
                    THERMO_AI(1, 1, 2) + qvqd*(THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3)) - THERMO_AI(1, 1, 3)

            dummy = (alpha/(ba_ratio*t) - THERMO_AI(1, 1, 3))/dummy

            qt = C_1_R - dummy

            if (qt > qsat) t_1 = t
            if (qt <= qsat) t_2 = t

        end do

        write (21, *) t*TREF - 273.15, t, qt*1000

! ###################################################################
    else if (iopt == 4) then
        write (21, *) '# T (C), T/298K, qt (g/kg)'

        do n = 1, nmax
            t = t_1 + (t_2 - t_1)*M_REAL(n - 1)/M_REAL(nmax - 1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t
            heat2 = heat1*(C_1_R + qvqd) - &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t

            alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

            dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))* &
                    qvqd*(C_1_R + qvqd/rd_ov_rv)
            dummy = dummy + &
                    THERMO_AI(1, 1, 2) + qvqd*(THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3)) - THERMO_AI(1, 1, 3)

            c2 = ba_ratio*dummy
            c1 = -(dummy*(C_1_R + ba_ratio) + (dummy + THERMO_AI(1, 1, 3))*ba_ratio - alpha*heat2/t)
            c0 = (C_1_R + ba_ratio)*(dummy + THERMO_AI(1, 1, 3)) - alpha*heat2/t

            qt = (-c1 + sqrt(c1*c1 - C_4_R*c0*c2))/(C_2_R*c2)

            if ((-c1 - sqrt(c1*c1 - C_4_R*c0*c2))/(C_2_R*c2) < C_1_R) then
                print *, 'error, second root also less than 1'
            end if

            if (qt < qsat) exit

            write (21, *) t*TREF - 273.15, t, qt*1000

        end do

! look for crossing point with saturation curve using bracketing
        dummy = (t_2 - t_1)/M_REAL(nmax - 1)
        t_2 = t
        t_1 = t - dummy
        do while ((t_2 - t_1)/t_1 > C_1EM6_R)

            t = C_05_R*(t_2 + t_1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t
            heat2 = heat1*(C_1_R + qvqd) - &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t

            alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

            dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))* &
                    qvqd*(C_1_R + qvqd/rd_ov_rv)
            dummy = dummy + &
                    THERMO_AI(1, 1, 2) + qvqd*(THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3)) - THERMO_AI(1, 1, 3)

            c2 = ba_ratio*dummy
            c1 = -(dummy*(C_1_R + ba_ratio) + (dummy + THERMO_AI(1, 1, 3))*ba_ratio - alpha*heat2/t)
            c0 = (C_1_R + ba_ratio)*(dummy + THERMO_AI(1, 1, 3)) - alpha*heat2/t

            qt = (-c1 + sqrt(c1*c1 - C_4_R*c0*c2))/(C_2_R*c2)

            if ((-c1 - sqrt(c1*c1 - C_4_R*c0*c2))/(C_2_R*c2) < C_1_R) then
                print *, 'error, second root also less than 1'
            end if

            if (qt > qsat) t_1 = t
            if (qt <= qsat) t_2 = t

        end do

        write (21, *) t*TREF - 273.15, t, qt*1000

! ###################################################################
    else if (iopt == 5) then
        write (21, *) '# T (C), T/298K, qt (g/kg)'

        do n = 1, nmax
            t = t_1 + (t_2 - t_1)*M_REAL(n - 1)/M_REAL(nmax - 1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t
            heat2 = heat1*(C_1_R + qvqd) - &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t

            alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

            dummy = (heat2 - ba_ratio)/t*alpha
            dummy = dummy - &
                    heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))* &
                    qvqd*(C_1_R + qvqd/rd_ov_rv)
            dummy = dummy - &
                    THERMO_AI(1, 1, 2) - qvqd*(THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))
            dummy = dummy/THERMO_AI(1, 1, 3)

            qt = dummy/(C_1_R + dummy)

            if (qt < qsat) exit

            write (21, *) t*TREF - 273.15, t, qt*1000

        end do

! look for crossing point with saturation curve using bracketing
        dummy = (t_2 - t_1)/M_REAL(nmax - 1)
        t_2 = t
        t_1 = t - dummy
        do while ((t_2 - t_1)/t_1 > C_1EM6_R)

            t = C_05_R*(t_2 + t_1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t
            heat2 = heat1*(C_1_R + qvqd) - &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t

            alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

            dummy = (heat2 - ba_ratio)/t*alpha
            dummy = dummy - &
                    heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))* &
                    qvqd*(C_1_R + qvqd/rd_ov_rv)
            dummy = dummy - &
                    THERMO_AI(1, 1, 2) - qvqd*(THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))
            dummy = dummy/THERMO_AI(1, 1, 3)

            qt = dummy/(C_1_R + dummy)

            if (qt > qsat) t_1 = t
            if (qt <= qsat) t_2 = t

        end do

        write (21, *) t*TREF - 273.15, t, qt*1000

! ###################################################################
    else if (iopt == 6) then
        write (21, *) '# T (C), T/298K, qt (g/kg)'

        do n = 1, nmax
            t = t_1 + (t_2 - t_1)*M_REAL(n - 1)/M_REAL(nmax - 1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t

            qt = (qsat*heat1 + ba_ratio)/(heat1 - (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t)

            if (qt > qsat) write (21, *) t*TREF - 273.15, t, qt*1000

        end do

! ###################################################################
    else if (iopt == 7) then
        write (21, *) '# T (C), T/298K, qt (g/kg)'

        do n = 1, nmax
            t = t_1 + (t_2 - t_1)*M_REAL(n - 1)/M_REAL(nmax - 1)

            call THERMO_POLYNOMIAL_PSAT(1, t, qvqd)
            qvqd = C_1_R/(MRATIO*p/qvqd - C_1_R)*rd_ov_rv
            qsat = qvqd/(C_1_R + qvqd)

            heat1 = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3) + &
                    (THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t

            qt = -ba_ratio/((THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 2))*t)

            if (qt < qsat) write (21, *) t*TREF - 273.15, t, qt*1000

        end do
    end if

    close (21)

    call TLAB_STOP

    stop
1010 format(10(1x, G_FORMAT_R))
end program REVERSAL
