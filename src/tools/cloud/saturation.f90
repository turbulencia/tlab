#include "dns_const.h"

program SATURATION

    use TLab_Constants, only: wp, wi
    use TLAB_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use Thermodynamics

    implicit none

    real(wp) t_min, t_max, t_del, t, psat(1), qsat(1), dummy(1), t_loc(1), p(1), dpsat(1)!, dpsat2
    integer(wi) iopt

! ###################################################################
    call TLab_Start

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call Thermodynamics_Initialize_Parameters()

    write (*, *) '1 - Saturation pressure as a function of T'
    write (*, *) '2 - Saturation specific humidity as a function of T-p'
    read (*, *) iopt

    write (*, *) 'Minimum T (cetigrade) ?'
    read (*, *) t_min
    write (*, *) 'Maximum T (centigrade) ?'
    read (*, *) t_max
    write (*, *) 'Increment T ?'
    read (*, *) t_del

    if (iopt == 2) then
        write (*, *) 'Pressure (hPa) ?'
        read (*, *) p
        p = p*100.0_wp
    end if

! ###################################################################
    open (21, file='vapor.dat')
    if (iopt == 1) then
        write (21, *) '# T (C), T (K), psat (Pa), L-ps (J/kg), L-cp (J/kg)'
    else if (iopt == 2) then
        write (21, *) '# T (C), T (K), qsat (g/kg)'
    end if

    t = t_min
    do while (t <= t_max)

        t_loc = t + 273.15
        call Thermo_Psat_Polynomial(1, t_loc, psat)
        call Thermo_dPsat_Polynomial(1, t_loc, dpsat)
        dummy = 1.0_wp/(p/psat - 1.0_wp)*rd_ov_rv
        qsat = dummy/(1.0_wp + dummy)
        if (iopt == 1) then
            write (21, 1000) t, t_loc, psat, dpsat*t_loc**2/psat*Rv, -(Cvl*t_loc + Lvl)
        else if (iopt == 2) then
            write (21, 2000) t, t_loc, qsat*1.0e3_wp
        end if

        t = t + t_del
    end do

    close (21)

    stop

1000 format(5(1x, G_FORMAT_R))
2000 format(3(1x, G_FORMAT_R))

end program SATURATION
