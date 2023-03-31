program SATURATION

#include "types.h"
#include "dns_const.h"

    use TLAB_VARS
    use TLAB_PROCS
    use THERMO_VARS

    implicit none

    TREAL t_min, t_max, t_del, t, psat, qsat, dummy, t_loc, p, dpsat!, dpsat2
    TINTEGER iopt

    integer, parameter :: i1 = 1

! ###################################################################
    call TLAB_START

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call THERMO_INITIALIZE()

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
        write (*, *) 'Pressure (bar) ?'
        read (*, *) p
    end if

! ###################################################################
    open (21, file='vapor.dat')
    if (iopt == 1) then
        write (21, *) '# T (C), T (K), psat (bar), L-ps (J/kg), L-cp (J/kg)'
    else if (iopt == 2) then
        write (21, *) '# T (C), T (K), qsat (g/kg)'
    end if

    t = t_min
    do while (t <= t_max)

        t_loc = (t + 273.15)/TREF
        call THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t_loc, psat)
        call THERMO_POLYNOMIAL_DPSAT(i1, i1, i1, t_loc, dpsat)
        dummy = C_1_R/(MRATIO*p/psat - C_1_R)*rd_ov_rv
        qsat = dummy/(C_1_R + dummy)
        if (iopt == 1) then
            ! dpsat2 = C_0_R
            ! DO ipsat = NPSAT,2,-1
            !    dpsat2 = dpsat2 *t_loc + THERMO_PSAT(ipsat)*M_REAL(ipsat-1)
            ! ENDDO
            ! PRINT*,dpsat-dpsat2
            write (21, 1000) t, t_loc*TREF, psat, dpsat*t_loc**2/psat*WGHT_INV(1)*(RGAS*TREF/WREF), &
                ((THERMO_AI(1, 1, 1) - THERMO_AI(1, 1, 3))*t_loc + THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3))*(RGAS*TREF/WREF)/GRATIO
        else if (iopt == 2) then
            write (21, 2000) t, t_loc, qsat*1.d3
        end if

        t = t + t_del
    end do

    close (21)

    stop

1000 format(5(1x, G_FORMAT_R))
2000 format(3(1x, G_FORMAT_R))

end program SATURATION
