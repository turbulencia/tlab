#include "dns_error.h"

subroutine THERMO_READ_CHEMKIN(name)
    use TLAB_CONSTANTS, only: wp, wi
    use THERMO_VARS
    use TLAB_CONSTANTS, only: efile, lfile
    use TLAB_PROCS

    implicit none

    character*(*) name

! -----------------------------------------------------------------------
    real(wp) T1, T2, T3
    integer(wi) i, il, is
    logical frun
    character*15 token
    character*225 wline
    character*80 line, line1, line2, line3
    integer(wi) THERMO_FLAG(MAX_NSP)

    integer, parameter :: i23 = 23

! #######################################################################
! Initialize thermodynamic data structure
    do is = 1, NSP
        THERMO_FLAG(is) = 0
    end do

! Read Thermodynamic file
    open (i23, file=name, status='old')

    rewind (i23)

! Read Header
    read (i23, *) line
    call TLAB_WRITE_ASCII(lfile, line)

    if (trim(adjustl(line)) /= 'THERMO') then
        call TLAB_WRITE_ASCII(efile, 'THERMO_READ_CHEMKIN. Thermodynamic file format error')
        call TLAB_STOP(DNS_ERROR_THERMOFORMAT)
    end if

! Read Temperature ranges
    read (i23, *) T1, T2, T3
    write (wline, *) T1, T2, T3
    call TLAB_WRITE_ASCII(lfile, wline(1:80))

! Remove comments
    frun = .true.
    do while (frun)
        read (i23, '(A80)', end=50) line
        if (line(1:1) /= '!') frun = .false.
    end do

    frun = .true.
    do while (frun)
!    Check for end of file
        read (line, '(A15)', end=50) token
        if (trim(adjustl(token)) == 'END') then
            frun = .false.
            goto 50
        end if

!    Read all relevant information
        read (i23, '(A80)', end=50) line1
        read (i23, '(A80)', end=50) line2
        read (i23, '(A80)', end=50) line3

!    Process lines
        do is = 1, NSP
            if (trim(adjustl(token)) == THERMO_SPNAME(is)) then
                call TLAB_WRITE_ASCII(lfile, line)
                call TLAB_WRITE_ASCII(lfile, line1)
                call TLAB_WRITE_ASCII(lfile, line2)
                call TLAB_WRITE_ASCII(lfile, line3)

!          Required species found, process information
!          Get limit temperatures
                do i = 1, 225
                    wline(i:i) = ' '
                end do
                wline = line(46:75)
                read (wline, *) (THERMO_TLIM(i, is), i=1, 3)

!          Concatenate lines so read is simpler
                wline = line1(1:75)//line2(1:75)//line3(1:75)

                do i = 1, 14
                    il = (i - 1)*15 + 1
                    read (wline(il:il + 14), *) THERMO_AI(i, 1, is)
                end do

                THERMO_FLAG(is) = 1
            end if
        end do

!    Read next line
        read (i23, '(A80)', end=50) line

    end do

50  close (i23)

    do is = 1, NSP
        if (THERMO_FLAG(is) == 0) then
            call TLAB_WRITE_ASCII(efile, 'THERMO_READ_CHEMKIN. Not all thermodynamic data contained in thermo file')
            call TLAB_STOP(DNS_ERROR_THERMOCONT)
        end if
    end do

    return
end subroutine THERMO_READ_CHEMKIN
