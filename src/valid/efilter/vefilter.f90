#include "types.h"
#include "dns_const.h"

program VEFILTER
    use TLAB_TYPES, only : filter_dt, grid_dt

    implicit none

#include "integers.h"

    TINTEGER imax, i, ik, i1bc, idummy
    parameter(imax=256)
    TREAL scalex
    TREAL x(imax), u(imax), uf(imax)
    TREAL wrk1d(imax, 5), wrk2d(imax), wrk3d(imax)
    type(filter_dt) filter
    type(grid_dt) g

! ###################################################################
    filter%type = DNS_FILTER_COMPACT_CUTOFF
    filter%size = imax
    filter%inb_filter = 7
    
    scalex = C_2_R*C_PI_R

    write (*, *) 'Periodic (0) or nonperiodic (1) case ?'
    read (*, *) i1bc

    if (i1bc == 0) then
        filter%periodic = .true.
        do i = 1, imax
            x(i) = M_REAL(i - 1)/M_REAL(imax)*scalex
        end do
    else
        filter%periodic = .false.
        open (21, file='y.dat')
        do i = 1, imax
!        x(i) = M_REAL(i-1)/M_REAL(imax-1)*scalex
            read (21, *) idummy, x(i)
        end do
        close (21)
        scalex = x(imax) - x(1)
    end if

! ###################################################################
! Define the function
    if (i1bc == 0) then
        write (*, *) 'Wavenumber ?'
        read (*, *) ik
        do i = 1, imax
            u(i) = SIN(C_2_R*C_PI_R/scalex*M_REAL(ik)*x(i) + C_PI_R*C_05_R)
            uf(i) = u(i)
        end do
    else
        open (21, file='f.dat')
        do i = 1, imax
            read (21, *) u(i)
            uf(i) = u(i)
        end do
        close (21)
    end if

! ###################################################################
    call OPR_FILTER_INITIALIZE( g, filter, wrk1d)

    call OPR_FILTER_1D(1, filter, u, uf, wrk1d, wrk2d, wrk3d)

    open (20, file='filter.dat')
    do i = 1, imax
        write (20, *) x(i), u(i), uf(i)
    end do
    close (20)

    open (20, file='transfer.dat')
    do ik = 1, imax /2
        u = SIN(C_2_R*C_PI_R/scalex*M_REAL(ik)*x)
        call OPR_FILTER_1D(1, filter, u, uf, wrk1d, wrk2d, wrk3d)
        write (20, *) ik, maxval(uf)
    end do
    close (20)

    stop
end program VEFILTER
