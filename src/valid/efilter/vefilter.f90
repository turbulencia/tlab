#include "types.h"
#include "dns_const.h"

program VEFILTER
    use TLAB_TYPES, only: filter_dt, grid_dt
    use TLAB_VARS, only: reynolds, schmidt
    USE OPR_FILTERS

    implicit none

    TINTEGER imax, i, ik, i1bc
    parameter(imax=257)
    TREAL x(imax, 50), u(imax), uf(imax)
    TREAL wrk1d(imax, 10)! , wrk2d(imax), wrk3d(imax) YOU NEED TO USE NEW MEM MANAGEMENT
    type(filter_dt) filter
    type(grid_dt) g

! ###################################################################
    g%size = imax
    g%scale = 2*C_PI_R
    g%mode_fdm = FDM_COM6_JACOBIAN
    reynolds = C_1_R
    schmidt = C_1_R

    write (*, *) 'Periodic (0) or nonperiodic (1) case ?'
    read (*, *) i1bc

    if (i1bc == 0) then
        g%periodic = .true.
        g%uniform = .true.
        filter%bcsmin = DNS_FILTER_BCS_PERIODIC
        filter%bcsmax = DNS_FILTER_BCS_PERIODIC
        do i = 1, imax
            x(i, 1) = M_REAL(i - 1)/M_REAL(imax)*g%scale
        end do
    else
        g%periodic = .false.
        g%uniform = .false.
        filter%bcsmin = DNS_FILTER_BCS_ZERO
        filter%bcsmax = DNS_FILTER_BCS_ZERO
        do i = 1, imax
            x(i, 1) = M_REAL(i - 1)/M_REAL(imax-1)*g%scale
        end do
        ! open (21, file='y.dat')
        ! do i = 1, imax
        !     read (21, *) x(i, 1)
        ! end do
        ! close (21)
        ! g%scale = x(imax, 1) - x(1, 1)
    end if

    call FDM_INITIALIZE(x, g, wrk1d)

    ! ###################################################################

    ! Define the function
    ! if (i1bc == 0) then
    write (*, *) 'Wavenumber ?'
    read (*, *) ik
    u(:) = SIN(C_2_R*C_PI_R/g%scale*M_REAL(ik)*x(:, 1))
    ! else
    !     open (21, file='f.dat')
    !     do i = 1, imax
    !         read (21, *) u(i)
    !         uf(i) = u(i)
    !     end do
    !     close (21)
    ! end if

! ###################################################################
    filter%size = g%size
    filter%periodic = g%periodic
    ! filter%type = DNS_FILTER_COMPACT_CUTOFF
    ! filter%inb_filter = 7
    filter%type = DNS_FILTER_COMPACT
    filter%inb_filter = 10
    filter%parameters(1) = 0.49

    call OPR_FILTER_INITIALIZE(g, filter)

    call OPR_FILTER_1D(1, filter, u, uf)
    ! call OPR_PARTIAL1(1, [0,0], g, u, uf, wrk2d)

    open (20, file='filter.dat')
    do i = 1, imax
        write (20, *) x(i,1), u(i), uf(i)
    end do
    close (20)

    open (20, file='transfer.dat')
    ! do ik = 1, imax/2
    do ik = 1, (imax-1)/2
            u = SIN(C_2_R*C_PI_R/g%scale*M_REAL(ik)*x(:,1))
        call OPR_FILTER_1D(1, filter, u, uf)
        ! call OPR_PARTIAL1(1, [0,0], g, u, uf, wrk2d)
        write (20, *) ik, maxval(uf)
    end do
    close (20)

    stop
end program VEFILTER
