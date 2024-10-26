#include "types.h"
#include "dns_const.h"

program VDIFFUSION

    use TLAB_VARS
    use IO_FIELDS

    implicit none

    TREAL, dimension(:, :), allocatable, save, target :: x, y, z
    TREAL, dimension(:, :), allocatable :: q, s, s_r
    TREAL, dimension(:), allocatable :: wrk1d, wrk2d, wrk3d

    TINTEGER i, j, ij, iopt
    TINTEGER isize_wrk3d
    TREAL dummy, error, pi_loc, factor, wavenumber, x_loc
    character*(32) fname

! ###################################################################
    call DNS_START

    call IO_READ_GLOBAL(ifile)

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (x(g(1)%size, g(1)%inb_grid))
    allocate (y(g(2)%size, g(2)%inb_grid))
    allocate (z(g(3)%size, g(3)%inb_grid))

    allocate (wrk1d(isize_wrk1d*inb_wrk1d))
    allocate (wrk2d(isize_wrk2d*inb_wrk2d))
    allocate (wrk3d(isize_wrk3d))
    allocate (q(isize_field, 3))
    allocate (s(isize_field, 1))
    allocate (s_r(isize_field, 1))

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

! ###################################################################
    wavenumber = C_1_R
    pi_loc = acos(-C_1_R)

    write (*, *) '1-ICs / 2-Error ?'; read (*, *) iopt

    if (iopt == 1) then
    do j = 1, jmax; do i = 1, imax
            ij = i + imax*(j - 1)
            s(ij, 1) = sin(C_2_R*pi_loc/scalex*wavenumber*x(i))
        end do; end do
    fname = trim(adjustl(tag_scal))//'0'
    call IO_WRITE_FIELDS(fname, IO_SCAL, imax, jmax, kmax, 1, s)

    q(:, 1) = mean_u; q(:, 2) = C_0_R; q(:, 3) = C_0_R
    fname = trim(adjustl(tag_flow))//'0'
    call IO_WRITE_FIELDS(fname, IO_FLOW, imax, jmax, kmax, 3, q)

! ###################################################################
    else if (iopt == 2) then
    write (*, *) 'Iteration?'; read (*, *) itime

    write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
    call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, 1, 0, s)

! Theoretical
    factor = exp(-visc*rtime*(C_2_R*pi_loc/scalex*wavenumber)**2)
    do j = 1, jmax; do i = 1, imax
            ij = i + imax*(j - 1)
            x_loc = x(i) - mean_u*rtime; 
            s_r(ij, 1) = factor*sin(C_2_R*pi_loc/scalex*wavenumber*x_loc)
        end do; end do

! Error
    error = C_0_R
    dummy = C_0_R
    do ij = 1, isize_field
        wrk3d(ij) = s(ij, 1) - s_r(ij, 1)
        error = error + wrk3d(ij)*wrk3d(ij)
        dummy = dummy + s_r(ij, 1)*s_r(ij, 1)
    end do

    if (dummy > C_0_R) then
        write (*, *) 'L-infinity .................: ', maxval(abs(wrk3d))
        write (*, *) 'Absolute error .............: ', sqrt(error/M_REAL(imax*jmax))
        write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
        fname = 'error'
        call IO_WRITE_FIELDS(fname, IO_SCAL, imax, jmax, kmax, 1, wrk3d)
    end if

    end if

    stop
end program VDIFFUSION
