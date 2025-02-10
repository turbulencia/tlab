#include "dns_const.h"

program VDIFFUSION
    use TLab_Constants, only: wp, wi
    use TLAB_VARS
    use IO_FIELDS
    use NavierStokes, only: NavierStokes_Initialize_Parameters

    implicit none

    real(wp), dimension(:, :), allocatable, save, target :: x, y, z
    real(wp), dimension(:, :), allocatable :: q, s, s_r
    real(wp), dimension(:), allocatable :: wrk1d, wrk2d, wrk3d

    integer(wi) i, j, ij, iopt
    integer(wi) isize_wrk3d
    real(wp) dummy, error, pi_loc, factor, wavenumber, x_loc, params(0)
    character*(32) fname

! ###################################################################
    call DNS_START

    call TLab_Initialize_Parameters(ifile)
    call NavierStokes_Initialize_Parameters(ifile)

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (wrk1d(isize_wrk1d*inb_wrk1d))
    allocate (wrk2d(isize_wrk2d*inb_wrk2d))
    allocate (wrk3d(isize_wrk3d))
    allocate (q(isize_field, 3))
    allocate (s(isize_field, 1))
    allocate (s_r(isize_field, 1))

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:,1), wrk1d(:,2), wrk1d(:,3))
    call FDM_Initialize(x, g(1), wrk1d(:,1), wrk1d(:,4))
    call FDM_Initialize(y, g(2), wrk1d(:,2), wrk1d(:,4))
    call FDM_Initialize(z, g(3), wrk1d(:,3), wrk1d(:,4))

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
    call IO_Write_Fields(fname, imax, jmax, kmax, itime, 1, s, io_header_s(1))

    q(:, 1) = mean_u; q(:, 2) = C_0_R; q(:, 3) = C_0_R
    fname = trim(adjustl(tag_flow))//'0'
    call IO_Write_Fields(fname, imax, jmax, kmax, itime, 3, q, io_header_q(1))

! ###################################################################
    else if (iopt == 2) then
    write (*, *) 'Iteration?'; read (*, *) itime

    write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
    call IO_Read_Fields(fname, imax, jmax, kmax, itime, 1, 0, s, params)

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
        call IO_Write_Fields(fname, imax, jmax, kmax, itime, 1, wrk3d, io_header_s(1))
    end if

    end if

    stop
end program VDIFFUSION
