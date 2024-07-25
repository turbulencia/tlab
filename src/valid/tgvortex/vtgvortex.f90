#include "types.h"
#include "dns_const.h"

program VTGVORTEX

    use TLAB_VARS
    use IO_FIELDS
    use OPR_FOURIER
    use OPR_ELLIPTIC
    use FI_SOURCES

    implicit none

    TREAL, dimension(:, :), allocatable, save, target :: x, y, z
    TREAL, dimension(:, :), allocatable :: txc, q
    TREAL, dimension(:), allocatable :: wrk1d, wrk2d, wrk3d

    TINTEGER ij, iv, iopt
    TREAL dummy, error
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
    allocate (q(isize_field, 4))
    allocate (txc(isize_txc_field, 4))

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call OPR_ELLIPTIC_INITIALIZE()

! ###################################################################
    call OPR_FOURIER_INITIALIZE()

    write (*, *) '1-ICs / 2-Error ?'; read (*, *) iopt

    if (iopt == 1) then
        rtime = C_0_R
        call FLOW_TAYLORGREEN(imax, jmax, kmax, rtime, visc, x, y, z, q(1, 1), q(1, 2), q(1, 3), q(1, 4))
        q(:, 1) = q(:, 1) + mean_u
        fname = 'pt0'
        call IO_WRITE_FIELDS(fname, IO_FLOW, imax, jmax, kmax, 3, q)

        q(:, 4) = C_0_R
        fname = 'sc0'
        call IO_WRITE_FIELDS(fname, IO_SCAL, imax, jmax, kmax, 1, q(1, 4))

! ###################################################################
    else if (iopt == 2) then
        write (*, *) 'Iteration?'
        read (*, *) itime

        write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
        call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, 3, 0, q)
        txc(:, 1) = C_0_R; txc(:, 4) = C_0_R
!  CALL FI_FORCING_1(imax,jmax,kmax,  &
!       rtime,visc, txc(1,1),txc(1,4), q(1,1),q(1,2),q(1,3),q(1,4))
!  CALL FI_FORCING_0(imax,jmax,kmax, rtime,visc, q(1,1),q(1,2), txc(1,1),txc(1,4))
!  CALL IO_READ_FIELDS(fname, IO_FLOW, imax,jmax,kmax, i3,i0, q, wrk3d)

        call FI_PRESSURE_BOUSSINESQ(q(1, 1), q(1, 2), q(1, 3), txc(1, 4), q(1, 4), &
                                    txc(1, 1), txc(1, 2), txc(1, 3))

! Theoretical Taylor-Green in array txc
        x = x - mean_u*rtime
        call FLOW_TAYLORGREEN(imax, jmax, kmax, rtime, visc, x, y, z, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
        txc(:, 1) = txc(:, 1) + mean_u

! Error
        do iv = 1, 4
            error = C_0_R
            dummy = C_0_R
            do ij = 1, isize_field !-imax
                dummy = dummy + txc(ij, iv)*txc(ij, iv)
                txc(ij, iv) = q(ij, iv) - txc(ij, iv)
                error = error + txc(ij, iv)*txc(ij, iv)
            end do
!     print*,dummy*dx(1)*dy(1)

            if (dummy > C_0_R) then
                write (*, *) 'L-infinity error ...: ', maxval(abs(txc(1:isize_field - imax, iv)))
                write (*, *) 'L-2 error ..........: ', sqrt(error*g(1)%jac(1, 1)*g(2)%jac(1, 1))
                write (*, *) 'Relative error .....: ', sqrt(error)/sqrt(dummy)
                write (fname, *) iv; fname = 'error'//trim(adjustl(fname))
                call IO_WRITE_FIELDS(fname, IO_SCAL, imax, jmax, kmax, 1, txc(1, iv))
            end if

        end do

    end if

    stop
end program VTGVORTEX

!########################################################################
!########################################################################

subroutine FLOW_TAYLORGREEN(nx, ny, nz, rtime, visc, x, y, z, u, v, w, p)

    implicit none

    TINTEGER nx, ny, nz
    TREAL rtime, visc
    TREAL, dimension(*) :: x, y, z
    TREAL, dimension(nx, ny, nz) :: u, v, w, p

! -----------------------------------------------------------------------
    TINTEGER i, j, k
    TREAL pi_loc, omega, factor, sigma

! #######################################################################
    pi_loc = acos(-C_1_R); omega = C_2_R*pi_loc

    do k = 1, nz
        do j = 1, ny; do i = 1, nx
!        u(i,j,k) = SIN(x(i)*omega)*COS(y(j)*omega)
!        v(i,j,k) =-COS(x(i)*omega)*SIN(y(j)*omega)
                u(i, j, k) = C_05_R*(sin((x(i) + y(j))*omega) + sin((x(i) - y(j))*omega))
                v(i, j, k) = -C_05_R*(sin((x(i) + y(j))*omega) - sin((x(i) - y(j))*omega))
                p(i, j, k) = (cos(C_2_R*omega*x(i)) + cos(C_2_R*omega*y(j)) - C_1_R)*C_025_R

!        u(i,j,k) = sin(omega*x(i))*       sin(C_2_R*omega*y(j))
!        v(i,j,k) =-cos(omega*x(i))*(C_1_R-cos(C_2_R*omega*y(j)))*C_05_R
!        p(i,j,k) = cos(C_2_R*omega*x(i))*(C_2_R-cos(C_2_R*omega*y(j)))/C_8_R &
!             - C_05_R*(sin(omega*y(j)))**4

!        u(i,j,k) = sin(pi_loc*x(i))*sin(pi_loc*x(i))*sin(C_2_R*pi_loc*y(j))
!        v(i,j,k) =-sin(C_2_R*pi_loc*x(i))*sin(pi_loc*y(j))*sin(pi_loc*y(j))
!        p(i,j,k) = sin(C_2_R*pi_loc*x(i))*sin(C_2_R*pi_loc*y(j))

            end do; end do
    end do
    w = C_0_R

! -----------------------------------------------------------------------
    sigma = omega*omega*C_2_R*visc
    factor = exp(-sigma*rtime)
!  factor = cos(omega*rtime) - sigma/omega*sin(omega*rtime)
!  factor = EXP( -sigma*rtime + C_1_R/omega*(C_1_R-cos(omega*rtime)) )
!  factor = cos(omega*rtime)

    u = u*factor; v = v*factor; w = w*factor
    p = p*factor*factor

    return
end subroutine FLOW_TAYLORGREEN
