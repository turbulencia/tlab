!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/23 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate compact interpolation schemes.
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
#include "dns_const.h"

program INTERPOL
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Types, only: grid_dt
    use TLab_WorkFlow

    implicit none

    type(grid_dt) :: g, g_pre
    integer(wi) :: jmax, kmax, i, l, test_type, periodic
    real(wp) :: lambda, error, sol

    integer(wi), parameter :: imax = 32, len = 10, inb_grid = 57
    integer(wi), parameter :: imaxp = imax - 1

    real(wp), dimension(imax, inb_grid) :: x
    real(wp), dimension(imaxp, inb_grid) :: x_pre ! pressure grid (for non-periodic case)

    real(wp), dimension(imax) :: x_int, x_aux
    real(wp), dimension(len, imax) :: u, u_int, u_aux, u_a, u_b
    real(wp), dimension(len, imax) :: dudx, dudx_int, dudx_aux
    real(wp), dimension(imax, 5) :: wrk1d
    real(wp), dimension(len) :: wrk2d

! ###################################################################
! Initialize
    g%size = imax
    g%scale = 1.0_wp
    g%uniform = .true.
    jmax = 1
    kmax = 1
    g%mode_fdm1 = FDM_COM6_JACOBIAN
    lambda = 1

! Input
    write (*, *) '........... Test type (1/2/3/4)? ...............'
    write (*, *) 'Interpolation from vel. to pre. grid        = 1'
    write (*, *) 'Interpolation from pre. to vel. grid        = 2'
    write (*, *) '1st Interpol. deriv. from vel. to pre. grid = 3'
    write (*, *) '1st Interpol. deriv. from pre. to vel. grid = 4'
    write (*, *) '................................................'
    read (*, *) test_type
    if (test_type < 1 .or. test_type > 4) then
        write (*, *) 'ERROR, STOP. Chose value between 1 and 4.'
        call TLab_Stop(0)
    end if
    write (*, *) 'Testing     periodic schemes = 0'
    write (*, *) 'Testing non-periodic schemes = 1'
    write (*, *) '................................................'
    read (*, *) periodic
    if (periodic /= 0 .and. periodic /= 1) then
        write (*, *) 'ERROR, STOP. Chose value between 0 and 1.'
        call TLab_Stop(0)
    end if
    if (periodic == 0) then
        g%periodic = .true.
        write (*, *) '.......... Testing periodic schemes ............'
    elseif (periodic == 1) then
        g%periodic = .false.
        write (*, *) '........ Testing non-periodic schemes ..........'
    end if
    write (*, *) '................................................'
! ###################################################################
! Initialize grid
    if (g%periodic) then
        do i = 1, imax
            x(i, 1) = real(i - 1)/real(imax)*g%scale
        end do
    else
        do i = 1, imax
            x(i, 1) = real(i - 1)/real(imax - 1)*g%scale
        end do
    end if

! Velocity grid
    call FDM_Initialize(x, g, wrk1d, wrk1d(:,2))

! Initialize grids (interpolation grid on midpoints)
    if (g%periodic) then
        do i = 1, imax
            x_int(i) = g%nodes(i) + 0.5*g%jac(i, 1)
        end do
    else
        do i = 1, imaxp
            x_int(i) = g%nodes(i) + 0.5*g%jac(i, 1)
        end do
        x_int(imax) = 0.0_wp
    end if
    x_aux(:) = x_int(:)

! Initialize pressure grid (only needed for non-periodic case)
! (here: periodic case implies the usage of uniform grids!)
    do i = 1, imaxp; 
        ! x_pre(i, 1) = x_int(i); 
        wrk1d(i, 1) = x_int(i); 
    end do
    g_pre%size = imaxp; g_pre%scale = x_int(imaxp); g_pre%uniform = g%uniform
    g_pre%mode_fdm1 = g%mode_fdm1; g_pre%periodic = .false.
    call FDM_Initialize(x_pre, g_pre, wrk1d, wrk1d(:,2))

! Define the function + deriv. on both grids
    do i = 1, imax
        do l = 1, len
            u(l, i) = &
                sin(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))
            u_int(l, i) = &
                sin(2.0_wp*pi_wp/g%scale*lambda*x_int(i))
            dudx(l, i) = (2.0_wp*pi_wp/g%scale*lambda) &
                         *cos(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))
            dudx_int(l, i) = (2.0_wp*pi_wp/g%scale*lambda) &
                             *cos(2.0_wp*pi_wp/g%scale*lambda*x_int(i))
            u_aux(l, i) = u_int(l, i)
            dudx_aux(l, i) = dudx_int(l, i)
        end do
    end do

! Switch grids and functions (according to interpolation direction [vp <--> pv])
    if (test_type == 2 .or. test_type == 4) then
        do i = 1, imax
            x_int(i) = x(i, 1)
            x(i, 1) = x_aux(i)
            do l = 1, len
                u_int(l, i) = u(l, i)
                u(l, i) = u_aux(l, i)
                !
                dudx_int(l, i) = dudx(l, i)
                dudx(l, i) = dudx_aux(l, i)
            end do
        end do
    end if
! ###################################################################
! Testing interpolation
! -------------------------------------------------------------------
    if (test_type == 1) then
        write (*, *) '..... Interpolation from vel. to pre. grid .....'
        if (g%periodic .eqv. .true.) then
            call FDM_C0INT6P_LHS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C0INTVP6P_RHS(imax, len, u, u_a)
            call TRIDPFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call TRIDPSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), u_a(1, 1), wrk2d(1))
        elseif (g%periodic .neqv. .true.) then
            call FDM_C0INTVP6_LHS(imaxp, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C0INTVP6_RHS(imax, imaxp, len, u, u_a)
            call TRIDFS(imaxp, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call TRIDSS(imaxp, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), u_a(1, 1))
        end if
! -------------------------------------------------------------------
    elseif (test_type == 2) then
        write (*, *) '..... Interpolation from pre. to vel. grid .....'
        if (g%periodic .eqv. .true.) then
            call FDM_C0INT6P_LHS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C0INTPV6P_RHS(imax, len, u, u_a)
            call TRIDPFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call TRIDPSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), u_a(1, 1), wrk2d(1))
        elseif (g%periodic .neqv. .true.) then
            call FDM_C0INTPV6_LHS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C0INTPV6_RHS(imax, imaxp, len, u, u_a)
            call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), u_a(1, 1))
        end if
! -------------------------------------------------------------------
! Testing interpolatory 1st derivative
! -------------------------------------------------------------------
    elseif (test_type == 3) then
        write (*, *) '1st Interpol. deriv. from vel. to pre. grid ....'
        if (g%periodic .eqv. .true.) then
            call FDM_C1INT6P_LHS(imax, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1INTVP6P_RHS(imax, len, u, u_a)
            call TRIDPFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call TRIDPSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), u_a(1, 1), wrk2d(1))
            do i = 1, imax
                do l = 1, len
                    u_int(l, i) = dudx_int(l, i)
                end do
            end do
        elseif (g%periodic .neqv. .true.) then
            call FDM_C1INTVP6_LHS(imaxp, g_pre%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1INTVP6_RHS(imax, imaxp, len, u, u_a)
            call TRIDFS(imaxp, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call TRIDSS(imaxp, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), u_a(1, 1))
            do i = 1, imax
                do l = 1, len
                    u_int(l, i) = dudx_int(l, i)
                end do
            end do
        end if
! -------------------------------------------------------------------
    elseif (test_type == 4) then
        write (*, *) '1st Interpol. deriv. from pre. to vel. grid ....'
        if (g%periodic .eqv. .true.) then
            call FDM_C1INT6P_LHS(imax, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1INTPV6P_RHS(imax, len, u, u_a)
            call TRIDPFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
            call TRIDPSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), u_a(1, 1), wrk2d(1))
            do i = 1, imax
                do l = 1, len
                    u_int(l, i) = dudx_int(l, i)
                end do
            end do
        elseif (g%periodic .neqv. .true.) then
            call FDM_C1INTPV6_LHS(imax, g%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call FDM_C1INTPV6_RHS(imax, imaxp, len, u, u_a)
            call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
            call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), u_a(1, 1))
            do i = 1, imax
                do l = 1, len
                    u_int(l, i) = dudx_int(l, i)
                end do
            end do
        end if
    end if
! ###################################################################
! IO - Error and function values
    if (g%periodic .neqv. .true.) then
        if (test_type == 1 .or. test_type == 3) then
            u_a(:, imax) = 0.0_wp
            u_int(:, imax) = 0.0_wp
            u(:, imax) = 0.0_wp
        end if
    end if

    open (20, file='interpol.dat')
    error = 0.0_wp
    sol = 0.0_wp
    do i = 1, imax
        do l = 1, len
            write (20, 1000) x(i, 1), x_int(i), u(l, i), u_int(l, i), u_a(l, i), u_a(l, i) - u_int(l, i)
            u_b(l, i) = abs(u_a(l, i) - u_int(l, i))
            error = error + u_b(l, i)*u_b(l, i)
            sol = sol + u_int(l, i)*u_int(l, i)
        end do
    end do
    close (20)

    write (*, 2000) 'Solution L2-norm ...........:', sqrt(g%jac(1, 1)*sol/real(len))
    if (sol == 0.0_wp) stop
    write (*, 2000) 'Error L2-norm ..............:', sqrt(g%jac(1, 1)*error/real(len))
    write (*, 2000) 'Error Linf-norm ............:', maxval(u_b(1, 1:imax))
    write (*, 2000) 'Relative error .............:', sqrt(error)/sqrt(sol)

    stop

1000 format(6(1x, e16.10))
2000 format(1(1x, A, 3x, e16.10))

end program INTERPOL
