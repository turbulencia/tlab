#include "dns_error.h"

!########################################################################
! Midpoint interpolation using compact schemes and
! the corresponding 1. order derivative, used in staggered grid for pressure.
! It assumes periodic boundary conditions
!########################################################################

module FDM_Interpolate
    use TLab_Constants, only: wp, wi, pi_wp, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use FDM_Com0_Jacobian
    implicit none
    private

    ! I wonder if the modified wavenumber for stagger_on case should be here...
    type, public :: fdm_interpol_dt
        sequence
        integer mode_fdm                            ! finite-difference method
        integer(wi) size
        real(wp), allocatable :: lu0i(:, :)         ! LU decomposition for interpolation
        real(wp), allocatable :: lu1i(:, :)         ! LU decomposition for 1. derivative inc. interp.
    end type fdm_interpol_dt

    public :: FDM_Interpol_Initialize
    public :: FDM_Interpol
    public :: FDM_Interpol_Der1

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Interpol_Initialize(x, dx, var, wn)
        real(wp), intent(in) :: x(:), dx(:)
        type(fdm_interpol_dt), intent(inout) :: var
        real(wp), intent(inout) :: wn(:)

        ! -------------------------------------------------------------------
        integer(wi) i, nx, nd
        real(wp) :: coef(5)

        !########################################################################
        nx = size(x)
        nd = 3 + 2         ! so far, only tridiagonal systems with periodic bcs are implemented

        var%size = nx

        !########################################################################
        ! interpolation
        if (allocated(var%lu0i)) deallocate (var%lu0i)
        allocate (var%lu0i(nx, nd))

        select case (var%mode_fdm)
        case DEFAULT
            call FDM_C0INT6P_LHS(nx, var%lu0i(:, 1), var%lu0i(:, 2), var%lu0i(:, 3))
        end select

        ! LU decomposition
        call TRIDPFS(nx, var%lu0i(1, 1), var%lu0i(1, 2), var%lu0i(1, 3), var%lu0i(1, 4), var%lu0i(1, 5))

        !########################################################################
        ! first interp. derivative
        if (allocated(var%lu1i)) deallocate (var%lu1i)
        allocate (var%lu1i(nx, nd))

        select case (var%mode_fdm)
        case DEFAULT
            call FDM_C1INT6P_LHS(nx, dx(:), var%lu1i(:, 1), var%lu1i(:, 2), var%lu1i(:, 3))
        end select

        ! LU decomposition
        call TRIDPFS(nx, var%lu1i(1, 1), var%lu1i(1, 2), var%lu1i(1, 3), var%lu1i(1, 4), var%lu1i(1, 5))

        ! -------------------------------------------------------------------
        ! modified wavenumbers; staggered case has different modified wavenumbers!
        do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
            if (i <= nx/2 + 1) then
                wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
            else
                wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
            end if
        end do

        select case (var%mode_fdm)
        case DEFAULT
            coef = [9.0_wp/62.0_wp, 0.0_wp, 63.0_wp/62.0_wp, 17.0_wp/62.0_wp, 0.0_wp]

        end select

        wn(:) = 2.0_wp*(coef(3)*sin(1.0_wp/2.0_wp*wn(:)) + coef(4)/3.0_wp*sin(3.0_wp/2.0_wp*wn(:))) &
                /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)))

        wn(:) = wn(:)/dx(1)           ! normalized by dx

        return
    end subroutine FDM_Interpol_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Interpol(dir, nlines, g, u, result, wrk2d)
        integer(wi), intent(in) :: dir      ! scalar direction flag
        !                                   0 'vp' --> vel. to pre. grid
        !                                   1 'pv' --> pre. to vel. grid
        integer(wi), intent(in) :: nlines   ! number of lines to be solved
        type(fdm_interpol_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: wrk2d(*)

        ! ###################################################################
        ! Interpolation, direction 'vp': vel. --> pre. grid
        if (dir == 0) then
            select case (g%mode_fdm)
            case DEFAULT
                call FDM_C0INTVP6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)

            ! Interpolation, direction 'pv': pre. --> vel. grid
        else if (dir == 1) then
            select case (g%mode_fdm)
            case DEFAULT
                call FDM_C0INTPV6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)
        end if

        return
    end subroutine FDM_Interpol

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Interpol_Der1(dir, nlines, g, u, result, wrk2d)
        integer(wi), intent(in) :: dir      ! scalar direction flag
        !                                   0 'vp' --> vel. to pre. grid
        !                                   1 'pv' --> pre. to vel. grid
        integer(wi), intent(in) :: nlines   ! number of lines to be solved
        type(fdm_interpol_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: wrk2d(*)

        ! ###################################################################
        ! 1st interpolatory derivative, direction 'vp': vel. --> pre. grid
        if (dir == 0) then
            select case (g%mode_fdm)
            case default
                call FDM_C1INTVP6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)

            ! 1st interpolatory derivative, direction 'pv': pre. --> vel. grid
        else if (dir == 1) then
            select case (g%mode_fdm)
            case default
                call FDM_C1INTPV6P_RHS(g%size, nlines, u, result)
            end select
            call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)
        end if

        return
    end subroutine FDM_Interpol_Der1

end module FDM_Interpolate
