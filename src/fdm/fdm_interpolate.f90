#include "dns_error.h"

module FDM_Interpolate
    use TLab_Constants, only: wp, wi, efile, wfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use FDM_Com0_Jacobian
    use FDM_Derivative, only: fdm_dt
    implicit none
    private

    ! I should create a new derived type with interpol data, which is used only here.
    type, public :: fdm_interpol_dt
        sequence
        integer mode_fdm                            ! finite-difference method
        real(wp), allocatable :: lu0i(:, :)         ! LU decomposition for interpolation
        real(wp), allocatable :: lu1i(:, :)         ! LU decomposition for 1. derivative inc. interp.
    end type fdm_interpol_dt

    public :: FDM_Interpol_Initialize
    public :: FDM_Interpol
    public :: FDM_Interpol_Der1

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Interpol_Initialize(g)!, var)
        type(fdm_dt), intent(inout) :: g
        ! type(fdm_interpol_dt), intent(inout) :: var

        ! -------------------------------------------------------------------
        integer(wi) nx, nd

        !########################################################################
        nx = g%size
        nd = 3 + 2         ! so far, only tridiagonal systems with periodic bcs are implemented

        !########################################################################
        ! interpolation
        ! if (allocated(var%lu0i)) deallocate (var%lu0i)
        ! allocate (var%lu0i(nx, nd))

        select case (g%mode_fdm1)
        case DEFAULT
            call FDM_C0INT6P_LHS(g%size, g%lu0i(:, 1), g%lu0i(:, 2), g%lu0i(:, 3))
        end select

        ! LU decomposition
        call TRIDPFS(g%size, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5))

        !########################################################################
        ! first interp. derivative
        ! if (allocated(var%lu1i)) deallocate (var%lu1i)
        ! allocate (var%lu1i(nx, nd))

        select case (g%mode_fdm1)
        case DEFAULT
            call FDM_C1INT6P_LHS(g%size, g%jac, g%lu1i(:, 1), g%lu1i(:, 2), g%lu1i(:, 3))
        end select

        ! LU decomposition
        call TRIDPFS(g%size, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5))

        return
    end subroutine FDM_Interpol_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Interpol(dir, nlines, g, u, result, wrk2d)
        integer(wi), intent(in) :: dir      ! scalar direction flag
        !                                   0 'vp' --> vel. to pre. grid
        !                                   1 'pv' --> pre. to vel. grid
        integer(wi), intent(in) :: nlines   ! number of lines to be solved
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: wrk2d(*)

        ! ###################################################################
        ! Interpolation, direction 'vp': vel. --> pre. grid
        if (dir == 0) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case DEFAULT
                    call FDM_C0INTVP6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'FDM_Interpol. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
            ! Interpolation, direction 'pv': pre. --> vel. grid
        else if (dir == 1) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case DEFAULT
                    call FDM_C0INTPV6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu0i(1, 1), g%lu0i(1, 2), g%lu0i(1, 3), g%lu0i(1, 4), g%lu0i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'FDM_Interpol. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
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
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: u(nlines, g%size)
        real(wp), intent(out) :: result(nlines, g%size)
        real(wp), intent(inout) :: wrk2d(*)

        ! ###################################################################
        ! 1st interpolatory derivative, direction 'vp': vel. --> pre. grid
        if (dir == 0) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case default
                    call FDM_C1INTVP6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'FDM_Interpol_Der1. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
            ! 1st interpolatory derivative, direction 'pv': pre. --> vel. grid
        else if (dir == 1) then
            if (g%periodic) then
                select case (g%mode_fdm1)
                case default
                    call FDM_C1INTPV6P_RHS(g%size, nlines, u, result)
                end select
                call TRIDPSS(g%size, nlines, g%lu1i(1, 1), g%lu1i(1, 2), g%lu1i(1, 3), g%lu1i(1, 4), g%lu1i(1, 5), result, wrk2d)
            else
                call TLab_Write_ASCII(efile, 'FDM_Interpol_Der1. Non-periodic case not implemented.')
                call TLab_Stop(DNS_ERROR_NOTIMPL)
            end if
        end if

        return
    end subroutine FDM_Interpol_Der1

end module FDM_Interpolate
