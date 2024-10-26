#include "dns_const.h"

!########################################################################
!#
!# Top-hat filter
!# Trapezoidal rule
!# Free boundary conditions
!#  (Ghost cells w/ linear extrapolation of function from the two nodes next to boundary)
!#
!########################################################################

module Filters_Tophat
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: FLT_T1_COEFFS
    public :: FLT_T1, FLT_T1ND, FLT_T1NDD2, FLT_T1NDD4, FLT_T1NDD6
    public :: FLT_T1P, FLT_T1PD2, FLT_T1PD4, FLT_T1P_ND                 ! periodic cases

contains
!########################################################################
!########################################################################
    subroutine FLT_T1_COEFFS(scalex, x, f, wrk1d)
        use TLab_Types, only: filter_dt
        real(wp), intent(IN) :: scalex
        real(wp), dimension(*), intent(IN) :: x
        type(filter_dt), intent(INOUT) :: f
        real(wp), dimension(f%size, 3), intent(INOUT) :: wrk1d

! -----------------------------------------------------------------------
        real(wp) dum
        integer(wi) i, ii, im, ic, ip, nx

! #######################################################################
        wrk1d = 0.0_wp
        nx = int(f%parameters(1)) ! delta

! #######################################################################
        if (f%uniform) then ! Only valid for free bcs
            if (.not. f%periodic) then ! I only need info for the two nodes next to the boundary
                do i = 1, nx/2
                    im = nx/2 - i + 1
                    dum = 0.0_wp
                    do ii = 1, im
                        dum = dum + real(ii, wp)
                    end do
                    ip = (i - 1)*2 + 1
                    f%coeffs(ip, 1) = dum + 0.5_wp*real(im + 1, wp)
                    if (nx == 2) then
                        f%coeffs(ip + 1, 1) = 0.5_wp - dum + 0.5_wp*real(im, wp)
                    else
                        f%coeffs(ip + 1, 1) = 1.0_wp - dum + 0.5_wp*real(im, wp)
                    end if
                end do
            end if

! #######################################################################
        else
            ! calculate delta(i)
            do i = 1, f%size - 1
                wrk1d(i, 1) = x(i + 1) - x(i)
            end do
            if (f%periodic) then
                wrk1d(f%size, 1) = scalex - (x(f%size) - x(1))
            else
                wrk1d(f%size, 1) = wrk1d(f%size - 1, 1)
            end if
            ! calculate deltasum(i) = delta(i-1)+delta(i)
            do i = 2, f%size
                wrk1d(i, 2) = wrk1d(i - 1, 1) + wrk1d(i, 1)
            end do
            if (f%periodic) then
                wrk1d(1, 2) = wrk1d(f%size, 1) + wrk1d(1, 1)
            else
                wrk1d(1, 2) = wrk1d(1, 1)*2.0_wp
            end if
            ! calculate deltaf(i) = delta(i-nx/2) + delta(i-nx/2+1) + ... + delta(i+nx/2-1)
            if (f%periodic) then
                do i = 1, f%size
                    wrk1d(i, 3) = 0.0_wp
                    do ii = i - nx/2, i + nx/2 - 1
                        im = ii + f%size - 1
                        im = mod(im, f%size) + 1
                        wrk1d(i, 3) = wrk1d(i, 3) + wrk1d(im, 1)
                    end do
                end do
            else
                do i = 1, nx/2
                    wrk1d(i, 3) = wrk1d(1, 1)*real(nx/2 - i + 1, wp)
                    do ii = 1, i + nx/2 - 1
                        wrk1d(i, 3) = wrk1d(i, 3) + wrk1d(ii, 1)
                    end do
                end do
                do i = 1 + nx/2, f%size - nx/2
                    do ii = i - nx/2, i + nx/2 - 1
                        wrk1d(i, 3) = wrk1d(i, 3) + wrk1d(ii, 1)
                    end do
                end do
                do i = f%size - nx/2 + 1, f%size
                    wrk1d(i, 3) = wrk1d(f%size, 1)*real(nx/2 - (f%size - i), wp)
                    do ii = i - nx/2, f%size - 1
                        wrk1d(i, 3) = wrk1d(i, 3) + wrk1d(ii, 1)
                    end do
                end do
            end if

! -----------------------------------------------------------------------
!    construct the coefficients array as if periodic; corrections for nonperiodic below
            do i = 1, f%size
                ii = i - nx/2               ! I need to use delta(ii)
                im = ii + f%size - 1          ! The index im deals with periodicity
                im = mod(im, f%size) + 1
                ic = ii - i + nx/2 + 1
                ip = (i - 1)*(nx + 1) + ic
                f%coeffs(ip, 1) = 0.5_wp*wrk1d(im, 1)/wrk1d(i, 3)

                do ii = i - nx/2 + 1, i + nx/2 - 1 ! I need to use deltasum(ii)
                    im = ii + f%size - 1       ! The index im deals with periodicity
                    im = mod(im, f%size) + 1
                    ic = ii - i + nx/2 + 1
                    ip = (i - 1)*(nx + 1) + ic
                    f%coeffs(ip, 1) = 0.5_wp*wrk1d(im, 2)/wrk1d(i, 3)
                end do

                ii = i + nx/2               ! I need to use delta(ii-1)
                im = (ii - 1) + f%size - 1      ! The index im deals with periodicity
                im = mod(im, f%size) + 1
                ic = ii - i + nx/2 + 1
                ip = (i - 1)*(nx + 1) + ic
                f%coeffs(ip, 1) = 0.5_wp*wrk1d(im, 1)/wrk1d(i, 3)

            end do

! -----------------------------------------------------------------------
! boundary treatment
            select case (f%BcsMin)

            case (DNS_FILTER_BCS_FREE)
                do i = 1, nx/2
                    im = nx/2 - i + 1
                    dum = 0.0_wp
                    do ii = 1, im
                        dum = dum + real(ii, wp)
                    end do
                    ic = nx/2 - i + 2
                    ip = (i - 1)*(nx + 1) + ic
                    f%coeffs(ip, 1) = f%coeffs(ip, 1) + &
                                      0.5_wp*(wrk1d(1, 2)*(dum - 1.0_wp) + wrk1d(1, 1)*real(im + 1, wp))/wrk1d(i, 3)
                    ic = nx/2 - i + 3
                    ip = (i - 1)*(nx + 1) + ic
                    f%coeffs(ip, 1) = f%coeffs(ip, 1) - &
                                      0.5_wp*(wrk1d(1, 2)*(dum - real(im, wp)) + wrk1d(1, 1)*real(im, wp))/wrk1d(i, 3)
! pad with ceros
                    do ic = 1, nx/2 - i + 1
                        ip = (i - 1)*(nx + 1) + ic
                        f%coeffs(ip, 1) = 0.0_wp
                    end do
                end do

            case (DNS_FILTER_BCS_SOLID)
                do i = 1, nx/2
                    im = nx/2 - i + 1
                    ic = nx/2 - i + 2
                    ip = (i - 1)*(nx + 1) + ic
                    f%coeffs(ip, 1) = 0.5_wp*real(2*im + 1, wp)*wrk1d(1, 1)/wrk1d(i, 3)
! pad with ceros
                    do ic = 1, nx/2 - i + 1
                        ip = (i - 1)*(nx + 1) + ic
                        f%coeffs(ip, 1) = 0.0_wp
                    end do
                end do

            end select

! -----------------------------------------------------------------------
            select case (f%BcsMax)

            case (DNS_FILTER_BCS_FREE)
                do i = f%size - nx/2 + 1, f%size
                    im = i - f%size + nx/2
                    dum = 0.0_wp
                    do ii = 1, im
                        dum = dum + real(ii, wp)
                    end do
                    ic = nx + 1 - im - 1
                    ip = (i - 1)*(nx + 1) + ic
                    f%coeffs(ip, 1) = f%coeffs(ip, 1) - &
                                      0.5_wp*(wrk1d(f%size, 2)*(dum - real(im, wp)) + wrk1d(f%size, 1)*real(im, wp))/wrk1d(i, 3)
                    ic = nx + 1 - im
                    ip = (i - 1)*(nx + 1) + ic
                    f%coeffs(ip, 1) = f%coeffs(ip, 1) + &
                                      0.5_wp*(wrk1d(f%size, 2)*(dum - 1.0_wp) + wrk1d(f%size, 1)*real(im + 1, wp))/wrk1d(i, 3)
! pad with ceros
                    do ic = nx + 1 - im + 1, nx + 1
                        ip = (i - 1)*(nx + 1) + ic
                        f%coeffs(ip, 1) = 0.0_wp
                    end do
                end do

            case (DNS_FILTER_BCS_SOLID)
                do i = f%size - nx/2 + 1, f%size
                    im = i - f%size + nx/2
                    ic = nx + 1 - im
                    ip = (i - 1)*(nx + 1) + ic
                    f%coeffs(ip, 1) = 0.5_wp*real(2*im + 1, wp)*wrk1d(f%size, 1)/wrk1d(i, 3)
! pad with ceros
                    do ic = nx + 1 - im + 1, nx + 1
                        ip = (i - 1)*(nx + 1) + ic
                        f%coeffs(ip, 1) = 0.0_wp
                    end do
                end do

            end select

        end if

        return
    end subroutine FLT_T1_COEFFS

!########################################################################
! Uniform grid
!########################################################################
    subroutine FLT_T1(kmax, ijmax, nx, cf, z1, zf1)
        integer(wi), intent(IN) :: kmax, ijmax     ! size of line, number of lines (or size of chunk)
        integer(wi), intent(IN) :: nx              ! filter size
        real(wp), intent(IN) :: cf(2, nx/2)      ! coefficients for BCs; only affect the two nodes next to boundary
        real(wp), intent(IN) :: z1(ijmax, kmax)  ! Field to filter, arranged as kmax chunks
        real(wp), intent(OUT) :: zf1(ijmax, kmax) ! Filtered field

! -------------------------------------------------------------------
! The implementation is based on local index ii varying around global index k
! The index im is just used at the upper boundary to express coefficients in terms of lower boundary
! The index ij is just the dummy index to apply filter over all ijmax lines
        integer(wi) k, ii, im
        integer(wi) ij
! In a uniform grid, I just need to divide by the size of the stencil nx
! There is a factor 1/2, however, that is factorized out from the trapezoidal rule
        real(wp) dmul05, dmul1

        dmul1 = 1.0_wp/real(nx, wp)
        dmul05 = 0.5_wp/real(nx, wp)

! ###########################################
! # First nx/2 points where bcs are imposed #
! ###########################################
        do k = 1, nx/2
! first 2 points reflect bc
            do ij = 1, ijmax
                zf1(ij, k) = dmul1*(z1(ij, 1)*cf(1, k) + z1(ij, 2)*cf(2, k))
            end do
! inner points
            do ii = 3, k + nx/2 - 1
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + dmul1*z1(ij, ii)
                end do
            end do
! last point
            if (nx > 2) then
                ii = k + nx/2
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + dmul05*z1(ij, ii)
                end do
            end if
        end do

! #################################################
! # Interior points where boundary does not enter #
! #################################################
        do k = 1 + nx/2, kmax - nx/2
! first point
            ii = k - nx/2
            do ij = 1, ijmax
                zf1(ij, k) = dmul05*z1(ij, ii)
            end do
! inner points
            do ii = k - nx/2 + 1, k + nx/2 - 1
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + dmul1*z1(ij, ii)
                end do
            end do
! last point
            ii = k + nx/2
            do ij = 1, ijmax
                zf1(ij, k) = zf1(ij, k) + dmul05*z1(ij, ii)
            end do
        end do

! ##########################################
! # Last nx/2 points where bcs are imposed #
! ##########################################
        do k = kmax - nx/2 + 1, kmax
! last 2 points account for the bc
            im = kmax - k + 1
            do ij = 1, ijmax
                zf1(ij, k) = dmul1*(z1(ij, kmax)*cf(1, im) + z1(ij, kmax - 1)*cf(2, im))
            end do
! inner points
            do ii = k - nx/2 + 1, kmax - 2
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + dmul1*z1(ij, ii)
                end do
            end do
! first point
            if (nx > 2) then
                ii = k - nx/2
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + dmul05*z1(ij, ii)
                end do
            end if

        end do

        return
    end subroutine FLT_T1

!########################################################################
! Non-uniform grid
!########################################################################
    subroutine FLT_T1ND(kmax, ijmax, nx, cf, z1, zf1)
        integer(wi), intent(IN) :: kmax, ijmax     ! size of line, number of lines (or size of chunk)
        integer(wi), intent(IN) :: nx              ! filter size
        real(wp), intent(IN) :: cf(nx + 1, kmax)   ! coefficients
        real(wp), intent(IN) :: z1(ijmax, kmax)  ! Field to filter, arranged as kmax chunks
        real(wp), intent(OUT) :: zf1(ijmax, kmax) ! Filtered field

! -------------------------------------------------------------------
        real(wp) dum
! The implementation is based on local index ii varying around global index k
! The offset iiorg is set to get counter ic equal to 1,2,3...
! The index ij is just the dummy index to apply filter over all ijmax lines
        integer(wi) k, ii, iiorg, ic
        integer(wi) ij

! ###########################################
! # First nx/2 points where bcs are imposed #
! ###########################################
        do k = 1, nx/2
            iiorg = nx/2 - k + 1

            ii = 1
            ic = iiorg + ii
            dum = cf(ic, k)
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, ii)*dum
            end do
            do ii = 2, k + nx/2
                ic = iiorg + ii
                dum = cf(ic, k)
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + z1(ij, ii)*dum
                end do
            end do

        end do

! #################################################
! # Interior points where boundary does not enter #
! #################################################
        do k = 1 + nx/2, kmax - nx/2
            iiorg = nx/2 - k + 1

            ii = k - nx/2
            ic = ii + iiorg
            dum = cf(ic, k)
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, ii)*dum
            end do
            do ii = k - nx/2 + 1, k + nx/2
                ic = ii + iiorg
                dum = cf(ic, k)
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + z1(ij, ii)*dum
                end do
            end do

        end do

! ##########################################
! # Last nx/2 points where bcs are imposed #
! ##########################################
        do k = kmax - nx/2 + 1, kmax
            iiorg = nx/2 - k + 1

            ii = k - nx/2
            ic = ii + iiorg
            dum = cf(ic, k)
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, ii)*dum
            end do
            do ii = k - nx/2 + 1, kmax
                ic = ii + iiorg
                dum = cf(ic, k)
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + z1(ij, ii)*dum
                end do
            end do

        end do

        return
    end subroutine FLT_T1ND

!########################################################################
! Non-uniform grid
! Particularized for a stencil of 3 points, nx=2
!########################################################################
    subroutine FLT_T1NDD2(kmax, ijmax, cf, z1, zf1)
        integer(wi) kmax, ijmax
        real(wp) cf(3, kmax)
        real(wp) z1(ijmax, *)
        real(wp) zf1(ijmax, *)

! -------------------------------------------------------------------
        integer(wi) ij, k
        real(wp) a1, a2, a3

! ######################################
! # First point where bc are imposed #
! ######################################
        a2 = cf(2, 1)
        a3 = cf(3, 1)
        do ij = 1, ijmax
            zf1(ij, 1) = z1(ij, 1)*a2 + z1(ij, 2)*a3
        end do

! #################################################
! # Interior points where boundary does not enter #
! #################################################
        do k = 2, kmax - 1
            a1 = cf(1, k)
            a2 = cf(2, k)
            a3 = cf(3, k)
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, k - 1)*a1 + z1(ij, k)*a2 + z1(ij, k + 1)*a3
            end do
        end do

! #####################################
! # Last point where bc are imposed #
! #####################################
        a1 = cf(1, kmax)
        a2 = cf(2, kmax)
        do ij = 1, ijmax
            zf1(ij, kmax) = z1(ij, kmax - 1)*a1 + z1(ij, kmax)*a2
        end do

        return
    end subroutine FLT_T1NDD2

!########################################################################
! Non-uniform grid
! Particularized for a stencil of 5 points, nx=4
!########################################################################
    subroutine FLT_T1NDD4(kmax, ijmax, cf, z1, zf1)
        integer(wi) kmax, ijmax
        real(wp) cf(5, kmax)
        real(wp) z1(ijmax, *)
        real(wp) zf1(ijmax, *)

! -------------------------------------------------------------------
        integer(wi) ij, k
        real(wp) a1, a2, a3, a4, a5

! #########################################
! # First 2 points where bc are imposed #
! #########################################
        a3 = cf(3, 1)
        a4 = cf(4, 1)
        a5 = cf(5, 1)
        do ij = 1, ijmax
            zf1(ij, 1) = z1(ij, 1)*a3 + z1(ij, 2)*a4 + z1(ij, 3)*a5
        end do

        a2 = cf(2, 2)
        a3 = cf(3, 2)
        a4 = cf(4, 2)
        a5 = cf(5, 2)
        do ij = 1, ijmax
            zf1(ij, 2) = z1(ij, 1)*a2 + z1(ij, 2)*a3 + z1(ij, 3)*a4 + z1(ij, 4)*a5
        end do

! #################################################
! # Interior points where boundary does not enter #
! #################################################
        do k = 3, kmax - 2
            a1 = cf(1, k)
            a2 = cf(2, k)
            a3 = cf(3, k)
            a4 = cf(4, k)
            a5 = cf(5, k)
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, k - 2)*a1 + z1(ij, k - 1)*a2 + z1(ij, k)*a3 + z1(ij, k + 1)*a4 + z1(ij, k + 2)*a5
            end do
        end do

! ########################################
! # Last 2 points where bc are imposed #
! ########################################
        a1 = cf(1, kmax - 1)
        a2 = cf(2, kmax - 1)
        a3 = cf(3, kmax - 1)
        a4 = cf(4, kmax - 1)
        do ij = 1, ijmax
            zf1(ij, kmax - 1) = z1(ij, kmax - 3)*a1 + z1(ij, kmax - 2)*a2 + z1(ij, kmax - 1)*a3 + z1(ij, kmax)*a4
        end do

        a1 = cf(1, kmax)
        a2 = cf(2, kmax)
        a3 = cf(3, kmax)
        do ij = 1, ijmax
            zf1(ij, kmax) = z1(ij, kmax - 2)*a1 + z1(ij, kmax - 1)*a2 + z1(ij, kmax)*a3
        end do

        return
    end subroutine FLT_T1NDD4

!########################################################################
! Non-uniform grid
! Particularized for a stencil of 7 points, nx=6
!########################################################################
    subroutine FLT_T1NDD6(kmax, ijmax, cf, z1, zf1)
        integer(wi) kmax, ijmax
        real(wp) cf(7, kmax)
        real(wp) z1(ijmax, *)
        real(wp) zf1(ijmax, *)

        integer(wi) ij, k
        real(wp) a1, a2, a3, a4, a5, a6, a7

! #########################################
! # First 3 points where bc are imposed #
! #########################################
        a4 = cf(4, 1)
        a5 = cf(5, 1)
        a6 = cf(6, 1)
        a7 = cf(7, 1)
        do ij = 1, ijmax
            zf1(ij, 1) = z1(ij, 1)*a4 + z1(ij, 2)*a5 + z1(ij, 3)*a6 + z1(ij, 4)*a7
        end do

        a3 = cf(3, 2)
        a4 = cf(4, 2)
        a5 = cf(5, 2)
        a6 = cf(6, 2)
        a7 = cf(7, 2)
        do ij = 1, ijmax
            zf1(ij, 2) = z1(ij, 1)*a3 + z1(ij, 2)*a4 + z1(ij, 3)*a5 + z1(ij, 4)*a6 + z1(ij, 5)*a7
        end do

        a2 = cf(2, 3)
        a3 = cf(3, 3)
        a4 = cf(4, 3)
        a5 = cf(5, 3)
        a6 = cf(6, 3)
        a7 = cf(7, 3)
        do ij = 1, ijmax
            zf1(ij, 3) = z1(ij, 1)*a2 + z1(ij, 2)*a3 + z1(ij, 3)*a4 + z1(ij, 4)*a5 + z1(ij, 5)*a6 + z1(ij, 6)*a7
        end do

! #################################################
! # Interior points where boundary does not enter #
! #################################################
        do k = 4, kmax - 3
            a1 = cf(1, k)
            a2 = cf(2, k)
            a3 = cf(3, k)
            a4 = cf(4, k)
            a5 = cf(5, k)
            a6 = cf(6, k)
            a7 = cf(7, k)
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, k - 3)*a1 + z1(ij, k - 2)*a2 + z1(ij, k - 1)*a3 + &
                             z1(ij, k)*a4 + z1(ij, k + 1)*a5 + z1(ij, k + 2)*a6 + z1(ij, k + 3)*a7
            end do
        end do

! ########################################
! # Last 3 points where bc are imposed #
! ########################################
        a1 = cf(1, kmax - 2)
        a2 = cf(2, kmax - 2)
        a3 = cf(3, kmax - 2)
        a4 = cf(4, kmax - 2)
        a5 = cf(5, kmax - 2)
        a6 = cf(6, kmax - 2)
        do ij = 1, ijmax
            zf1(ij, kmax - 2) = z1(ij, kmax - 5)*a1 + z1(ij, kmax - 4)*a2 + &
                                z1(ij, kmax - 3)*a3 + z1(ij, kmax - 2)*a4 + z1(ij, kmax - 1)*a5 + z1(ij, kmax)*a6
        end do

        a1 = cf(1, kmax - 1)
        a2 = cf(2, kmax - 1)
        a3 = cf(3, kmax - 1)
        a4 = cf(4, kmax - 1)
        a5 = cf(5, kmax - 1)
        do ij = 1, ijmax
            zf1(ij, kmax - 1) = z1(ij, kmax - 4)*a1 + z1(ij, kmax - 3)*a2 + &
                                z1(ij, kmax - 2)*a3 + z1(ij, kmax - 1)*a4 + z1(ij, kmax)*a5
        end do

        a1 = cf(1, kmax)
        a2 = cf(2, kmax)
        a3 = cf(3, kmax)
        a4 = cf(4, kmax)
        do ij = 1, ijmax
            zf1(ij, kmax) = z1(ij, kmax - 3)*a1 + z1(ij, kmax - 2)*a2 + z1(ij, kmax - 1)*a3 + z1(ij, kmax)*a4
        end do

        return
    end subroutine FLT_T1NDD6

!########################################################################
!########################################################################
    ! Periodic boundary conditions

!########################################################################
! Uniform grid
!########################################################################
    subroutine FLT_T1P(kmax, ijmax, nx, z1, zf1)
        integer(wi) kmax, ijmax
        integer(wi) nx
        real(wp) z1(ijmax, *)
        real(wp) zf1(ijmax, *)

        ! -------------------------------------------------------------------
        real(wp) dmul05
        real(wp) dmul1
        integer(wi) k, ij
        integer(wi) ii, im

        dmul1 = 1.0_wp/real(nx, wp)
        dmul05 = 0.5_wp/real(nx, wp)

        ! ##################################################
        ! # First nx/2 points where periodicity is imposed #
        ! ##################################################
        do k = 1, nx/2
            ! first point and last point
            ii = k - nx/2
            ii = ii + kmax - 1
            ii = mod(ii, kmax) + 1
            im = k + nx/2
            do ij = 1, ijmax
                !        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
                zf1(ij, k) = dmul05*(z1(ij, ii) + z1(ij, im))
            end do
            ! inner points
            do ii = k - nx/2 + 1, k + nx/2 - 1
                im = ii + kmax - 1
                im = mod(im, kmax) + 1
                do ij = 1, ijmax
                    !           zf1(ij,k) = dmul05*(zf1(ij,k)+2.0_wp*z1(ij,im))
                    zf1(ij, k) = zf1(ij, k) + dmul1*z1(ij, im)
                end do
            end do
        end do

        ! #################################################
        ! # Interior points where boundary does not enter #
        ! #################################################
        do k = 1 + nx/2, kmax - nx/2
            ! first point and last point
            ii = k - nx/2
            im = k + nx/2
            do ij = 1, ijmax
                !        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
                zf1(ij, k) = dmul05*(z1(ij, ii) + z1(ij, im))
            end do
            ! inner points
            do ii = k - nx/2 + 1, k + nx/2 - 1
                do ij = 1, ijmax
                    !           zf1(ij,k) = dmul05*(zf1(ij,k)+2.0_wp*z1(ij,ii))
                    zf1(ij, k) = zf1(ij, k) + dmul1*z1(ij, ii)
                end do
            end do
        end do

        ! #################################################
        ! # Last nx/2 points where periodicity is imposed #
        ! #################################################
        do k = kmax - nx/2 + 1, kmax
            ! first point and last point
            ii = k - nx/2
            im = k + nx/2
            im = im + kmax - 1
            im = mod(im, kmax) + 1
            do ij = 1, ijmax
                !        zf1(ij,k) = z1(ij,ii)+z1(ij,im)
                zf1(ij, k) = dmul05*(z1(ij, ii) + z1(ij, im))
            end do
            ! inner points
            do ii = k - nx/2 + 1, k + nx/2 - 1
                im = ii + kmax - 1
                im = mod(im, kmax) + 1
                do ij = 1, ijmax
                    !           zf1(ij,k) = dmul05*(zf1(ij,k)+2.0_wp*z1(ij,im))
                    zf1(ij, k) = zf1(ij, k) + dmul1*z1(ij, im)
                end do
            end do
        end do

        return
    end subroutine FLT_T1P

    !########################################################################
    ! Particularized for a stencil of 3 points, nx=2
    !########################################################################
    subroutine FLT_T1PD2(kmax, ijmax, z1, zf1)
        integer(wi) kmax, ijmax
        real(wp) z1(ijmax, *)
        real(wp) zf1(ijmax, *)

        ! -------------------------------------------------------------------
        integer(wi) k, ij

        ! ############################################
        ! # First point where periodicity is imposed #
        ! ############################################
        do ij = 1, ijmax
            zf1(ij, 1) = 0.25_wp*(z1(ij, kmax) + z1(ij, 2)) + 0.5_wp*z1(ij, 1)
        end do

        ! #################################################
        ! # Interior points where boundary does not enter #
        ! #################################################
        do k = 2, kmax - 1
            do ij = 1, ijmax
                zf1(ij, k) = 0.25_wp*(z1(ij, k - 1) + z1(ij, k + 1)) + 0.5_wp*z1(ij, k)
            end do
        end do

        ! ###########################################
        ! # Last point where periodicity is imposed #
        ! ###########################################
        do ij = 1, ijmax
            zf1(ij, kmax) = 0.25_wp*(z1(ij, kmax - 1) + z1(ij, 1)) + 0.5_wp*z1(ij, kmax)
        end do

        return
    end subroutine FLT_T1PD2

    !########################################################################
    ! Particularized for a stencil of 5 points, nx=4
    !########################################################################
    subroutine FLT_T1PD4(kmax, ijmax, z1, zf1)
        integer(wi) kmax, ijmax
        real(wp) z1(ijmax, *)
        real(wp) zf1(ijmax, *)

        ! -------------------------------------------------------------------
        integer(wi) k, ij

        ! ###############################################
        ! # First 2 points where periodicity is imposed #
        ! ###############################################
        do ij = 1, ijmax
            zf1(ij, 1) = 0.125_wp*(z1(ij, kmax - 1) + z1(ij, 3)) + 0.25_wp*(z1(ij, kmax) + z1(ij, 1) + z1(ij, 2))
        end do

        do ij = 1, ijmax
            zf1(ij, 2) = 0.125_wp*(z1(ij, kmax) + z1(ij, 4)) + 0.25_wp*(z1(ij, 1) + z1(ij, 2) + z1(ij, 3))
        end do

        ! #################################################
        ! # Interior points where boundary does not enter #
        ! #################################################
        do k = 3, kmax - 2
            do ij = 1, ijmax
                zf1(ij, k) = 0.125_wp*(z1(ij, k - 2) + z1(ij, k + 2)) + 0.25_wp*(z1(ij, k - 1) + z1(ij, k) + z1(ij, k + 1))
            end do
        end do

        ! ##############################################
        ! # Last 2 points where periodicity is imposed #
        ! ##############################################
        do ij = 1, ijmax
            zf1(ij, kmax - 1) = 0.125_wp*(z1(ij, kmax - 3) + z1(ij, 1)) + &
                                0.25_wp*(z1(ij, kmax - 2) + z1(ij, kmax - 1) + z1(ij, kmax))
        end do

        do ij = 1, ijmax
            zf1(ij, kmax) = 0.125_wp*(z1(ij, kmax - 2) + z1(ij, 2)) + &
                            0.25_wp*(z1(ij, kmax - 1) + z1(ij, kmax) + z1(ij, 1))
        end do

        return
    end subroutine FLT_T1PD4

    !########################################################################
    ! Non-uniform grid
    !########################################################################
    subroutine FLT_T1P_ND(kmax, ijmax, nx, cf, z1, zf1)

        implicit none

        integer(wi) kmax, ijmax
        integer(wi) nx
        real(wp) cf(nx + 1, kmax)
        real(wp) z1(ijmax, *)
        real(wp) zf1(ijmax, *)

        ! -------------------------------------------------------------------
        integer(wi) k, ij, ic
        integer(wi) ii, im

        ! ##################################################
        ! # First nx/2 points where periodicity is imposed #
        ! ##################################################
        do k = 1, nx/2
            ii = k - nx/2
            im = ii + kmax - 1
            im = mod(im, kmax) + 1
            ic = ii - k + nx/2 + 1
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, im)*cf(ic, k)
            end do
            !     DO ii = k-nx/2,k+nx/2
            do ii = k - nx/2 + 1, k + nx/2
                im = ii + kmax - 1
                im = mod(im, kmax) + 1
                ic = ii - k + nx/2 + 1
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + z1(ij, im)*cf(ic, k)
                end do
            end do
        end do

        ! #################################################
        ! # Interior points where boundary does not enter #
        ! #################################################
        do k = 1 + nx/2, kmax - nx/2
            ii = k - nx/2
            ic = ii - k + nx/2 + 1
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, ii)*cf(ic, k)
            end do
            !     DO ii = k-nx/2,k+nx/2
            do ii = k - nx/2 + 1, k + nx/2
                ic = ii - k + nx/2 + 1
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + z1(ij, ii)*cf(ic, k)
                end do
            end do
        end do

        ! #################################################
        ! # Last nx/2 points where periodicity is imposed #
        ! #################################################
        do k = kmax - nx/2 + 1, kmax
            ii = k - nx/2
            im = ii + kmax - 1
            im = mod(im, kmax) + 1
            ic = ii - k + nx/2 + 1
            do ij = 1, ijmax
                zf1(ij, k) = z1(ij, im)*cf(ic, k)
            end do
            !     DO ii = k-nx/2,k+nx/2
            do ii = k - nx/2 + 1, k + nx/2
                im = ii + kmax - 1
                im = mod(im, kmax) + 1
                ic = ii - k + nx/2 + 1
                do ij = 1, ijmax
                    zf1(ij, k) = zf1(ij, k) + z1(ij, im)*cf(ic, k)
                end do
            end do
        end do

        return
    end subroutine FLT_T1P_ND

end module Filters_Tophat
