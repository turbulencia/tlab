module Integration
    use TLAB_CONSTANTS, only: wp, wi, BCS_MIN, BCS_MAX !, pi_wp
    implicit none
    private

    public :: Int_Trapezoidal, Int_Trapezoidal_InPlace
    public :: Int_Simpson, Int_Simpson_v

contains
!########################################################################
!########################################################################
    ! Second-order integral
    subroutine Int_Trapezoidal(f, x, result, ibc)
        real(wp), intent(in) :: f(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout) :: result(:, :)   ! contains bcs
        integer, intent(in), optional :: ibc

! -------------------------------------------------------------------
        integer(wi) j, nmax
        integer ibc_loc

! ###################################################################
        nmax = size(f, 2)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)
            do j = 2, nmax
                result(:, j) = result(:, j - 1) + 0.5_wp*(f(:, j) + f(:, j - 1))*(x(j) - x(j - 1))
            end do

        case (BCS_MAX)
            do j = nmax - 1, 1, -1
                result(:, j) = result(:, j + 1) - 0.5_wp*(f(:, j) + f(:, j + 1))*(x(j + 1) - x(j))
            end do

        end select

    end subroutine Int_Trapezoidal

!########################################################################
!########################################################################
! Second-order integral; assumes zero boundary condition
    subroutine Int_Trapezoidal_InPlace(f, x, ibc)
        real(wp), intent(inout) :: f(:, :)
        real(wp), intent(in) :: x(:)
        integer, intent(in), optional :: ibc

! -------------------------------------------------------------------
        integer(wi) j, nmax
        integer ibc_loc

! ###################################################################
        nmax = size(f, 2)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)
            j = 1
            f(:, j) = 0.5_wp*(f(:, j) + f(:, j + 1))*(x(j + 1) - x(j))
            do j = 2, nmax - 1
                f(:, j) = f(:, j - 1) + 0.5_wp*(f(:, j) + f(:, j + 1))*(x(j + 1) - x(j))
            end do
            do j = nmax, 2, -1
                f(:, j) = f(:, j - 1)
            end do
            f(:, 1) = 0.0_wp

        case (BCS_MAX)
            j = nmax
            f(:, j) = -0.5_wp*(f(:, j) + f(:, j - 1))*(x(j) - (x(j - 1)))
            do j = nmax - 1, 2, -1
                f(:, j) = f(:, j + 1) - 0.5_wp*(f(:, j) + f(:, j - 1))*(x(j) - x(j - 1))
            end do
            do j = 1, nmax - 1
                f(:, j) = f(:, j + 1)
            end do
            f(:, nmax) = 0.0_wp

        end select

    end subroutine Int_Trapezoidal_InPlace

! ###################################################################
! ###################################################################
    function Int_Simpson(nmax, u, x) result(integral)
        integer(wi), intent(IN) :: nmax
        real(wp), intent(IN) :: u(nmax), x(nmax)
        real(wp) integral

        ! -------------------------------------------------------------------
        integer(wi) i, nn
        real(wp) dx21, dx20, dx10, du20, du10, du21, b, c
        real(wp) c13

! ###################################################################
        c13 = 1.0_wp/3.0_wp

! Correct the last element contribution
        if (MOD(nmax, 2) == 0) then
            dx21 = x(nmax) - x(nmax - 1)
            dx20 = x(nmax) - x(nmax - 2)
            dx10 = x(nmax - 1) - x(nmax - 2)
            du20 = u(nmax) - u(nmax - 2)
            du10 = u(nmax - 1) - u(nmax - 2)
            du21 = u(nmax) - u(nmax - 1)
            c = (du21/dx21 - du10/dx10)/dx20
            b = (du21/dx21 - c*dx21)*0.5_wp
            integral = dx21*(u(nmax - 1) + dx21*(b + c*dx21*c13))
            nn = nmax - 1
        else
            integral = 0.0_wp
            nn = nmax
        end if

        do i = 2, nn - 1, 2
            dx21 = x(i + 1) - x(i)
            dx20 = x(i + 1) - x(i - 1)
            dx10 = x(i) - x(i - 1)
            du20 = u(i + 1) - u(i - 1)
            du10 = u(i) - u(i - 1)
            c = (du20/dx20 - du10/dx10)/dx21
            b = (du20/dx20 - c*dx20)*0.5_wp
            integral = integral + dx20*(u(i - 1) + dx20*(b + c*dx20*c13))
        end do

        return
    end function Int_Simpson

! ###################################################################
! ###################################################################
    subroutine Int_Simpson_v(u, x, result, wrk2d)
        real(wp), intent(IN) :: u(:, :)
        real(wp), intent(IN) :: x(:)
        real(wp), intent(OUT) :: result(:)
        real(wp), intent(inout) :: wrk2d(:, :)

        target wrk2d

        ! -------------------------------------------------------------------
        integer(wi) i, nn, nmax
        real(wp) dx21, dx20, dx10
        real(wp), dimension(:), pointer :: du20, du10, du21, b, c
        real(wp) c13

! ###################################################################
        nmax = size(u, 2)

        du20 => wrk2d(:, 1)
        du10 => wrk2d(:, 2)
        du21 => wrk2d(:, 3)
        b => wrk2d(:, 4)
        c => wrk2d(:, 5)

        if (nmax == 2) then
            result = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
            return
        end if

        c13 = 1.0_wp/3.0_wp

! Correct the last element contribution
        if (MOD(nmax, 2) == 0) then
            dx21 = x(nmax) - x(nmax - 1)
            dx20 = x(nmax) - x(nmax - 2)
            dx10 = x(nmax - 1) - x(nmax - 2)
            du20 = u(:, nmax) - u(:, nmax - 2)
            du10 = u(:, nmax - 1) - u(:, nmax - 2)
            du21 = u(:, nmax) - u(:, nmax - 1)
            c = (du21/dx21 - du10/dx10)/dx20
            b = (du21/dx21 - c*dx21)*0.5_wp
            result = dx21*(u(:, nmax - 1) + dx21*(b + c*dx21*c13))
            nn = nmax - 1
        else
            result(:) = 0.0_wp
            nn = nmax
        end if

        do i = 2, nn - 1, 2
            dx21 = x(i + 1) - x(i)
            dx20 = x(i + 1) - x(i - 1)
            dx10 = x(i) - x(i - 1)
            du20 = u(:, i + 1) - u(:, i - 1)
            du10 = u(:, i) - u(:, i - 1)
            c = (du20/dx20 - du10/dx10)/dx21
            b = (du20/dx20 - c*dx20)*0.5_wp
            result = result + dx20*(u(:, i - 1) + dx20*(b + c*dx20*c13))
        end do

        return
    end subroutine Int_Simpson_v

!########################################################################
!########################################################################
    ! function SIMPSON_D(nmax, u, dx) result(integral)
    !     integer(wi), intent(IN) :: nmax
    !     real(wp), intent(IN) :: u(nmax), dx(nmax)
    !     real(wp) integral

    !     ! -------------------------------------------------------------------
    !     integer(wi) i, nn
    !     real(wp) slast
    !     real(wp) dx21, dx20, dx10, du20, du10, du21, b, c

    !     !########################################################################
    !     if (MOD(nmax, 2) == 0) then
    !         dx21 = dx(nmax)
    !         dx20 = dx(nmax) + dx(nmax - 1)
    !         dx10 = dx(nmax - 1)
    !         du20 = u(nmax) - u(nmax - 2)
    !         du10 = u(nmax - 1) - u(nmax - 2)
    !         du21 = u(nmax) - u(nmax - 1)
    !         c = (du21/dx21 - du10/dx10)/dx20
    !         b = (du21/dx21 - c*dx21)/2.0_wp
    !         slast = dx21*(u(nmax - 1) + dx21*(b + c*dx21/3.0_wp))
    !         nn = nmax - 1
    !     else
    !         slast = 0.0_wp
    !         nn = nmax
    !     end if

    !     integral = u(1)*dx(1) + u(nn)*dx(nn)

    !     do i = 2, nn - 1, 2
    !         integral = integral + 4.0_wp*u(i)*dx(i)
    !     end do

    !     do i = 3, nn - 2, 2
    !         integral = integral + 2.0_wp*u(i)*dx(i)
    !     end do

    !     integral = integral/3.0_wp + slast

    !     return
    ! end function SIMPSON

end module Integration
