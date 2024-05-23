module Integration
    use TLAB_CONSTANTS, only: wp, wi, BCS_MIN, BCS_MAX !, pi_wp
    implicit none
    private

    public :: Int_Trapezoidal, Int_Trapezoidal_v                            ! the last one is a vector version
    public :: Int_Trapezoidal_f                                             ! antiderivatives (function of x)
    public :: Int_Trapezoidal_Increments_InPlace
    public :: Int_Simpson, Int_Simpson_v, Int_Simpson_Biased_v
    public :: Int_Simpson_Biased_f
    public :: Int_Simpson_Biased_Increments, Int_Simpson_Biased_Increments_InPlace

contains
!########################################################################
!########################################################################
    function Int_Trapezoidal(f, x) result(integral)
        real(wp), intent(in) :: f(:)
        real(wp), intent(in) :: x(:)
        real(wp) integral

! -------------------------------------------------------------------
        integer(wi) n

! ###################################################################
        integral = 0.0
        do n = 2, size(x)
            integral = integral + 0.5_wp*(f(n) + f(n - 1))*(x(n) - x(n - 1))
        end do

        return
    end function Int_Trapezoidal

!########################################################################
!########################################################################
    subroutine Int_Trapezoidal_v(f, x, integral)
        ! function Int_Trapezoidal(f, x) result(integral(:))
        real(wp), intent(in) :: f(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(out) :: integral(:)

! -------------------------------------------------------------------
        integer(wi) n

! ###################################################################
        integral = 0.0
        do n = 2, size(x)
            integral(:) = integral(:) + 0.5_wp*(f(:, n) + f(:, n - 1))*(x(n) - x(n - 1))
        end do

        return
        ! end function Int_Trapezoidal
    end subroutine Int_Trapezoidal_v

!########################################################################
!########################################################################
    ! Second-order integral, antiderivative (function of x)
    subroutine Int_Trapezoidal_f(f, x, result, ibc)
        real(wp), intent(in) :: f(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout) :: result(:, :)   ! contains bcs
        integer, intent(in), optional :: ibc

! -------------------------------------------------------------------
        integer(wi) j, nmax
        integer ibc_loc

! ###################################################################
        nmax = size(x)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)              ! Backwards stencil. Result is I_n = int_{x_{1}}^{x_n} u dx
            do j = 2, nmax
                result(:, j) = result(:, j - 1) + 0.5_wp*(f(:, j) + f(:, j - 1))*(x(j) - x(j - 1))
            end do

        case (BCS_MAX)              ! Forward stencil. Result is I_n = int_{x_{n}}^{x_{nmax}} u dx
            do j = nmax - 1, 1, -1
                result(:, j) = result(:, j + 1) - 0.5_wp*(f(:, j) + f(:, j + 1))*(x(j + 1) - x(j))
            end do

        end select

    end subroutine Int_Trapezoidal_f

!########################################################################
!########################################################################
    subroutine Int_Trapezoidal_Increments_InPlace(f, x, ibc)
        real(wp), intent(inout) :: f(:, :)
        real(wp), intent(in) :: x(:)
        integer, intent(in), optional :: ibc

! -------------------------------------------------------------------
        integer(wi) j, nmax
        integer ibc_loc

! ###################################################################
        nmax = size(x)

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)              ! Backwards stencil. Result is I_n = int_{x_{n-1}}^{x_n} u dx
            do j = nmax, 2, -1      ! Starting from the top, f(:,n) is no longer used with the backwards stencil
                f(:, j) = 0.5_wp*(f(:, j) + f(:, j - 1))*(x(j) - x(j - 1))
            end do

        case (BCS_MAX)              ! Forward stencil. Result is I_n = int_{x_{n}}^{x_{n+1}} u dx
            do j = 1, nmax - 1      ! Starting from the bottom, f(:,n) is no longer used with the forward stencil
                f(:, j) = 0.5_wp*(f(:, j) + f(:, j + 1))*(x(j + 1) - x(j))
            end do

        end select

    end subroutine Int_Trapezoidal_Increments_InPlace

! ###################################################################
! ###################################################################
    function Int_Simpson(u, x) result(integral)
        real(wp), intent(IN) :: u(:), x(:)
        real(wp) integral

        ! -------------------------------------------------------------------
        integer(wi) i, nn, nmax
        real(wp) dx21, dx20, dx10, du20, du10, du21, b, c
        real(wp) c13

! ###################################################################
        nmax = size(x)

        c13 = 1.0_wp/3.0_wp

! Correct the last element contribution
        if (mod(nmax, 2) == 0) then
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
! to be updated to allow for the correction to be at the last or the first element
! see biased routine below
    subroutine Int_Simpson_v(u, x, result)
        real(wp), intent(IN) :: u(:, :)
        real(wp), intent(IN) :: x(:)
        real(wp), intent(OUT) :: result(:)

        ! -------------------------------------------------------------------
        integer(wi) n, nmax
        real(wp) a, b, c, dxm2, dxm1, dxp1, c16

! ###################################################################
        nmax = size(x)

        if (nmax == 2) then
            result = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
            return
        end if

        c16 = 1.0_wp/6.0_wp

        result = 0.0_wp

        do n = 2, nmax - 1, 2
            dxm1 = x(n) - x(n - 1)
            dxp1 = x(n + 1) - x(n)
            a = (2.0_wp - dxp1/dxm1)*(dxm1 + dxp1)*c16
            b = (dxm1 + dxp1)**2.0_wp/(dxm1*dxp1)*(dxm1 + dxp1)*c16
            c = (2.0_wp - dxm1/dxp1)*(dxm1 + dxp1)*c16

            result = result + a*u(:, n - 1) + b*u(:, n) + c*u(:, n + 1)

        end do

        if (mod(nmax, 2) == 0) then ! Correct the last element contribution
            n = nmax
            ! this might lead to oscillations near sharp gradients
            dxm1 = x(n) - x(n - 1)
            dxm2 = x(n - 1) - x(n - 2)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            result = result + a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            ! fall back to first order, monotone
            ! result = result + 0.5_wp*(u(:, n) + u(:, n - 1))*(x(n) - x(n - 1))

        end if

        return
    end subroutine Int_Simpson_v

! ###################################################################
! ###################################################################
    ! Calculate integral step by step instead every 2 steps using biased, 2 point formula
    ! Bias direction is set by ibc=[BCS_MIN,BCS_MAX] parameter
    ! ..... x x ..... calculate integral in this interval between to consecutive grid points
    ! ... + + + ..... using this stencil (BCS_MIN)
    ! ..... + + + ... or this stencil (BCS_MAX)
    subroutine Int_Simpson_Biased_v(u, x, result, ibc)
        real(wp), intent(IN) :: u(:, :)
        real(wp), intent(IN) :: x(:)
        real(wp), intent(OUT) :: result(:)
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n, nmax
        real(wp) a, b, c, dxm2, dxm1, c16
        integer ibc_loc

! ###################################################################
        nmax = size(x)

        if (nmax == 2) then
            result = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
            return
        end if

        c16 = 1.0_wp/6.0_wp

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)              ! Backwards stencil
            n = 1                   ! At the lower boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n + 1) - x(n)
            dxm2 = x(n + 2) - x(n + 1)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            result = a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2)

            do n = 3, nmax
                dxm1 = x(n) - x(n - 1)
                dxm2 = x(n - 1) - x(n - 2)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                result = result + a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            end do

        case (BCS_MAX)              ! Forward stencil
            n = nmax                ! At the upper boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n) - x(n - 1)
            dxm2 = x(n - 1) - x(n - 2)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            result = a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            do n = nmax - 2, 1, -1
                dxm1 = x(n + 1) - x(n)
                dxm2 = x(n + 2) - x(n + 1)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                result = result + a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2)

            end do

        end select

        return
    end subroutine Int_Simpson_Biased_v

! ###################################################################
! ###################################################################
    ! Calculate integral step by step instead every 2 steps using biased, 2 point formula
    ! Bias direction is set by ibc=[BCS_MIN,BCS_MAX] parameter
    ! ..... x x ..... calculate integral in this interval between to consecutive grid points
    ! ... + + + ..... using this stencil (BCS_MIN). Result is I_n = int_{x_1}^{x_{n}}
    ! ..... + + + ... or this stencil (BCS_MAX). Result is I_n = int_{x_n}^{x_{nmax}}
    subroutine Int_Simpson_Biased_f(u, x, result, ibc)
        real(wp), intent(IN) :: u(:, :)
        real(wp), intent(IN) :: x(:)
        real(wp), intent(OUT) :: result(:, :)
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n, nmax
        real(wp) a, b, c, dxm2, dxm1, c16
        integer ibc_loc

! ###################################################################
        nmax = size(x)

        c16 = 1.0_wp/6.0_wp

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)              ! Backwards stencil. Result is I_n = int_{x_{1}}^{x_n} u dx
            if (nmax == 2) then
                result(:, 2) = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
                return
            end if

            n = 1                   ! At the lower boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n + 1) - x(n)
            dxm2 = x(n + 2) - x(n + 1)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            result(:, 2) = a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2)

            do n = 3, nmax
                dxm1 = x(n) - x(n - 1)
                dxm2 = x(n - 1) - x(n - 2)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                result(:, n) = result(:, n - 1) + a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            end do

        case (BCS_MAX)              ! Forward stencil. Result is I_n = int_{x_n}^{x_{nmax}}
            if (nmax == 2) then
                result(:, 1) = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
                return
            end if

            n = nmax                ! At the upper boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n) - x(n - 1)
            dxm2 = x(n - 1) - x(n - 2)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            result(:, nmax - 1) = a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            do n = nmax - 2, 1, -1
                dxm1 = x(n + 1) - x(n)
                dxm2 = x(n + 2) - x(n + 1)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                result(:, n) = result(:, n + 1) + a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2)

            end do

        end select

        return
    end subroutine Int_Simpson_Biased_f

! ###################################################################
! ###################################################################
    ! Calculate integral step by step instead every 2 steps using biased, 2 point formula
    ! Bias direction is set by ibc=[BCS_MIN,BCS_MAX] parameter
    ! ..... x x ..... calculate integral in this interval between to consecutive grid points
    ! ... + + + ..... using this stencil (BCS_MIN). Result is I_n = int_{x_{n-1}}^{x_n} u dx
    ! ..... + + + ... or this stencil (BCS_MAX). Result is I_n = int_{x_{n}}^{x_{n+1}} u dx
    subroutine Int_Simpson_Biased_Increments(u, x, result, ibc)
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(out) :: result(:, :)
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n, nmax
        real(wp) a, b, c, dxm2, dxm1, c16
        integer ibc_loc

! ###################################################################
        nmax = size(x)

        c16 = 1.0_wp/6.0_wp

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)              ! Backwards stencil. Result is I_n = int_{x_{n-1}}^{x_n} u dx
            if (nmax == 2) then
                result(:, nmax) = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
                return
            end if

            n = 1                   ! At the lower boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n + 1) - x(n)
            dxm2 = x(n + 2) - x(n + 1)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            result(:, n + 1) = a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2)
            do n = nmax, 3, -1      ! Starting from the top, u(:,n) is no longer used with the backwards stencil
                dxm1 = x(n) - x(n - 1)
                dxm2 = x(n - 1) - x(n - 2)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                result(:, n) = a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            end do

        case (BCS_MAX)              ! Forward stencil. Result is I_n = int_{x_{n}}^{x_{n+1}} u dx
            if (nmax == 2) then
                result(:, 1) = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
                return
            end if

            n = nmax                ! At the upper boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n) - x(n - 1)
            dxm2 = x(n - 1) - x(n - 2)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            result(:, n - 1) = a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            do n = 1, nmax - 2      ! Starting from the bottom, u(:,n) is no longer used with the forward stencil
                dxm1 = x(n + 1) - x(n)
                dxm2 = x(n + 2) - x(n + 1)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                result(:, n) = a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2)

            end do

        end select

        return
    end subroutine Int_Simpson_Biased_Increments

! ###################################################################
! ###################################################################
    subroutine Int_Simpson_Biased_Increments_Inplace(u, x, wrk2d, ibc)
        real(wp), intent(inout) :: u(:, :)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout) :: wrk2d(:)
        integer, intent(in), optional :: ibc

        ! -------------------------------------------------------------------
        integer(wi) n, nmax
        real(wp) a, b, c, dxm2, dxm1, c16
        integer ibc_loc

! ###################################################################
        nmax = size(x)

        c16 = 1.0_wp/6.0_wp

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_MIN
        end if

        select case (ibc_loc)
        case (BCS_MIN)              ! Backwards stencil. Result is I_n = int_{x_{n-1}}^{x_n} u dx
            if (nmax == 2) then
                u(:, nmax) = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
                return
            end if

            n = 1                   ! At the lower boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n + 1) - x(n)
            dxm2 = x(n + 2) - x(n + 1)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            wrk2d(1:size(u, 1)) = a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2) ! save for later

            do n = nmax, 3, -1      ! Starting from the top, u(:,n) is no longer used with the backwards stencil
                dxm1 = x(n) - x(n - 1)
                dxm2 = x(n - 1) - x(n - 2)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                u(:, n) = a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2)

            end do

            u(:, 2) = wrk2d(1:size(u, 1))

        case (BCS_MAX)              ! Forward stencil. Result is I_n = int_{x_{n}}^{x_{n+1}} u dx
            if (nmax == 2) then
                u(:, 1) = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
                return
            end if

            n = nmax                ! At the upper boundary, I need to reverse the stencil and use the other biased formula
            dxm1 = x(n) - x(n - 1)
            dxm2 = x(n - 1) - x(n - 2)
            a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
            b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
            c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

            wrk2d(1:size(u, 1)) = a*u(:, n) + b*u(:, n - 1) - c*u(:, n - 2) ! save for later

            do n = 1, nmax - 2      ! Starting from the bottom, u(:,n) is no longer used with the forward stencil
                dxm1 = x(n + 1) - x(n)
                dxm2 = x(n + 2) - x(n + 1)
                a = c16*(2.0_wp*dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/(dxm2 + dxm1)
                b = c16*(dxm1*dxm1 + 3.0_wp*dxm1*dxm2)/dxm2
                c = c16*dxm1*dxm1*dxm1/dxm2/(dxm2 + dxm1)

                u(:, n) = a*u(:, n) + b*u(:, n + 1) - c*u(:, n + 2)

            end do

            u(:, nmax - 1) = wrk2d(1:size(u, 1))

        end select

        return
    end subroutine Int_Simpson_Biased_Increments_Inplace

! ! ###################################################################
! ! ###################################################################
!     subroutine Int_Simpson_v_old(u, x, result, wrk2d)
!         real(wp), intent(IN) :: u(:, :)
!         real(wp), intent(IN) :: x(:)
!         real(wp), intent(OUT) :: result(:)
!         real(wp), intent(inout) :: wrk2d(:, :)

!         target wrk2d

!         ! -------------------------------------------------------------------
!         integer(wi) i, nn, nmax
!         real(wp) dx21, dx20, dx10
!         real(wp), dimension(:), pointer :: du20, du10, du21, b, c
!         real(wp) c13

! ! ###################################################################
!         nmax = size(u, 2)

!         du20 => wrk2d(:, 1)
!         du10 => wrk2d(:, 2)
!         du21 => wrk2d(:, 3)
!         b => wrk2d(:, 4)
!         c => wrk2d(:, 5)

!         if (nmax == 2) then
!             result = 0.5_wp*(u(:, 1) + u(:, 2))*(x(2) - x(1))
!             return
!         end if

!         c13 = 1.0_wp/3.0_wp

! ! Correct the last element contribution
!         if (MOD(nmax, 2) == 0) then
!             dx21 = x(nmax) - x(nmax - 1)
!             dx20 = x(nmax) - x(nmax - 2)
!             dx10 = x(nmax - 1) - x(nmax - 2)
!             du20 = u(:, nmax) - u(:, nmax - 2)
!             du10 = u(:, nmax - 1) - u(:, nmax - 2)
!             du21 = u(:, nmax) - u(:, nmax - 1)
!             c = (du21/dx21 - du10/dx10)/dx20
!             b = (du21/dx21 - c*dx21)*0.5_wp
!             result = dx21*(u(:, nmax - 1) + dx21*(b + c*dx21*c13))
!             nn = nmax - 1
!         else
!             result(:) = 0.0_wp
!             nn = nmax
!         end if

!         do i = 2, nn - 1, 2
!             dx21 = x(i + 1) - x(i)
!             dx20 = x(i + 1) - x(i - 1)
!             dx10 = x(i) - x(i - 1)
!             du20 = u(:, i + 1) - u(:, i - 1)
!             du10 = u(:, i) - u(:, i - 1)
!             c = (du20/dx20 - du10/dx10)/dx21
!             b = (du20/dx20 - c*dx20)*0.5_wp
!             result = result + dx20*(u(:, i - 1) + dx20*(b + c*dx20*c13))
!         end do

!         return
!     end subroutine Int_Simpson_v_old

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
