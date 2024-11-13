program vintegration
    use TLab_Constants, only: wp, wi, BCS_MAX, pi_wp
    use Integration
    implicit none

    integer(wi), parameter :: nmax = 768
    integer(wi) l, n
    real(wp) x(nmax), scale, f(1, nmax), I_a(1, nmax), I_n(1, nmax)
    real(wp) kappa, parameter1, xref

    open (21, file='y.dat')
    do n = 1, nmax
        read (21, *) x(n)
        scale = x(nmax) - x(1)
    end do
    close (21)

    write (*, *) 'Parameter ?'; read (*, *) parameter1

! ###################################################################
! Define the function f and analytic solution
    do n = 1, nmax
! linear
        ! f(:, n) = x(n)
        ! I_a(:, n) = 0.5*(x(n)**2.0_wp - x(1)**2.0_wp)

! sinusoidal
        ! kappa = 2.0_wp*pi_wp/scale*parameter1
        ! f(:, n) = sin(kappa*x(n))
        ! I_a(:, n) = (cos(kappa*x(1)) - cos(kappa*x(n)))/kappa

! exponential
        ! f(:, n) = exp(parameter1*x(n))
        ! I_a(:, n) = (f(:, n) - f(:, 1))/parameter1

! step
        xref = 0.5_wp*(x(nmax) - x(1))
        f(:, n) = sign(1.0_wp, x(n) - xref)
        I_a(:, n) = max(x(1) - x(n), x(n) - x(nmax))

    end do

    I_n(:, 1) = 0.0_wp
    do n = 2, nmax
        call Int_Trapezoidal_v(f(:, 1:n), x(1:n), I_n(:, n))
    end do
    call check(I_a, I_n, 'integral1.dat')

    I_n(:, 1) = 0.0_wp
    do n = 2, nmax
        call Int_Simpson_Biased_v(f(:, 1:n), x(1:n), I_n(:, n))
        ! call Int_Simpson_v(f(:, 1:n), x(1:n), I_n(:, n))
    end do
    ! I_n(:, nmax) = 0.0_wp
    ! do n = nmax - 1, 1, -1
    !     call Int_Simpson_Biased_v(f(:, n:nmax), x(n:nmax), I_n(:, n), BCS_MAX)
    !     ! call Int_Simpson_v(f(:, n:nmax), x(n:nmax), I_n(:, n))
    ! end do
    ! do n = nmax, 1, -1
    !     I_n(:, n) = I_n(:, 1) - I_n(:, n) 
    ! end do
    call check(I_a, I_n, 'integral2.dat')

    stop

    ! ###################################################################
contains
    subroutine check(I_a, I_n, name)
        real(wp), intent(in) :: I_a(:, :), I_n(:, :)
        character(len=*), optional :: name

        real(wp) dummy, error_l2, error_max

        if (present(name)) then
            open (20, file=name)
        end if
        error_l2 = 0.0_wp
        error_max = 0.0_wp
        dummy = 0.0_wp
        do n = 1, size(I_a, 2)
            do l = 1, size(I_a, 1)
                if (present(name)) then
                    write (20, 1000) x(n), I_a(l, n), I_n(l, n), abs(I_a(l, n) - I_n(l, n))
                end if
                dummy = dummy + I_a(l, n)*I_a(l, n)
                error_l2 = error_l2 + (I_a(l, n) - I_n(l, n))**2.0_wp
                error_max = max(error_max, abs(I_a(l, n) - I_n(l, n)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Absolute Error Linf-norm ...:', error_max
        if (dummy == 0.0_wp) return
        ! write (*, *) 'Relative Error L2-norm .....:', sqrt(g%jac(1, 1)*error_l2)/maxval(abs(u))
        write (*, *) 'Relative Error Linf-norm ...:', error_max/maxval(abs(I_a))

        return
1000    format(5(1x, e12.5))
    end subroutine check

end program vintegration
