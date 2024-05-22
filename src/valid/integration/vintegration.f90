program vintegration
    use TLAB_CONSTANTS, only: wp, wi
    use Integration

    integer(wi), parameter :: nmax = 768
    integer(wi) n
    real(wp) x(nmax), f(1, nmax), I_a(1,nmax), I_n(1, nmax)

    open (21, file='y.dat')
    do n = 1, nmax
        read (21, *) x(n)
    end do
    close (21)

! ###################################################################
! Define the function f and analytic solution
    do n = 1, nmax
! linear
        f(:, n) = x(n)
        I_a(:,n) = 0.5*(x(n)**2.0_wp - x(1)**2.0_wp)

! step

! sinusoidal

! exponential
        
    end do

    call Int_Trapezoidal_v(f, x, I_n(:, nmax))
    call check(I_a, I_n, 'integral1.dat')

    call Int_Simpson_v(f, x, I_n(:, nmax))
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
        do i = 1, size(I_a, 2)
            do l = 1, size(I_a, 1)
                if (present(name)) then
                    write (20, 1000) x(i), I_a(l, i), I_n(l, i), I_a(l, i) - I_n(l, i)
                end if
                dummy = dummy + I_a(l, i)*I_a(l, i)
                error_l2 = error_l2 + (I_a(l, i) - I_n(l, i))**2.0_wp
                error_max = max(error_max, abs(I_a(l, i) - I_n(l, i)))
            end do
        end do
        if (present(name)) then
            close (20)
        end if

        if (dummy == 0.0_wp) return
        ! write (*, *) 'Relative Error L2-norm .....:', sqrt(g%jac(1, 1)*error_l2)/maxval(abs(u))
        write (*, *) 'Relative Error Linf-norm ...:', error_max/maxval(abs(I_a))

        return
1000    format(5(1x, e12.5))
    end subroutine check

end program vintegration
