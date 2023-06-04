!########################################################################
!# Compact FDMs for non-uniform grids from Shukla and Zhong (2005), JCP, 204, 404–429
!#
!# We set up mmax linear system of size nmax of the form
!#                         A up = B u
!# in this routine. The matrix A is tridiagonal, and B is at most pentadiagonal.
!#
!# We normalize system to get a diagonal 1 in B and we only need to store 4 RHS diagonals
!# The RHS diagonals are then O(1), the LHS diagonal O(h^2)
!########################################################################
module FDM_COM_DIRECT
    use TLAB_CONSTANTS
    implicit none
    private

    public FDM_C2N4ND_INITIALIZE
    public FDM_C2N6ND_INITIALIZE
    public FDM_C2NXND_RHS

contains
    ! #######################################################################
    ! Calcualte f = B u, assuming B is penta-diagonal with center diagonal is 1
    ! #######################################################################
    subroutine FDM_C2NXND_RHS(nmax, mmax, rhs, u, f)
        integer(wi), intent(in) :: nmax, mmax       ! m linear systems or size n
        real(wp), intent(in) :: rhs(nmax, 4)        ! RHS diagonals (#=5-1 because of normalization)
        real(wp), intent(in) :: u(mmax, nmax)       ! function u
        real(wp), intent(out) :: f(mmax, nmax)      ! RHS, f = B u
        
        ! -------------------------------------------------------------------
        integer(wi) n

        ! -------------------------------------------------------------------
        ! Boundary
        n = 1 ! rhs(1,1) contains 3. superdiagonal
        f(:, n) = &
            +u(:, n) &
            + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4) + u(:, n + 3)*rhs(n, 1)

        n = 2
        f(:, n) = u(:, n - 1)*rhs(n, 2) &
                  + u(:, n) &
                  + u(:, n + 1)*rhs(n, 3)

        ! Interior points
        do n = 3, nmax - 2
            f(:, n) = u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
                      + u(:, n) &
                      + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4)
        end do

        ! Boundary
        n = nmax - 1
        f(:, n) = u(:, n - 1)*rhs(n, 2) &
                  + u(:, n) &
                  + u(:, n + 1)*rhs(n, 3)

        n = nmax ! rhs(1,4) contains 3. subdiagonal
        f(:, n) = u(:, n - 3)*rhs(n, 4) + u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
                  + u(:, n)

        return
    end subroutine FDM_C2NXND_RHS

    !########################################################################
    !# 6th-order approximation to 2nd-order derivative:
    !########################################################################
    subroutine FDM_C2N6ND_INITIALIZE(nmax, x, lhs, rhs)
        integer(wi), intent(in) :: nmax
        real(wp), intent(in) :: x(nmax)
        real(wp), intent(out) :: lhs(nmax, 3) ! LHS diagonals (#=3)
        real(wp), intent(out) :: rhs(nmax, 4) ! RHS diagonals (#=5-1 because of normalization)

        ! -------------------------------------------------------------------
        real(wp) am1, a, ap1                ! Left-hand side; for clarity below
        real(wp) bm2, bm1, b, bp1, bp2      ! Right-hand side
        real(wp) dx, dxp, dxm               ! basic increments
        real(wp) D                          ! discriminant of linear systems
        real(wp) dummy
        real(wp) coef(6)
        integer(wi) n, idx(2)

        ! #######################################################################
        ! Equations (16) for the first/last points.
        ! Notation in paper for the interpolation:
        !
        !       +                    I_n: set of points where the function and derivatives are given
        !   +---+---+---+---...
        !   +       +   +            I_m: set of points where only the function is given.
        ! #######################################################################
        n = 1

        coef = coef_c2n3_biased(x, n)               ! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2

        dummy = 1.0_wp/coef(3)                      ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals

        lhs(n, 1) = 0.0_wp                          ! not used
        lhs(n, 2:3) = coef(1:2)*dummy               ! a, ap1

        rhs(n, 2) = 0.0_wp                          ! not used
        rhs(n, [3, 4, 1]) = coef([4, 5, 6])*dummy   ! bp2, bp3, bp4; bp4 is saved into rhs(1)

        ! -------------------------------------------------------------------
        n = nmax

        coef = coef_c2n3_biased(x, n, backwards=.true.) ! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2

        dummy = 1.0_wp/coef(3)                          ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals

        lhs(n, 3) = 0.0_wp                              ! not used
        lhs(n, [2, 1]) = coef(1:2)*dummy                ! am1, a

        rhs(n, 3) = 0.0_wp                              ! not used
        rhs(n, [2, 1, 4]) = coef([4, 5, 6])*dummy       ! bm2, bm3, bm4; bm4 is saved into rhs(4)

        ! #######################################################################
        ! Table B.2 for the second/second-to-last points.
        ! Notation in paper for the interpolation:
        !
        !         +       +         I_n: set of points where the function and derivatives are given
        !   ...---+---+---+---...
        !             +             I_m: set of points where only the function is given.
        ! #######################################################################
        idx = [2, nmax - 1]

        do n = 1, size(idx)
            coef = coef_c2n4(x, idx(n))                 ! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2

            dummy = 1.0_wp/coef(5)                      ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals

            lhs(idx(n), 1:3) = coef(1:3)*dummy          ! am1, a, ap1

            rhs(idx(n), [1, 4]) = 0.0_wp                ! not used
            rhs(idx(n), 2:3) = coef([4, 6])*dummy       ! bm1, b, bp1

        end do

        ! #######################################################################
        ! Equations (15) + Table 2 for the interior points (beware of 2 typos in paper, as indicated below).
        !
        ! Notation in paper for the interpolation:
        !
        !             +       +             I_n: set of points where the function and derivatives are given
        !   ...---+---+---+---+---+---...
        !         +       +       +         I_m: set of points where only the function is given.
        ! #######################################################################
        do n = 3, nmax - 2

            D = D_coef(x, n)

! left-hand side
            a = 1.0_wp
            ap1 = a2n6_coef(x, n - 1, n + 1, n)/D
            am1 = a2n6_coef(x, n + 1, n - 1, n)/D

! right-hand side
            bp1 = b2n6_coef(x, n - 1, n + 1, n)
            bm1 = b2n6_coef(x, n + 1, n - 1, n)

            dx = x(n + 1) - x(n - 1)
            dxp = x(n) - x(n + 1)
            dxm = x(n) - x(n - 1)

            b = 2.0_wp*C2D(x, n, n)/D + 2.0_wp*C1D(x, n, n)/D*((dxm + dxp)/(dxp*dxm) + Lag_p(x, n, n, [n - 2, n, n + 2])) &
                + (2.0_wp + 2.0_wp*Lag_p(x, n, n, [n - 2, n, n + 2])*(dxm + dxp))/(dxp*dxm) + Lag_pp_3(x, n, n, [n - 2, n, n + 2])
            bp2 = c2n6_coef(x, n + 2, n)
            bm2 = c2n6_coef(x, n - 2, n)

            ! if uniform, we should have ( 2/11 1 2/11 ) and ( 3/44 12/11 -51/22 12/11 3/44 )/h^2
            ! print*, [am1, a, ap1]
            ! print*, [bm2, bm1, b, bp1, bp2]*(x(2)-x(1))**2.0

            ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
            dummy = 1.0_wp/b
            lhs(n, 1:3) = [am1, a, ap1]*dummy
            rhs(n, 1:4) = [bm2, bm1, bp1, bp2]*dummy

        end do

        return
    end subroutine FDM_C2N6ND_INITIALIZE

    !########################################################################
    !# 4th-order approximation to 2nd-order derivative; extracted from above
    !########################################################################
    subroutine FDM_C2N4ND_INITIALIZE(nmax, x, lhs, rhs)
        integer(wi), intent(in) :: nmax
        real(wp), intent(in) :: x(nmax)
        real(wp), intent(out) :: lhs(nmax, 3) ! LHS diagonals (#=3)
        real(wp), intent(out) :: rhs(nmax, 4) ! RHS diagonals (#=5-1 because of normalization)

        real(wp) dummy
        real(wp) coef(6)
        integer(wi) n

        ! #######################################################################
        ! Equations (16) for the first/last points.
        ! #######################################################################
        n = 1

        coef = coef_c2n3_biased(x, n)               ! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2

        dummy = 1.0_wp/coef(3)                      ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals

        lhs(n, 1) = 0.0_wp                          ! not used
        lhs(n, 2:3) = coef(1:2)*dummy               ! a, ap1

        rhs(n, 2) = 0.0_wp                          ! not used
        rhs(n, [3, 4, 1]) = coef([4, 5, 6])*dummy   ! bp1, bp2, bp3; bp3 is saved into rhs(1)

        ! -------------------------------------------------------------------
        n = nmax

        coef = coef_c2n3_biased(x, n, backwards=.true.) ! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2

        dummy = 1.0_wp/coef(3)                          ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals

        lhs(n, 3) = 0.0_wp                              ! not used
        lhs(n, [2, 1]) = coef(1:2)*dummy                ! am1, a

        rhs(n, 3) = 0.0_wp                              ! not used
        rhs(n, [2, 1, 4]) = coef([4, 5, 6])*dummy       ! am1, am2, am3; am3 is saved into rhs(4)

        ! #######################################################################
        ! Table B.2 for the interior points.
        ! #######################################################################
        do n = 2, nmax - 1
            coef = coef_c2n4(x, n)                 ! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2

            dummy = 1.0_wp/coef(5)                 ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals

            lhs(n, 1:3) = coef(1:3)*dummy          ! am1, a, ap1

            rhs(n, [1, 4]) = 0.0_wp                ! not used
            rhs(n, 2:3) = coef([4, 6])*dummy       ! bm1, b, bp1

        end do

        return
    end subroutine FDM_C2N4ND_INITIALIZE

    !########################################################################
    ! First line in Table B.2
    !########################################################################
    function coef_c2n4(x, i) result(coef)           ! Interval around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        real(wp) coef(6)

        real(wp) am1, a, ap1                ! Left-hand side; for clarity below
        real(wp) bm1, b, bp1                ! Right-hand side
        real(wp) dx, dxp, dxm               ! basic increments
        real(wp) D                          ! discriminant of linear systems

        dx = x(i + 1) - x(i - 1)
        dxp = x(i + 1) - x(i)
        dxm = x(i) - x(i - 1)

        D = dxp*dxm + dx**2.0_wp

        am1 = (dxm**2.0_wp - dxp**2.0_wp + dxp*dxm)*dxp/dx/D
        a = 1.0_wp
        ap1 = (dxp**2.0_wp - dxm**2.0_wp + dxp*dxm)*dxm/dx/D

        bm1 = dxp/dx*12.0_wp/D
        b = -12.0_wp/D
        bp1 = dxm/dx*12.0_wp/D

        coef = [am1, a, ap1, bm1, b, bp1]

        return
    end function

    !########################################################################
    ! Equations (16)
    !########################################################################
    function coef_c2n3_biased(x, i, backwards) result(coef)
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        logical, intent(in), optional :: backwards
        real(wp) coef(6)

        real(wp) a1, a2, b1, b2, b3, b4
        real(wp) dx1, dx3, dx4
        real(wp) D
        integer(wi) set_m(3), i1, i2, i3, i4

        i1 = i
        if (present(backwards)) then
            ! same as fowards, but changing the signs of the increments w.r.t. i
            ! To understand it, e.g., define a new variable k = -j, where k is the
            ! discrete variable moving around i
            i2 = i - 1
            i3 = i - 2
            i4 = i - 3
        else
            i2 = i + 1
            i3 = i + 2
            i4 = i + 3
        end if
        dx1 = x(i2) - x(i1)
        dx3 = x(i2) - x(i3)
        dx4 = x(i2) - x(i4)
        set_m = [i1, i3, i4]

        ! -------------------------------------------------------------------
        a1 = 1.0_wp

        a2 = (0.5_wp*dx1*Pi_pp_3(x, i1, set_m) - Pi_p(x, i1, set_m))/Pi_p(x, i2, set_m)

        b2 = Pi_pp_3(x, i1, set_m) &
             + 0.5_wp*dx1*Pi_pp_3(x, i1, set_m)*Pi_pp_3(x, i2, set_m)/Pi_p(x, i2, set_m) &
             - Pi_p(x, i1, set_m)/Pi_p(x, i2, set_m)*Pi_pp_3(x, i2, set_m)
        b2 = b2/Pi(x, i2, set_m)

        ! -------------------------------------------------------------------
        D = Lag(x, i2, i1, set_m) + dx1*Lag_p(x, i2, i1, set_m)
        b1 = -2.0_wp*Lag_p(x, i1, i1, set_m)*(Lag(x, i2, i1, set_m) + 2.0_wp*dx1*Lag_p(x, i2, i1, set_m)) &
             + 2.0_wp*Lag_p(x, i2, i1, set_m)
        b1 = b1/D/dx1 + Lag_pp_3(x, i1, i1, set_m)

        D = Lag(x, i2, i3, set_m) + dx3*Lag_p(x, i2, i3, set_m)
        b3 = (Lag(x, i2, i3, set_m) + dx1*Lag_p(x, i2, i3, set_m))*dx1*Lag_pp_3(x, i1, i3, set_m) &
             - 2.0_wp*(Lag(x, i2, i3, set_m) + 2.0_wp*dx1*Lag_p(x, i2, i3, set_m))*Lag_p(x, i1, i3, set_m)
        b3 = b3/D/dx3

        D = Lag(x, i2, i4, set_m) + dx4*Lag_p(x, i2, i4, set_m)
        b4 = (Lag(x, i2, i4, set_m) + dx1*Lag_p(x, i2, i4, set_m))*dx1*Lag_pp_3(x, i1, i4, set_m) &
             - 2.0_wp*(Lag(x, i2, i4, set_m) + 2.0_wp*dx1*Lag_p(x, i2, i4, set_m))*Lag_p(x, i1, i4, set_m)
        b4 = b4/D/dx4

        coef = [a1, a2, b1, b2, b3, b4]

        return
    end function

    !########################################################################
    ! Equations (15) + Table 2 for the interior points (beware of 2 typos in paper, as indicated below).
    !########################################################################
    function a2n6_coef(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f

        real(wp) dx, dxp, dxm               ! basic increments
        real(wp) f1, f2

        dx = x(ip) - x(im)
        dxp = x(i) - x(ip)
        dxm = x(i) - x(im)

        f1 = B1D(x, ip, im, i)*(dxm + dxp) + B2D(x, im, ip, i)*dxp*(dxp + 2.0_wp*dxm)
        f1 = f1*2.0_wp*PI_p(x(:), i, [i - 2, i, i + 2])
        f2 = B1D(x, ip, im, i) + B2D(x, im, ip, i)*dxp
        f2 = f2*PI_pp_3(x(:), i, [i - 2, i, i + 2])*dxp*dxm
        f = -(f1 + f2)/dx/Pi(x, ip, [i - 2, i, i + 2])

        return
    end function

    function b2n6_coef(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f

        real(wp) dx, dxp, dxm               ! basic increments
        real(wp) D                          ! discriminant of linear systems
        real(wp) f1, f2

        dx = x(ip) - x(im)
        dxp = x(i) - x(ip)
        dxm = x(i) - x(im)
        D = D_coef(x, i)

        f1 = 1.0_wp + A1D(x, im, ip, i)/D*(dxm + dxp) + A2D(x, im, ip, i)/D*dxp*(dxp + 2.0_wp*dxm)
        f1 = f1*2.0_wp*PI_p(x(:), i, [i - 2, i, i + 2])
        f2 = 1.0_wp + A1D(x, im, ip, i)/D*dxp + A2D(x, im, ip, i)/D*dxp**2.0_wp
        f2 = f2*PI_pp_3(x(:), i, [i - 2, i, i + 2])*dxm
        f = (f1 + f2)/dx/Pi(x, ip, [i - 2, i, i + 2])

        return
    end function

    function c2n6_coef(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i
        real(wp) f

        real(wp) dx, dxp, dxm               ! basic increments
        real(wp) D                          ! discriminant of linear systems
        real(wp) f1, f2, dxm2, dxp2

        dx = x(i) - x(j)
        dxp = x(i) - x(i + 1)
        dxm = x(i) - x(i - 1)
        dxp2 = x(j) - x(i + 1)
        dxm2 = x(j) - x(i - 1)
        D = D_coef(x, i)

        f1 = C1D(x, j, i)/D*(1.0_wp + dx/dxp + dx/dxm) &
             + C2D(x, j, i)/D*(2.0_wp + dx/dxp + dx/dxm)*dx + 1.0_wp/dxp + 1.0_wp/dxm
        f1 = f1*2.0_wp*Lag_p(x, i, j, [i - 2, i, i + 2])*dxp*dxm/(dxp2*dxm2) !dxp*dxp/(dxp2*dxp2) this was a mistake in the paper
        f2 = 1.0_wp + C1D(x, j, i)/D*dx + C2D(x, j, i)/D*dx**2.0_wp
        f2 = f2*Lag_pp_3(x, j, j, [i - 2, i, i + 2])*dxp*dxm/(dxp2*dxm2)
        f = f1 + f2

        return
    end function

    function PIp_o_PI(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        f = (x(j) - x(i + 2))*(x(j) - x(i - 2)) &
            + (x(j) - x(i))*(x(j) - x(i - 2)) &
            + (x(j) - x(i))*(x(j) - x(i + 2))
        f = f/Pi(x, j, [i - 2, i, i + 2])

        return
    end function

    function PIpp_o_PI(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        f = x(j) - x(i + 2) + x(j) - x(i - 2) + x(j) - x(i)
        f = 2.0_wp*f/Pi(x, j, [i - 2, i, i + 2]) 

        return
    end function

    function D_coef(x, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        real(wp) f

        real(wp) dx

        dx = x(i + 1) - x(i - 1)

        f = 6.0_wp + 4.0_wp*dx*(PIp_o_PI(x(:), i + 1, i) - PIp_o_PI(x(:), i - 1, i)) &
            - 2.0_wp*dx**2.0_wp*PIp_o_PI(x(:), i + 1, i)*PIp_o_PI(x(:), i - 1, i)

        return
    end function

    function A1D(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f

        real(wp) dx

        dx = x(ip) - x(im)

        f = -4.0_wp*PIp_o_PI(x(:), ip, i) &
            - 2.0_wp*PIp_o_PI(x(:), im, i) &
            + 2.0_wp*dx*(PIp_o_PI(x(:), ip, i)*PIp_o_PI(x(:), im, i) - PIpp_o_PI(x(:), ip, i)) &
            + dx**2.0_wp*PIpp_o_PI(x(:), ip, i)*PIp_o_PI(x(:), im, i)

        return
    end function

    function A2D(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f

        real(wp) dx

        dx = x(ip) - x(im)

        f = 4.0_wp*PIp_o_PI(x(:), ip, i)*PIp_o_PI(x(:), im, i) &
            - PIpp_o_PI(x(:), ip, i) &
            - 2.0_wp/dx*(PIp_o_PI(x(:), ip, i) - PIp_o_PI(x(:), im, i)) &
            + dx*PIpp_o_PI(x(:), ip, i)*PIp_o_PI(x(:), im, i)

        return
    end function

    function B1D(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f

        real(wp) dx

        dx = x(ip) - x(im)

        f = -(2.0_wp/dx + PIp_o_PI(x(:), ip, i))*dx**2.0_wp

        return
    end function

    function B2D(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f

        real(wp) dx

        dx = x(ip) - x(im)

        f = 1.0_wp - dx*PIp_o_PI(x(:), im, i)

        return
    end function

    function C1D(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i
        real(wp) f

        real(wp) dx, dxp, dxm               ! basic increments

        dx = x(i + 1) - x(i - 1)
        dxp = x(i + 1) - x(j)
        dxm = x(j) - x(i - 1)

!        f = (dxp - dxm)/(dxp*dxm)*(10.0_wp - 4.0_wp*dx**2.0_wp/(dxp*dxm))  this was a mistake in the paper
        f = (dxp - dxm)/(dxp*dxm)*(6.0_wp - 4.0_wp*dx**2.0_wp/(dxp*dxm)) &
            + 2.0_wp*dx*(dxm/dxp - dxp/dxm)*PIp_o_PI(x(:), i - 1, i)*PIp_o_PI(x(:), i + 1, i) &
            + PIp_o_PI(x(:), i - 1, i)*(4.0_wp*dx/dxp - 4.0_wp*dx/dxm - 2.0_wp*dx**2.0/dxp**2.0) &
            - PIp_o_PI(x(:), i + 1, i)*(4.0_wp*dx/dxp - 4.0_wp*dx/dxm + 2.0_wp*dx**2.0/dxm**2.0)

        return
    end function

    function C2D(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i
        real(wp) f

        real(wp) dx, dxp, dxm               ! basic increments

        dx = x(i + 1) - x(i - 1)
        dxp = x(i + 1) - x(j)
        dxm = x(j) - x(i - 1)

        f = 2.0_wp*(1.0_wp/dxp**2.0 + 1.0_wp/dxm**2.0 - 1.0_wp/(dxp*dxm)) &
            + 2.0_wp*dx**2.0/(dxp*dxm)*PIp_o_PI(x(:), i + 1, i)*PIp_o_PI(x(:), i - 1, i) &
            - 2.0_wp*PIp_o_PI(x(:), i + 1, i)*dx/dxm*(1.0_wp/dxp - 1.0_wp/dxm) &
            - 2.0_wp*PIp_o_PI(x(:), i - 1, i)*dx/dxp*(1.0_wp/dxp - 1.0_wp/dxm)

        return
    end function

    !########################################################################
    ! Product functions and Lagrange polynomials over interval given by idx(:)
    !########################################################################
    function Pi(x, j, idx) result(f)    ! Product function defined over interval given by idx(:)
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            f = f*(x(j) - x(idx(k)))
        end do

        return
    end function

    function Pi_p(x, j, idx) result(f)
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        real(wp) dummy
        integer(wi) k, m

        f = 0.0_wp
        do k = 1, size(idx)
            dummy = 1.0_wp
            do m = 1, size(idx)
                if (m /= k) then
                    dummy = dummy*(x(j) - x(idx(m)))
                end if
            end do
            f = f + dummy
        end do

        return
    end function

    function Pi_pp_3(x, j, idx) result(f)       ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, idx(:)
        real(wp) f

        f = 2.0_wp*(x(j) - x(idx(1)) + x(j) - x(idx(2)) + x(j) - x(idx(3)))

        return
    end function

    function Lag(x, j, i, idx) result(f)        ! Lagrange polynomials on idx(:) around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k

        f = 1.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f*(x(j) - x(idx(k)))/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

    function Lag_p(x, j, i, idx) result(f)        ! 1. derivative of Lagrange polynomials on idx(:) around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j, idx(:)
        real(wp) f

        integer(wi) k, m
        real(wp) den, dummy

        den = 1.0_wp
        f = 0.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                dummy = 1.0_wp
                do m = 1, size(idx)
                    if (idx(m) /= i .and. m /= k) then
                        dummy = dummy*(x(j) - x(idx(m)))
                    end if
                end do
                f = f + dummy
                den = den*(x(i) - x(idx(k)))
            end if
        end do
        f = f/den

        return
    end function

    function Lag_pp_3(x, j, i, idx) result(f)    ! It assumes idx has only 3 points
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i, idx(:)
        real(wp) f

        integer(wi) k

        f = 2.0_wp
        do k = 1, size(idx)
            if (idx(k) /= i) then
                f = f/(x(i) - x(idx(k)))
            end if
        end do

        return
    end function

end module FDM_COM_DIRECT
