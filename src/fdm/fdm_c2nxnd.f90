!########################################################################
!#
!# Compact FDMs for non-uniform grids from Shukla and Zhong (2005), JCP, 204, 404–429
!#
!# We set up mmax linear system of size nmax of the form
!#                         A up = B u
!# in this routine. The matrix A is tridiagonal, and B is pentadiagonal,
!# except for the boundary points.
!#
!# We normalize system to get a diagonal 1 in B.
!#
!# 6th-order approximation to 2nd-order derivative:
!#
!# Equations (15) + Table 2 for the interior points (beware of 2 typos in paper, as indicated below).
!# Equations (16) for the first/last points.
!# Table B.2 for the second/second-to-last points.
!#
!# Notation in paper for the interpolation:
!#
!#             +       +             I_n: set of points where the function and derivatives are given
!#   ...---+---+---+---+---+---...
!#         +       +       +         I_m: set of points where only the function is given.
!#
!########################################################################
module FDM_COM_DIRECT
    use TLAB_CONSTANTS
    implicit none
    private

    public FDM_C2N4ND_INITIALIZE
    public FDM_C2N6ND_INITIALIZE
    public FDM_C2NXND_RHS

contains
!########################################################################
! Create diagonals
!########################################################################
    subroutine FDM_C2N6ND_INITIALIZE(nmax, x, lhs, rhs)

        integer(wi), intent(IN) :: nmax
        real(wp), dimension(nmax), intent(IN) :: x
        real(wp), dimension(nmax, 3), intent(OUT) :: lhs ! LHS diagonals (#=3)
        real(wp), dimension(nmax, 4), intent(OUT) :: rhs ! RHS diagonals (#=5-1 because of normalization)

! -------------------------------------------------------------------
        integer(wi) n, idx(2)
        real(wp) am1, a, ap1                            ! Left-hand side
        real(wp) bm2, bm1, b, bp1, bp2, bp3             ! Right-hand side
        real(wp) tmp1, coef(6)
        real(wp) D, dx, dxp, dxm

! #######################################################################
! n = 1
! #######################################################################
        n = 1

! left-hand side
        am1 = 0.0_wp
        a = 1.0_wp
        ap1 = (x(n + 1) - x(n))*(x(n) - x(n + 2)) + (x(n + 1) - x(n))*(x(n) - x(n + 3)) - (x(n) - x(n + 2))*(x(n) - x(n + 3))
        ap1 = ap1/ &
   ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) + (x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

! right-hand side
        bm2 = 0.0_wp
        bm1 = 0.0_wp
        b = (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + 2.0_wp*(x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
        b = b/ &
            ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
        b = b*((x(n) - x(n + 2)) + (x(n) - x(n + 3)))/ &
            (x(n) - x(n + 1)) + 1.0_wp
        b = b/ &
            ((x(n) - x(n + 2))*(x(n) - x(n + 3)))
        tmp1 = ((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))/ &
               ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
        b = (b + tmp1/(x(n + 1) - x(n)))*2.0_wp

        bp1 = (x(n) - x(n + 2)) + (x(n) - x(n + 3)) + ap1*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
        bp1 = bp1*2.0_wp/ &
              ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

        bp2 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 3))) + &
              (x(n) - x(n + 3))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 3)))
        bp2 = bp2/ &
              ((x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) - (x(n + 2) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 3))))
        bp2 = bp2*2.0_wp*(x(n + 1) - x(n))/ &
              ((x(n + 2) - x(n + 1))*(x(n + 2) - x(n))*(x(n + 2) - x(n + 3)))

        bp3 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 2))) + &
              (x(n) - x(n + 2))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 2)))
        bp3 = bp3/ &
              ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) - (x(n + 3) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2))))
        bp3 = bp3*2.0_wp*(x(n + 1) - x(n))/ &
              ((x(n + 3) - x(n + 1))*(x(n + 3) - x(n))*(x(n + 3) - x(n + 2)))

! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
        tmp1 = 1.0_wp/b

        lhs(n, 1) = 0.0_wp
        lhs(n, 2) = a*tmp1
        lhs(n, 3) = ap1*tmp1

        rhs(n, 1) = bp3*tmp1 ! saving here the last term that goes into the 3. superdiagonal
        rhs(n, 2) = 0.0_wp
        rhs(n, 3) = bp1*tmp1
        rhs(n, 4) = bp2*tmp1

! #######################################################################
! n = nmax; same as n = 1, but changing the signs of the increments w.r.t. n
! To understand it, e.g., define a new variable i = -j, where i is the
! discrete variable moving around n
! #######################################################################
        n = nmax

! left-hand side
        am1 = 0.0_wp
        a = 1.0_wp
        ap1 = (x(n - 1) - x(n))*(x(n) - x(n - 2)) + (x(n - 1) - x(n))*(x(n) - x(n - 3)) - (x(n) - x(n - 2))*(x(n) - x(n - 3))
        ap1 = ap1/ &
   ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) + (x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

! right-hand side
        bm2 = 0.0_wp
        bm1 = 0.0_wp
        b = (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + 2.0_wp*(x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
        b = b/ &
            ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
        b = b*((x(n) - x(n - 2)) + (x(n) - x(n - 3)))/ &
            (x(n) - x(n - 1)) + 1.0_wp
        b = b/ &
            ((x(n) - x(n - 2))*(x(n) - x(n - 3)))
        tmp1 = ((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))/ &
               ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
        b = (b + tmp1/(x(n - 1) - x(n)))*2.0_wp

        bp1 = (x(n) - x(n - 2)) + (x(n) - x(n - 3)) + ap1*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
        bp1 = bp1*2.0_wp/ &
              ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

        bp2 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 3))) + &
              (x(n) - x(n - 3))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 3)))
        bp2 = bp2/ &
              ((x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) - (x(n - 2) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 3))))
        bp2 = bp2*2.0_wp*(x(n - 1) - x(n))/ &
              ((x(n - 2) - x(n - 1))*(x(n - 2) - x(n))*(x(n - 2) - x(n - 3)))

        bp3 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 2))) + &
              (x(n) - x(n - 2))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 2)))
        bp3 = bp3/ &
              ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) - (x(n - 3) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2))))
        bp3 = bp3*2.0_wp*(x(n - 1) - x(n))/ &
              ((x(n - 3) - x(n - 1))*(x(n - 3) - x(n))*(x(n - 3) - x(n - 2)))

! if uniform, we should have ( 11 1 0 ) and ( -1 15 -27 13 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
        tmp1 = 1.0_wp/b

        lhs(n, 3) = 0.0_wp
        lhs(n, 2) = a*tmp1
        lhs(n, 1) = ap1*tmp1

        rhs(n, 4) = bp3*tmp1 ! saving here the last term that goes into the 3. subdiagonal
        rhs(n, 3) = 0.0_wp
        rhs(n, 2) = bp1*tmp1
        rhs(n, 1) = bp2*tmp1

! #######################################################################
! second and second-to-last points
! #######################################################################
        idx = [2, nmax - 1]

        do n = 1, size(idx)
            coef = coef_c2n4(x, idx(n))

! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
            tmp1 = 1.0_wp/coef(5)

            lhs(idx(n), 1) = coef(1)*tmp1
            lhs(idx(n), 2) = coef(2)*tmp1
            lhs(idx(n), 3) = coef(3)*tmp1

            rhs(idx(n), 1) = 0.0_wp
            rhs(idx(n), 2) = coef(4)*tmp1
            rhs(idx(n), 3) = coef(6)*tmp1
            rhs(idx(n), 4) = 0.0_wp

        end do

! #######################################################################
! 2 < n < nmax-1
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

            b = 2.0_wp*C2D(x, n, n)/D + 2.0_wp*C1D(x, n, n)/D*((dxm + dxp)/(dxp*dxm) + lp(x, n, n)) &
                + (2.0_wp + 2.0_wp*lp(x, n, n)*(dxm + dxp))/(dxp*dxm) + lpp(x, n, n)
            bp2 = c2n6_coef(x, n + 2, n)
            bm2 = c2n6_coef(x, n - 2, n)

! -------------------------------------------------------------------
! if uniform, we should have ( 2/11 1 2/11 ) and ( 3/44 12/11 -51/22 12/11 3/44 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
            ! print *, am1, a, ap1
            ! print *, bm2*dxp**2., bm1*dxp**2., b*dxp**2., bp1*dxp**2., bp2*dxp**2.

            tmp1 = 1.0_wp/b

            lhs(n, 1) = am1*tmp1
            lhs(n, 2) = a*tmp1
            lhs(n, 3) = ap1*tmp1

            rhs(n, 1) = bm2*tmp1
            rhs(n, 2) = bm1*tmp1
            rhs(n, 3) = bp1*tmp1
            rhs(n, 4) = bp2*tmp1

        end do

        return
    end subroutine FDM_C2N6ND_INITIALIZE

!########################################################################
! First line in Table B.2
!########################################################################
    function coef_c2n4(x, i) result(coef)           ! Interval around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i
        real(wp) coef(6)
        real(wp) am1, a, ap1, bm1, b, bp1           ! for clarity below
        real(wp) D, dx, dxm, dxp

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
!########################################################################
    function a2n6_coef(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f, f1, f2, dx, dxm, dxp

        dx = x(ip) - x(im)
        dxp = x(i) - x(ip)
        dxm = x(i) - x(im)

        f1 = B1D(x, ip, im, i)*(dxm + dxp) + B2D(x, im, ip, i)*dxp*(dxp + 2.0_wp*dxm)
        f1 = f1*2.0_wp*PIp(x(:), i, i)
        f2 = B1D(x, ip, im, i) + B2D(x, im, ip, i)*dxp
        f2 = f2*PIpp(x(:), i, i)*dxp*dxm
        f = -(f1 + f2)/dx/PI(x, ip, i)

        return
    end function

    function b2n6_coef(x, im, ip, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: ip, im, i
        real(wp) f, f1, f2, dxm, dxp, D

        real(wp) dx

        dx = x(ip) - x(im)
        dxp = x(i) - x(ip)
        dxm = x(i) - x(im)
        D = D_coef(x, i)

        f1 = 1.0_wp + A1D(x, im, ip, i)/D*(dxm + dxp) + A2D(x, im, ip, i)/D*dxp*(dxp + 2.0_wp*dxm)
        f1 = f1*2.0_wp*PIp(x(:), i, i)
        f2 = 1.0_wp + A1D(x, im, ip, i)/D*dxp + A2D(x, im, ip, i)/D*dxp**2.0_wp
        f2 = f2*PIpp(x(:), i, i)*dxm
        f = (f1 + f2)/dx/PI(x, ip, i)

        return
    end function

    function c2n6_coef(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: j, i
        real(wp) f, f1, f2, dxm, dxp, dxm2, dxp2, D

        real(wp) dx

        dx = x(i) - x(j)
        dxp = x(i) - x(i + 1)
        dxm = x(i) - x(i - 1)
        dxp2 = x(j) - x(i + 1)
        dxm2 = x(j) - x(i - 1)
        D = D_coef(x, i)

        f1 = C1D(x, j, i)/D*(1.0_wp + dx/dxp + dx/dxm) &
             + C2D(x, j, i)/D*(2.0_wp + dx/dxp + dx/dxm)*dx + 1.0_wp/dxp + 1.0_wp/dxm
        f1 = f1*2.0_wp*lp(x, j, i)*dxp*dxm/(dxp2*dxm2) !dxp*dxp/(dxp2*dxp2) this was a mistake in the paper
        f2 = 1.0_wp + C1D(x, j, i)/D*dx + C2D(x, j, i)/D*dx**2.0_wp
        f2 = f2*lpp(x, j, i)*dxp*dxm/(dxp2*dxm2)
        f = f1 + f2

        return
    end function

    function PI(x, j, i) result(f)        ! Interval m around i, i.e, {i-2,i,i+2}
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        f = (x(j) - x(i))*(x(j) - x(i + 2))*(x(j) - x(i - 2))

        return
    end function

    function PIp_o_PI(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        f = (x(j) - x(i + 2))*(x(j) - x(i - 2)) &
            + (x(j) - x(i))*(x(j) - x(i - 2)) &
            + (x(j) - x(i))*(x(j) - x(i + 2))
        f = f/PI(x, j, i)

        return
    end function

    function PIpp_o_PI(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        f = x(j) - x(i + 2) + x(j) - x(i - 2) + x(j) - x(i)
        f = 2.0_wp*f/PI(x, j, i)

        return
    end function

    function PIp(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        f = (x(j) - x(i + 2))*(x(j) - x(i - 2)) &
            + (x(j) - x(i))*(x(j) - x(i - 2)) &
            + (x(j) - x(i))*(x(j) - x(i + 2))

        return
    end function

    function PIpp(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        f = 2.0_wp*(x(j) - x(i + 2) + x(j) - x(i - 2) + x(j) - x(i))

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

    function lp(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        real(wp) dx, dxp, dxm

        dx = x(i + 2) - x(i - 2)
        dxp = x(i + 2) - x(i)
        dxm = x(i) - x(i - 2)

        select case (j - i)
        case (-2)
            f = -dxp/(dxm*dx)

        case (0)
            f = (dxp - dxm)/(dxp*dxm)

        case (+2)
            f = dxm/(dxp*dx)

        end select

        return
    end function

    function lpp(x, j, i) result(f)        ! Interval m around i
        real(wp), intent(in) :: x(:)
        integer(wi), intent(in) :: i, j
        real(wp) f

        real(wp) dx, dxp, dxm

        dx = x(i + 2) - x(i - 2)
        dxp = x(i + 2) - x(i)
        dxm = x(i) - x(i - 2)

        select case (j - i)
        case (-2)
            f = 2.0_wp/(dxm*dx)

        case (0)
            f = -2.0_wp/(dxp*dxm)

        case (+2)
            f = 2.0_wp/(dxp*dx)

        end select

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

        real(wp) dx, dxp, dxm

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

        real(wp) dx, dxp, dxm

        dx = x(i + 1) - x(i - 1)
        dxp = x(i + 1) - x(j)
        dxm = x(j) - x(i - 1)

        f = 2.0_wp*(1.0_wp/dxp**2.0 + 1.0_wp/dxm**2.0 - 1.0_wp/(dxp*dxm)) &
            + 2.0_wp*dx**2.0/(dxp*dxm)*PIp_o_PI(x(:), i + 1, i)*PIp_o_PI(x(:), i - 1, i) &
            - 2.0_wp*PIp_o_PI(x(:), i + 1, i)*dx/dxm*(1.0_wp/dxp - 1.0_wp/dxm) &
            - 2.0_wp*PIp_o_PI(x(:), i - 1, i)*dx/dxp*(1.0_wp/dxp - 1.0_wp/dxm)

        return
    end function

! #######################################################################
! Constructing forcing term
! #######################################################################
    subroutine FDM_C2NXND_RHS(nmax, mmax, rhs, u, d)

        integer(wi), intent(in) :: nmax, mmax ! m linear systems or size n
        real(wp), intent(in) :: rhs(nmax, 4)       ! RHS diagonals (#=5-1 because of normalization)
        real(wp), intent(in) :: u(mmax, nmax)         ! function
        real(wp), intent(out) :: d(mmax, nmax)         ! RHS

! -------------------------------------------------------------------
        integer(wi) n

! #######################################################################
! Boundary
        n = 1 ! rhs(1,1) contains 3. superdiagonal
        d(:, n) = &
            +u(:, n) &
            + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4) + u(:, n + 3)*rhs(n, 1)

        n = 2
        d(:, n) = u(:, n - 1)*rhs(n, 2) &
                  + u(:, n) &
                  + u(:, n + 1)*rhs(n, 3)

! Interior points
        do n = 3, nmax - 2
            d(:, n) = u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
                      + u(:, n) &
                      + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4)
        end do

! Boundary
        n = nmax - 1
        d(:, n) = u(:, n - 1)*rhs(n, 2) &
                  + u(:, n) &
                  + u(:, n + 1)*rhs(n, 3)

        n = nmax ! rhs(1,4) contains 3. subdiagonal
        d(:, n) = u(:, n - 3)*rhs(n, 4) + u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
                  + u(:, n)

        return
    end subroutine FDM_C2NXND_RHS

!########################################################################
!# 4th-order approximation to 2nd-order derivative:
!########################################################################
    subroutine FDM_C2N4ND_INITIALIZE(nmax, x, lhs, rhs)

        integer(wi), intent(in) :: nmax
        real(wp), intent(in) :: x(nmax)
        real(wp), intent(out) :: lhs(nmax, 3) ! LHS diagonals (#=3)
        real(wp), intent(out) :: rhs(nmax, 4) ! RHS diagonals (#=5-1 because of normalization)

! -------------------------------------------------------------------
        integer(wi) n
        real(wp) am1, a, ap1                            ! Left-hand side
        real(wp) bm2, bm1, b, bp1, bp2, bp3             ! Right-hand side
        real(wp) tmp1                                   ! Intermediate ops
        real(wp) coef(6)

! #######################################################################
! n = 1
! #######################################################################
        n = 1

! left-hand side
        am1 = 0.0_wp
        a = 1.0_wp
        ap1 = (x(n + 1) - x(n))*(x(n) - x(n + 2)) + (x(n + 1) - x(n))*(x(n) - x(n + 3)) - (x(n) - x(n + 2))*(x(n) - x(n + 3))
        ap1 = ap1/ &
   ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) + (x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

! right-hand side
        bm2 = 0.0_wp
        bm1 = 0.0_wp
        b = (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + 2.0_wp*(x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
        b = b/ &
            ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
        b = b*((x(n) - x(n + 2)) + (x(n) - x(n + 3)))/ &
            (x(n) - x(n + 1)) + 1.0_wp
        b = b/ &
            ((x(n) - x(n + 2))*(x(n) - x(n + 3)))
        tmp1 = ((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))/ &
               ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
        b = (b + tmp1/(x(n + 1) - x(n)))*2.0_wp

        bp1 = (x(n) - x(n + 2)) + (x(n) - x(n + 3)) + ap1*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
        bp1 = bp1*2.0_wp/ &
              ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

        bp2 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 3))) + &
              (x(n) - x(n + 3))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 3)))
        bp2 = bp2/ &
              ((x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) - (x(n + 2) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 3))))
        bp2 = bp2*2.0_wp*(x(n + 1) - x(n))/ &
              ((x(n + 2) - x(n + 1))*(x(n + 2) - x(n))*(x(n + 2) - x(n + 3)))

        bp3 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 2))) + &
              (x(n) - x(n + 2))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 2)))
        bp3 = bp3/ &
              ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) - (x(n + 3) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2))))
        bp3 = bp3*2.0_wp*(x(n + 1) - x(n))/ &
              ((x(n + 3) - x(n + 1))*(x(n + 3) - x(n))*(x(n + 3) - x(n + 2)))

! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
        tmp1 = 1.0_wp/b

        lhs(n, 1) = 0.0_wp
        lhs(n, 2) = a*tmp1
        lhs(n, 3) = ap1*tmp1

        rhs(n, 1) = bp3*tmp1 ! saving here the last term that goes into the 3. superdiagonal
        rhs(n, 2) = 0.0_wp
        rhs(n, 3) = bp1*tmp1
        rhs(n, 4) = bp2*tmp1

! #######################################################################
! n = 2
! #######################################################################
        do n = 2, nmax - 1
            coef = coef_c2n4(x, n)

            ! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2
            ! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
            ! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
            tmp1 = 1.0_wp/coef(5)

            lhs(n, 1) = coef(1)*tmp1
            lhs(n, 2) = coef(2)*tmp1
            lhs(n, 3) = coef(3)*tmp1

            rhs(n, 1) = 0.0_wp
            rhs(n, 2) = coef(4)*tmp1
            rhs(n, 3) = coef(6)*tmp1
            rhs(n, 4) = 0.0_wp

        end do
        
! #######################################################################
! n = nmax; same as n = 1, but changing the signs of the increments w.r.t. n
! To understand it, e.g., define a new variable i = -j, where i is the
! discrete variable moving around n
! #######################################################################
        n = nmax

! left-hand side
        am1 = 0.0_wp
        a = 1.0_wp
        ap1 = (x(n - 1) - x(n))*(x(n) - x(n - 2)) + (x(n - 1) - x(n))*(x(n) - x(n - 3)) - (x(n) - x(n - 2))*(x(n) - x(n - 3))
        ap1 = ap1/ &
   ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) + (x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

! right-hand side
        bm2 = 0.0_wp
        bm1 = 0.0_wp
        b = (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + 2.0_wp*(x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
        b = b/ &
            ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
        b = b*((x(n) - x(n - 2)) + (x(n) - x(n - 3)))/ &
            (x(n) - x(n - 1)) + 1.0_wp
        b = b/ &
            ((x(n) - x(n - 2))*(x(n) - x(n - 3)))
        tmp1 = ((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))/ &
               ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
        b = (b + tmp1/(x(n - 1) - x(n)))*2.0_wp

        bp1 = (x(n) - x(n - 2)) + (x(n) - x(n - 3)) + ap1*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
        bp1 = bp1*2.0_wp/ &
              ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

        bp2 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 3))) + &
              (x(n) - x(n - 3))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 3)))
        bp2 = bp2/ &
              ((x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) - (x(n - 2) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 3))))
        bp2 = bp2*2.0_wp*(x(n - 1) - x(n))/ &
              ((x(n - 2) - x(n - 1))*(x(n - 2) - x(n))*(x(n - 2) - x(n - 3)))

        bp3 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 2))) + &
              (x(n) - x(n - 2))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 2)))
        bp3 = bp3/ &
              ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) - (x(n - 3) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2))))
        bp3 = bp3*2.0_wp*(x(n - 1) - x(n))/ &
              ((x(n - 3) - x(n - 1))*(x(n - 3) - x(n))*(x(n - 3) - x(n - 2)))

! if uniform, we should have ( 11 1 0 ) and ( -1 15 -27 13 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
        tmp1 = 1.0_wp/b

        lhs(n, 3) = 0.0_wp
        lhs(n, 2) = a*tmp1
        lhs(n, 1) = ap1*tmp1

        rhs(n, 4) = bp3*tmp1 ! saving here the last term that goes into the 3. subdiagonal
        rhs(n, 3) = 0.0_wp
        rhs(n, 2) = bp1*tmp1
        rhs(n, 1) = bp2*tmp1

        return
    end subroutine FDM_C2N4ND_INITIALIZE

end module FDM_COM_DIRECT
