module GRID_LOCAL
    use TLAB_TYPES, only: cp, ci
    implicit none
    save

    integer(ci), parameter :: GTYPE_UNIFORM = 0
    integer(ci), parameter :: GTYPE_TANH = 5
    integer(ci), parameter :: GTYPE_EXP = 6

    integer(ci), parameter :: MAX_PARAMES = 9
    integer(ci), parameter :: MAX_OPTIONS = 5
    integer(ci), parameter :: MAX_SEGMENT = 10

    type grid_build_dt                              ! Information to construct grid in one direction
        sequence
        integer(ci) nseg                            ! number of segments in this direction
        logical mirrored                            ! It true, mirror the grid
        real(cp) fixed_scale                        ! If positive, rescale grid to this value
        integer(ci) SIZE(MAX_SEGMENT)               ! number of points in each segment
        real(cp) end(MAX_SEGMENT)                   ! physical end of each segment
        integer(ci) opts(MAX_OPTIONS, MAX_SEGMENT)  ! 0   uniform segment
        ! 1   Colonius, Lele and Moin stretching
        ! 2   Second order polynomial stretching
        ! 3   Third order polynomial stretching
        ! 4   Geometric progression
        ! 5   Hyperbolic tangent
        ! 6   Exponential
        real(cp) vals(MAX_PARAMES, MAX_SEGMENT)
    end type grid_build_dt

    type(grid_build_dt) g_build(3)

contains
    !########################################################################
    !# The grid spacing dx/ds is given by a hyperbolic tangent function
    !#
    !# dx/ds = f_0 + [(f_1-f_0)/2]*[ 1 + TANH[(s-s_1)/(2*delta_1)] ]
    !#             + [(f_2-f_0)/2]*[ 1 + TANH[(s-s_2)/(2*delta_2)] ]
    !#             + ...
    !########################################################################
    subroutine BLD_TANH(idir, iseg, x, nmax, work)
        implicit none

        integer(ci), intent(IN) :: idir, iseg, nmax
        real(cp), intent(INOUT) :: x(nmax), work(nmax)

        ! -----------------------------------------------------------------------
        real(cp) st(3), f(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
        integer(ci) im

        ! #######################################################################
        ! define local variables for readability below; strides of 3
        st = g_build(idir)%vals(1 :: 3, iseg)     ! transition point in uniform grid
        f = g_build(idir)%vals(2 :: 3, iseg)     ! ratio dx(i)/dx_0
        delta = g_build(idir)%vals(3 :: 3, iseg)

        ! create mapping from grid s to grid x
        work = 0.0_cp
        do im = 1, 3                           ! 3 modes at most
            if (ABS(delta(im)) > 0.0_cp) then
                work(:) = work(:) + (f(im) - 1.0_cp)*delta(im)*LOG(EXP((x(:) - st(im))/delta(im)) + 1.0_cp)
            end if
        end do
        work = work - work(1) ! Integration constant

        x = x + work

        return
    end subroutine BLD_TANH

    !########################################################################
    !# The stretching factor is given by a hyperbolic tangent function, which means (see manual)
    !#
    !# dx/ds = [ 1 + EXP[(s-s_1)/delta_1] ] **(delta_1*f1/h_0)
    !#        *[ 1 + EXP[(s-s_2)/delta_2] ] **(delta_2*f2/h_0)
    !#        *...
    !#
    !# Using compact schemes to integrate this equation
    !########################################################################
    subroutine BLD_EXP(idir, iseg, x, nmax, w)
        implicit none

        integer(ci), intent(IN) :: idir, iseg, nmax
        real(cp), intent(INOUT) :: x(nmax), w(nmax, 8)

        ! -----------------------------------------------------------------------
        real(cp) st(3), df(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
        integer(ci) im, ibc, i1
        real(cp) ds

        ! #######################################################################
        ds = (x(nmax) - x(1))/real(nmax - 1,cp)

        ! define local variables for readability below; strides of 3
        st = g_build(idir)%vals(1 :: 3, iseg)     ! transition point in uniform grid
        df = g_build(idir)%vals(2 :: 3, iseg)/ds ! normalized maximum stretching factor
        delta = g_build(idir)%vals(3 :: 3, iseg)

        ! create mapping from grid s to grid x
#define jac(j)    w(j,6)
#define rhs(j)    w(j,7)
#define result(j) w(j,8)
        rhs(:) = 1.0_cp
        do im = 1, 3               ! 3 modes at most
            if (ABS(delta(im)) > 0.0_cp) then
                rhs(:) = rhs(:)*(EXP((x(:) - st(im))/delta(im)) + 1.0_cp)**(df(im)*delta(im))
            end if
        end do
        jac(:) = ds
        ibc = 1                           ! Boundary condition at the first node
        i1 = 1                            ! One equation to be solved
        call INT_C1N6_LHS(nmax, ibc, w(1, 1), w(1, 2), w(1, 3), w(1, 4), w(1, 5))
        call INT_C1N6_RHS(nmax, i1, ibc, jac(1), rhs(1), result(1))
        call PENTADFS(nmax, w(2, 1), w(2, 2), w(2, 3), w(2, 4), w(2, 5))
        call PENTADSS(nmax, i1, w(2, 1), w(2, 2), w(2, 3), w(2, 4), w(2, 5), result(2))
        x(2:nmax) = x(1) + result(2:nmax)  ! Boundary condition
#undef jac
#undef rhs
#undef result

        return
    end subroutine BLD_EXP

    !########################################################################
    ! Geometric progression and other options
    !########################################################################
    subroutine BLD_THEREST(idir, iseg, x, nmax)
        implicit none

        integer(ci), intent(IN) :: idir, iseg, nmax
        real(cp), intent(INOUT) :: x(nmax)

        ! -----------------------------------------------------------------------
        real(cp) a, b, c, d, e
        real(cp) eta, deta, dx
        integer(ci) n

        ! #######################################################################
        call BLD_CONSTANTS(nmax, x(1), x(nmax), &
                           g_build(idir)%opts(1, iseg), g_build(idir)%opts(2, iseg), &
                           g_build(idir)%vals(1, iseg), g_build(idir)%vals(2, iseg), &
                           g_build(idir)%vals(3, iseg), g_build(idir)%vals(4, iseg), a, b, c, d, e)

        ! Go from computational eta in [0,1] to physical domains
        dx = 1.0_cp
        if (nmax == 1) then
            deta = 1.0_cp
        else; 
            deta = 1.0_cp/real(nmax - 1,cp)
        end if

        do n = 2, nmax
            eta = real(n - 1,cp)*deta

            select case (g_build(idir)%opts(1, iseg))
            case (1)     ! Colonius, Lele and Moin
                x(n) = e + a*eta + d*LOG(EXP(c*(a*eta - b)) + 1.0_cp - EXP(-b*c))
            case (2, 3)   ! 2nd- 3rd-Order Polynomial Stretching
                x(n) = a + b*eta + c*eta*eta + d*eta*eta*eta
            case (4)     ! Geometric prograssion
                dx = dx*g_build(idir)%vals(1, iseg)
                x(n) = x(n - 1) + dx

            end select

        end do

        return
    end subroutine BLD_THEREST

    ! #######################################################################
    subroutine BLD_CONSTANTS(imax, vbeg, vend, iopt1_loc, iopt2_loc, val_1, val_2, val_3, val_4, a, b, c, d, e)

        implicit none

        integer(ci) imax, iopt1_loc, iopt2_loc
        real(cp) vbeg, vend, val_1, val_2, val_3, val_4
        real(cp) a, b, c, d, e

        real(cp) x1, x2, x3, x4, valmx
        real(cp) z1, z2, z3, z4
        integer(ci) indx1, indx2, indx3, indx4

        ! ######################################
        ! # Colonius, Lele and Moin Stretching #
        ! ######################################

        if (iopt1_loc == 1) then
            x2 = val_4 - vbeg
            x3 = vend - vbeg

            !
            ! Calculate Constants
            ! -------------------

            a = FLOAT(imax - 1)*val_1
            b = (a*(1.0 + val_2/val_1) - x3)/(val_2/val_1)
            c = val_3/val_1 - 1.0
            c = LOG(val_2/(c*val_1))/(b - x2)
            d = val_2/(c*val_1)
            e = vbeg

            !
            ! Force Values to be Exact at Domain Edges (x=0 at Xi=0)
            ! ------------------------------------------------------

            valmx = a + d*LOG(EXP(c*(a - b)) + 1.0 - EXP(-b*c))
            a = a*(x3/valmx)
            b = b*(x3/valmx)
            c = c/(x3/valmx)
            d = d*(x3/valmx)

            ! ######################################
            ! # Second order polynomial stretching #
            ! ######################################

        elseif (iopt1_loc == 2) then

            !
            ! Points to be Used in Calculating Stretching Constants
            ! -----------------------------------------------------

            ! Clustering points at i=1

            if (iopt2_loc == 1) then
                x1 = vbeg
                indx1 = 1
                x2 = x1 + val_1
                indx2 = 2
                x3 = vend
                indx3 = imax

                ! Clustering points at i=i_max

            elseif (iopt2_loc == 2) then
                x1 = vbeg
                indx1 = 1
                x2 = vend - val_1
                indx2 = imax - 1
                x3 = vend
                indx3 = imax

                ! Unknown stretching

            else
                write (*, *) 'Unknown stretching option iopt2_loc =', iopt2_loc
                stop
            end if

            !
            ! Calculate Stretching Constants
            ! ------------------------------

            z1 = FLOAT(indx1 - 1)/FLOAT(imax - 1)
            z2 = FLOAT(indx2 - 1)/FLOAT(imax - 1)
            z3 = FLOAT(indx3 - 1)/FLOAT(imax - 1)

            a = (-(x3*z1**2*z2) + x3*z1*z2**2 + x2*z1**2*z3 - &
                 x1*z2**2*z3 - x2*z1*z3**2 + x1*z2*z3**2)/ &
                ((-z1 + z2)*(-z1 + z3)*(-z2 + z3))
            b = (-(x2*z1**2) + x3*z1**2 + x1*z2**2 - x3*z2**2 - x1*z3**2 + &
                 x2*z3**2)/((-z1 + z2)*(-z1 + z3)*(-z2 + z3))
            c = (x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)/ &
                ((-z1 + z2)*(-z1 + z3)*(-z2 + z3))
            d = 0.0

            !
            ! Force Values to be Exact at Indx=1
            ! ----------------------------------

            a = a - (a + b*z1 + c*z1*z1 + d*z1*z1*z1 - x1)

            ! #####################################
            ! # Third order polynomial stretching #
            ! #####################################

        elseif (iopt1_loc == 3) then

            !
            ! Points to be Used in Calculating Stretching Constants
            ! -----------------------------------------------------

            ! Clustering Points Both Ends

            if (iopt2_loc == 1) then
                x1 = vbeg
                indx1 = 1
                x2 = x1 + val_1
                indx2 = 2
                x3 = vend - val_2
                indx3 = imax - 1
                x4 = vend
                indx4 = imax

                ! Clustering Points at an Internal Point

            elseif (iopt2_loc == 2) then
                x1 = vbeg
                indx1 = 1
                x2 = val_2 - val_1/2.0
                indx2 = INT(val_3*FLOAT(imax))
                x3 = val_2 + val_1/2.0
                indx3 = indx2 + 1
                x4 = vend
                indx4 = imax

                !  Unknown Stretching Option

            else
                write (*, *) 'Unknown stretching option iopt2_loc =', iopt2_loc
                stop
            end if

            !
            ! Calculate Stretching Constants
            ! ------------------------------

            z1 = FLOAT(indx1 - 1)/FLOAT(imax - 1)
            z2 = FLOAT(indx2 - 1)/FLOAT(imax - 1)
            z3 = FLOAT(indx3 - 1)/FLOAT(imax - 1)
            z4 = FLOAT(indx4 - 1)/FLOAT(imax - 1)

            a = (x4*z1**3*z2**2*z3 - x4*z1**2*z2**3*z3 - &
                 x4*z1**3*z2*z3**2 + x4*z1*z2**3*z3**2 + x4*z1**2*z2*z3**3 - &
                 x4*z1*z2**2*z3**3 - x3*z1**3*z2**2*z4 + x3*z1**2*z2**3*z4 + &
                 x2*z1**3*z3**2*z4 - x1*z2**3*z3**2*z4 - x2*z1**2*z3**3*z4 + &
                 x1*z2**2*z3**3*z4 + x3*z1**3*z2*z4**2 - x3*z1*z2**3*z4**2 - &
                 x2*z1**3*z3*z4**2 + x1*z2**3*z3*z4**2 + x2*z1*z3**3*z4**2 - &
                 x1*z2*z3**3*z4**2 - x3*z1**2*z2*z4**3 + x3*z1*z2**2*z4**3 + &
                 x2*z1**2*z3*z4**3 - x1*z2**2*z3*z4**3 - x2*z1*z3**2*z4**3 + &
                 x1*z2*z3**2*z4**3)/ &
                ((-z1 + z2)*(-z1 + z3)*(-z2 + z3)*(-z1 + z4)*(-z2 + z4)* &
                 (-z3 + z4))
            b = (x3*z1**3*z2**2 - x4*z1**3*z2**2 - x3*z1**2*z2**3 + &
                 x4*z1**2*z2**3 - x2*z1**3*z3**2 + x4*z1**3*z3**2 + &
                 x1*z2**3*z3**2 - x4*z2**3*z3**2 + x2*z1**2*z3**3 - &
                 x4*z1**2*z3**3 - x1*z2**2*z3**3 + x4*z2**2*z3**3 + &
                 x2*z1**3*z4**2 - x3*z1**3*z4**2 - x1*z2**3*z4**2 + &
                 x3*z2**3*z4**2 + x1*z3**3*z4**2 - x2*z3**3*z4**2 - &
                 x2*z1**2*z4**3 + x3*z1**2*z4**3 + x1*z2**2*z4**3 - &
                 x3*z2**2*z4**3 - x1*z3**2*z4**3 + x2*z3**2*z4**3)/ &
                ((-z1 + z2)*(-z1 + z3)*(-z2 + z3)*(-z1 + z4)*(-z2 + z4)* &
                 (-z3 + z4))
            c = (-x3*z1**3*z2 + x4*z1**3*z2 + x3*z1*z2**3 - x4*z1*z2**3 + &
                 x2*z1**3*z3 - x4*z1**3*z3 - x1*z2**3*z3 + x4*z2**3*z3 - &
                 x2*z1*z3**3 + x4*z1*z3**3 + x1*z2*z3**3 - x4*z2*z3**3 - &
                 x2*z1**3*z4 + x3*z1**3*z4 + x1*z2**3*z4 - x3*z2**3*z4 - &
                 x1*z3**3*z4 + x2*z3**3*z4 + x2*z1*z4**3 - x3*z1*z4**3 - &
                 x1*z2*z4**3 + x3*z2*z4**3 + x1*z3*z4**3 - x2*z3*z4**3)/ &
                ((-z1 + z2)*(-z1 + z3)*(-z2 + z3)*(-z1 + z4)*(-z2 + z4)* &
                 (-z3 + z4))
            d = -(((-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 &
                               + x1*z3 - x2*z3)* &
                   (-z1 + z4)*(-z2 + z4) + &
                   (-z1 + z2)*(z1 - z3)*(z2 - z3)* &
                   (-(x2*z1) + x4*z1 + x1*z2 - x4*z2 - x1*z4 + x2*z4))/ &
                  ((-z1 + z2)**2*(z1 - z3)*(z2 - z3)*(z1 - z4)*(z2 - z4)* &
                   (-z3 + z4)))

            !
            ! Force Values to be Exact at Indx=1
            ! ----------------------------------

            a = a - (a + b*z1 + c*z1*z1 + d*z1*z1*z1 - x1)

        end if

        return
    end subroutine BLD_CONSTANTS

end module GRID_LOCAL
