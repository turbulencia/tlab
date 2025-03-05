module GRID_LOCAL
    use TLab_Constants, only: wp, wi
    implicit none
    save

    integer(wi), parameter :: GTYPE_UNIFORM = 0
    integer(wi), parameter :: GTYPE_TANH = 5
    integer(wi), parameter :: GTYPE_EXP = 6

    integer(wi), parameter :: MAX_PARAMES = 9
    integer(wi), parameter :: MAX_OPTIONS = 5
    integer(wi), parameter :: MAX_SEGMENT = 10

    type grid_build_dt                              ! Information to construct grid in one direction
        sequence
        integer(wi) nseg                            ! number of segments in this direction
        logical mirrored                            ! It true, mirror the grid
        real(wp) fixed_scale                        ! If positive, rescale grid to this value
        integer(wi) size(MAX_SEGMENT)               ! number of points in each segment
        real(wp) end(MAX_SEGMENT)                   ! physical end of each segment
        integer(wi) opts(MAX_OPTIONS, MAX_SEGMENT)  ! 0   uniform segment
        ! 1   Colonius, Lele and Moin stretching
        ! 2   Second order polynomial stretching
        ! 3   Third order polynomial stretching
        ! 4   Geometric progression
        ! 5   Hyperbolic tangent
        ! 6   Exponential
        real(wp) vals(MAX_PARAMES, MAX_SEGMENT)
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

        integer(wi), intent(IN) :: idir, iseg, nmax
        real(wp), intent(INOUT) :: x(nmax), work(nmax)

        ! -----------------------------------------------------------------------
        real(wp) st(3), f(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
        integer(wi) im

        ! #######################################################################
        ! define local variables for readability below; strides of 3
        st = g_build(idir)%vals(1 :: 3, iseg)     ! transition point in uniform grid
        f = g_build(idir)%vals(2 :: 3, iseg)     ! ratio dx(i)/dx_0
        delta = g_build(idir)%vals(3 :: 3, iseg)

        ! create mapping from grid s to grid x
        work = 0.0_wp
        do im = 1, 3                           ! 3 modes at most
            if (abs(delta(im)) > 0.0_wp) then
                work(:) = work(:) + (f(im) - 1.0_wp)*delta(im)*log(exp((x(:) - st(im))/delta(im)) + 1.0_wp)
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
        use TLab_Constants, only: BCS_MIN
        use FDM, only: fdm_dt, FDM_Initialize, FDM_COM6_JACOBIAN
        use FDM_Integral, only: FDM_Int1_Solve, fdm_integral_dt
        integer(wi), intent(IN) :: idir, iseg, nmax
        real(wp), intent(INOUT) :: x(nmax), w(nmax, 8)

        ! -----------------------------------------------------------------------
        type(fdm_dt) g
        type(fdm_integral_dt) fdmi(2)
        real(wp) st(3), df(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
        integer(wi) im
        real(wp) ds, aux(2)

        ! #######################################################################
        ds = (x(nmax) - x(1))/real(nmax - 1, wp)

        ! define local variables for readability below; strides of 3
        st = g_build(idir)%vals(1 :: 3, iseg)       ! transition point in uniform grid
        df = g_build(idir)%vals(2 :: 3, iseg)/ds    ! normalized maximum stretching factor
        delta = g_build(idir)%vals(3 :: 3, iseg)

        ! create mapping from grid s to grid x
#define rhs(j)    w(j,7)
#define result(j) w(j,8)
        rhs(:) = 1.0_wp
        do im = 1, 3               ! 3 modes at most
            if (abs(delta(im)) > 0.0_wp) then
                rhs(:) = rhs(:)*(exp((x(:) - st(im))/delta(im)) + 1.0_wp)**(df(im)*delta(im))
            end if
        end do
        !
        g%size = nmax
        g%uniform = .true.
        g%periodic = .false.
        g%mode_fdm1 = FDM_COM6_JACOBIAN
        g%mode_fdm2 = FDM_COM6_JACOBIAN
        call FDM_Initialize(x, g, fdmi)
        ! x(1) is already set
        call FDM_Int1_Solve(1, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, rhs(:), result(:), aux)
        x(:) = result(:)
#undef rhs
#undef result

        return
    end subroutine BLD_EXP

    !########################################################################
    ! Geometric progression and other options
    !########################################################################
    subroutine BLD_THEREST(idir, iseg, x, nmax)
        implicit none

        integer(wi), intent(IN) :: idir, iseg, nmax
        real(wp), intent(INOUT) :: x(nmax)

        ! -----------------------------------------------------------------------
        real(wp) a, b, c, d, e
        real(wp) eta, deta, dx
        integer(wi) n

        ! #######################################################################
        call BLD_CONSTANTS(nmax, x(1), x(nmax), &
                           g_build(idir)%opts(1, iseg), g_build(idir)%opts(2, iseg), &
                           g_build(idir)%vals(1, iseg), g_build(idir)%vals(2, iseg), &
                           g_build(idir)%vals(3, iseg), g_build(idir)%vals(4, iseg), a, b, c, d, e)

        ! Go from computational eta in [0,1] to physical domains
        dx = 1.0_wp
        if (nmax == 1) then
            deta = 1.0_wp
        else; 
            deta = 1.0_wp/real(nmax - 1, wp)
        end if

        do n = 2, nmax
            eta = real(n - 1, wp)*deta

            select case (g_build(idir)%opts(1, iseg))
            case (1)     ! Colonius, Lele and Moin
                x(n) = e + a*eta + d*log(exp(c*(a*eta - b)) + 1.0_wp - exp(-b*c))
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

        integer(wi) imax, iopt1_loc, iopt2_loc
        real(wp) vbeg, vend, val_1, val_2, val_3, val_4
        real(wp) a, b, c, d, e

        real(wp) x1, x2, x3, x4, valmx
        real(wp) z1, z2, z3, z4
        integer(wi) indx1, indx2, indx3, indx4

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
            c = log(val_2/(c*val_1))/(b - x2)
            d = val_2/(c*val_1)
            e = vbeg

            !
            ! Force Values to be Exact at Domain Edges (x=0 at Xi=0)
            ! ------------------------------------------------------

            valmx = a + d*log(exp(c*(a - b)) + 1.0 - exp(-b*c))
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
                indx2 = int(val_3*FLOAT(imax))
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
