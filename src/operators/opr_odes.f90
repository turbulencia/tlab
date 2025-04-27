! # Solvers for boundary value problems of
! # linear ordinary differential equations with constant coefficients
! #
! # Mellado & Ansorge, 2012: Z. Angew. Math. Mech., 1-12, 10.1002/zamm.201100078
! #
module OPR_ODES
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_MIN, BCS_MAX
    use FDM_Integral, only: fdm_integral_dt, FDM_Int1_Solve
    implicit none
    private

    ! First-order ODEs. Shall I move here FDM_Integral module?

    ! Second-order ODEs. Reducing problem to a system of 2 first-order equations
    public :: OPR_ODE2_Factorize_DD                   ! Dirichlet/Dirichlet boundary conditions
    public :: OPR_ODE2_Factorize_NN                   ! Neumann/Neumann boundary conditions

    public :: OPR_ODE2_Factorize_DD_Sing
    public :: OPR_ODE2_Factorize_DN_Sing
    public :: OPR_ODE2_Factorize_ND_Sing
    public :: OPR_ODE2_Factorize_NN_Sing

contains
    !########################################################################
    !#
    !#     (u')'_i - \lambda^2 u_i = f_i        N - 2 eqns (singular case when \lambda = 0)
    !#     u_1, u'_1, u_N, or u'_N given        2   eqns
    !#     Au' = Bu                             N   eqns
    !#     A(u')' = Bu'                         N   eqns
    !#
    !########################################################################

    !########################################################################
    !########################################################################
    ! Dirichlet/Neumann boundary conditions
    subroutine OPR_ODE2_Factorize_DN_Sing(nlines, fdmi, u, f, bcs, v, wrk1d, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(out) :: u(nlines, size(fdmi(1)%lhs, 1)) ! solution
        real(wp), intent(out) :: v(nlines, size(fdmi(1)%lhs, 1)) ! derivative of solution
        real(wp), intent(inout) :: wrk1d(size(fdmi(1)%lhs), 3)
        real(wp), intent(inout) :: wrk2d(nlines, 3)

        ! -----------------------------------------------------------------------
        integer(wi) nx, i
        real(wp) du1_n(1), f1

        ! #######################################################################
        nx = size(fdmi(1)%lhs, 1)

#define f1(i) wrk1d(i,1)
#define v1(i) wrk1d(i,2)
#define u1(i) wrk1d(i,3)
#define du0_n(i) wrk2d(i,3)

        ! -----------------------------------------------------------------------
        ! solve for v^(0) in v' = f , v_n given
        f(:, 1) = 0.0_wp
        v(:, nx) = bcs(:, 2)
        call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, f, v, wrk2d)

        ! solve for v^(1)
        f1(:) = 0.0_wp; f1(1) = 1.0_wp
        v1(nx) = 0.0_wp
        call FDM_Int1_Solve(1, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, f1(:), v1(:), wrk2d)

        ! -----------------------------------------------------------------------
        ! solve for u^(0) in u' = v, u_1 given
        u(:, 1) = bcs(:, 1)
        call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, v, u, wrk2d, du0_n(:))

        ! solve for u^(1)
        u1(1) = 0.0_wp
        call FDM_Int1_Solve(1, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, v1(:), u1(:), wrk2d, du1_n(1))

        ! -----------------------------------------------------------------------
        ! Constraint
        f1 = 1.0_wp/(du1_n(1) - v1(1))
        du0_n(:) = (v(:, 1) - du0_n(:))*f1

        ! Result
        do i = 1, nx
            u(:, i) = u(:, i) + du0_n(:)*u1(i)
            v(:, i) = v(:, i) + du0_n(:)*v1(i)
        end do

#undef du0_n
#undef f1
#undef v1
#undef u1

        return
    end subroutine OPR_ODE2_Factorize_DN_Sing

    !########################################################################
    !########################################################################
    ! Neumann/Dirichlet boundary conditions
    subroutine OPR_ODE2_Factorize_ND_Sing(nlines, fdmi, u, f, bcs, v, wrk1d, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(out) :: u(nlines, size(fdmi(1)%lhs, 1)) ! solution
        real(wp), intent(out) :: v(nlines, size(fdmi(1)%lhs, 1)) ! derivative of solution
        real(wp), intent(inout) :: wrk1d(size(fdmi(1)%lhs), 3)
        real(wp), intent(inout) :: wrk2d(nlines, 3)

        ! -----------------------------------------------------------------------
        integer(wi) nx, i
        real(wp) du1_n(1), fn

        ! #######################################################################
        nx = size(fdmi(1)%lhs, 1)

#define f1(i) wrk1d(i,1)
#define v1(i) wrk1d(i,2)
#define u1(i) wrk1d(i,3)
#define du0_n(i) wrk2d(i,3)

        ! -----------------------------------------------------------------------
        ! solve for v^(0) in v' = f , v_1 given
        f(:, nx) = 0.0_wp
        v(:, 1) = bcs(:, 1)
        call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, f, v, wrk2d)

        ! solve for v^(1)
        f1(:) = 0.0_wp; f1(nx) = 1.0_wp
        v1(1) = 0.0_wp
        call FDM_Int1_Solve(1, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, f1(:), v1(:), wrk2d)

        ! -----------------------------------------------------------------------
        ! solve for u^(0) in u' = v, u_n given
        u(:, nx) = bcs(:, 2)
        call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, v, u, wrk2d, du0_n(:))

        ! solve for u^(1)
        u1(nx) = 0.0_wp
        call FDM_Int1_Solve(1, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, v1(:), u1(:), wrk2d, du1_n(1))

        ! -----------------------------------------------------------------------
        ! Constraint
        fn = 1.0_wp/(du1_n(1) - v1(nx))
        du0_n(:) = (v(:, nx) - du0_n(:))*fn

        ! Result
        do i = 1, nx
            u(:, i) = u(:, i) + du0_n(:)*u1(i)
            v(:, i) = v(:, i) + du0_n(:)*v1(i)
        end do

#undef du0_n
#undef f1
#undef v1
#undef u1

        return
    end subroutine OPR_ODE2_Factorize_ND_Sing

    !########################################################################
    !########################################################################
    ! Neumann/Neumann boundary conditions; must be compatible!
    subroutine OPR_ODE2_Factorize_NN_Sing(nlines, fdmi, u, f, bcs, v, wrk1d, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(inout) :: bcs(nlines, 2)
        real(wp), intent(out) :: u(nlines, size(fdmi(1)%lhs, 1)) ! solution
        real(wp), intent(out) :: v(nlines, size(fdmi(1)%lhs, 1)) ! derivative of solution
        real(wp), intent(inout) :: wrk1d(size(fdmi(1)%lhs), 3)
        real(wp), intent(inout) :: wrk2d(nlines, 3)

        ! #######################################################################
        ! Assumes compatible problem, i.e., bcs_n -bcs_1 = int f
        ! We write it in terms of the Dirichlet/Neumann problem
        ! (We could have written it in terms of the Neumann/Dirichlet problem as well)
        bcs(:, 1) = 0.0_wp
        call OPR_ODE2_Factorize_DN_Sing(nlines, fdmi, u, f, bcs, v, wrk1d, wrk2d)

        return
    end subroutine OPR_ODE2_Factorize_NN_Sing

    !########################################################################
    !########################################################################
    ! Dirichlet/Dirichlet boundary conditions
    subroutine OPR_ODE2_Factorize_DD_Sing(nlines, fdmi, u, f, bcs, v, wrk1d, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(out) :: u(nlines, size(fdmi(1)%lhs, 1)) ! solution
        real(wp), intent(out) :: v(nlines, size(fdmi(1)%lhs, 1)) ! derivative of solution
        real(wp), intent(inout) :: wrk1d(size(fdmi(1)%lhs, 1), 4)
        real(wp), intent(inout) :: wrk2d(nlines, 3)

        ! -----------------------------------------------------------------------
        integer(wi) nx, i
        real(wp) dummy, du1_n(1), fn

        ! #######################################################################
        nx = size(fdmi(1)%lhs, 1)

#define f1(i) wrk1d(i,1)
#define v1(i) wrk1d(i,2)
#define u1(i) wrk1d(i,3)
#define sp(i) wrk1d(i,4)
#define du0_n(i) wrk2d(i,3)

        ! -----------------------------------------------------------------------
        ! solve for v^(0) in v' = f , v_1 given (0 for now, to be found later on)
        f(:, nx) = 0.0_wp
        v(:, 1) = 0.0_wp
        call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, f, v, wrk2d)

        ! solve for v^(1)
        f1(:) = 0.0_wp; f1(nx) = 1.0_wp
        v1(1) = 0.0_wp
        call FDM_Int1_Solve(1, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, f1(:), v1(:), wrk2d)

        ! -----------------------------------------------------------------------
        ! solve for u^(0) in u' = v, u_n given
        u(:, nx) = bcs(:, 2)
        call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, v, u, wrk2d, du0_n(:))

        ! solve for u^(1)
        u1(nx) = 0.0_wp
        call FDM_Int1_Solve(1, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, v1(:), u1(:), wrk2d, du1_n(1))

        ! solve for s^(+); this is the grid node displacements x(:) - x_n,
        ! but we calculate it here to avoid having to pass the array with node positions
        f1(:) = 1.0_wp
        sp(nx) = 0.0_wp
        call FDM_Int1_Solve(1, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, f1(:), sp(:), wrk2d)

        ! -----------------------------------------------------------------------
        ! Constraint
        fn = 1.0_wp/(du1_n(1) - v1(nx))
        du0_n(:) = (v(:, nx) - du0_n(:))*fn

        ! Contribution from v_1 to satisfy bc at the bottom
        dummy = 1.0_wp/sp(1)
        v(:, 1) = (bcs(:, 1) - (u(:, 1) + du0_n(:)*u1(1)))*dummy
        u(:, 1) = bcs(:, 1)

        ! Result
        do i = 2, nx
            u(:, i) = u(:, i) + du0_n(:)*u1(i) + v(:, 1)*sp(i)
            v(:, i) = v(:, i) + du0_n(:)*v1(i) + v(:, 1)
        end do

#undef du0_n
#undef f1
#undef v1
#undef u1
#undef sp

        return
    end subroutine OPR_ODE2_Factorize_DD_Sing

    !########################################################################
    !########################################################################
    ! Neumann/Neumann boundary conditions
    subroutine OPR_ODE2_Factorize_NN(nlines, fdmi, rhsi_b, rhsi_t, u, f, bcs, v, wrk1d, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(inout) :: fdmi(2)
        real(wp), intent(in) :: rhsi_b(:, :), rhsi_t(:, :)
        real(wp), intent(inout) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(out) :: u(nlines, size(fdmi(1)%lhs, 1)) ! solution
        real(wp), intent(out) :: v(nlines, size(fdmi(1)%lhs, 1)) ! derivative of solution
        real(wp), intent(inout) :: wrk1d(3, size(fdmi(1)%lhs, 1), 2)
        real(wp), intent(inout) :: wrk2d(max(nlines, 3), 3)

        ! -----------------------------------------------------------------------
        integer(wi) nx, i
        real(wp) lambda, a(3, 3), der_bcs(3)

        ! #######################################################################
        lambda = fdmi(BCS_MIN)%lambda

        nx = size(fdmi(1)%lhs, 1)

#define f1(j,i) wrk1d(j,i,1)

#define v1(i) wrk1d(1,i,2)
#define em(i) wrk1d(2,i,2)
#define dd(i) wrk1d(3,i,2)

#define u1(i) wrk1d(1,i,1)
#define sp(i) wrk1d(2,i,1)
#define ep(i) wrk1d(3,i,1)

#define du1_n der_bcs(1)
#define dsp_n der_bcs(2)
#define dep_n der_bcs(3)

#define du0_n(i) wrk2d(i,3)
#define fn(i) wrk2d(i,3)

        ! -----------------------------------------------------------------------
        ! solve for v^(0) in v' + lambda v= f , v_1 given (0 for now, to be found later on)
        f(:, nx) = 0.0_wp
        v(:, 1) = 0.0_wp
        call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), rhsi_b, f, v, wrk2d)

        ! solve for v^(1) and e^(-); 1 line is not used but we need a 3 x nx array
        f1(:, :) = 0.0_wp
        f1(1, nx) = 1.0_wp; v1(1) = 0.0_wp
        em(1) = 1.0_wp
        dd(1) = 0.0_wp
        call FDM_Int1_Solve(3, fdmi(BCS_MIN), rhsi_b, f1(1, 1), v1(1), wrk2d)

        ! -----------------------------------------------------------------------
        ! solve for u^(0) in u' - lambda u = v, u_n given (0 for now, to be found later on)
        u(:, nx) = 0.0_wp
        call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), rhsi_t, v, u, wrk2d, du0_n(:))

        ! solve for u^(1) and s^(+) and e^(+)
        u1(nx) = 0.0_wp
        sp(nx) = 0.0_wp
        dd(:) = 0.0_wp; ep(nx) = 1.0_wp
        call FDM_Int1_Solve(3, fdmi(BCS_MAX), rhsi_t, v1(1), u1(1), wrk2d, der_bcs)

        ! -----------------------------------------------------------------------
        ! Constraint and boundary conditions
        ! System
        a(1, 1) = 1.0_wp + lambda*sp(1)
        a(2, 1) = em(nx)
        a(3, 1) = dsp_n

        a(1, 2) = lambda*ep(1)
        a(2, 2) = lambda
        a(3, 2) = dep_n

        a(1, 3) = lambda*u1(1)
        a(2, 3) = v1(nx)
        a(3, 3) = du1_n

        ! LU decomposition
        a(1, 2) = a(1, 2)/a(1, 1)
        a(2, 2) = a(2, 2) - a(2, 1)*a(1, 2)
        a(3, 2) = a(3, 2) - a(3, 1)*a(1, 2)

        a(1, 3) = a(1, 3)/a(1, 1)
        a(2, 3) = (a(2, 3) - a(2, 1)*a(1, 3))/a(2, 2)
        a(3, 3) = a(3, 3) - a(3, 1)*a(1, 3) - a(3, 2)*a(2, 3)

        ! Solution
        v(:, 1) = (bcs(:, 1) - lambda*u(:, 1))/a(1, 1)
        u(:, nx) = (bcs(:, 2) - v(:, nx) - a(2, 1)*v(:, 1))/a(2, 2)
        fn(:) = (bcs(:, 2) - du0_n(:) - a(3, 1)*v(:, 1) - a(3, 2)*u(:, nx))/a(3, 3)

        u(:, nx) = u(:, nx) - a(2, 3)*fn(:)
        v(:, 1) = v(:, 1) - a(1, 2)*u(:, nx) - a(1, 3)*fn(:)

        ! Result
        i = nx
        v(:, i) = v(:, i) + fn(:)*v1(i) + v(:, 1)*em(i) + lambda*u(:, i)
        do i = nx - 1, 2, -1
            u(:, i) = u(:, i) + fn(:)*u1(i) + v(:, 1)*sp(i) + u(:, nx)*ep(i)
            v(:, i) = v(:, i) + fn(:)*v1(i) + v(:, 1)*em(i) + lambda*u(:, i)
        end do
        i = 1
        u(:, i) = u(:, i) + fn(:)*u1(i) + v(:, 1)*sp(i) + u(:, nx)*ep(i)
        v(:, i) = v(:, i) + lambda*u(:, i)      ! Final correction to get u' instead of v = u'- lambda u.

#undef fn
#undef du0_n

#undef f1
#undef v1
#undef em
#undef dd

#undef u1
#undef ep
#undef sp

#undef du1_n
#undef dsp_n
#undef dep_n

        return
    end subroutine OPR_ODE2_Factorize_NN

    !########################################################################
    !########################################################################
    ! Dirichlet/Dirichlet boundary conditions
    subroutine OPR_ODE2_Factorize_DD(nlines, fdmi, rhsi_b, rhsi_t, u, f, bcs, v, wrk1d, wrk2d)
        integer(wi) nlines
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(in) :: rhsi_b(:, :), rhsi_t(:, :)
        real(wp), intent(inout) :: f(nlines, size(fdmi(1)%lhs, 1))
        real(wp), intent(in) :: bcs(nlines, 2)
        real(wp), intent(out) :: u(nlines, size(fdmi(1)%lhs, 1)) ! solution
        real(wp), intent(out) :: v(nlines, size(fdmi(1)%lhs, 1)) ! derivative of solution
        real(wp), intent(inout) :: wrk1d(2, size(fdmi(1)%lhs, 1), 2)
        real(wp), intent(inout) :: wrk2d(max(nlines, 2), 3)

        ! -----------------------------------------------------------------------
        integer(wi) nx, i
        real(wp) lambda, dummy, aa, bb, der_bcs(2)

        ! #######################################################################
        lambda = fdmi(BCS_MIN)%lambda

        nx = size(fdmi(1)%lhs, 1)

#define f1(j,i) wrk1d(j,i,1)

#define v1(i) wrk1d(1,i,2)
#define em(i) wrk1d(2,i,2)

#define u1(i) wrk1d(1,i,1)
#define sp(i) wrk1d(2,i,1)

#define du1_n der_bcs(1)
#define dsp_n der_bcs(2)

#define du0_n(i) wrk2d(i,3)
#define fn(i) wrk2d(i,3)

        ! -----------------------------------------------------------------------
        ! solve for v^(0) in v' + lambda v= f , v_1 given (0 for now, to be found later on)
        f(:, nx) = 0.0_wp
        v(:, 1) = 0.0_wp
        call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), rhsi_b, f, v, wrk2d)

        ! solve for v^(1) and e^(-); 1 line is not used but we need a 2 x nx array
        f1(:, :) = 0.0_wp
        f1(1, nx) = 1.0_wp; v1(1) = 0.0_wp
        em(1) = 1.0_wp
        call FDM_Int1_Solve(2, fdmi(BCS_MIN), rhsi_b, f1(1, 1), v1(1), wrk2d)

        ! -----------------------------------------------------------------------
        ! solve for u^(0) in u' - lambda u = v, u_n given
        u(:, nx) = bcs(:, 2)
        call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), rhsi_t, v, u, wrk2d, du0_n(:))

        ! solve for u^(1) and s^(+)
        u1(nx) = 0.0_wp
        sp(nx) = 0.0_wp
        call FDM_Int1_Solve(2, fdmi(BCS_MAX), rhsi_t, v1(1), u1(1), wrk2d, der_bcs)

        ! -----------------------------------------------------------------------
        ! Constraint and bottom boundary condition
        aa = du1_n - v1(nx)
        bb = dsp_n - em(nx)
        dummy = 1.0_wp/(aa*sp(1) - bb*u1(1))
        v(:, 1) = (aa*(bcs(:, 1) - u(:, 1)) - u1(1)*(lambda*bcs(:, 2) - du0_n(:) + v(:, nx)))*dummy   ! q1
        fn(:) = (sp(1)*(lambda*bcs(:, 2) - du0_n(:) + v(:, nx)) - bb*(bcs(:, 1) - u(:, 1)))*dummy

        ! Result
        do i = nx, 2, -1
            u(:, i) = u(:, i) + fn(:)*u1(i) + v(:, 1)*sp(i)
            v(:, i) = v(:, i) + fn(:)*v1(i) + v(:, 1)*em(i) + lambda*u(:, i)
        end do
        i = 1
        u(:, i) = bcs(:, 1)
        v(:, i) = v(:, i) + lambda*u(:, i)

#undef fn
#undef du0_n

#undef f1
#undef v1
#undef em

#undef u1
#undef sp

#undef du1_n
#undef dsp_n

        return
    end subroutine OPR_ODE2_Factorize_DD

end module OPR_ODES
