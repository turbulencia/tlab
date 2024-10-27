#include "dns_const.h"

program VPARTIAL
    use TLab_Constants
    use TLab_Types, only:  grid_dt
    use TLAB_VARS, only: g, imax, jmax, kmax, isize_field, isize_wrk1d, inb_wrk1d, isize_wrk2d, inb_wrk2d, isize_wrk3d, inb_txc, isize_txc_field
    use TLAB_VARS, only: visc, schmidt, area
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays, only: wrk1d, wrk2d, txc, x, y, z, wrk3d
    use FDM_ComX_Direct
    use FDM_PROCS
    use FDM_Com1_Jacobian
    use FDM_Com2_Jacobian
    use OPR_PARTIAL
    use IBM_VARS
    use IO_FIELDS
    use OPR_PARTIAL

    implicit none

    ! type(grid_dt) :: g

    integer(wi) :: i, l, len, j

    real(wp), dimension(:), pointer :: u
    real(wp), dimension(:), pointer :: du1_a, du1_b, du1_c, du1_n
    real(wp), dimension(:), pointer :: du2_a, du2_n1, du2_n2, du2_n3
    real(wp), dimension(:, :), pointer :: bcs
    real(wp), dimension(:), pointer    :: epsi, epsj, epsk
    integer(wi) bcs_aux(2, 2)
    real(wp) :: lambda, coef(5), dummy
    integer(wi) :: test_type, ibc, ip
    integer(wi) :: nmin, nmax, nsize
    integer, parameter :: i1 = 1, cases(4) = [BCS_DD, BCS_ND, BCS_DN, BCS_NN]

! ###################################################################
    call TLab_Start()
    call IO_READ_GLOBAL(ifile)
    ! call Particle_Initialize_Parameters(ifile)
    ! call DNS_READ_LOCAL(ifile)
    call IBM_READ_INI(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
! Initialize
    
    len = jmax*kmax

    visc = 1.0_wp   ! Needed in FDM_INITIALIZE
    schmidt = 1.0_wp

    ! g%inb_grid = 16
    ! g%size = imax
    ! g%scale = 1.0_wp
    ! g%uniform = .false.
    ! g(1)%size = imax

    ! isize_wrk1d = max(imax,jmax,kmax)
    ! isize_wrk2d = len
    ! inb_wrk1d = 20
    ! inb_wrk2d = 2
    inb_txc = 12

    call TLab_Initialize_Memory(__FILE__)
    call IBM_ALLOCATE(__FILE__)
    u(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 1)
    du1_a(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 2)
    du1_b(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 3)
    du1_c(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 4)
    du1_n(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 5)
    du2_a(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 6)
    du2_n1(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 7)
    du2_n2(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 8)
    du2_n3(1:imax*jmax*kmax) => txc(1:imax*jmax*kmax, 9)
    epsi => txc(:,10)
    epsj => txc(:,11)
    epsk => txc(:,12)
    allocate (bcs(len, 2))
    ! Valid settings
    test_type = 1

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    lambda = 1 ! WRITE(*,*) 'Eigenvalue ?'; READ(*,*) lambda
    g%mode_fdm1 = FDM_COM6_JACOBIAN ! FDM_COM6_JACOBIAN_PENTA
    ! g%mode_fdm1 = FDM_COM6_DIRECT
    g%mode_fdm2 = g%mode_fdm1
    
    ibm_partial = .true.
    ims_pro_ibm_x = .true.
!  ###################################################################
    call IBM_GENERATE_GEOMETRY_BOX(wrk3d)
    ! write(*,*) 'eps feild', eps
    call IBM_GEOMETRY_TRANSPOSE(epsi, epsj, epsk, wrk3d)
    call IBM_GENERATE_GEOMETRY(epsi, epsj, epsk)

    ! if (g%mode_fdm1 == FDM_COM6_JACOBIAN) C1N6M_ALPHA = 0.56_wp

!  ###################################################################

    ! if (g(1)%periodic) then
    !     do i = 1, imax
    !         ! do j = 1, jmax
    !             ! x(i, j) = real(i - 1, wp)/real(imax, wp)*g(1)%scale
    !             x(i, 1) = real(i - 1, wp)/real(imax, wp)*g(1)%scale
    !         ! end do
    !     end do
    ! else
    !     ! do i = 1, imax
    !     !     x(i, 1) = real(i - 1, wp)/real(imax - 1, wp)*g%scale
    !     ! end do
    !     open (21, file='y.dat')
    !     do i = 1, imax
    !         read (21, *) x(i, 1)
    !     end do
    !     close (21)
    !     g%scale = x(imax, 1) - x(1, 1)
    ! end if

! Bcs
    ! bcs(1, 1) = 0.0_wp
    ! bcs(1, 2) = 0.0_wp
    ! bcs_aux = 0

! ###################################################################
! Define the function and analytic derivatives
    do i = 1, imax
        do l = 1, len
! single-mode
            ! u(l, i) = 1.0_wp + &
            !           sin(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))!+pi_wp/C_4_R)
            ! du1_a(l, i) = (2.0_wp*pi_wp/g%scale*lambda) &
            !               *cos(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))!+pi_wp/C_4_R)
            ! ! u(l, i) = 1.0_wp + &
            ! !           cos(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))!+pi_wp/C_4_R)
            ! ! du1_a(l, i) = -(2.0_wp*pi_wp/g%scale*lambda) &
            ! !               *sin(2.0_wp*pi_wp/g%scale*lambda*g%nodes(i))!+pi_wp/C_4_R)

            ! du2_a(l, i) = -(2.0_wp*pi_wp/g%scale*lambda)**2 &
            !               *u(l, i)

            ! Gaussian
            ! dummy = 1.0_wp / ( 2.0_wp*(g%scale/M_REAL(lambda*l))**2 )
            ! u(l,i)     = EXP(-dummy*(g%nodes(i)-x_0*g%scale)**2)
            ! du1_a(l,i) =-2.0_wp *dummy *(g%nodes(i)-x_0*g%scale) *u(l,i)
            ! du2_a(l,i) =-2.0_wp *dummy *(g%nodes(i)-x_0*g%scale) *du1_a(l,i) - 2.0_wp *dummy *u(l,i)
! Exponential
            ! u(l,i)     = EXP(g%nodes(i)/(g%scale/lambda))
            ! du1_a(l,i) = lambda/g%scale*u(l,i)
            ! du2_a(l,i) = lambda/g%scale*du1_a(l,i)
! delta-function
            ! u(i)     = MAX(0.0_wp,2.0_wp-M_REAL(i))
            ! du1_a(i) = 0.0_wp
            ! du2_a(i) = 0.0_wp
! hyperboic tangent
            ! u(l,i)     = lambda*LOG(1.0_wp+EXP(g%nodes(i)/lambda))
            ! du1_a(l,i) = C_05_R*(1.0_wp+TANH(C_05_R*g%nodes(i)/lambda))
            ! du2_a(l,i) = C_025_R/lambda/(COSH(C_05_R*g%nodes(i)/lambda))**2
! Polynomialx(i, j) = real(i - 1, wp)/real(imax, wp)*g(1)%scale
            ! dummy = 4.0_wp
            ! u(l, i) = ((g%scale - g%nodes(i))/lambda)**dummy
            ! du1_a(l, i) = -dummy*((g%scale - g%nodes(i))/lambda)**(dummy - 1.0_wp)
            ! du2_a(l, i) = dummy*(dummy - 1.0_wp)*((g%scale - g%nodes(i))/lambda)**(dummy - 2.0_wp)
! Trignometrical
            u((l-1)*imax + i)     =                             sin(2*pi_wp*g(1)%nodes(i)/g(1)%scale)
            du1_a((l-1)*imax + i) = (2*pi_wp/g(1)%scale)      * cos(2*pi_wp*g(1)%nodes(i)/g(1)%scale) ! (2*pi_wp/7.5)
            du2_a((l-1)*imax + i) = ((2*pi_wp/g(1)%scale)**2) * sin(2*pi_wp*g(1)%nodes(i)/g(1)%scale)
        end do
    end do

! ###################################################################
    if (test_type == 1) then
        bcs_aux = 0
    ! PRINT *,nobi,nobj,nobk  
! -------------------------------------------------------------------
! Testing first-order derivatives
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g(1), u, du1_b)
        print *, '1. order, Jacobian 4'
        ! call FDM_C1N4_Jacobian(imax, g%jac, g%lu1(:, :), g%rhs1(:, :), coef, g%periodic)
        ! call MatMul_3d_antisym(imax, len, g%rhs1(:, 1), g%rhs1(:, 2), g%rhs1(:, 3), u, du1_n, g%periodic)
        ! call TRIDFS(g%size, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3))
        ! call TRIDSS(g%size, len, g%lu1(:, 1), g%lu1(:, 2), g%lu1(:, 3), du1_n)
        call check(u, du1_a, du1_b, 'partial_Case7.csv')
        ! call check(u, du1_a, du1_b, 'partial_Case4.csv')
        ! call check(u, du1_a, du1_b, 'partial.csv')

        print *, '1. order, Jacobian 6'
        call FDM_C1N6_Jacobian(imax, g(1)%jac, g(1)%lu1(:, :), g(1)%rhs1(:, :), coef, g(1)%periodic)
        call MatMul_5d_antisym(imax, len, g(1)%rhs1(:, 1), g(1)%rhs1(:, 2), g(1)%rhs1(:, 3), g(1)%rhs1(:, 4), g(1)%rhs1(:, 5), u, du1_n, g(1)%periodic)
        call TRIDFS(g(1)%size, g(1)%lu1(:, 1), g(1)%lu1(:, 2), g(1)%lu1(:, 3))
        call TRIDSS(g(1)%size, len, g(1)%lu1(:, 1), g(1)%lu1(:, 2), g(1)%lu1(:, 3), du1_n)
        ! call TRIDPFS(g(1)%size, g(1)%lu1(:, 1), g(1)%lu1(:, 2), g(1)%lu1(:, 3), g(1)%lu1(:, 4), g(1)%lu1(:, 5))
        ! call TRIDPSS(g(1)%size, len, g(1)%lu1(:, 1), g(1)%lu1(:, 2), g(1)%lu1(:, 3), g(1)%lu1(:, 4), g(1)%lu1(:, 5), du1_n, wrk3d)
        ! call check(u, du1_a, du1_n, 'partial.csv')

        ! Direct metrics

! -------------------------------------------------------------------
! Testing second-order derivatives
        ! Jacobian based
        ! call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs_aux, g(1), u, du2_n2, du1_n)
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_aux, g(1), du1_n, du2_n1)
        !

        print *, '2. order, Jacobian 4'
        call FDM_C2N4_Jacobian(imax, g(1)%jac, g(1)%lu2(:, :), g(1)%rhs2(:, :), coef, g(1)%periodic)
        ! call FDM_Bcs(g(1)%lu2(:, 1:3), BCS_DD)
        call TRIDFS(g(1)%size, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3))
        call MatMul_5d_sym(imax, len, g(1)%rhs2(:, 1), g(1)%rhs2(:, 2), g(1)%rhs2(:, 3), g(1)%rhs2(:, 4), g(1)%rhs2(:, 5), u, du2_n2, g(1)%periodic)
        ip = 5
        call MatMul_3d_add(imax, len, g(1)%rhs2(:, ip + 1), g(1)%rhs2(:, ip + 2), g(1)%rhs2(:, ip + 3), du1_a, du2_n2)
        call TRIDSS(g(1)%size, len, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3), du2_n2)
        call check(u, du2_a, du2_n2, 'partial.dat')

        print *, '2. order, Jacobian 6'
        call FDM_C2N6_Jacobian(imax, g(1)%jac, g(1)%lu2(:, :), g(1)%rhs2(:, :), coef, g(1)%periodic)
        ! call FDM_Bcs(g(1)%lu2(:, 1:3), BCS_DD)
        call TRIDFS(g(1)%size, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3))
        call MatMul_5d_sym(imax, len, g(1)%rhs2(:, 1), g(1)%rhs2(:, 2), g(1)%rhs2(:, 3), g(1)%rhs2(:, 4), g(1)%rhs2(:, 5), u, du2_n2, g(1)%periodic)
        ip = 5
        call MatMul_3d_add(imax, len, g(1)%rhs2(:, ip + 1), g(1)%rhs2(:, ip + 2), g(1)%rhs2(:, ip + 3), du1_a, du2_n2)
        call TRIDSS(g(1)%size, len, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3), du2_n2)
        call check(u, du2_a, du2_n2, 'partial.dat')

        print *, '2. order, Jacobian 6 hyper'
        call FDM_C2N6_Hyper_Jacobian(imax, g(1)%jac, g(1)%lu2(:, :), g(1)%rhs2(:, :), coef, g(1)%periodic)
        ! do i = 1,imax
        !     print*,g(1)%lu2(i, 1:3)/g(1)%jac(1,1)/g(1)%jac(1,1)
        !     print*,'rhs',g(1)%rhs2(i, 1:7)
        ! end do
        ! call FDM_Bcs(g(1)%lu2(:, 1:3), BCS_DD)
        call TRIDFS(g(1)%size, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3))
        call MatMul_7d_sym(imax, len, g(1)%rhs2(:, 1), g(1)%rhs2(:, 2), g(1)%rhs2(:, 3), g(1)%rhs2(:, 4), g(1)%rhs2(:, 5), g(1)%rhs2(:, 6), g(1)%rhs2(:, 7), u, du2_n2, g(1)%periodic)
        ip = 7
        call MatMul_3d_add(imax, len, g(1)%rhs2(:, ip + 1), g(1)%rhs2(:, ip + 2), g(1)%rhs2(:, ip + 3), du1_a, du2_n2)
        call TRIDSS(g(1)%size, len, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3), du2_n2)
        ! call TRIDPFS(g(1)%size, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3), g(1)%lu2(:, 4), g(1)%lu2(:, 5))
        ! call TRIDPSS(g(1)%size, len, g(1)%lu2(:, 1), g(1)%lu2(:, 2), g(1)%lu2(:, 3), g(1)%lu2(:, 4), g(1)%lu2(:, 5), du2_n2, wrk3d)
        call check(u, du2_a, du2_n2, 'partial.dat')

        ! Direct metrics
        print *, '2. order, Direct 4'
        call FDM_C2N4_Direct(imax, x, wrk1d(:, 1), wrk1d(:, 4))
        call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
        call MatMul_5d(imax, len, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), u, du2_n2)
        call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du2_n2)
        call check(u, du2_a, du2_n2, 'partial.dat')

        print *, '2. order, Direct 6'
        call FDM_C2N6_Direct(imax, x, wrk1d(:, 1), wrk1d(:, 4))
        call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
        call MatMul_5d(imax, len, wrk1d(:, 4), wrk1d(:, 5), wrk1d(:, 6), wrk1d(:, 7), u, du2_n2)
        call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du2_n2)
        call check(u, du2_a, du2_n2, 'partial.dat')

! ###################################################################
    ! elseif (test_type == 3) then ! Testing new BCs routines
    !     do ip = 1, size(cases)
    !         ibc = cases(ip)
    !         print *, new_line('a'), 'Case ', ibc

    !         nmin = 1
    !         nmax = g(1)%size
    !         if (any([BCS_ND, BCS_NN] == ibc)) then
    !             nmin = nmin + 1
    !         end if
    !         if (any([BCS_DN, BCS_NN] == ibc)) then
    !             nmax = nmax - 1
    !         end if
    !         nsize = nmax - nmin + 1

    !         du1_n(:, 1) = du1_a(:, 1)
    !         du1_n(:, imax) = du1_a(:, imax)

    !         print *, '1. order, Jacobian 4'
    !         call FDM_C1N4_Jacobian(imax, g(1)%jac, g(1)%lu1(:, :), g(1)%rhs1(:, :), coef, g(1)%periodic)
    !         call FDM_Bcs_Neumann(ibc, g(1)%lu1(:, 1:3), g(1)%rhs1(:, 1:3), g(1)%rhs1_b, g(1)%rhs1_t)
    !         call TRIDFS(nsize, g(1)%lu1(nmin:nmax, 1), g(1)%lu1(nmin:nmax, 2), g(1)%lu1(nmin:nmax, 3))
    !     call MatMul_3d_antisym(imax, len, g(1)%rhs1(:, 1), g(1)%rhs1(:, 2), g(1)%rhs1(:, 3), u(:, :), du1_n(:, :), g(1)%periodic, ibc, g(1)%rhs1_b, g(1)%rhs1_t, wrk2d(:,1), wrk2d(:,2))
    !         call TRIDSS(nsize, len, g(1)%lu1(nmin:nmax, 1), g(1)%lu1(nmin:nmax, 2), g(1)%lu1(nmin:nmax, 3), du1_n(:, nmin:nmax))
    !         call check(u(:, nmin:nmax), du1_a(:, nmin:nmax), du1_n(:, nmin:nmax), 'partial.dat')
    !         if (any([BCS_ND, BCS_NN] == ibc)) then
    !             print *, u(:, 1), g(1)%lu1(1, 2)*du1_n(:, 1) + g(1)%lu1(1, 3)*du1_n(:, 2) + wrk2d(1:len, 1)
    !         end if
    !         if (any([BCS_DN, BCS_NN] == ibc)) then
    !             print *, u(:, imax), g(1)%lu1(imax, 1)*du1_n(:, imax - 1) + g(1)%lu1(imax, 2)*du1_n(:, imax) + wrk2d(1:len, 2)
    !         end if

    !         print *, '1. order, Jacobian 6'
    !         call FDM_C1N6_Jacobian(imax, g(1)%jac, g(1)%lu1(:, :), g(1)%rhs1(:, :), coef, g(1)%periodic)
    !         call FDM_Bcs_Neumann(ibc, g(1)%lu1(:, 1:3), g(1)%rhs1(:, 1:5), g(1)%rhs1_b, g(1)%rhs1_t)
    !         call TRIDFS(nsize, g(1)%lu1(nmin:nmax, 1), g(1)%lu1(nmin:nmax, 2), g(1)%lu1(nmin:nmax, 3))
    !     call MatMul_5d_antisym(imax, len, g(1)%rhs1(:, 1), g(1)%rhs1(:, 2), g(1)%rhs1(:, 3), g(1)%rhs1(:, 4), g(1)%rhs1(:, 5), u(:, :), du1_n(:, :), g(1)%periodic, ibc, g(1)%rhs1_b, g(1)%rhs1_t, wrk2d(:,1), wrk2d(:,2))
    !         call TRIDSS(nsize, len, g(1)%lu1(nmin:nmax, 1), g(1)%lu1(nmin:nmax, 2), g(1)%lu1(nmin:nmax, 3), du1_n(:, nmin:nmax))
    !         call check(u(:, nmin:nmax), du1_a(:, nmin:nmax), du1_n(:, nmin:nmax), 'partial.dat')
    !         if (any([BCS_ND, BCS_NN] == ibc)) then
    !             print *, u(:, 1), g(1)%lu1(1, 2)*du1_a(:, 1) + g(1)%lu1(1, 3)*du1_n(:, 2) + wrk2d(1:len, 1)
    !         end if
    !         if (any([BCS_DN, BCS_NN] == ibc)) then
    !             print *, u(:, imax), g(1)%lu1(imax, 1)*du1_n(:, imax - 1) + g(1)%lu1(imax, 2)*du1_a(:, imax) + wrk2d(1:len, 2)
    !         end if
    !     end do
    end if

! ! ###################################################################
!     elseif (test_type == 2) then ! Testing new BCs routines

!         if (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN) then
!             ! call FDM_C1N6_BCS_LHS(imax, ibc, g(1)%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
!             ! call FDM_C1N6_BCS_RHS(imax, len, ibc, u, du1_b)
!         elseif (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
!             call FDM_C1N6M_BCS_LHS(imax, ibc, g(1)%jac, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
!             call FDM_C1N6M_BCS_RHS(imax, len, ibc, u, du1_b)
!         end if

!         if (ibc == 0) then
!             if (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN) then
!                 call TRIDFS(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
!                 call TRIDSS(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du1_b)
!             elseif (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
!                 call PENTADFS2(imax, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
!                 call PENTADSS2(imax, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), du1_b(1, 1))
!             end if

!         else if (ibc == 1) then
!             if (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN) then
!                 wrk1d(:, 4) = 0.0_wp
!                 wrk1d(1, 4) = 1.0_wp
!                 wrk1d(2, 4) = wrk1d(1, 1)
!                 wrk1d(3, 4) = wrk1d(2, 1)
!                 call TRIDFS(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3))
!                 call TRIDSS(imax - 1, len, wrk1d(2:, 1), wrk1d(2:, 2), wrk1d(2:, 3), du1_b(:, 2:))
!                 call TRIDSS(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4))
!                 bcs(:, 1) = du1_b(:, 1)
!                 du1_b(:, 1) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 4) ! BCs
!                 end do
!                 bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 2)*du1_b(:, 1) + wrk1d(1, 3)*du1_b(:, 2))
!                 write (*, *) bcs(:, 1), u(:, 1)
!             elseif (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
!                 wrk1d(:, 6) = 0.0_wp
!                 wrk1d(1, 6) = 1.0_wp
!                 wrk1d(2, 6) = wrk1d(1, 2)
!                 wrk1d(3, 6) = wrk1d(2, 2)
!                 call PENTADFS2(imax - 1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
!                 call PENTADSS2(imax - 1, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), du1_b(1, 2))
!                 call PENTADSS2(imax - 1, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
!                 bcs(:, 1) = du1_b(:, 1)
!                 du1_b(:, 1) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 6) ! BCs
!                 end do
!                 bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 3)*du1_b(:, 1) + wrk1d(1, 4)*du1_b(:, 2))
!                 write (*, *) bcs(:, 1), u(:, 1)
!             end if

!         else if (ibc == 2) then
!             if (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN) then
!                 wrk1d(:, 5) = 0.0_wp
!                 wrk1d(imax, 5) = 1.0_wp
!                 wrk1d(imax - 2, 5) = wrk1d(imax - 1, 3)
!                 wrk1d(imax - 1, 5) = wrk1d(imax, 3)
!                 call TRIDFS(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3))
!                 call TRIDSS(imax - 1, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), du1_b)
!                 call TRIDSS(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 5))
!                 bcs(:, 2) = du1_b(:, imax)
!                 du1_b(:, imax) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 5) ! BCs
!                 end do
!                 bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 1)*du1_b(:, imax - 1) + wrk1d(imax, 2)*du1_b(:, imax))
!                 write (*, *) bcs(:, 2), u(:, imax)
!             elseif (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
!                 wrk1d(:, 6) = 0.0_wp
!                 wrk1d(imax, 6) = 1.0_wp
!                 wrk1d(imax - 2, 6) = wrk1d(imax - 1, 4)
!                 wrk1d(imax - 1, 6) = wrk1d(imax, 4)
!                 call PENTADFS2(imax - 1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5))
!                 call PENTADSS2(imax - 1, len, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), du1_b(1, 1))
!                 call PENTADSS2(imax - 1, i1, wrk1d(1, 1), wrk1d(1, 2), wrk1d(1, 3), wrk1d(1, 4), wrk1d(1, 5), wrk1d(1, 6))
!                 bcs(:, 2) = du1_b(:, imax)
!                 du1_b(:, imax) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 6) ! BCs
!                 end do
!                 bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 2)*du1_b(:, imax - 1) + wrk1d(imax, 3)*du1_b(:, imax))
!                 write (*, *) bcs(:, 2), u(:, imax)
!             end if

!         else if (ibc == 3) then
!             if (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN) then
!                 wrk1d(:, 4) = 0.0_wp
!                 wrk1d(1, 4) = 1.0_wp
!                 wrk1d(2, 4) = wrk1d(1, 1)
!                 wrk1d(3, 4) = wrk1d(2, 1)
!                 wrk1d(:, 5) = 0.0_wp
!                 wrk1d(imax, 5) = 1.0_wp
!                 wrk1d(imax - 2, 5) = wrk1d(imax - 1, 3)
!                 wrk1d(imax - 1, 5) = wrk1d(imax, 3)
!                 call TRIDFS(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3))
!                 call TRIDSS(imax - 2, len, wrk1d(2:, 1), wrk1d(2:, 2), wrk1d(2:, 3), du1_b(:, 2:))
!                 call TRIDSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4))
!                 call TRIDSS(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 5))
!                 bcs(:, 1) = du1_b(:, 1)
!                 du1_b(:, 1) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 4) ! BCs
!                 end do
!                 bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 2)*du1_b(:, 1) + wrk1d(1, 3)*du1_b(:, 2))
!                 write (*, *) bcs(:, 1), u(:, 1)
!                 bcs(:, 2) = du1_b(:, imax)
!                 du1_b(:, imax) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 5) ! BCs
!                 end do
!                 bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 1)*du1_b(:, imax - 1) + wrk1d(imax, 2)*du1_b(:, imax))
!                 write (*, *) bcs(:, 2), u(:, imax)
!             elseif (g(1)%mode_fdm1 == FDM_COM6_JACOBIAN_PENTA) then
!                 wrk1d(:, 6) = 0.0_wp
!                 wrk1d(1, 6) = 1.0_wp
!                 wrk1d(2, 6) = wrk1d(1, 2)
!                 wrk1d(3, 6) = wrk1d(2, 2)
!                 wrk1d(:, 7) = 0.0_wp
!                 wrk1d(imax, 7) = 1.0_wp
!                 wrk1d(imax - 2, 7) = wrk1d(imax - 1, 4)
!                 wrk1d(imax - 1, 7) = wrk1d(imax, 4)
!                 call PENTADFS2(imax - 2, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5))
!                 call PENTADSS2(imax - 2, len, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), du1_b(1, 2))
!                 call PENTADSS2(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 6))
!                 call PENTADSS2(imax - 2, i1, wrk1d(2, 1), wrk1d(2, 2), wrk1d(2, 3), wrk1d(2, 4), wrk1d(2, 5), wrk1d(2, 7))
!                 bcs(:, 1) = du1_b(:, 1)
!                 du1_b(:, 1) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, 1)*wrk1d(i, 6) ! BCs
!                 end do
!                 bcs(:, 1) = bcs(:, 1) + (wrk1d(1, 3)*du1_b(:, 1) + wrk1d(1, 4)*du1_b(:, 2))
!                 write (*, *) bcs(:, 1), u(:, 1)
!                 bcs(:, 2) = du1_b(:, imax)
!                 du1_b(:, imax) = 0.0_wp
!                 do i = 1, imax
!                     du1_b(:, i) = du1_b(:, i) + du1_a(:, imax)*wrk1d(i, 7) ! BCs
!                 end do
!                 bcs(:, 2) = bcs(:, 2) + (wrk1d(imax, 2)*du1_b(:, imax - 1) + wrk1d(imax, 3)*du1_b(:, imax))
!                 write (*, *) bcs(:, 2), u(:, imax)
!             end if
!         end if
!     end if

    stop

    ! ###################################################################
contains
    subroutine check(u, du_a, du_n, name)
        real(wp), intent(in) :: u(:), du_a(:), du_n(:)
        character(len=*), optional :: name
        integer(wi) axis
        real(wp) dummy, error_l2, error_max

        if (present(name)) then
            open (20, file=name)
        end if
        error_l2 = 0.0_wp
        error_max = 0.0_wp
        dummy = 0.0_wp
        axis = 1
        do i = 1, size(u, 1)
            if (present(name)) then
                write (20, 1000) g(axis)%nodes(i), u(i), du_a(i), du_n(i), du_a(i) - du_n(i)
            end if
            dummy = dummy + du_a(i)*du_a(i)
            error_l2 = error_l2 + (du_a(i) - du_n(i))**2.0_wp
            error_max = max(error_max, abs(du_a(i) - du_n(i)))
            write(*,*) ' u: ', u(i) ,  ' du_a: ', du_a(i), ' du_n: ', du_n(i)
        end do
        if (present(name)) then
            close (20)
        end if

        write (*, *) 'Solution L2-norm ...........:', sqrt(g(axis)%jac(1, 1)*dummy)/real(len, wp)
        if (dummy == 0.0_wp) return
        write (*, *) 'Relative Error L2-norm .....:', sqrt(g(axis)%jac(1, 1)*error_l2)/maxval(abs(du1_a))
        write (*, *) 'Relative Error Linf-norm ...:', error_max/maxval(abs(du1_a))

        ! write(*,*) ' u: ', u !,  ' du_a ', du_a, ' du_n ', du_n 
        ! write(*,*) 'du_a:', imax*ims_npro_i,' x ', jmax, ' x ', kmax*ims_npro_k
        ! write(*,*) 'du_n:', imax*ims_npro_i,' x ', jmax, ' x ', kmax*ims_npro_k
        return
1000    format(5(1x, e12.5))
    end subroutine check

end program VPARTIAL
