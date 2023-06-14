#include "dns_const.h"

#define C_FILE_LOC "VPOISSON"

! We can check for
!       lap a = f, given f,
! or for
!       lap a = f, given a (constructing first f and then solving)

program VPOISSON
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_PROCS
    use TLAB_ARRAYS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use IO_FIELDS
    use OPR_FILTERS
    use OPR_FOURIER
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use AVGS

    implicit none

    real(wp), dimension(:, :), allocatable :: bcs_hb, bcs_ht
    real(wp), dimension(:, :, :), pointer :: a, b, c, d, e, f
    real(wp) mean, lambda
    ! real(wp) SIMPSON_NU, delta

    integer(wi) i, j, k, ig, bcs(2, 2)
    integer(wi) type_of_operator, type_of_problem
    integer ibc

! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    isize_wrk3d = isize_txc_field
    inb_txc = 8

    call TLAB_ALLOCATE(__FILE__)

    allocate (bcs_ht(imax, kmax), bcs_hb(imax, kmax))

    a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    d(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    e(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 7)
    f(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 8)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call OPR_ELLIPTIC_INITIALIZE()

! Staggering of the pressure grid not implemented here
    if (stagger_on) then
        call TLAB_WRITE_ASCII(wfile, C_FILE_LOC//'. Staggering of the pressure grid not implemented here.')
        stagger_on = .false. ! turn staggering off for OPR_POISSON_FXZ(...)
    end if
    if (any(PressureFilter%type /= DNS_FILTER_NONE)) then
        call TLAB_WRITE_ASCII(wfile, C_FILE_LOC//'. Pressure and dpdy Filter not implemented here.')
    end if

    bcs = 0

    call OPR_FOURIER_INITIALIZE()
    call OPR_CHECK()

    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
    end do

    type_of_operator = 1   ! Poisson routines
    ! type_of_operator = 2   ! Helmholtz routines
    if (type_of_operator == 2) then
        write (*, *) 'Eigenvalue ?'
        read (*, *) lambda
    end if

    ! type_of_problem = 1     ! the forcing in the rhs is given
    type_of_problem = 2     ! the field in the lhs is given

    select case (type_of_problem)
! ###################################################################
    case (1) ! The input field f is used as rhs in lap a = f and we solve for a

        call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, f)
        ! ! remove 2\Delta x wave
        ! call OPR_FILTER(imax, jmax, kmax, Dealiasing, f, txc)

        a = f
        bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
        ! For Neumann conditions, we need to satisfy the compatibility constraint dpdy_top-dpdy_bottom=int f
        ! mean = AVG_IK(imax, 1, kmax, 1, bcs_hb, g(1)%jac, g(3)%jac, area)
        ! call AVG_IK_V(imax, jmax, kmax, jmax, a, g(1)%jac, g(3)%jac, wrk1d(:, 1), wrk1d(:, 2), area)
        ! delta = mean + SIMPSON_NU(jmax, wrk1d, g(2)%nodes)
        ! mean = AVG_IK(imax, 1, kmax, 1, bcs_ht, g(1)%jac, g(3)%jac, area)
        ! bcs_ht = bcs_ht - mean + delta

        if (type_of_operator == 1) then
            call OPR_POISSON_FXZ(imax, jmax, kmax, g, 3, a, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht, d)
            call OPR_POISSON_FXZ_D(imax, jmax, kmax, g, 3, a, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht, d)

        else if (type_of_operator == 2) then
            ! call OPR_HELMHOLTZ_FXZ(imax, jmax, kmax, g, 0, lambda, a, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht)
            call OPR_HELMHOLTZ_FXZ_D(imax, jmax, kmax, g, 0, lambda, a, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht)

        end if

        ! -------------------------------------------------------------------
        ! With the calculated a, we calculate the b = lap a
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, c)
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), c, b)
        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, txc(:, 1))

        ! call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, c)
        ! call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), c, e)
        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), a, e, txc(:, 1))
        b = b + e

        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, c)
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), c, e)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, e, txc(:, 1))
        b = b + e

        if (type_of_operator == 2) then
            b = b + lambda*a
        end if

        ! -------------------------------------------------------------------
        b(:, 1, :) = f(:, 1, :); b(:, jmax, :) = f(:, jmax, :)  ! The boundary points do not satisfy the PDE
        call IO_WRITE_FIELDS('field.out', IO_SCAL, imax, jmax, kmax, 1, b)
        call check(f, b, txc(:, 1), 'field.dif')

! ###################################################################
    case (2) ! The input field a is used to construct the forcing term as lap a = f
        call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, a)
        ! ! remove 2\Delta x wave
        ! call OPR_FILTER(imax, jmax, kmax, Dealiasing, f, txc)
        ! lambda = 4.0
        ! do j = 1, jmax
        !     ! a(:,j,:) = sin(2.0_wp*pi_wp/g(2)%scale*lambda*g(2)%nodes(j))!+pi_wp/C_4_R)
        !     a(:,j,:) = exp(lambda*g(2)%nodes(j))
        ! end do
        

        ! -------------------------------------------------------------------
        ! DO j = 1,jmax
        !    mean = AVG_IK(imax,jmax,kmax, j, a, dx,dz, area)
        !    a(:,j,:) = a(:,j,:) - mean
        ! ENDDO

        ! DC level at lower boundary set to zero
        mean = AVG_IK(imax, jmax, kmax, 1, a, g(1)%jac, g(3)%jac, area)
        a = a - mean

        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), a, c)
        ! call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), c, b)
        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c)

        ! call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, c)
        ! call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), c, d)
        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), a, d, c)
        b = b + d

        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), a, c)
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), c, d)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, d, c)
        b = b + d

        if (type_of_operator == 2) then
            b = b + lambda*a
        end if

        ibc = BCS_NN
        select case (ibc)
        case (BCS_DD)
            bcs_hb(:, :) = a(:, 1, :); bcs_ht(:, :) = a(:, jmax, :)
        case (BCS_DN)
            bcs_hb(:, :) = a(:, 1, :); bcs_ht(:, :) = c(:, jmax, :)
            ! bcs_ht(:, :) = (11.0_wp*a(:, jmax, :)-18.0_wp*a(:, jmax-1, :)+9.0_wp*a(:, jmax-2, :)-2.0_wp*a(:, jmax-3, :))/6.0_wp/g(2)%jac(jmax,1)
        case (BCS_ND)
            bcs_hb(:, :) = c(:, 1, :); bcs_ht(:, :) = a(:, jmax, :)
            ! bcs_hb(:, :) = (-3.0_wp*a(:, 1, :)+4.0_wp*a(:, 2, :)-a(:, 3, :))/2.0_wp/g(2)%jac(1,1)
            ! bcs_hb(:, :) = (-11.0_wp*a(:, 1, :) + 18.0_wp*a(:, 2, :) - 9.0_wp*a(:, 3, :) + 2.0_wp*a(:, 4, :))/6.0_wp/g(2)%jac(1, 1)
        case (BCS_NN)
            bcs_hb(:, :) = c(:, 1, :); bcs_ht(:, :) = c(:, jmax, :)
            ! bcs_hb(:, :) = (-3.0_wp*a(:, 1, :)+4.0_wp*a(:, 2, :)-a(:, 3, :))/2.0_wp/g(2)%jac(1,1)
            ! bcs_hb(:, :) = (-11.0_wp*a(:, 1, :) + 18.0_wp*a(:, 2, :) - 9.0_wp*a(:, 3, :) + 2.0_wp*a(:, 4, :))/6.0_wp/g(2)%jac(1, 1)
            ! bcs_ht(:, :) = (11.0_wp*a(:, jmax, :)-18.0_wp*a(:, jmax-1, :)+9.0_wp*a(:, jmax-2, :)-2.0_wp*a(:, jmax-3, :))/6.0_wp/g(2)%jac(jmax,1)
        end select

        if (type_of_operator == 1) then
            ! call OPR_POISSON_FXZ(imax, jmax, kmax, g, ibc, b, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht, d)
            ! call OPR_POISSON_FXZ_D_TRANSPOSE(imax, jmax, kmax, g, ibc, b, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht, d)
            call OPR_POISSON_FXZ_D(imax, jmax, kmax, g, ibc, b, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht, d)

        else if (type_of_operator == 2) then
            ! call OPR_HELMHOLTZ_FXZ(imax, jmax, kmax, g, ibc, lambda, b, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht)
            call OPR_HELMHOLTZ_FXZ_D(imax, jmax, kmax, g, ibc, lambda, b, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), b, d)

        end if

        ! -------------------------------------------------------------------
        call IO_WRITE_FIELDS('field.out', IO_SCAL, imax, jmax, kmax, 1, b)
        call check(a, b, txc(:, 1), 'field.dif')

        call check(c, d, txc(:, 1))

    end select

    stop

! ###################################################################
contains
    subroutine check(a1, a2, dif, name)
        real(wp), intent(in) :: a1(imax, jmax, kmax), a2(imax, jmax, kmax)
        real(wp), intent(inout) :: dif(imax, jmax, kmax)
        character(len=*), optional :: name

        real(wp) dummy, error

        error = 0.0_wp
        dummy = 0.0_wp
        do k = 1, kmax
            do j = 1, jmax
                ! do j = 2, jmax - 1
                do i = 1, imax
                    dif(i, j, k) = a2(i, j, k) - a1(i, j, k)
                    error = error + dif(i, j, k)*dif(i, j, k)
                    dummy = dummy + a1(i, j, k)*a1(i, j, k)
                end do
            end do
        end do
        write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
        if (present(name)) then
            call IO_WRITE_FIELDS(name, IO_SCAL, imax, jmax, kmax, 1, dif)
        end if

        return
    end subroutine check

end program VPOISSON
