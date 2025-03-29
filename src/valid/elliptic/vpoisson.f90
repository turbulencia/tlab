#include "dns_const.h"

#define C_FILE_LOC "VPOISSON"

! We can check for
!       lap a = f, given f,
! or for
!       lap a = f, given a (constructing first f and then solving)

program VPOISSON
    use TLab_Constants, only: wp, wi, pi_wp, BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: ifile, gfile, wfile
    use TLab_Time, only: itime
    use TLab_Memory, only: imax, jmax, kmax, inb_txc
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, stagger_on
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use FDM, only: g, FDM_Initialize
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use Tlab_Background, only: TLab_Initialize_Background
    use TLab_Grid
    use IO_Fields
    use OPR_PARTIAL
    use OPR_FOURIER
    use OPR_FILTERS
    use OPR_ELLIPTIC
    use Averages

    implicit none

    real(wp), dimension(:, :), allocatable :: bcs_hb, bcs_ht
    real(wp), dimension(:, :, :), pointer :: a, b, c, d, e, f
    real(wp) mean, lambda, params(0)
    ! real(wp) Int_Simpson, delta

    integer(wi) i, j, k, ig, bcs(2, 2)
    integer(wi) type_of_operator, type_of_problem
    integer ibc, ib, bcs_cases(4)

! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Transpose_Initialize(ifile)
#endif
    call NavierStokes_Initialize_Parameters(ifile)

    inb_txc = 8

    call TLab_Initialize_Memory(__FILE__)

    allocate (bcs_ht(imax, kmax), bcs_hb(imax, kmax))

    a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)
    d(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 6)
    e(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 7)
    f(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 8)

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call TLab_Initialize_Background(ifile)

    call OPR_Elliptic_Initialize(ifile)

! Staggering of the pressure grid not implemented here
    if (stagger_on) then
        call TLab_Write_ASCII(wfile, C_FILE_LOC//'. Staggering of the pressure grid not implemented here.')
        stagger_on = .false. ! turn staggering off for OPR_Poisson_FourierXZ_Factorize(...)
    end if
    if (any(PressureFilter%type /= DNS_FILTER_NONE)) then
        call TLab_Write_ASCII(wfile, C_FILE_LOC//'. Pressure and dpdy Filter not implemented here.')
    end if

    bcs = 0

    call OPR_FOURIER_INITIALIZE()
    call OPR_CHECK()

    print*, '1. Poisson routines'
    print*, '2. Helmholtz routines'
    read(*,*)  type_of_operator

    if (type_of_operator == 2) then
        write (*, *) 'Eigenvalue (negative)?'
        read (*, *) lambda
        if (lambda > 0.0_wp) then
            print *, 'Eigenvalue for Helmholtz operator needs to be negative'
            stop
        end if
    end if

    ! type_of_problem = 1     ! the forcing in the rhs is given
    type_of_problem = 2     ! the field in the lhs is given

    select case (type_of_problem)
! ###################################################################
    case (1) ! The input field f is used as rhs in lap a = f and we solve for a

        call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, f, params)
        ! ! remove 2\Delta x wave
        ! call OPR_FILTER(imax, jmax, kmax, Dealiasing, f, txc)

        a = f
        ! For Neumann conditions, we need to satisfy the compatibility constraint dpdy_top-dpdy_bottom=int f
        ! mean = AVG_IK(imax, 1, kmax, 1, bcs_hb)
        ! call AVG_IK_V(imax, jmax, kmax, jmax, a, wrk1d(:, 1), wrk1d(:, 2))
        ! delta = mean + Int_Simpson(wrk1d(1:jmax,1), g(2)%nodes(1:jmax))
        ! mean = AVG_IK(imax, 1, kmax, 1, bcs_ht)
        ! bcs_ht = bcs_ht - mean + delta
        ibc = BCS_NN
        select case (ibc)
        case (BCS_DD)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
        case (BCS_DN)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
        case (BCS_ND)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
        case (BCS_NN)
            bcs_hb(:, :) = 0.0_wp; bcs_ht(:, :) = 0.0_wp
        end select

        if (type_of_operator == 1) then
            call OPR_Poisson(imax, jmax, kmax, g, ibc, a, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht, d)

        else if (type_of_operator == 2) then
            call OPR_Helmholtz(imax, jmax, kmax, g, ibc, lambda, a, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht)

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
        call IO_Write_Fields('field.out', imax, jmax, kmax, itime, 1, b, io_header_s(1:1))
        call check(f, b, txc(:, 1), 'field.dif')

! ###################################################################
    case (2) ! The input field a is used to construct the forcing term as lap a = f
        ! Reading field
        call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, a, params)
        ! call random_seed()
        ! call random_number(a)

        ! remove 2\Delta x wave
        do ig = 1, 3
            call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig))
        end do
        call OPR_FILTER(imax, jmax, kmax, FilterDomain, a, txc)

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

        ! ! Creating a field
        ! lambda = 1.0
        ! do j = 1, jmax
        !     a(:, j, :) = sin(2.0_wp*pi_wp/g(2)%scale*lambda*g(2)%nodes(j))
        !     ! a(:,j,:) = exp(lambda*g(2)%nodes(j))
        ! end do
        ! b = -(2.0_wp*pi_wp/g(2)%scale*lambda)**2.0*a
        ! c = (2.0_wp*pi_wp/g(2)%scale*lambda)*cos(2.0_wp*pi_wp/g(2)%scale*lambda*g(2)%nodes(j))

        ! -------------------------------------------------------------------
        ! DC level at lower boundary set to zero
        mean = AVG_IK(imax, jmax, kmax, 1, a)
        a = a - mean

        if (type_of_operator == 2) then
            b = b + lambda*a
        end if

        bcs_cases(1:4) = [BCS_DD, BCS_NN, BCS_DN, BCS_ND]

        do ib = 1, 2
            ibc = bcs_cases(ib)
            print *, new_line('a')

            select case (ibc)
            case (BCS_DD)
                print *, 'Dirichlet/Dirichlet'
                bcs_hb(:, :) = a(:, 1, :); bcs_ht(:, :) = a(:, jmax, :)
            case (BCS_DN)
                print *, 'Dirichlet/Neumann'
                bcs_hb(:, :) = a(:, 1, :); bcs_ht(:, :) = c(:, jmax, :)
            case (BCS_ND)
                print *, 'Neumann/Dirichlet'
                bcs_hb(:, :) = c(:, 1, :); bcs_ht(:, :) = a(:, jmax, :)
            case (BCS_NN)
                print *, 'Neumann/Neumann'
                bcs_hb(:, :) = c(:, 1, :); bcs_ht(:, :) = c(:, jmax, :)
            end select

            e = b       ! to save b for other cases in the loop
            if (type_of_operator == 1) then
                call OPR_Poisson(imax, jmax, kmax, g, ibc, e, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht, d)

            else if (type_of_operator == 2) then
                call OPR_Helmholtz(imax, jmax, kmax, g, ibc, lambda, e, txc(1, 1), txc(1, 2), bcs_hb, bcs_ht)
                call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), e, d)

            end if

            ! -------------------------------------------------------------------
            io_datatype = IO_TYPE_SINGLE
            call IO_Write_Fields('field.out', imax, jmax, kmax, itime, 1, e)!, io_header_s(1:1))
            call check(a, e, txc(:, 1), 'field.dif')

            call check(c, d, txc(:, 1))

        end do

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
            call IO_Write_Fields(name, imax, jmax, kmax, itime, 1, dif) !, io_header_s(1:1))
        end if

        return
    end subroutine check

end program VPOISSON
