#include "dns_const.h"

program VPARTIAL3D
    use TLab_Constants
    use TLAB_VARS
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
#endif
    use FDM, only: g,  FDM_Initialize
    use IO_FIELDS
    use OPR_FILTERS
    use OPR_FOURIER
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use Averages

    implicit none

    real(wp), dimension(:, :), allocatable :: bcs_hb, bcs_ht
    real(wp), dimension(:, :, :), pointer :: a, b, c, d, e, f

    integer(wi) i, j, k, bcs(2, 2)
    integer(wi) type_of_problem

! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
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

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:,1), wrk1d(:,2), wrk1d(:,3))
    call FDM_Initialize(x, g(1), wrk1d(:,1), wrk1d(:,4))
    call FDM_Initialize(y, g(2), wrk1d(:,2), wrk1d(:,4))
    call FDM_Initialize(z, g(3), wrk1d(:,3), wrk1d(:,4))

    bcs = 0

    ! type_of_problem = 1     ! 1. order derivative
    type_of_problem = 2     ! 2. order derivative

    call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, f)

    select case (type_of_problem)
! ###################################################################
    case (1) ! to be done

! ###################################################################
    case (2)
        g(2)%mode_fdm1 = FDM_COM6_JACOBIAN
        call FDM_Initialize(y, g(2), wrk1d(:,2), wrk1d(:,4))
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), f, c)
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), c, a)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), f, a, c)
        call IO_WRITE_FIELDS('field.out1', IO_SCAL, imax, jmax, kmax, 1, a)

        g(2)%mode_fdm1 = FDM_COM4_DIRECT
        call FDM_Initialize(y, g(2), wrk1d(:,2), wrk1d(:,4))
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), f, d)
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), d, b)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), f, b, d)
        call IO_WRITE_FIELDS('field.out2', IO_SCAL, imax, jmax, kmax, 1, b)

        ! -------------------------------------------------------------------
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

end program VPARTIAL3D
