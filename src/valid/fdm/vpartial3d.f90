#include "dns_const.h"

program VPARTIAL3D
    use TLab_Constants, only: wp, wi, gfile, ifile
    use TLab_Time, only: itime
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: imax, jmax, kmax, inb_txc
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use FDM, only: g, FDM_Initialize
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use IO_Grid
    use IO_FIELDS
    use OPR_FOURIER
    use OPR_PARTIAL
    use OPR_ELLIPTIC
    use Averages

    implicit none

    real(wp), dimension(:, :), allocatable :: bcs_hb, bcs_ht
    real(wp), dimension(:, :, :), pointer :: a, b, c, d, e, f

    integer(wi) i, j, k, bcs(2, 2)
    integer(wi) type_of_problem
    real(wp) params(0)

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

    call IO_READ_GRID(gfile, x, y, z, [g(1)%size, g(2)%size, g(3)%size])
    call FDM_Initialize(g(1), x)
    call FDM_Initialize(g(2), y)
    call FDM_Initialize(g(3), z)

    bcs = 0

    ! type_of_problem = 1     ! 1. order derivative
    type_of_problem = 2     ! 2. order derivative

    call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, f, params)

    select case (type_of_problem)
! ###################################################################
    case (1) ! to be done

! ###################################################################
    case (2)
        g(2)%mode_fdm1 = FDM_COM6_JACOBIAN
        call FDM_Initialize(g(2), wrk1d(:, 2))
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), f, c)
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), c, a)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), f, a, c)
        call IO_Write_Fields('field.out1', imax, jmax, kmax, itime, 1, a)

        g(2)%mode_fdm1 = FDM_COM4_DIRECT
        call FDM_Initialize(g(2), wrk1d(:, 2))
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), f, d)
        ! call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), d, b)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), f, b, d)
        call IO_Write_Fields('field.out2', imax, jmax, kmax, itime, 1, b)

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
            call IO_Write_Fields(name, imax, jmax, kmax, itime, 1, dif)
        end if

        return
    end subroutine check

end program VPARTIAL3D
