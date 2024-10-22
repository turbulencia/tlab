#include "dns_const.h"

program VBURGERS

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLAB_ARRAYS
    use TLAB_POINTERS_3D, only: tmp1
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
    use TLabMPI_VARS
#endif
    use IO_FIELDS
    use OPR_PARTIAL
    use OPR_BURGERS
    use OPR_FILTERS
    implicit none

#ifdef USE_MPI
    real(wp) error2, dummy2
#else
    integer(wi), parameter :: ims_pro = 0
#endif

    real(wp), dimension(:, :, :), pointer :: a, b, c

    integer(wi) i, j, k, ig, bcs(2, 2)
    real(wp) dummy, error

! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif

    inb_txc = 4

    call TLab_Initialize_Memory(__FILE__)

    a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 2)
    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)

    visc = 1.0_wp/big_wp    ! inviscid

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE()

    bcs = 0

    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
    end do

! ###################################################################
! Define forcing term
! ###################################################################
    call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, 1, 0, a)

    visc = 1.0_wp/big_wp

! ###################################################################
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c)
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
                ! b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
            end do
        end do
    end do
    call IO_WRITE_FIELDS('fieldXdirect.out', IO_SCAL, imax, jmax, kmax, 1, b)

    call OPR_BURGERS_X(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(1), a, a, c, tmp1)
    call IO_WRITE_FIELDS('fieldXburgers.out', IO_SCAL, imax, jmax, kmax, 1, c)

    c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Relative error X ...........: ', sqrt(error)/sqrt(dummy)
    end if
    call IO_WRITE_FIELDS('fieldX.dif', IO_SCAL, imax, jmax, kmax, 1, c)

! ###################################################################
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, b, c)
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
                ! b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
            end do
        end do
    end do
    call IO_WRITE_FIELDS('fieldYdirect.out', IO_SCAL, imax, jmax, kmax, 1, b)

    call OPR_BURGERS_Y(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(2), a, a, c, tmp1)
    call IO_WRITE_FIELDS('fieldYburgers.out', IO_SCAL, imax, jmax, kmax, 1, c)

    c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Relative error Y ...........: ', sqrt(error)/sqrt(dummy)
    end if
    call IO_WRITE_FIELDS('fieldY.dif', IO_SCAL, imax, jmax, kmax, 1, c)

! ###################################################################
    if (g(3)%size > 1) then

        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), a, b, c)
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    b(i, j, k) = b(i, j, k)*visc - a(i, j, k)*c(i, j, k)
                    ! b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
                end do
            end do
        end do
        call IO_WRITE_FIELDS('fieldZdirect.out', IO_SCAL, imax, jmax, kmax, 1, b)

        call OPR_BURGERS_Z(OPR_B_SELF, 0, imax, jmax, kmax, bcs, g(3), a, a, c, tmp1)
        call IO_WRITE_FIELDS('fieldZburgers.out', IO_SCAL, imax, jmax, kmax, 1, c)

        c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
        error2 = error; dummy2 = dummy
        call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
        if (ims_pro == 0) then
            write (*, *) 'Relative error Z ...........: ', sqrt(error)/sqrt(dummy)
        end if
        call IO_WRITE_FIELDS('fieldZ.dif', IO_SCAL, imax, jmax, kmax, 1, c)

    end if

    call TLAB_STOP(0)
end program VBURGERS
