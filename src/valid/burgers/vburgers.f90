#include "types.h"
#include "dns_const.h"

program VBURGERS

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
    use TLAB_MPI_VARS
#endif
    use IO_FIELDS
    use OPR_PARTIAL
    use OPR_BURGERS
    implicit none

#include "integers.h"
#ifdef USE_MPI
    TREAL error2, dummy2
#else
    TINTEGER, parameter :: ims_pro = 0
#endif

    TREAL, dimension(:, :), allocatable, save, target :: x, y, z
    TREAL, dimension(:, :, :), allocatable :: a, b, c
    TREAL, dimension(:, :), allocatable :: wrk1d, wrk2d
    TREAL, dimension(:), allocatable :: wrk3d, tmp1

    TINTEGER i, j, k, bcs(2, 2)
    TREAL dummy, error

! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL('dns.ini')
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    isize_wrk3d = isize_txc_field

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (x(g(1)%size, g(1)%inb_grid))
    allocate (y(g(2)%size, g(2)%inb_grid))
    allocate (z(g(3)%size, g(3)%inb_grid))

    allocate (wrk1d(isize_wrk1d, inb_wrk1d + 1))
    allocate (wrk2d(isize_wrk2d, inb_wrk2d))
    allocate (a(imax, jmax, kmax), b(imax, jmax, kmax), c(imax, jmax, kmax))
    allocate (tmp1(isize_txc_field), wrk3d(isize_wrk3d))

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE()

    bcs = 0

! ###################################################################
! Define forcing term
! ###################################################################
    call IO_READ_FIELDS('field.inp', IO_SCAL, imax, jmax, kmax, i1, i0, a, wrk3d)

! ###################################################################
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), a, b, c, wrk2d, wrk3d)
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
!           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
                b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
            end do
        end do
    end do

    call OPR_BURGERS_X(i0, i0, imax, jmax, kmax, bcs, g(1), a, a, a, c, tmp1)

    c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
    end if
! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)

! ###################################################################
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), a, b, c, wrk2d, wrk3d)
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
!           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
                b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
            end do
        end do
    end do

    call OPR_BURGERS_Y(i0, i0, imax, jmax, kmax, bcs, g(2), a, a, a, c, tmp1)

    c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
    error2 = error; dummy2 = dummy
    call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
    call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
    if (ims_pro == 0) then
        write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
    end if
! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)

! ###################################################################
    if (g(3)%size > 1) then

        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), a, b, c, wrk2d, wrk3d)
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    !           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
                    b(i, j, k) = b(i, j, k)*visc*ribackground(j) - a(i, j, k)*c(i, j, k)
                end do
            end do
        end do

        call OPR_BURGERS_Z(i0, i0, imax, jmax, kmax, bcs, g(3), a, a, a, c, tmp1)

        c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
        error2 = error; dummy2 = dummy
        call MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        call MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
        if (ims_pro == 0) then
            write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
        end if
!    CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)

    end if

    call TLAB_STOP(0)
end program VBURGERS
