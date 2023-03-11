#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif
program VHELMHOLTZ_FXZ

    use TLAB_TYPES, only: pointers_dt
    use TLAB_VARS, only: imax, jmax, kmax, inb_wrk1d, inb_wrk2d, isize_wrk1d, isize_wrk2d, gfile, isize_txc_field
    use TLAB_PROCS
    use IO_FIELDS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use OPR_FOURIER
    use OPR_PARTIAL
    use OPR_ELLIPTIC

    implicit none

    TREAL, dimension(:, :), allocatable, save, target :: x, y, z
    TREAL, dimension(:, :, :), allocatable :: b, c, d, h
    TREAL, dimension(:, :, :, :), allocatable :: a, f, e
    TREAL, dimension(:, :), allocatable :: txc
    TREAL, dimension(:, :), allocatable :: wrk1d, wrk2d
    TREAL, dimension(:, :, :), allocatable :: bcs_hb, bcs_ht
    TREAL, dimension(:), allocatable :: cx, cy, cz, wrk3d

    TINTEGER i, j, k, ibc_x(4), ibc_y(4), ibc_z(4), ifield, imeasure
    TINTEGER nfield, nmeasure
    parameter(nfield=4, nmeasure=1)
    TINTEGER opt
    TREAL dummy, error, rms, max_error, beta, t_new, t_old
    character*8 :: date
    character*10 :: time1, time2
    TINTEGER :: h1, m1, s1, n1, h2, m2, s2, n2
    type(pointers_dt), dimension(nfield) :: data

    target a

#ifdef USE_MPI
#else
    TINTEGER :: ims_pro
    parameter(ims_pro=0)
#endif

! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL('dns.ini')
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    isize_wrk3d = imax*jmax*kmax
    isize_wrk3d = max(isize_wrk3d, isize_txc_field)

! --------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (x(g(1)%size, g(1)%inb_grid))
    allocate (y(g(2)%size, g(2)%inb_grid))
    allocate (z(g(3)%size, g(3)%inb_grid))

    allocate (wrk1d(isize_wrk1d, 4*nfield + 7))
    allocate (wrk2d(isize_wrk2d, inb_wrk2d))
    allocate (bcs_ht(imax, kmax, nfield), bcs_hb(imax, kmax, nfield))
    allocate (a(imax, jmax, kmax, nfield), b(imax, jmax, kmax), c(imax, jmax, kmax))
    allocate (d(imax, jmax, kmax), e(imax, jmax, kmax, nfield), f(imax, jmax, kmax, nfield), h(imax, jmax, kmax))
    allocate (txc(isize_txc_field, nfield + 1), wrk3d(isize_wrk3d))
    allocate (cx(6*imax), cy(6*jmax), cz(6*g(3)%size))

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call OPR_FOURIER_INITIALIZE()
    call OPR_FOURIER_INITIALIZE()
! ###################################################################
! Define forcing term
! ###################################################################

    beta = -C_1_R/0.05
    ! Velocity BCs
    bcs_ht(:, :, 1) = C_1_R
    bcs_ht(:, :, 2:3) = C_0_R
    bcs_hb(:, :, :) = C_0_R

    ! Scalar BCs
    if (nfield > 3) then
        bcs_hb(:, :, 4:nfield) = C_1_R
        bcs_ht(:, :, 4:nfield) = C_0_R
    end if

    t_new = 0.0
    t_old = 0.0
    do opt = 1, 2
        call IO_READ_FIELDS('field.inp', IO_FLOW, imax, jmax, kmax, nfield, 0, a)
        f = a
        if (opt == 1) then
            do imeasure = 1, nmeasure
                a = f

                do ifield = 1, nfield
                    data(ifield)%field => a(:, 1, 1, ifield)
                end do

                call date_and_time(date, time1)
                call OPR_HELMHOLTZ_FXZ_D_N(imax, jmax, kmax, h, nfield, i0, beta, &
                                           data, txc(1, 1), txc(1, nfield + 1), &
                                           bcs_hb(1, 1, 1), bcs_ht(1, 1, 1))
                call date_and_time(date, time2)
                read (time1(1:10), '(i2,i2,i2,1x,i3)') h1, m1, s1, n1
                read (time2(1:10), '(i2,i2,i2,1x,i3)') h2, m2, s2, n2
                t_new = t_new + (3600*(h2 - h1) + 60*(m2 - m1) + (s2 - s1)) + (n2 - n1)/1.0e+03
            end do
        else if (opt == 2) then
            do imeasure = 1, nmeasure
                a = f
                call date_and_time(date, time1)
                do ifield = 1, nfield
                    call OPR_HELMHOLTZ_FXZ_D(imax, jmax, kmax, h, 0, beta, &
                                             a(1, 1, 1, ifield), txc(1, 1), txc(1, 2), &
                                             bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield))
                end do
                call date_and_time(date, time2)
                read (time1(1:10), '(i2,i2,i2,1x,i3)') h1, m1, s1, n1
                read (time2(1:10), '(i2,i2,i2,1x,i3)') h2, m2, s2, n2
                t_old = t_old + (3600*(h2 - h1) + 60*(m2 - m1) + (s2 - s1)) + (n2 - n1)/1.0e+03
            end do
        end if
#ifdef USE_MPI
        dummy = t_new
        call MPI_Reduce(dummy, t_new, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        dummy = t_old
        call MPI_Reduce(dummy, t_old, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        if (ims_pro == 0) then
#endif
            if (opt == 1) then
                write (*, *) 'NEW SOLVER', t_new/nmeasure, 'ms'
            else if (opt == 2) then
                write (*, *) 'OLD SOLVER', t_old/nmeasure, 'ms'
            end if
#ifdef USE_MPI
        end if
#endif
        ! normalize solution
        a(:, :, :, :) = a(:, :, :, :)*beta

        ! ###################################################################
        ! Error
        ! ###################################################################
        do ifield = 1, nfield
            call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), a(1, 1, 1, ifield), b, h)
            call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), a(1, 1, 1, ifield), c, h)
            call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), a(1, 1, 1, ifield), d, h)

            error = C_0_R
            max_error = C_0_R
            rms = C_0_R

            do k = 1, kmax
                do j = 2, jmax - 1
                    do i = 1, imax
                        b(i, j, k) = a(i, j, k, ifield) + (b(i, j, k) + c(i, j, k) + d(i, j, k))/beta
                        e(i, j, k, ifield) = b(i, j, k) - f(i, j, k, ifield)
                        error = error + e(i, j, k, ifield)*e(i, j, k, ifield)
                        rms = rms + a(i, j, k, ifield)*a(i, j, k, ifield)

                        if (error > max_error) max_error = error
                    end do
                end do
            end do
#ifdef USE_MPI
            dummy = error
            call MPI_Reduce(dummy, error, i1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
            dummy = rms
            call MPI_Reduce(dummy, rms, i1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
            dummy = max_error
            call MPI_Reduce(dummy, max_error, i1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
            if (ims_pro == 0) then
#endif
                write (*, *) ims_pro, ifield, 'Relative error .............: ', sqrt(error)/sqrt(rms)
                write (*, *) ims_pro, ifield, 'Maximum error ..............: ', max_error
#ifdef USE_MPI
            end if
#endif
        end do
    end do

    call TLAB_STOP(0)
end program VHELMHOLTZ_FXZ
