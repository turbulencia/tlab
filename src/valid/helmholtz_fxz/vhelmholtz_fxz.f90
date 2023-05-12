#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif
program VHELMHOLTZ_FXZ
    use TLAB_CONSTANTS
    use TLAB_TYPES, only: pointers_dt
    use TLAB_VARS
    use TLAB_PROCS
    use TLAB_ARRAYS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS
#endif
    use IO_FIELDS
    use OPR_FOURIER
    use OPR_PARTIAL
    use OPR_ELLIPTIC

    implicit none

    real(wp), dimension(:, :, :), pointer :: b, c, d, e
    real(wp), dimension(:, :, :, :), allocatable :: a, f
    real(wp), dimension(:, :, :), allocatable :: bcs_hb, bcs_ht

    integer(wi), parameter :: nfield = 1, nmeasure = 1
    integer(wi) i, j, k, ifield, imeasure, bcs(2, 2)
    integer(wi) opt
    real(wp) error, rms, max_error, beta, t_new, t_old
    character*8 :: date
    character*10 :: time1, time2
    integer(wi) :: h1, m1, s1, n1, h2, m2, s2, n2
    type(pointers_dt), dimension(nfield) :: data

    target a

#ifdef USE_MPI
    dummy
#else
    integer(wi), parameter :: ims_pro = 0
#endif

! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    isize_wrk3d = imax*jmax*kmax
    isize_wrk3d = max(isize_wrk3d, isize_txc_field)
    inb_txc = max(5, nfield + 1)

    call TLAB_ALLOCATE(__FILE__)

    allocate (a(imax, jmax, kmax, nfield), f(imax, jmax, kmax, nfield))
    allocate (bcs_ht(imax, kmax, nfield), bcs_hb(imax, kmax, nfield))

    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 2)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    d(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    e(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call OPR_FOURIER_INITIALIZE()

    bcs = 0

! ###################################################################
! Define forcing term
! ###################################################################

    beta = -1.0_wp/0.05_wp
    ! Velocity BCs
    bcs_hb(:, :, :) = 0.0_wp
    bcs_ht(:, :, :) = 1.0_wp
    bcs_ht(:, :, 2:) = 0.0_wp

    ! Scalar BCs
    if (nfield > 3) then
        bcs_hb(:, :, 4:nfield) = 1.0_wp
        bcs_ht(:, :, 4:nfield) = 0.0_wp
    end if

    t_new = 0.0_wp
    t_old = 0.0_wp
    do opt = 1, 2
        call IO_READ_FIELDS('field.inp', IO_FLOW, imax, jmax, kmax, nfield, 0, f)
        if (opt == 1) then
            do imeasure = 1, nmeasure
                a = f

                do ifield = 1, nfield
                    data(ifield)%field => a(:, 1, 1, ifield)
                end do

                call date_and_time(date, time1)
                call OPR_HELMHOLTZ_FXZ_D_N(imax, jmax, kmax, nfield, g, 0, beta, &
                                           data, txc(1, 1), txc(1, nfield + 1), &
                                           bcs_hb(1, 1, 1), bcs_ht(1, 1, 1))
                call date_and_time(date, time2)
                read (time1(1:10), '(i2,i2,i2,1x,i3)') h1, m1, s1, n1
                read (time2(1:10), '(i2,i2,i2,1x,i3)') h2, m2, s2, n2
                t_new = t_new + (3600.0_wp*(h2 - h1) + 60.0_wp*(m2 - m1) + (s2 - s1)) + (n2 - n1)/1.0e+03_wp
            end do
        else if (opt == 2) then
            do imeasure = 1, nmeasure
                a = f
                call date_and_time(date, time1)
                do ifield = 1, nfield
                    call OPR_HELMHOLTZ_FXZ_D(imax, jmax, kmax, g, 0, beta, &
                                             a(1, 1, 1, ifield), txc(1, 1), txc(1, 2), &
                                             bcs_hb(1, 1, ifield), bcs_ht(1, 1, ifield))
                end do
                call date_and_time(date, time2)
                read (time1(1:10), '(i2,i2,i2,1x,i3)') h1, m1, s1, n1
                read (time2(1:10), '(i2,i2,i2,1x,i3)') h2, m2, s2, n2
                t_old = t_old + (3600.0_wp*(h2 - h1) + 60.0_wp*(m2 - m1) + (s2 - s1)) + (n2 - n1)/1.0e+03_wp
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
            call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), a(1, 1, 1, ifield), b, txc(:, 1))
            call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), a(1, 1, 1, ifield), c, txc(:, 1))
            call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), a(1, 1, 1, ifield), d, txc(:, 1))

            error = 0.0_wp
            max_error = 0.0_wp
            rms = 0.0_wp

            do k = 1, kmax
                do j = 2, jmax - 1
                    do i = 1, imax
                        b(i, j, k) = a(i, j, k, ifield) + (b(i, j, k) + c(i, j, k) + d(i, j, k))/beta
                        e(i, j, k) = b(i, j, k) - f(i, j, k, ifield)
                        error = error + e(i, j, k)*e(i, j, k)
                        rms = rms + a(i, j, k, ifield)*a(i, j, k, ifield)

                        if (error > max_error) max_error = error
                    end do
                end do
            end do
#ifdef USE_MPI
            dummy = error
            call MPI_Reduce(dummy, error, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
            dummy = rms
            call MPI_Reduce(dummy, rms, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
            dummy = max_error
            call MPI_Reduce(dummy, max_error, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
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
