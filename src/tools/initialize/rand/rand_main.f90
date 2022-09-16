#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "INIRAND"

program INIRAND
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_PROCS
#endif
    use RAND_LOCAL
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_pro
#endif
    use IO_FIELDS

    implicit none

    ! -------------------------------------------------------------------
    integer(wi) iq, is

    ! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call RAND_READ_LOCAL(ifile)
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    isize_wrk3d = isize_txc_field

    inb_txc = 3

    call TLAB_ALLOCATE(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE(wrk1d)

    ! ###################################################################
    call TLAB_WRITE_ASCII(lfile, 'Initializing random fiels.')

#ifdef USE_MPI
    seed = seed + ims_pro         ! seed for random generator
#endif
    seed = -ABS(seed)

    if (ifourier == 1) then
        call OPR_FOURIER_INITIALIZE(txc, wrk1d, wrk2d, wrk3d)
    end if

    call OPR_CHECK(imax, jmax, kmax, q, txc, wrk2d, wrk3d)

    itime = 0; rtime = 0.0_wp

    do iq = 1, inb_flow
        call RAND_FIELD(ucov(iq), q(1, iq), txc(1, 1), txc(1, 2), txc(1, 3), wrk2d, wrk3d)
    end do
    if (ipdf == 2) then ! Gaussian PDF
        call RAND_COVARIANCE(ucov, q(:, 1), q(:, 2), q(:, 3))
    end if
    call IO_WRITE_FIELDS('flow.rand', IO_FLOW, imax, jmax, kmax, inb_flow, q, txc)

    do is = 1, inb_scal
        call RAND_FIELD(ucov(is), s(1, is), txc(1, 1), txc(1, 2), txc(1, 3), wrk2d, wrk3d)
    end do
    call IO_WRITE_FIELDS('scal.rand', IO_SCAL, imax, jmax, kmax, inb_scal, s, txc)

    call TLAB_STOP(0)
end program INIRAND
