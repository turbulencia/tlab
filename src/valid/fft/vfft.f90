#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

program VFFT

    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: isize_txc_dimz
    use TLab_WorkFlow
    use OPR_FOURIER
#ifdef USE_MPI
    use MPI
    use TLabMPI_PROCS
#endif

    implicit none

    TREAL, dimension(:), allocatable :: trans, trans2, trans_ref
    TREAL, dimension(:), allocatable :: tmp1, tmp2, tmp3, tmp4, wrk2d, wrk3d

    TINTEGER :: i, j, k, ip, ip_ref, bad_count, good_count, check_mode, bad
    TINTEGER :: isize_fft3d, isize_trn3d
    TINTEGER :: err_count, case_count
    TREAL :: norm

    err_count = i0
    case_count = i0

    call TLab_Start()
    call TLab_Initialize_Parameters('tlab.ini')
#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif
    call NavierStokes_Initialize_Parameters(ifile)

    isize_fft3d = isize_txc_dimz*kmax
    isize_trn3d = (imax/2)*jmax*(2*kmax)

    allocate ( &
        trans(isize_trn3d), &
        trans_ref(isize_trn3d), &
        trans2(isize_trn3d), &
        tmp1(isize_fft3d), &
        tmp2(isize_fft3d), &
        tmp3(isize_fft3d), &
        tmp4(isize_fft3d), &
        wrk3d(isize_fft3d), &
        wrk2d(isize_fft3d))

    call OPR_FOURIER_INITIALIZE()

    call FFT_CHECK(2, err_count, case_count, &
                   trans, &
                   trans_ref, &
                   tmp1, &
                   tmp2, &
                   tmp3, &
                   tmp4, &
                   wrk3d, &
                   wrk2d)

    call FFT_CHECK(1, err_count, case_count, &
                   trans, &
                   trans_ref, &
                   tmp1, &
                   tmp2, &
                   tmp3, &
                   tmp4, &
                   wrk3d, &
                   wrk2d)

    call FFT_CHECK(3, err_count, case_count, &
                   trans, &
                   trans_ref, &
                   tmp1, &
                   tmp2, &
                   tmp3, &
                   tmp4, &
                   wrk3d, &
                   wrk2d)

#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        write (*, 1000) err_count, case_count
1000    format('fft-check completed. ', I3, ' Errors in ', I3, ' Checks. For details see file tlab.log')
#ifdef USE_MPI
    end if
#endif

    call TLab_Stop(0)
end program VFFT
