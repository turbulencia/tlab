#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI

#endif

program VFFT
    use TLab_Constants, only: wp, wi
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_Memory, only: isize_txc_dimz
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use OPR_FOURIER
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_PROCS, only: TLabMPI_Initialize
use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif

    implicit none

    real(wp), dimension(:), allocatable :: trans, trans2, trans_ref
    real(wp), dimension(:), allocatable :: tmp1, tmp2, tmp3, tmp4, wrk2d, wrk3d

    integer(wi) :: i, j, k, ip, ip_ref, bad_count, good_count, check_mode, bad
    integer(wi) :: isize_fft3d, isize_trn3d
    integer(wi) :: err_count, case_count
    real(wp) :: norm

    err_count = i0
    case_count = i0

    call TLab_Start()
    call TLab_Initialize_Parameters('tlab.ini')
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
call TLabMPI_Transpose_Initialize(ifile)
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
