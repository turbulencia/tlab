#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

PROGRAM VFFT

USE DNS_GLOBAL, ONLY:  imax,jmax,kmax
USE DNS_GLOBAL, ONLY:  isize_txc_dimz
USE TLAB_CORE
#ifdef USE_MPI
USE TLAB_MPI_PROCS
#endif

IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif


TREAL, DIMENSION(:), ALLOCATABLE    :: trans, trans2, trans_ref
TREAL, DIMENSION(:), ALLOCATABLE    :: tmp1,tmp2,tmp3,tmp4,wrk2d,wrk3d

TINTEGER                            :: i,j,k,ip,ip_ref,bad_count,good_count,check_mode,bad
TINTEGER                            :: isize_fft3d, isize_trn3d
TINTEGER                            :: err_count, case_count
TREAL                               :: norm

err_count = i0
case_count= i0

CALL TLAB_START()
CALL DNS_READ_GLOBAL('dns.ini')
#ifdef USE_MPI
CALL DNS_MPI_INITIALIZE
#endif

isize_fft3d = isize_txc_dimz*kmax
isize_trn3d = (imax/2)* jmax   *(2*kmax)

ALLOCATE( &
     trans (isize_trn3d),  &
     trans_ref(isize_trn3d), &
     trans2(isize_trn3d), &
     tmp1 (isize_fft3d),  &
     tmp2 (isize_fft3d),  &
     tmp3 (isize_fft3d),  &
     tmp4 (isize_fft3d),  &
     wrk3d(isize_fft3d),  &
     wrk2d(isize_fft3d)  )

CALL OPR_FOURIER_INITIALIZE(tmp1,tmp2,wrk2d,wrk3d)


CALL FFT_CHECK(2,err_count,case_count,&
     trans,     &
     trans_ref, &
     tmp1 ,     &
     tmp2 ,     &
     tmp3 ,     &
     tmp4 ,     &
     wrk3d,     &
     wrk2d )

CALL FFT_CHECK(1,err_count,case_count, &
     trans,     &
     trans_ref, &
     tmp1 ,     &
     tmp2 ,     &
     tmp3 ,     &
     tmp4 ,     &
     wrk3d,     &
     wrk2d  )

CALL FFT_CHECK(3,err_count,case_count,&
     trans,     &
     trans_ref, &
     tmp1 ,     &
     tmp2 ,     &
     tmp3 ,     &
     tmp4 ,     &
     wrk3d,     &
     wrk2d )


#ifdef USE_MPI
IF (ims_pro .EQ. 0 ) THEN
#endif
   WRITE(*,1000) err_count, case_count
1000 FORMAT('fft-check completed. ', I3, ' Errors in ', I3, ' Checks. For details see file dns.log')
#ifdef USE_MPI
ENDIF
#endif

CALL TLAB_STOP(0)
END PROGRAM VFFT
