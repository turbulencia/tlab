#include "types.h"
#include "dns_error.h"
!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2022/03/07 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Spectral filtering of 1d-array with erf-function.
!# In case of R2C dfftw-transform comment out second loop,
!# here designed for C2C dfftw.
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
subroutine OPR_STAGGERING_INITIALIZE(aux)

    use TLAB_CONSTANTS, only: efile
    use TLAB_VARS, only: g
    use OPR_FOURIER, only: fft_plan_fy1d, fft_plan_by1d
    use TLAB_PROCS

    implicit none

#include "integers.h"
#ifdef USE_FFTW
#include "fftw3.f"
#endif

    TCOMPLEX, dimension(g(2)%size, 2), intent(INOUT) :: aux
! ###################################################################

! FFTW library
#ifdef USE_FFTW

! Build 1d-plans
#ifdef _DEBUG
    call dfftw_plan_dft_1d(fft_plan_fy1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(fft_plan_by1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_BACKWARD, FFTW_ESTIMATE)
#else
    call dfftw_plan_dft_1d(fft_plan_fy1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_FORWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(fft_plan_by1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_BACKWARD, FFTW_MEASURE)
#endif

#else
    call TLAB_WRITE_ASCII(efile, 'OPR_STAGGERING_INITIALIZE. FFTW needed for vertical pressure filtering.')
    call TLAB_STOP(DNS_ERROR_UNDEVELOP)
#endif

    return
end subroutine OPR_STAGGERING_INITIALIZE
! ###################################################################
! ###################################################################
subroutine FILTER_VERTICAL_PRESSURE(a, b, n, lcut, wrk1)

    use OPR_FOURIER, only: fft_plan_fy1d, fft_plan_by1d

    implicit none

#include "integers.h"
#ifdef USE_FFTW
#include "fftw3.f"
#endif

    TCOMPLEX, dimension(n), intent(INOUT) :: a, b     ! 1d-arrays to be filtered
    ! here: p and dpdy
    TINTEGER, intent(INOUT) :: n       ! ny
    TREAL, intent(INOUT) :: lcut    ! filter parameter
    TCOMPLEX, dimension(n, 2), intent(INOUT) :: wrk1    ! aux 1d-arrays

! -------------------------------------------------------------------

    TREAL :: norm

! ###################################################################
    norm = C_1_R/M_REAL(n)

! Execute + Filtering
    call dfftw_execute_dft(fft_plan_fy1d, a, wrk1(1, 1))
    call dfftw_execute_dft(fft_plan_fy1d, b, wrk1(1, 2))

    call FILTER_ERF_1D(wrk1(1, 1), n, lcut)
    call FILTER_ERF_1D(wrk1(1, 2), n, lcut)

    call dfftw_execute_dft(fft_plan_by1d, wrk1(1, 1), a)
    call dfftw_execute_dft(fft_plan_by1d, wrk1(1, 2), b)

! Normalize
    a = a*norm
    b = b*norm

    return
end subroutine FILTER_VERTICAL_PRESSURE
