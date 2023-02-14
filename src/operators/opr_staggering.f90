#include "types.h"
#include "dns_error.h"
!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2023/10/02 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Vertical midpoint pressure filter (p and dpdy, boundaries untouched)
!#
!########################################################################
subroutine FILTER_VERTICAL_PRESSURE(ny, vfilter_param, a, b)

    use TLAB_CONSTANTS, only: wp, wi

    implicit none
    
    integer(wi),                intent(in   ) :: ny
    real(wp),                   intent(in   ) :: vfilter_param
    complex(wp), dimension(ny), intent(inout) :: a, b ! here: p and dpdy      

    integer(wi)                               :: j
    real(wp)                                  :: w1, w2
    complex(wp), dimension(2,2)               :: buf
    ! -------------------------------------------------------------------

    ! filter parameters (weights)
    if      (vfilter_param == 0) then ! max filter strength (potentially unstable)
        w2 = 0.0_wp; w1 = 0.5_wp
    else if (vfilter_param == 2) then ! filter is off
        w2 = 1.0_wp; w1 = 0.0_wp
    else                              ! vfilter_param = [0,2]
        w2 =  vfilter_param / 2.0_wp
        w1 = (1.0_wp - w2)  / 2.0_wp
    end if

    ! initialize
    buf(1,1) = a(1); buf(1,2) = b(1)   
    buf(2,1) = a(2); buf(2,2) = b(2)  
    
    ! buffered filter, boundaries untouched
    do j = 2, ny-1
        a(j) = w1*buf(1,1) + w2*buf(2,1) + w1*a(j+1) 
        b(j) = w1*buf(1,2) + w2*buf(2,2) + w1*b(j+1) 
        
        ! reorder
        buf(1,1) = buf(2,1); buf(1,2) = buf(2,2)
        buf(2,1) = a(j+1);   buf(2,2) = b(j+1)
    end do

    return
end subroutine FILTER_VERTICAL_PRESSURE
!########################################################################
!# Vertical spectral pressure filter with erf-cut-off function
!#
!# only for quasi periodic BCs in the vertical, e.g. closed channel flow 
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
! subroutine OPR_STAGGERING_INITIALIZE(aux)

!     use TLAB_CONSTANTS, only: efile
!     use TLAB_VARS, only: g
!     use OPR_FOURIER, only: fft_plan_fy1d, fft_plan_by1d
!     use TLAB_PROCS

!     implicit none

! #include "integers.h"
! #ifdef USE_FFTW
! #include "fftw3.f"
! #endif

!     TCOMPLEX, dimension(g(2)%size, 2), intent(INOUT) :: aux
! ! ###################################################################

! ! FFTW library
! #ifdef USE_FFTW

! ! Build 1d-plans
! #ifdef _DEBUG
!     call dfftw_plan_dft_1d(fft_plan_fy1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_FORWARD, FFTW_ESTIMATE)
!     call dfftw_plan_dft_1d(fft_plan_by1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_BACKWARD, FFTW_ESTIMATE)
! #else
!     call dfftw_plan_dft_1d(fft_plan_fy1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_FORWARD, FFTW_MEASURE)
!     call dfftw_plan_dft_1d(fft_plan_by1d, g(2)%size, aux(1, 1), aux(1, 2), FFTW_BACKWARD, FFTW_MEASURE)
! #endif

! #else
!     call TLAB_WRITE_ASCII(efile, 'OPR_STAGGERING_INITIALIZE. FFTW needed for vertical pressure filtering.')
!     call TLAB_STOP(DNS_ERROR_UNDEVELOP)
! #endif

!     return
! end subroutine OPR_STAGGERING_INITIALIZE
! ! ###################################################################
! subroutine FILTER_VERTICAL_PRESSURE(a, b, n, lcut, wrk1)

!     use OPR_FOURIER, only: fft_plan_fy1d, fft_plan_by1d

!     implicit none

! #include "integers.h"
! #ifdef USE_FFTW
! #include "fftw3.f"
! #endif

!     TCOMPLEX, dimension(n), intent(INOUT) :: a, b     ! 1d-arrays to be filtered
!     ! here: p and dpdy
!     TINTEGER, intent(INOUT) :: n       ! ny
!     TREAL, intent(INOUT) :: lcut    ! filter parameter
!     TCOMPLEX, dimension(n, 2), intent(INOUT) :: wrk1    ! aux 1d-arrays

! ! -------------------------------------------------------------------

!     TREAL :: norm

! ! ###################################################################
!     norm = C_1_R/M_REAL(n)

! ! Execute + Filtering
!     call dfftw_execute_dft(fft_plan_fy1d, a, wrk1(1, 1))
!     call dfftw_execute_dft(fft_plan_fy1d, b, wrk1(1, 2))

!     call FILTER_ERF_1D(wrk1(1, 1), n, lcut)
!     call FILTER_ERF_1D(wrk1(1, 2), n, lcut)

!     call dfftw_execute_dft(fft_plan_by1d, wrk1(1, 1), a)
!     call dfftw_execute_dft(fft_plan_by1d, wrk1(1, 2), b)

! ! Normalize
!     a = a*norm
!     b = b*norm

!     return
! end subroutine FILTER_VERTICAL_PRESSURE