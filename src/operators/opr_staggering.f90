#include "types.h"
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
SUBROUTINE OPR_STAGGERING_INITIALIZE(aux)

  USE TLAB_CONSTANTS, ONLY : efile
  ! USE TLAB_TYPES,     ONLY : grid_dt
  USE TLAB_VARS,      ONLY : g
  USE TLAB_VARS,      ONLY : fft_plan_fy1d, fft_plan_by1d
  USE TLAB_PROCS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_FFTW
#include "fftw3.f"
#endif

  ! TYPE(grid_dt),                    INTENT(IN)    :: g(3)
  TCOMPLEX, DIMENSION(g(2)%size,2), INTENT(INOUT) :: aux  
! ###################################################################

! FFTW library
#ifdef USE_FFTW

! Build 1d-plans 
#ifdef _DEBUG
  CALL dfftw_plan_dft_1d(fft_plan_fy1d, g(2)%size, aux(1,1), aux(1,2), FFTW_FORWARD,  FFTW_ESTIMATE)
  CALL dfftw_plan_dft_1d(fft_plan_by1d, g(2)%size, aux(1,1), aux(1,2), FFTW_BACKWARD, FFTW_ESTIMATE)
#else
  CALL dfftw_plan_dft_1d(fft_plan_fy1d, g(2)%size, aux(1,1), aux(1,2), FFTW_FORWARD,  FFTW_MEASURE)
  CALL dfftw_plan_dft_1d(fft_plan_by1d, g(2)%size, aux(1,1), aux(1,2), FFTW_BACKWARD, FFTW_MEASURE)
#endif

#else
  CALL TLAB_WRITE_ASCII(efile,'OPR_STAGGERING_INITIALIZE. FFTW needed for vertical pressure filtering.')
  CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
#endif

  RETURN
END SUBROUTINE OPR_STAGGERING_INITIALIZE
! ###################################################################
! ###################################################################
SUBROUTINE FILTER_VERTICAL_PRESSURE(a, b, n, lcut, wrk1)

  USE TLAB_VARS, ONLY : fft_plan_fy1d, fft_plan_by1d

  IMPLICIT NONE
  
#include "integers.h"
#ifdef USE_FFTW
#include "fftw3.f"
#endif 
  
  TCOMPLEX, DIMENSION(n  ), INTENT(INOUT) :: a,b     ! 1d-arrays to be filtered
                                                     ! here: p and dpdy
  TINTEGER,                 INTENT(INOUT) :: n       ! ny
  TREAL,                    INTENT(INOUT) :: lcut    ! filter parameter
  TCOMPLEX, DIMENSION(n,2), INTENT(INOUT) :: wrk1    ! aux 1d-arrays
  
! -------------------------------------------------------------------

  TREAL                                   :: norm 
  
! ###################################################################
  norm = C_1_R / M_REAL(n)
  
! Execute + Filtering
  CALL dfftw_execute_dft(fft_plan_fy1d, a, wrk1(1,1))
  CALL dfftw_execute_dft(fft_plan_fy1d, b, wrk1(1,2))

  CALL FILTER_ERF_1D(wrk1(1,1), n, lcut)
  CALL FILTER_ERF_1D(wrk1(1,2), n, lcut)

  CALL dfftw_execute_dft(fft_plan_by1d, wrk1(1,1), a)
  CALL dfftw_execute_dft(fft_plan_by1d, wrk1(1,2), b)

! Normalize
  a = a * norm
  b = b * norm

  RETURN
END SUBROUTINE FILTER_VERTICAL_PRESSURE