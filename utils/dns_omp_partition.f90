#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2011/11/01 - C. Ansorge
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Defining loops for individual threads
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE DNS_OMP_PARTITION(len, omp_srt, omp_end, omp_siz) 
#ifdef USE_OPENMP
  USE OMP_LIB
#endif 

#ifdef USE_OPENMP 
  USE DNS_GLOBAL, ONLY : dns_omp_numThreads 
#endif

  IMPLICIT NONE 

  TINTEGER, INTENT(IN)    :: len 
  TINTEGER, INTENT(INOUT) :: omp_srt, omp_end, omp_siz

! -----------------------------------------------------------------------
#ifdef USE_OPENMP 
  TINTEGER                :: residual 
  TINTEGER                :: dns_omp_threadID
#endif

! #######################################################################
#ifdef USE_OPENMP 
  dns_omp_threadID = omp_get_thread_num() 

  IF ( len .GT. dns_omp_numThreads ) THEN
     residual = MOD(len, dns_omp_numThreads)
     IF ( dns_omp_threadID .GE. residual ) THEN
        omp_siz = len / dns_omp_numThreads 
        omp_srt = dns_omp_threadID * (len / dns_omp_numThreads) + residual + 1
        omp_end = omp_srt + omp_siz -1
     ELSE
        omp_siz = len / dns_omp_numThreads + 1
        omp_srt = dns_omp_threadID * (len / dns_omp_numThreads + 1 ) + 1 
        omp_end = omp_srt + omp_siz -1
     ENDIF
  ELSE 
     IF ( dns_omp_threadID .EQ. 0 ) THEN 
        omp_siz = len
        omp_srt = 1
        omp_end = len
     ELSE
        omp_srt = 1
        omp_end = -1
        omp_siz = -1
     ENDIF
  END IF

#else 
  omp_siz = len 
  omp_srt = 1 
  omp_end = len

#endif

  RETURN 

END SUBROUTINE DNS_OMP_PARTITION
