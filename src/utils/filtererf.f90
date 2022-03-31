#include "types.h"
!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2022/02/06 - C. Ansorge
!#              Created           
!#
!# 2022/03/07 - J. Kostelecky
!#              Modified        
!#
!########################################################################
!# DESCRIPTION
!#
!# Spectral filtering of 1d-array with erf-function.
!# In case of R2C dfftw-transform comment out second loop, 
!# here designed for C2C dfftw-transform. 
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE FILTER_ERF_1D(a, n, lcut)

  IMPLICIT NONE
 
#include "integers.h"
 
  TCOMPLEX, DIMENSION(n), INTENT(INOUT) :: a       ! 1d-array to be filtered
  TINTEGER,               INTENT(IN   ) :: n       ! size
  TREAL,                  INTENT(IN   ) :: lcut    ! filter parameter
 
 ! -------------------------------------------------------------------
   
  TREAL                                 :: k_ref, k_rel
  TINTEGER                              :: i, j, n_ny
 
 ! ###################################################################
  n_ny    = n/2+1
  k_ref   = C_2_R * n / lcut
 
  DO i=1,n_ny
    k_rel = (i-1) / k_ref
    a(i)  = a(i) * (ERF(C_8_R * (C_1_R - k_rel)) + C_1_R) / C_2_R
  ENDDO

  i = n_ny+1

  DO j=n_ny-1,n-(n-2),-1 
    k_rel = (j-1) / k_ref
    a(i)  = a(i) * (ERF(C_8_R * (C_1_R - k_rel)) + C_1_R) / C_2_R
    i     = i + 1
  ENDDO

  RETURN
END SUBROUTINE FILTER_ERF_1D