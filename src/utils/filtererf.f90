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
subroutine FILTER_ERF_1D(a, n, lcut)

    implicit none

    TCOMPLEX, dimension(n), intent(INOUT) :: a       ! 1d-array to be filtered
    TINTEGER, intent(IN) :: n       ! size
    TREAL, intent(IN) :: lcut    ! filter parameter

    ! -------------------------------------------------------------------

    TREAL :: k_ref, k_rel
    TINTEGER :: i, j, n_ny

    ! ###################################################################
    n_ny = n/2 + 1
    k_ref = C_2_R*n/lcut

    do i = 1, n_ny
        k_rel = (i - 1)/k_ref
        a(i) = a(i)*(erf(C_8_R*(C_1_R - k_rel)) + C_1_R)/C_2_R
    end do

    i = n_ny + 1

    do j = n_ny - 1, n - (n - 2), -1
        k_rel = (j - 1)/k_ref
        a(i) = a(i)*(erf(C_8_R*(C_1_R - k_rel)) + C_1_R)/C_2_R
        i = i + 1
    end do

    return
end subroutine FILTER_ERF_1D
