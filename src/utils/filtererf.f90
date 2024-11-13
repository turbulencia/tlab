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
    use TLab_Constants, only: wp, wi
    implicit none

    complex(wp), dimension(n), intent(INOUT) :: a       ! 1d-array to be filtered
    integer(wi), intent(IN) :: n       ! size
    real(wp), intent(IN) :: lcut    ! filter parameter

    ! -------------------------------------------------------------------

    real(wp) :: k_ref, k_rel
    integer(wi) :: i, j, n_ny

    ! ###################################################################
    n_ny = n/2 + 1
    k_ref = 2.0_wp*n/lcut

    do i = 1, n_ny
        k_rel = (i - 1)/k_ref
        a(i) = a(i)*(erf(8.0_wp*(1.0_wp - k_rel)) + 1.0_wp)/2.0_wp
    end do

    i = n_ny + 1

    do j = n_ny - 1, n - (n - 2), -1
        k_rel = (j - 1)/k_ref
        a(i) = a(i)*(erf(8.0_wp*(1.0_wp - k_rel)) + 1.0_wp)/2.0_wp
        i = i + 1
    end do

    return
end subroutine FILTER_ERF_1D
