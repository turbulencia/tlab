!########################################################################
!# Tool/Library SUPERLAYER
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/14 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Sampling the fields given in array b along the surface given by sl
!# Data stored in array c
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
subroutine SL_BOUNDARY_SAMPLE(imax, jmax, kmax, nfield_loc, nfield, y, sl, b, c)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, nfield_loc, nfield
    real(wp) y(jmax)
    real(wp) sl(imax, kmax)
    real(wp) b(imax, jmax, kmax, nfield_loc)
    real(wp) c(nfield, imax, kmax)

! -------------------------------------------------------------------
    real(wp) dy_u, dy_loc
    integer(wi) i, k, ifield, jm

! ###################################################################
    dy_u = y(2) - y(1)

! ###################################################################
! Loop on the points
! ###################################################################
    do k = 1, kmax
        do i = 1, imax
            jm = INT((sl(i, k) - y(1))/dy_u) + 1
            dy_loc = sl(i, k) - y(jm)
            do ifield = 1, nfield_loc
                c(ifield, i, k) = b(i, jm, k, ifield) + &
                                  (b(i, jm + 1, k, ifield) - b(i, jm, k, ifield))/(y(jm + 1) - y(jm))*dy_loc
            end do
        end do
    end do

    return
end subroutine SL_BOUNDARY_SAMPLE
