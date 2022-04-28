#include "types.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/05 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   compute gamma (conditional/intrinsic averaging) [cf. Pope. p.170]
!#   eps field needs to be adjusted because of the treatment of interface
!#   points (Riemann midpoint rule), these points are weighted with 0.5*dz
!# 
!#
!########################################################################
!# ARGUMENTS 
!#                            
!#                           
!########################################################################
!# REQUIREMENTS                 
!#
!#
!########################################################################

subroutine IBM_AVG_GAMMA(gamma, eps, wrk3d, wrk1d)

  use TLAB_VARS, only : imax, jmax, kmax, g, area

  implicit none

  TREAL, dimension(     jmax     ), intent(out  ) ::  gamma
  TREAL, dimension(imax,jmax,kmax), intent(in   ) ::  eps
  TREAL, dimension(imax,jmax,kmax), intent(inout) ::  wrk3d
  TREAL, dimension(     jmax     ), intent(inout) ::  wrk1d

  TINTEGER                                        :: i,j,k

  ! ================================================================== !
  ! exclude body points (filled with zeros)
  wrk3d(:,:,:) = (C_1_R - eps(:,:,:))

  ! in x-direction - inner points - boundary points
  do k = 1, kmax   
    do j = 1, jmax 
      i = 1
      if (        eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
        wrk3d(i+1,j,k) = C_05_R 
      end if
      do i = 2, imax-1
        if (      eps(i,j,k) == 0 .and. eps(i+1,j,k) == 1 ) then
          wrk3d(i+1,j,k) = C_05_R 
        else if ( eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
          wrk3d(i-1,j,k) = C_05_R 
        end if
      end do 
      i = imax
      if (        eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
        wrk3d(i-1,j,k) = C_05_R 
      end if
    end do
  end do

  ! in z-direction - inner points - boundary points
  do i = 1, imax   
    do j = 1, jmax 
      k = 1        
      if (        eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
        wrk3d(i,j,k+1) = C_05_R 
      end if 
      do k = 2, kmax-1
        if (      eps(i,j,k) == 0 .and. eps(i,j,k+1) == 1 ) then
          wrk3d(i,j,k+1) = C_05_R 
        else if ( eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
          wrk3d(i,j,k-1) = C_05_R 
        end if
      end do 
      k = kmax
      if (        eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
        wrk3d(i,j,k-1) = C_05_R 
      end if
    end do
  end do

  ! horizontal average
  CALL AVG_IK_V(imax,jmax,kmax, jmax, wrk3d, g(1)%jac, g(3)%jac, gamma, wrk1d, area)

  return
end subroutine IBM_AVG_GAMMA

!########################################################################