#include "types.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/05 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   horizontal staggering of geometry field eps to epsp
!#    
!# 
!########################################################################
!# ARGUMENTS 
!#  epsilon field 'eps' is an indicator field:
!#    epsp(i,j,k) = 0 for fluid domain
!#    epsp(i,j,k) = 1 for solid domain
!#                            
!#                           
!########################################################################
!# REQUIREMENTS                 
!#
!#
!########################################################################

subroutine IBM_STAGGER_GEOMETRY(eps, epsp)

  use TLAB_VARS, only : imax, jmax, kmax

  implicit none

  TREAL, dimension(imax,jmax,kmax), intent(inout) ::  eps, epsp

  TINTEGER                                        :: i,j,k

  ! ================================================================== !
   
  ! horizontal staggering - move right interfaces one step to the left,
  ! don't touch vertical direction since no staggering is applied here
 
  epsp(:,:,:) = eps(:,:,:)
  
  do k = 1, kmax   
    do j = 1, jmax 
      do i = 2, imax
        if ( eps(i,j,k) == 0 .and. eps(i-1,j,k) == 1 ) then
          epsp(i-1,j,k) = C_0_R 
        end if
      end do 
    end do
  end do
  
  do i = 1, imax   
    do j = 1, jmax 
      do k = 2, kmax
        if ( eps(i,j,k) == 0 .and. eps(i,j,k-1) == 1 ) then
          epsp(i,j,k-1) = C_0_R 
        end if
      end do 
    end do
  end do

  return
end subroutine IBM_STAGGER_GEOMETRY

!########################################################################