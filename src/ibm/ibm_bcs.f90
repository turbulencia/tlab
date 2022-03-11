#include "types.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#  
!#  
!#    
!#    
!#
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!#                            
!#                           
!#
!########################################################################
subroutine IBM_INITIALIZE_SCAL(s)
  
  use DNS_IBM,   only : eps, xbars_geo
  use DNS_IBM,   only : ibmscaljmin, ibmscaljmax 
  use TLAB_VARS, only : isize_field,inb_scal
  use TLAB_VARS, only : imax,jmax,kmax

  implicit none

#include "integers.h"
  
  TREAL, dimension(isize_field,inb_scal), intent(inout) :: s

  TINTEGER                                              :: is, j, k, ip

! ================================================================== !
! get scalar dirichlet boundary values of ini scalar field
! (assuming always homogenous horizontal temperature on upper/lower boundaries)
  do is = 1, inb_scal
    ibmscaljmin(is) = s(1,                is)
    ibmscaljmax(is) = s(imax*(jmax-1) + 1,is)
  end do

! set scalar values in solid to zero
  call IBM_BCS_FIELD_COMBINED(i1,s)

! default, set only scalar value on lower boundary
! if objects also on upper boundary present, 
! values needs to be overwritten, cf. next loop
  do is = 1, inb_scal
    s(:,is) = s(:,is) + eps(:) * ibmscaljmin(is) 
  end do

! in case of objects on upper boundary, set different temperature here
  if ( xbars_geo%mirrored ) then
    do is = 1, inb_scal
      ip = imax*(jmax-xbars_geo%height) + 1
      do k = 1,kmax
        do j = 1, xbars_geo%height   
          s(ip:ip+imax,is) = (C_1_R - eps(ip:ip+imax)) * s(ip:ip+imax,is) + eps(ip:ip+imax)*ibmscaljmax(is)
          ip  = ip + imax
        end do
        ip = imax*(jmax-xbars_geo%height) + k*imax*jmax
      end do
    end do
  end if

  return
end subroutine IBM_INITIALIZE_SCAL
!########################################################################
subroutine IBM_BCS_FIELD_COMBINED(is,fld)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : isize_field, inb_flow, inb_scal

  implicit none
  
  TINTEGER,                        intent(in)    :: is
  TREAL, dimension(isize_field,*), intent(inout) :: fld

  TINTEGER                                       :: i

  ! ================================================================== !
  ! apply IBM BCs on many fields

  if ( is == 0 ) then
    do i = 1, inb_flow
      call IBM_BCS_FIELD(fld(1,i))
    end do
  elseif ( is == 1 ) then
    do i = 1, inb_scal
      call IBM_BCS_FIELD(fld(1,i))
    end do
  end if

  return
end subroutine IBM_BCS_FIELD_COMBINED
!########################################################################
subroutine IBM_BCS_FIELD(fld)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : isize_field

  implicit none
  
  TREAL, dimension(isize_field), intent(inout) :: fld

  TINTEGER                                     :: i

  ! ================================================================== !
  ! apply IBM BCs on scalar/flow fields
  do i = 1, isize_field
    fld(i) = (C_1_R - eps(i)) * fld(i)  
  end do

  return
end subroutine IBM_BCS_FIELD