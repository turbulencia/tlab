#include "types.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2021/XX/XX - J. Kostelecky
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

subroutine IBM_BCS_FIELD(fld,iq_fld)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : isize_field

  implicit none
  
  TINTEGER,                                intent(in)    :: iq_fld
  TREAL,    dimension(isize_field,iq_fld), intent(inout) :: fld
  TINTEGER                                               :: iq

  ! ================================================================== !
  ! apply IBM BCs on flow fields
  do iq = 1, iq_fld
    fld(:,iq) = (C_1_R - eps(:)) * fld(:,iq)  
  end do

  return
end subroutine IBM_BCS_FIELD