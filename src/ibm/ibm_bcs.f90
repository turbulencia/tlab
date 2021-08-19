#include "types.h"

!########################################################################
!# HISTORY / ATHORS
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

subroutine IBM_BCS_FLOW(q,iq_flow)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : isize_field

  implicit none
  
  TINTEGER,                                 intent(in)    :: iq_flow
  TREAL,    dimension(isize_field,iq_flow), intent(inout) :: q
  TINTEGER                                                :: iq

  ! ================================================================== !
  ! apply IBM BCs on flow fields
  do iq = 1, iq_flow
    q(:,iq) = (C_1_R - eps(:)) * q(:,iq)  
  end do

  return
end subroutine IBM_BCS_FLOW

!########################################################################

subroutine IBM_BCS_SCAL(s, iq_scal)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : isize_field

  
  TINTEGER,                                 intent(in)    :: iq_scal
  TREAL,    dimension(isize_field,iq_flow), intent(inout) :: s
  TINTEGER                                                :: is                                      

  ! ================================================================== !
  ! apply IBM BCs on scal fields
  do is = 1, iq_scal
    s(:,is) = (C_1_R - eps(:)) * s(:,is)  
  end do

  return
end subroutine IBM_BCS_SCAL

!########################################################################