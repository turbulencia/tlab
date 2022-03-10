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
subroutine IBM_IO_READ(wrk3d)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : imax,jmax,kmax
  use IO_FIELDS

  implicit none

#include "integers.h"
  
  TREAL, dimension(*), intent(inout) ::  wrk3d

  character(len=32)                  :: fname

  ! ================================================================== !
  ! read eps field (filename 'eps0.1')
  write(fname,*) i0; 
  fname = trim(adjustl('eps'))//trim(adjustl(fname))
  
  call IO_READ_FIELDS(fname, IO_FLOW, imax,jmax,kmax, i1, i1, eps, wrk3d)
  
  return
end subroutine IBM_IO_READ