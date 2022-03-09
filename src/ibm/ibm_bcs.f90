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
subroutine IBM_INITIALIZE_SCAL(s,iq_scal)
  
  use DNS_IBM,   only : eps
  use DNS_IBM,   only : ibmscaljmin, ibmscaljmax 
  use TLAB_VARS, only : isize_field
  use TLAB_VARS, only : imax,jmax,kmax

  implicit none

#include "integers.h"
  
  TREAL,    dimension(isize_field,iq_scal), intent(inout) :: s
  TINTEGER,                                 intent(in)    :: iq_scal

  TINTEGER                                                :: iq

  ! ================================================================== !
  ! get scalar dirichlet boundary values of scalar field
  ! (assuming always homogenous horizontal temperature on boundaries)

  do iq = 1, iq_scal
    ibmscaljmin(iq) = s(1,                iq)
    ibmscaljmax(iq) = s(imax*(jmax-1) + 1,iq)
  end do

  call IBM_BCS_FIELD(i1,s,iq_scal)

  return
end subroutine IBM_INITIALIZE_SCAL
!########################################################################
!########################################################################
subroutine IBM_BCS_FIELD(is,fld,iq_fld)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : isize_field
  use DNS_IBM,   only : ibmscaljmin, ibmscaljmax, xbars_geo
  use TLAB_VARS, only : imax,jmax,kmax

  implicit none

#include "integers.h"
  
  TINTEGER,                                intent(in)    :: is  ! flag, is==[0,1] flow/scal
  TREAL,    dimension(isize_field,iq_fld), intent(inout) :: fld
  TINTEGER,                                intent(in)    :: iq_fld

  TINTEGER                                               :: iq, ip_b, ip_b2, k, j

  ! ================================================================== !
  ! apply IBM BCs on fields
  do iq = 1, iq_fld
    fld(:,iq) = (C_1_R - eps(:)) * fld(:,iq)  
  end do

  ! default, set only scalar value on lower boundary
  ! if objects also on upper boundary present, 
  ! values needs to be overwritten, cf. next loop
  if ( is == 1) then
    do iq = 1, iq_fld
      fld(:,iq) = fld(:,iq) + eps(:) * ibmscaljmin(iq) 
    end do

  ! in case of objects on upper boundary, set different temperature here
    if ( xbars_geo%mirrored ) then
      do iq = 1, iq_fld
        ip_b = imax*(jmax-xbars_geo%height) + 1
        do k = 1,kmax
          do j = 1, xbars_geo%height   
            ip_b2 = ip_b + imax
            fld(ip_b:ip_b2,iq) = (C_1_R - eps(ip_b:ip_b2)) * fld(ip_b:ip_b2,iq) + eps(ip_b:ip_b2)*ibmscaljmax(iq)
            ip_b  = ip_b + imax
          end do
          ip_b = imax*(jmax-xbars_geo%height) + k*imax*jmax
        end do
      end do
    end if

  end if

  return
end subroutine IBM_BCS_FIELD