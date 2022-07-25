#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF MODLE
!#   Add new geometry types or ibm related data types in this module
!#                    
!#
!########################################################################

module IBM_TYPES

  implicit none
  
  save
  
  type ibm_geo_dt
      sequence
      CHARACTER(32) :: name
      TINTEGER      :: number, height, width
      LOGICAL       :: mirrored
  end type ibm_geo_dt

end module IBM_TYPES

!########################################################################