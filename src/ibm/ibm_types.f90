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

  use TLAB_CONSTANTS, only : wi

  implicit none
  
  save
  
  type ibm_geo_dt
      sequence
      character(32) :: name
      integer(wi)   :: number, height, width, hill_slope
      logical       :: mirrored
  end type ibm_geo_dt

end module IBM_TYPES

!########################################################################