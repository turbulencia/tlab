#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / ATHORS
!#
!# 2021/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   called once in dns_main.f90 before time integration starts,
!#   all relevant geometry informations needed for the IBM are 
!#   available (and maybe written to disk?)
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

subroutine IBM_INITIALIZE_GEOMETRY(txc, wrk3d)  
  
  use DNS_IBM
  
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  
  TREAL, dimension(*), intent(inout) :: txc, wrk3d

  ! ================================================================== !

  ! generate native 3d-geometry field (eps_aux) of immersed objects (define your own geomtry here)
  call IBM_GENERATE_GEOMETRY_XBARS(wrk3d) 

  ! transpose eps in epsi, epsj, epsk and allocate neccessary memory
  call IBM_GEOMETRY_TRANSPOSE(wrk3d,txc)

  ! generate relevant geometry fields for IBM routines (nobi, nobj, nobk)
  call IBM_GENERATE_GEOMETRY(wrk3d,txc) ! txc for DEBUG
  
  ! not coded yet
  ! read/write geometry fields from/to disk
  ! call IBM_READ_GEOMETRY() ! call IBM_WRITE_GEOMETRY()           
  
  return
end subroutine IBM_INITIALIZE_GEOMETRY

!########################################################################