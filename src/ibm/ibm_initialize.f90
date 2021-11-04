#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
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

  ! check idle procs
  if (ibm_procs_idle) then
    call IBM_CHECK_PROCS()
  else
    ims_pro_ibm_x = .true.; ims_pro_ibm_y = .true.; ims_pro_ibm_z = .true.
  end if
  
  ! not implemented yet
  ! read/write geometry fields from/to disk
  ! call IBM_READ_GEOMETRY() ! call IBM_WRITE_GEOMETRY()
  
  ! deallocate not needed arrays after ini
  call IBM_DEALLOCATE()
  
  return
end subroutine IBM_INITIALIZE_GEOMETRY

!########################################################################