#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

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

subroutine IBM_CHECK_PROCS()
  
  use DNS_IBM
  use DNS_MPI, only: ims_pro
  
  implicit none
  
#include "integers.h"
#ifdef USE_MPI 
#include "mpif.h"
#include "dns_const_mpi.h"     
#endif

  ! ================================================================== !

#ifdef USE_MPI 
  ! Check in X
  if (sum(epsi) == 0) then
    ims_pro_ibm_x = .false.
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in x'
  else 
    ims_pro_ibm_x = .true.
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in x'
  end if 

  ! Check in Y
  if (sum(epsj) == 0) then
    ims_pro_ibm_y = .false.
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in y' 
  else 
    ims_pro_ibm_y = .true.
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in y'
  end if 

  ! Check in Z
  if (sum(epsk) == 0) then
    ims_pro_ibm_z = .false.
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in z' 
  else 
    ims_pro_ibm_z = .true.
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in z'
  end if 
#else
  ims_pro_ibm_x = .true.; ims_pro_ibm_y = .true.; ims_pro_ibm_z = .true.
#endif

  return
end subroutine IBM_CHECK_PROCS

!########################################################################