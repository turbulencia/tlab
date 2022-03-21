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

subroutine IBM_CHECK_PROCS()
  
  use DNS_IBM

#ifdef USE_MPI 
  use MPI
  use TLAB_MPI_VARS, only: ims_pro
#endif
  
  implicit none
  
#include "integers.h"

  ! ================================================================== !

#ifdef USE_MPI 
  ! Check in X
  if (sum(epsi) == 0) then
    ims_pro_ibm_x = .false.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in x'
#endif
  else 
    ims_pro_ibm_x = .true.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in x'
#endif
  end if 

  ! Check in Y
  if (sum(epsj) == 0) then
    ims_pro_ibm_y = .false.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in y' 
#endif
  else 
    ims_pro_ibm_y = .true.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in y'
#endif
  end if 

  ! Check in Z
  if (sum(epsk) == 0) then
    ims_pro_ibm_z = .false.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' idle   for IBM spline generation in z' 
#endif
  else 
    ims_pro_ibm_z = .true.
#ifdef IBM_DEBUG
    write(*,*) 'Task: ', ims_pro, ' active for IBM spline generation in z'
#endif
  end if 
#else
  ims_pro_ibm_x = .true.; ims_pro_ibm_y = .true.; ims_pro_ibm_z = .true. ! one task with full domain
#endif

  return
end subroutine IBM_CHECK_PROCS

!########################################################################