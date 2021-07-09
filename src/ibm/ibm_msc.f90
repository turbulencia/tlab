SUBROUTINE DNS_IBM_TEST()
  ! USE DNS_IBM
  WRITE(*,*) 'TESTING IBM COMPILATION'
END SUBROUTINE DNS_IBM_TEST


subroutine BOUNDARY_BCS_IBM_FLOW()

  implicit none
  
#ifdef USE_MPI 
#include "mpif.h"
#endif 

  return
end subroutine BOUNDARY_BCS_IBM_FLOW
!########################################################################
subroutine WRITE_GEOMETRY()
  !
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  
  return
end subroutine WRITE_GEOMETRY
!########################################################################
subroutine READ_GEOMETRY() ! if needed: restart/run with already generated geometry
  !
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  !
  return
end subroutine READ_GEOMETRY
!########################################################################  
subroutine BOUNDARY_BCS_IBM_SCAL()
  !
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  !
  return
end subroutine BOUNDARY_BCS_IBM_SCAL
!########################################################################
subroutine IBM_FINALIZE() ! dealloc of arrays? maybe not needed...
  !
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  !
  return
end subroutine IBM_FINALIZE
!########################################################################