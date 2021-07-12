subroutine IBM_TEST()

  write(*,*) 'TESTING IBM COMPILATION'

end subroutine IBM_TEST
!########################################################################
subroutine IBM_WRITE_GEOMETRY()
  !
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  
  return
end subroutine IBM_WRITE_GEOMETRY
!########################################################################
subroutine IBM_READ_GEOMETRY() ! if needed: restart/run with already generated geometry
  !
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  !
  return
end subroutine IBM_READ_GEOMETRY
!########################################################################  
subroutine IBM_BOUNDARY_BCS_IBM_SCAL()
  !
  implicit none

#ifdef USE_MPI 
#include "mpif.h"
#endif 
  !
  return
end subroutine IBM_BOUNDARY_BCS_IBM_SCAL
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