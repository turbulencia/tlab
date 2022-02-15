#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Calculate the three invariants of the velocity gradient tensor, div u 
!#
!########################################################################

!########################################################################
! First invariant
! (caution: div(u)=0 condition only holds on pressure nodes)
!########################################################################
SUBROUTINE FI_INVARIANT_P(nx,ny,nz, u,v,w, result, tmp1, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imode_ibm
  USE DNS_IBM,   ONLY : ibm_partial
  
! ############################################# ! 
! DEBUG
#ifdef IBM_DEBUG
#ifdef USE_MPI
  use TLAB_MPI_VARS, only : ims_pro
  use IO_FIELDS
#endif
#endif
! ############################################# ! 

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)

  ! ############################################# ! 
  ! DEBUG
#ifdef IBM_DEBUG
#ifdef USE_MPI
#else
  TINTEGER, parameter                       :: ims_pro=0
#endif
#endif
  ! ############################################# ! 
  
! ###################################################################
  bcs = 0

  ! -------------------------------------------------------------------
  ! IBM   (if .true., OPR_PARTIAL_X/Y/Z uses modified fields for derivatives)
  IF ( imode_ibm == 1 ) ibm_partial = .true.
#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) '========================================================='
  if (ims_pro == 0) write(*,*) 'ibm_partial in dns_control "FI_INVARIANT_P"', ibm_partial
#endif
  ! -------------------------------------------------------------------
  
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, result, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp1,   wrk3d, wrk2d,wrk3d)
#ifdef IBM_DEBUG 
  call IO_WRITE_FIELDS('dudx', IO_FLOW, nx,ny,nz, i1, result, wrk3d)
  call IO_WRITE_FIELDS('dvdy', IO_FLOW, nx,ny,nz, i1, tmp1,   wrk3d)
#endif
  result =  result + tmp1
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp1,   wrk3d, wrk2d,wrk3d)
#ifdef IBM_DEBUG   
  call IO_WRITE_FIELDS('dwdz', IO_FLOW, nx,ny,nz, i1, tmp1,   wrk3d)
#endif
  result =-(result + tmp1)

  ! -------------------------------------------------------------------
  IF ( imode_ibm == 1 ) ibm_partial = .false.  

#ifdef IBM_DEBUG
  if (ims_pro == 0) write(*,*) '========================================================='
  if (ims_pro == 0) write(*,*) 'ibm_partial in dns_control "FI_INVARIANT_P"', ibm_partial
  call IO_WRITE_FIELDS('dil',  IO_FLOW, nx,ny,nz, i1, result, wrk3d)
#endif

  RETURN
END SUBROUTINE FI_INVARIANT_P

!########################################################################
! Second invariant of the velocity gradient tensor, Q, as defined
! by Chong et al. (1990). Regions with high positive values of Q are 
! vorticity dominated, regions with high negative values of Q are 
! strain dominated.
! If incompressible, P=0 and then Q=(w^2-2s^2)/4, where P is the first
! invariant, P=-div(v). 
!########################################################################
SUBROUTINE FI_INVARIANT_Q(nx,ny,nz, u,v,w, result, tmp1,tmp2,tmp3, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! -------------------------------------------------------------------
! diagonal terms
! -------------------------------------------------------------------
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
  result = tmp1 *tmp2 +tmp2 *tmp3 +tmp3 *tmp1

! -------------------------------------------------------------------
! off-diagonal terms
! -------------------------------------------------------------------
! UyVx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result -tmp1 *tmp2

! UzWx
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp2, wrk3d, wrk2d,wrk3d)
  result = result -tmp1 *tmp2

! WyVz
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp2, wrk3d, wrk2d,wrk3d)
  result = result -tmp1 *tmp2

  RETURN
END SUBROUTINE FI_INVARIANT_Q

!########################################################################
! Calculate third invariant of the velocity gradient tensor (-determinant)
! The derivatives v,y and v,z are repeated to avoid more arrays
!########################################################################
SUBROUTINE FI_INVARIANT_R(nx,ny,nz, u,v,w, result, tmp1,tmp2,tmp3,tmp4,tmp5, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : g
  
  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: u,v,w
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: result
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4,tmp5
  TREAL, DIMENSION(*),        INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TINTEGER bcs(2,2)
  
! ###################################################################
  bcs = 0

! ###################################################################
! Term u,x (v,y w,z-w,y v,z)
! ###################################################################
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), w, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), w, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), u, tmp5, wrk3d, wrk2d,wrk3d)
  result = tmp5 *( tmp1 *tmp3 -tmp2 *tmp4 )

! ###################################################################
! Term v,x (u,z w,y-w,z u,y)
! ###################################################################
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), u, tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), v, tmp5, wrk3d, wrk2d,wrk3d)
  result = result + tmp5 *( tmp2 *tmp3 -tmp1 *tmp4 )

! ###################################################################
! Term v,x (u,z w,y-w,z u,y)
! ###################################################################
  CALL OPR_PARTIAL_Z(OPR_P1, nx,ny,nz, bcs, g(3), v, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1, nx,ny,nz, bcs, g(1), w, tmp5, wrk3d, wrk2d,wrk3d)
  result = result + tmp5 *( tmp1 *tmp4 -tmp2 *tmp3 )

  result =-result ! set the right sign

  RETURN
END SUBROUTINE FI_INVARIANT_R
