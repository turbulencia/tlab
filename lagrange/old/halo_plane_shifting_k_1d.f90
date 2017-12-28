!########################################################################
!# Tool/Library
!#
!########################################################################
!# DESCRIPTION
!#
!# shifting of data to ghost plane
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

SUBROUTINE  halo_plane_shifting_k_1d(field,halo_field, buffer_send, buffer_recv, diagonal_point_send, upper_left_point)

USE DNS_GLOBAL, ONLY: imax,jmax,kmax
#ifdef USE_MPI
  USE DNS_MPI
#endif

IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"

  integer source
  integer dest , l
  integer  mpireq(ims_npro*2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)

#endif

  TREAL, DIMENSION(imax,jmax,kmax),INTENT(in)       ::field
  TREAL, DIMENSION(imax,jmax,2),INTENT(inout)       ::halo_field !ghost plane
  TREAL, DIMENSION(imax,jmax,1)                     ::buffer_send, buffer_recv
  TREAL, DIMENSION(jmax)                            ::diagonal_point_send, upper_left_point


! ######################################################################
! first array halo
! local last slice transfered to halo_field
! ######################################################################

  halo_field(1:imax,1:jmax,1)=field(1:imax,1:jmax,kmax) !back row of normal grid

! ######################################################################
! second array of the halo field (in genral from other proc)
! ######################################################################

  IF (ims_npro_k .EQ. 1) THEN ! Only one proc, no MPI needed
    halo_field(1:imax,1:jmax,2)=field(1:imax,1:jmax,1)
   
  ELSE !General case 
    
  ! ######################################################################
  ! start parallel computing
  ! ######################################################################
#ifdef USE_MPI
    
      mpireq(1:ims_npro*2)=MPI_REQUEST_NULL !need to be set for all mpireqs
    
    ! ######################################################################
    ! Transfer array data to halo_field_i
    ! ######################################################################
      IF (ims_pro_k .EQ. 0) THEN ! Fisrt row
         dest= ims_pro_i + ims_npro_i*(ims_npro_k-1)!Destination of the message
         source= ims_pro_i + ims_npro_i !Source of the message
      ELSE IF (ims_pro_k .EQ. (ims_npro_k-1))THEN !Last row
         dest=  ims_pro_i + ims_npro_i*(ims_npro_k-2)
         source=  ims_pro_i
      ELSE !Any case
         dest= ims_pro_i +  ims_npro_i*(ims_pro_k-1)  !Destination of the message
         source= ims_pro_i +  ims_npro_i*(ims_pro_k+1) !Source of the message
      ENDIF
      
      l = 2*ims_pro +1
      buffer_send(1:imax,1:jmax,1)=field(1:imax,1:jmax,1)
      CALL MPI_ISEND(buffer_send,imax*jmax,MPI_REAL8,dest,0,MPI_COMM_WORLD,mpireq(l), ims_err)
      CALL MPI_IRECV(buffer_recv,imax*jmax,MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
      
      CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

      halo_field(1:imax,1:jmax,2)=buffer_recv(1:imax,1:jmax,1)
      diagonal_point_send(1:jmax)=halo_field(1,1:jmax,2)
      upper_left_point(1:jmax)=halo_field(imax,1:jmax,2)


#endif
  END IF
  
  RETURN
END SUBROUTINE halo_plane_shifting_k_1d
