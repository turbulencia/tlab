!########################################################################
!# Tool/Library
!#
!########################################################################
!# DESCRIPTION  halo_plane_shifting
!#
!# Data from the main fields is transfered into the data to halo plane
!# In cae of MPI, this means transfer additional transfer of data from
!# neigthbour processors. Two MPI transfer are done, taking into account
!# the third halo zone.
!#
!# 1 refers to the creation of the first halo (only one MPI transfer)
!#
!# 12.2013 Created Lukas Muessle
!#
!#######################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

SUBROUTINE  halo_plane_shifting_i_1d(field,halo_field, halo_field_ik,buffer_send,buffer_recv, &
                                  diagonal_point_send, diagonal_point_recv, upper_left_point)

USE DNS_GLOBAL, ONLY: imax,jmax,kmax, imax_total, kmax_total
#ifdef USE_MPI
  USE DNS_MPI
#endif

IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"

  integer source
  integer dest
  integer l  ! Counter for messages
  integer  mpireq(ims_npro*2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)


#endif

  TREAL, DIMENSION(imax,jmax,kmax)        ::field ! Array of the main fields  
  TREAL, DIMENSION(2,jmax,kmax)           ::halo_field !halo plane in i direction
  TREAL, DIMENSION(1,jmax,kmax+1)         ::buffer_send, buffer_recv !Buffers needed for the mpi
  TREAL, DIMENSION(jmax)                  ::diagonal_point_send, diagonal_point_recv, upper_left_point
  TREAL, DIMENSION(2,jmax,2)              ::halo_field_ik !ghost plane
  TINTEGER  j 

! ######################################################################
! first array halo
! local last slice transfered to halo_field
! ######################################################################

  halo_field(1,1:jmax,1:kmax)=field(imax,1:jmax,1:kmax) !last row of imax in ghost plane x-direction

! ######################################################################
! second array of the halo field (in genral from other proc)
! ######################################################################
 
  IF (ims_npro_i .EQ. 1) THEN   ! Serial code in i-direction 
    halo_field(2,1:jmax,1:kmax)=field(1,1:jmax,1:kmax) !first row of imax in ghost plane
  ELSE  ! General

  ! ######################################################################
  ! start parallel computing
  ! ######################################################################
   
#ifdef USE_MPI
    ! ######################################################################
    ! Using ISend and IRecv (non blocking communication)
    ! ######################################################################
    
    !  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
      mpireq(1:ims_npro*2)=MPI_REQUEST_NULL !need to be set for all mpireqs
    
    
    ! ######################################################################
    ! Transfer array data to halo_field_i and halo_fiel_ik
    ! ######################################################################
      IF (ims_pro_i .EQ. 0) THEN ! Fisrt row
        dest= ims_npro_i -1 + ims_npro_i*ims_pro_k !Source of the message
        source= ims_pro_i+1 + ims_npro_i*ims_pro_k !Destination of the message

      ELSE IF (ims_pro_i .EQ. (ims_npro_i-1))THEN !Last row
        dest=  ims_pro_i-1 + ims_npro_i*ims_pro_k 
        source=  ims_npro_i*ims_pro_k

      ELSE !Any case
        dest= ims_pro_i-1 + ims_npro_i*ims_pro_k !Destination of the message
        source= ims_pro_i+1 + ims_npro_i*ims_pro_k !Source of the message
      ENDIF
      l = 2*ims_pro +1
      
      buffer_send(1,1:jmax,1:kmax)=field(1,1:jmax,1:kmax)

      ! Pack buffer with the additional point for diagonal halo_field_ik
        DO j=1,jmax
          buffer_send(1,j,kmax+1)=diagonal_point_send(j)
        END DO
    
      CALL MPI_ISEND(buffer_send,jmax*(kmax+1),MPI_REAL8,dest,0,MPI_COMM_WORLD,mpireq(l), ims_err)
      CALL MPI_IRECV(buffer_recv,jmax*(kmax+1),MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
       
      CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

      halo_field(2,1:jmax,1:kmax)=buffer_recv(1,1:jmax,1:kmax)
      ! Extract the additional diagonal point which was sended
        DO j=1,jmax
          diagonal_point_recv(j)=buffer_recv(1,j,kmax+1)
        END DO
    
   
    
    
    IF (ims_npro_k .EQ. 1) THEN   !MPI in i-direcition but serial in k-direction
      halo_field_ik(1,1:jmax,1)=field(imax,1:jmax,kmax) !point 1,1 is down left
      halo_field_ik(1,1:jmax,2)=halo_field(2,1:jmax,kmax) !1,2 is down right
      halo_field_ik(2,1:jmax,1)=field(imax,1:jmax,1) !2,1 is top left
      halo_field_ik(2,1:jmax,2)=halo_field(2,1:jmax,1) !2,2 is top right
    ELSE    ! MPI in both direction (i and k)
#endif
      halo_field_ik(1,1:jmax,1)=field(imax,1:jmax,kmax) !point 1,1 is down left
      halo_field_ik(1,1:jmax,2)=halo_field(2,1:jmax,kmax) !1,2 is down right
      halo_field_ik(2,1:jmax,1)=upper_left_point(1:jmax) !2,1 is top left
      halo_field_ik(2,1:jmax,2)=diagonal_point_recv(1:jmax) !2,2 is top right
#ifdef USE_MPI     
    END IF
#endif
  END IF    

  RETURN
END SUBROUTINE halo_plane_shifting_i_1d
