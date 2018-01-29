#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_HALO_SERIAL(nvar, data, halo_field_1, halo_field_2, halo_field_3)

  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers3d_dt
  USE DNS_GLOBAL,     ONLY : imax,jmax,kmax
  USE DNS_GLOBAL,     ONLY : g

  IMPLICIT NONE

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar) :: data
  TREAL, DIMENSION(2,jmax,kmax,nvar) :: halo_field_1
  TREAL, DIMENSION(imax,jmax,2,nvar) :: halo_field_2
  TREAL, DIMENSION(2,jmax,2,nvar)    :: halo_field_3

! -------------------------------------------------------------------
  TINTEGER i
  
!#######################################################################
  DO i = 1,nvar
     halo_field_1(1,1:jmax,1:kmax,i) = data(i)%field(g(1)%size,1:jmax,1:kmax)  
     halo_field_1(2,1:jmax,1:kmax,i) = data(i)%field(1,        1:jmax,1:kmax)  
     
     halo_field_2(1:imax,1:jmax,1,i) = data(i)%field(1:imax,1:jmax,g(3)%size)
     halo_field_2(1:imax,1:jmax,2,i) = data(i)%field(1:imax,1:jmax,1        )
     
     halo_field_3(1,1:jmax,1,i)      = data(i)%field(g(1)%size,1:jmax,g(3)%size)
     halo_field_3(1,1:jmax,2,i)      = data(i)%field(g(1)%size,1:jmax,1        )
     halo_field_3(2,1:jmax,1,i)      = data(i)%field(1,        1:jmax,g(3)%size)
     halo_field_3(2,1:jmax,2,i)      = data(i)%field(1,        1:jmax,1        )
  END DO

  RETURN
END SUBROUTINE PARTICLE_HALO_SERIAL

#ifdef USE_MPI

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_HALO_K(nvar, data, halo_field, buffer_send, buffer_recv, diagonal_point_send, upper_left_point)

  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers3d_dt
  USE DNS_GLOBAL,     ONLY : imax,jmax,kmax

  USE DNS_MPI, ONLY : ims_pro, ims_npro, ims_pro_k, ims_npro_k, ims_map_k
  USE DNS_MPI, ONLY : ims_err
  
  IMPLICIT NONE

#include "mpif.h"

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar)  :: data
  TREAL, DIMENSION(imax,jmax,2,nvar)  :: halo_field !ghost plane
  TREAL, DIMENSION(imax,jmax,1,nvar)  :: buffer_send, buffer_recv
  TREAL, DIMENSION(jmax,nvar)         :: diagonal_point_send, upper_left_point

! -------------------------------------------------------------------
  TINTEGER i
  
  integer source, dest, l
  integer mpireq(ims_npro*2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)

! ######################################################################
  IF ( ims_npro_k .EQ. 1 ) THEN  ! Serial code in k-direction
     DO i = 1,nvar
        halo_field(1:imax,1:jmax,1,i) = data(i)%field(1:imax,1:jmax,kmax)
        halo_field(1:imax,1:jmax,2,i) = data(i)%field(1:imax,1:jmax,1)
        
     ENDDO
     
! ######################################################################
  ELSE
     dest   = ims_map_k(MOD(ims_pro_k-1 +ims_npro_k,ims_npro_k) +1)
     source = ims_map_k(MOD(ims_pro_k+1 +ims_npro_k,ims_npro_k) +1)
     l = 2*ims_pro +1
     
     DO i = 1,nvar
        halo_field(1:imax,1:jmax,1,i) = data(i)%field(1:imax,1:jmax,kmax)
        buffer_send(1:imax,1:jmax,1,i)= data(i)%field(1:imax,1:jmax,1) ! Data to be transfered
     ENDDO
     
     mpireq(1:ims_npro*2) = MPI_REQUEST_NULL !need to be set for all mpireqs

     CALL MPI_ISEND(buffer_send,imax*jmax*nvar,MPI_REAL8,dest,  0,          MPI_COMM_WORLD,mpireq(l),   ims_err)
     CALL MPI_IRECV(buffer_recv,imax*jmax*nvar,MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
     
     CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)
     
     halo_field(1:imax,1:jmax,2,1:nvar) = buffer_recv(1:imax,1:jmax,1,1:nvar)
     diagonal_point_send(1:jmax,1:nvar) = halo_field(1,1:jmax,2,1:nvar)
     upper_left_point(1:jmax,1:nvar)    = halo_field(imax,1:jmax,2,1:nvar)
     
  END IF
  
  RETURN
END SUBROUTINE PARTICLE_HALO_K

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_HALO_I(nvar, data, halo_field, halo_field_ik, buffer_send, buffer_recv, &
     diagonal_point_send, diagonal_point_recv, upper_left_point)

  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers3d_dt
  USE DNS_GLOBAL,     ONLY : imax,jmax,kmax

  USE DNS_MPI, ONLY : ims_pro, ims_npro, ims_pro_i, ims_npro_i, ims_npro_k, ims_map_i
  USE DNS_MPI, ONLY : ims_err
  
  IMPLICIT NONE

#include "mpif.h"

  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar)   :: data
  TREAL, DIMENSION(2,jmax,kmax,nvar)   :: halo_field  !halo plane in i direction
  TREAL, DIMENSION(2,jmax,2,nvar)      :: halo_field_ik  !ghost plane
  TREAL, DIMENSION(1,jmax,kmax+1,nvar) :: buffer_send, buffer_recv !Buffers needed for the mpi
  TREAL, DIMENSION(jmax,nvar)          :: diagonal_point_send, diagonal_point_recv, upper_left_point

! -------------------------------------------------------------------
  TINTEGER i
  
  integer source
  integer dest
  integer l  ! Counter for messages
  integer mpireq(ims_npro*2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)

! ######################################################################
  IF ( ims_npro_i .EQ. 1 ) THEN   ! Serial code in i-direction 
     DO i = 1,nvar
        halo_field(1,1:jmax,1:kmax,i) = data(i)%field(imax,1:jmax,1:kmax)  
        halo_field(2,1:jmax,1:kmax,i) = data(i)%field(1,   1:jmax,1:kmax)
        
     END DO

! ######################################################################
  ELSE
     dest   = ims_map_i(MOD(ims_pro_i-1 +ims_npro_i,ims_npro_i) +1)
     source = ims_map_i(MOD(ims_pro_i+1 +ims_npro_i,ims_npro_i) +1)
     l = 2*ims_pro +1
           
     DO i = 1,nvar
        halo_field(1,1:jmax,1:kmax,i) = data(i)%field(imax,1:jmax,1:kmax)  

        buffer_send(1,1:jmax,1:kmax,i)= data(i)%field(1,   1:jmax,1:kmax) ! Data to be transfered
        buffer_send(1,1:jmax,kmax+1,i)= diagonal_point_send(1:jmax,i) ! diagonal halo_field_ik
        
     ENDDO
     
     mpireq(1:ims_npro*2)=MPI_REQUEST_NULL !need to be set for all mpireqs

     CALL MPI_ISEND(buffer_send,jmax*(kmax+1)*nvar,MPI_REAL8,dest,  0,          MPI_COMM_WORLD,mpireq(l),   ims_err)
     CALL MPI_IRECV(buffer_recv,jmax*(kmax+1)*nvar,MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
     
     CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

     halo_field(2,1:jmax,1:kmax,1:nvar) = buffer_recv(1,1:jmax,1:kmax,1:nvar)
      
! Extract the additional diagonal point which was sended
     ! DO i=1,nvar
     !    diagonal_point_recv(1:jmax,i)=buffer_recv(1,1:jmax,kmax+1,i)
     ! END DO
     diagonal_point_recv(1:jmax,1:nvar)=buffer_recv(1,1:jmax,kmax+1,1:nvar)
   
     IF ( ims_npro_k .EQ. 1 ) THEN   !MPI in i-direcition but serial in k-direction
        DO i=1,nvar
           halo_field_ik(1,1:jmax,1,i) = data(i)%field(imax,1:jmax,kmax) !point 1,1 is down left
           halo_field_ik(1,1:jmax,2,i) = halo_field(2,1:jmax,kmax,i) !1,2 is down right
           halo_field_ik(2,1:jmax,1,i) = data(i)%field(imax,1:jmax,1) !2,1 is top left
           halo_field_ik(2,1:jmax,2,i) = halo_field(2,1:jmax,1,i) !2,2 is top right
           
        ENDDO
        
     ELSE    ! MPI in both direction (i and k)
        DO i=1,nvar
           halo_field_ik(1,1:jmax,1,i) = data(i)%field(imax,1:jmax,kmax) !point 1,1 is down left
           halo_field_ik(1,1:jmax,2,i) = halo_field(2,1:jmax,kmax,i) !1,2 is down right
           halo_field_ik(2,1:jmax,1,i) = upper_left_point(1:jmax,i) !2,1 is top left
           halo_field_ik(2,1:jmax,2,i) = diagonal_point_recv(1:jmax,i) !2,2 is top right
           
        ENDDO
        
     END IF

  END IF

  RETURN
END SUBROUTINE PARTICLE_HALO_I

#endif
