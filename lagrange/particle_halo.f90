#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_HALO_K(nvar, data, halo_field_k, buffer_send, buffer_recv)

  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers3d_dt
  USE DNS_GLOBAL,     ONLY : imax,jmax,kmax

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro, ims_npro, ims_pro_k, ims_npro_k, ims_map_k
  USE DNS_MPI, ONLY : ims_err
#endif
  
  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif
  
  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar):: data
  TREAL, DIMENSION(imax,jmax,2,nvar)  :: halo_field_k
  TREAL, DIMENSION(imax,jmax,1,nvar)  :: buffer_send, buffer_recv

! -------------------------------------------------------------------
  TINTEGER i
  
#ifdef USE_MPI
  integer source, dest, l, size
  integer mpireq(ims_npro*2+2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)
#endif
  
! ######################################################################
#ifdef USE_MPI
  IF ( ims_npro_k .EQ. 1 ) THEN
#endif
     DO i = 1,nvar
        halo_field_k(1:imax,1:jmax,1,i) = data(i)%field(1:imax,1:jmax,kmax)
        halo_field_k(1:imax,1:jmax,2,i) = data(i)%field(1:imax,1:jmax,1   )        
     ENDDO
     
#ifdef USE_MPI
! ######################################################################
  ELSE
     DO i = 1,nvar
        halo_field_k(1:imax,1:jmax,1,i) = data(i)%field(1:imax,1:jmax,kmax)
        buffer_send(1:imax,1:jmax,1,i)  = data(i)%field(1:imax,1:jmax,1   ) ! data to be transfered
     ENDDO
     
     mpireq(1:ims_npro*2) = MPI_REQUEST_NULL
     l      = 2*ims_pro +1
     dest   = ims_map_k(MOD(ims_pro_k-1 +ims_npro_k,ims_npro_k) +1)
     source = ims_map_k(MOD(ims_pro_k+1 +ims_npro_k,ims_npro_k) +1)
     size   = imax*jmax*nvar
     CALL MPI_ISEND(buffer_send,size,MPI_REAL8,dest,  0,          MPI_COMM_WORLD,mpireq(l),   ims_err)
     CALL MPI_IRECV(buffer_recv,size,MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
     CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)
     
     halo_field_k(1:imax,1:jmax,2,1:nvar) = buffer_recv(1:imax,1:jmax,1,1:nvar)

  END IF
#endif
  
  RETURN
END SUBROUTINE PARTICLE_HALO_K

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_HALO_I(nvar, data, halo_field_i, halo_field_k, halo_field_ik, buffer_send, buffer_recv)

  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers3d_dt
  USE DNS_GLOBAL,     ONLY : imax,jmax,kmax

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro, ims_npro, ims_pro_i, ims_npro_i, ims_map_i
  USE DNS_MPI, ONLY : ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif
  
  TINTEGER nvar
  TYPE(pointers3d_dt), DIMENSION(nvar):: data
  TREAL, DIMENSION(2,jmax,kmax,nvar)  :: halo_field_i
  TREAL, DIMENSION(imax,jmax,2,nvar)  :: halo_field_k 
  TREAL, DIMENSION(2,jmax,2,nvar)     :: halo_field_ik
  TREAL, DIMENSION(1,jmax,kmax+1,nvar):: buffer_send, buffer_recv

! -------------------------------------------------------------------
  TINTEGER i
  
#ifdef USE_MPI
  integer source, dest, l, size
  integer mpireq(ims_npro*2+2)
  integer status(MPI_STATUS_SIZE,ims_npro*2)
#endif
  
! ######################################################################
#ifdef USE_MPI
  IF ( ims_npro_i .EQ. 1 ) THEN
#endif
     DO i = 1,nvar
        halo_field_i (1,1:jmax,1:kmax,i) = data(i)%field(imax,1:jmax,1:kmax) 
        halo_field_i (2,1:jmax,1:kmax,i) = data(i)%field(1,   1:jmax,1:kmax) 
        halo_field_ik(2,1:jmax,2     ,i) = halo_field_k (1,   1:jmax,2     ,i) ! top-right corner
     END DO

#ifdef USE_MPI
! ######################################################################
  ELSE
     DO i = 1,nvar
        halo_field_i(1,1:jmax,1:kmax,i) = data(i)%field(imax,1:jmax,1:kmax)  
        buffer_send(1,1:jmax,1:kmax,i)  = data(i)%field(1,   1:jmax,1:kmax)   ! data to be transfered
        buffer_send(1,1:jmax,kmax+1,i)  = halo_field_k (1,   1:jmax,2     ,i)
     ENDDO
     
     mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
     l      = 2*ims_pro +1
     dest   = ims_map_i(MOD(ims_pro_i-1 +ims_npro_i,ims_npro_i) +1)
     source = ims_map_i(MOD(ims_pro_i+1 +ims_npro_i,ims_npro_i) +1)
     size   = jmax*(kmax+1)*nvar
     CALL MPI_ISEND(buffer_send,size,MPI_REAL8,dest,  0,          MPI_COMM_WORLD,mpireq(l),   ims_err)
     CALL MPI_IRECV(buffer_recv,size,MPI_REAL8,source,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1), ims_err)
     CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

     halo_field_i (2,1:jmax,1:kmax,1:nvar) = buffer_recv(1,1:jmax,1:kmax,1:nvar)
     halo_field_ik(2,1:jmax,2     ,1:nvar) = buffer_recv(1,1:jmax,kmax+1,1:nvar) ! top-right corner
     
  END IF
#endif
  
  halo_field_ik(1,1:jmax,1,1:nvar) = halo_field_i(1,   1:jmax,kmax,1:nvar)
  halo_field_ik(2,1:jmax,1,1:nvar) = halo_field_i(2,   1:jmax,kmax,1:nvar)
  halo_field_ik(1,1:jmax,2,1:nvar) = halo_field_k(imax,1:jmax,2   ,1:nvar)
  
  RETURN
END SUBROUTINE PARTICLE_HALO_I
