#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY
!#
!# 2014/03 - L. Muessle
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Sends field_halo to neighbouring processors 
!# Field has size imax+1 and kmax+1 
!# Only halo columns are needed for interpolation to field
!#
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_EAST(f_buffer_1,f_buffer_2, field )

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_MPI

  IMPLICIT NONE

#include "mpif.h"

  TINTEGER source_east
  TINTEGER dest_east
  TINTEGER l
  TINTEGER mpireq(ims_npro*2)
  TINTEGER status(MPI_STATUS_SIZE, ims_npro*2)

  TREAL, DIMENSION(imax+1,jmax,kmax+1) :: field
  TREAL, DIMENSION(1,jmax,kmax+1) :: f_buffer_1, f_buffer_2 

  f_buffer_1=C_0_R
  f_buffer_2=C_0_R

!###################################################################
!Calculating the source/dest for east and north
!###################################################################
  IF (ims_pro_i .EQ. 0) THEN ! Fisrt row
     dest_east= ims_pro_i +1 + ims_npro_i*ims_pro_k !Destination of the message - to East
     source_east= ims_npro_i-1 + ims_npro_i*ims_pro_k !Source of the message - from West

  ELSE IF (ims_pro_i .EQ. (ims_npro_i-1))THEN !Last row
     dest_east=  ims_npro_i*ims_pro_k
     source_east=  ims_pro_i-1 + ims_npro_i*ims_pro_k 

  ELSE !Any case
     dest_east= ims_pro_i+1 + ims_npro_i*ims_pro_k !Dest of the message
     source_east= ims_pro_i-1 + ims_npro_i*ims_pro_k !source of the message
  ENDIF

!#################################################################
!Setting up the plane which needs to be send EAST
!Send this plane EAST and receive the plane from WEST
!#################################################################
!The buffer is one point larger because the diagonal point is first send EAST
!and then NORTH
  f_buffer_1(1,1:jmax,1:(kmax+1))=field(imax+1,1:jmax,1:(kmax+1))

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(f_buffer_1,jmax*(kmax+1),MPI_REAL8,dest_east,0,MPI_COMM_WORLD,mpireq(l),ims_err) ! Send particles to the west in buffer2 
  CALL MPI_IRECV(f_buffer_2,jmax*(kmax+1),MPI_REAL8,source_east,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err) !Get particles from the east in buffer1

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)  

  field(1,1:jmax,1:(kmax+1))=field(1,1:jmax,1:(kmax+1)) + f_buffer_2(1,1:jmax,1:(kmax+1))

!    field(1,1:jmax,1:(kmax-1))=field(1,1:jmax,1:(kmax-1)) + f_buffer_2(1,1:jmax,1:(kmax-1))
!    field(1,1:jmax,kmax)=field(1,1:jmax,kmax) + f_buffer_2(1,1:jmax,kmax)

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_EAST
