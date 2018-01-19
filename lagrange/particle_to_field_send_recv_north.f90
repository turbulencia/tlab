#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# DESCRIPTION
!#
!# Sends field_halo to neighbouring processors 
!# Field has size imax+1 and kmax+1 
!# Only halo columns are needed for interpolation to field
!#
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_NORTH(f_buffer_1,f_buffer_2, field )
  
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_MPI

  IMPLICIT NONE

#include "mpif.h"

  TINTEGER  source_north
  TINTEGER  dest_north
  TINTEGER l
  TINTEGER mpireq(ims_npro*2)
  TINTEGER status(MPI_STATUS_SIZE, ims_npro*2)

  TREAL, DIMENSION(imax+1,jmax,kmax+1) :: field
  TREAL, DIMENSION(imax+1,jmax,1) :: f_buffer_1, f_buffer_2 

  f_buffer_1=C_0_R
  f_buffer_2=C_0_R

!###################################################################
!Calculating the source/dest for east and north
!###################################################################
  IF (ims_pro_k .EQ. 0) THEN ! Fisrt row
     dest_north= ims_pro_i + ims_npro_i * (ims_pro_k+1) !dest of the message - to east
     source_north= ims_pro_i + ims_npro_i * (ims_npro_k-1) !source of the message - from wesst

  ELSE IF (ims_pro_k .EQ. (ims_npro_k-1))THEN !Last row
     dest_north=  ims_pro_i
     source_north=  ims_pro_i + ims_npro_i * (ims_pro_k-1)

  ELSE !Any case
     dest_north= ims_pro_i + ims_npro_i * (ims_pro_k+1) !Dest of the message
     source_north= ims_pro_i + ims_npro_i * (ims_pro_k-1)  !source of the message

  ENDIF

!#################################################################
!Setting up the plane which needs to be send NORTH
!Send this plane NORTH and receive the plane from SOUTH
!#################################################################
  f_buffer_1(1:(imax+1),1:jmax,1)=field(1:(imax+1),1:jmax,kmax+1)

!CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(f_buffer_1,(imax+1)*jmax,MPI_REAL8,dest_north,0,MPI_COMM_WORLD,mpireq(l),ims_err) !Send p_buffer_2 to the east
  CALL MPI_IRECV(f_buffer_2,(imax+1)*jmax,MPI_REAL8,source_north,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err) !Get particles from the west into p_buffer_1

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err) 

  field(1:(imax+1),1:jmax,1)=field(1:(imax+1),1:jmax,1) + f_buffer_2(1:(imax+1),1:jmax,1)

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_NORTH
