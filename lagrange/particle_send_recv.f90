#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# HISTORY
!#
!# 2014/02 - L. Muessle
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Send particles to neighbouring processors 
!# Subroutine is divided into 2 optional parts  
!# One for East - West direction or South - North
!# First amount of to be sended particles is communicated
!# Them the actual send is realized
!#
!########################################################################

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_SEND_RECV_I(nzone_grid, nzone_west, nzone_east, &
     p_buffer_1, p_buffer_2, l_q, l_hq, l_tags)
  
  USE DNS_GLOBAL, ONLY: isize_particle, inb_particle
  
  USE DNS_MPI

  IMPLICIT NONE
  
#include "mpif.h"

  TINTEGER nzone_grid, nzone_west, nzone_east

  TREAL, DIMENSION(isize_particle, inb_particle) :: l_q
  TREAL, DIMENSION(isize_particle, inb_particle) :: l_hq
  TREAL, DIMENSION(isize_particle)               :: l_tags !Attention. Chosen TREAL on purpose. 
  TREAL, DIMENSION(*) :: p_buffer_1, p_buffer_2 !allocation = isize_particle/4*7
  
! -------------------------------------------------------------------
  TINTEGER nzone_send_west, nzone_send_east
  TINTEGER source_west, source_east
  TINTEGER dest_west, dest_east
  TINTEGER l
  TINTEGER mpireq(ims_npro*4)
  TINTEGER status(MPI_STATUS_SIZE, ims_npro*4)
  TINTEGER i, k, k_loc, m

!#######################################################################
  m=(inb_particle*2)+1 !Sending size of the buffer_parts

  ims_size_p(ims_pro+1) = nzone_grid

!###################################################################
!Setting up destination and source + send/recv number of particles which will be send
!###################################################################
  IF (ims_pro_i .EQ. 0) THEN ! Fisrt row
     dest_west= ims_npro_i -1 + ims_npro_i*ims_pro_k !Destination of the message - to West
     source_west= ims_pro_i+1 + ims_npro_i*ims_pro_k !Source of the message - from East

  ELSE IF (ims_pro_i .EQ. (ims_npro_i-1))THEN !Last row
     dest_west=  ims_pro_i-1 + ims_npro_i*ims_pro_k 
     source_west=  ims_npro_i*ims_pro_k

  ELSE !Any case
     dest_west= ims_pro_i-1 + ims_npro_i*ims_pro_k 
     source_west= ims_pro_i+1 + ims_npro_i*ims_pro_k 
  ENDIF

  mpireq(1:ims_npro*4)=MPI_REQUEST_NULL
  l = 4*ims_pro +1
!Amount of particles which will be sent/received
  CALL MPI_ISEND(nzone_west,1,MPI_INTEGER4,dest_west,0,MPI_COMM_WORLD,mpireq(l),ims_err)
  CALL MPI_IRECV(nzone_send_west,1,MPI_INTEGER4,source_west,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err)

!    CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err) !ONLY ONE MPI_WAITALL FOR THE 4 COMMUNICATIONS

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

!    mpireq(1:ims_npro*2)=MPI_REQUEST_NULL  !ONLY ONE MPI_WAITALL NOW
!    l = 2*ims_pro +1       
!Amount of particles which will be sent/received
  CALL MPI_ISEND(nzone_east,1,MPI_INTEGER4,dest_east,0,MPI_COMM_WORLD,mpireq(l+2),ims_err)
  CALL MPI_IRECV(nzone_send_east,1,MPI_INTEGER4,source_east,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+3),ims_err)

  CALL MPI_Waitall(ims_npro*4,mpireq,status,ims_err)  

!#################################################################
! Send nzone_west particles to west and get nzone_send_west particles from east
!#################################################################
!I get nzone_send_west particles form the east  
  ims_size_p(ims_pro+1)=ims_size_p(ims_pro+1)+nzone_send_west

!Construct sending buffer to west in p_buffer_2
  IF (nzone_west .NE. 0) THEN !I have something to send

!Send particles and co. west
     k=0  ! Array index 1,3 for positions, 4,6 for rhs velocities and 7 for id
     DO k_loc=1,inb_particle !Position        
        k=k+1
        DO i=1,nzone_west
           p_buffer_2(i+((k-1)*nzone_west)) = l_q(nzone_grid+i,k_loc)
        END DO
     END DO

     DO k_loc=1,inb_particle !Right hand side
        k=k+1
        DO i=1,nzone_west
           p_buffer_2(i+((k-1)*nzone_west)) = l_hq(nzone_grid+i,k_loc)
        END DO

     END DO

     k=k+1
     DO i=1,nzone_west
        p_buffer_2(i+((k-1)*nzone_west)) = l_tags(nzone_grid+i)
     END DO

  END IF

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(p_buffer_2,nzone_west*m,MPI_REAL8,dest_west,0,MPI_COMM_WORLD,mpireq(l),ims_err) ! Send particles to the west in buffer2 
  CALL MPI_IRECV(p_buffer_1,nzone_send_west*m,MPI_REAL8,source_west,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err) !Get particles from the east in buffer1

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)  

!Construct sending buffer to east in p_buffer_2 
  IF (nzone_east .NE. 0) THEN  ! I have something to send to the east

!Send particles and co. east
     k=0
     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_east
           p_buffer_2(i+((k-1)*nzone_east)) = l_q(nzone_grid+nzone_west+i,k_loc)
        END DO

     END DO

     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_east
           p_buffer_2(i+((k-1)*nzone_east)) = l_hq(nzone_grid+nzone_west+i,k_loc)
        END DO

     END DO

     k=k+1
     DO i=1,nzone_east
        p_buffer_2(i+((k-1)*nzone_east)) = l_tags(nzone_grid+nzone_west+i)
     END DO

  END IF

!Write nzone_send_west recieved particles from east (step one) that are in p_buffer_1 to the particle array
!Afterwards, p_buffer_1 is free
  IF (nzone_send_west .NE. 0) THEN
!Received particles and co. from east
     k=0 
     DO k_loc=1,inb_particle     
        k=k+1
        DO i=1,nzone_send_west 
           l_q(nzone_grid+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_west))
        END DO
     END DO

     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_send_west 
           l_hq(nzone_grid+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_west))
        END DO

     END DO

     k=k+1
     DO i=1,nzone_send_west 
        l_tags(nzone_grid+i) = p_buffer_1(i+((k-1)*nzone_send_west))
     END DO

  END IF

!#################################################################
! Send to east and get from west
!#################################################################
! Increase the particle number by the particles recieved from west
  ims_size_p(ims_pro+1)=ims_size_p(ims_pro+1)+nzone_send_east

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(p_buffer_2,nzone_east*m,MPI_REAL8,dest_east,0,MPI_COMM_WORLD,mpireq(l),ims_err) !Send p_buffer_2 to the east
  CALL MPI_IRECV(p_buffer_1,nzone_send_east*m,MPI_REAL8,source_east,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err) !Get particles from the west into p_buffer_1

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err) 

!Write nzone_send_east recieved particles from west (step one) that are in p_buffer_1 to the particle array
  IF (nzone_send_east .NE. 0) THEN
!Received particles and co. from west
     k=0 
     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_send_east 
           l_q(nzone_grid+nzone_send_west+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_east))
        END DO

     END DO

     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_send_east
           l_hq(nzone_grid+nzone_send_west+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_east))
        END DO

     END DO

     k=k+1
     DO i=1,nzone_send_east 
        l_tags(nzone_grid+nzone_send_west+i) = p_buffer_1(i+((k-1)*nzone_send_east))

     END DO

  END IF

  RETURN
END SUBROUTINE PARTICLE_SEND_RECV_I

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_SEND_RECV_K(nzone_grid, nzone_south, nzone_north, &
     p_buffer_1, p_buffer_2, l_q, l_hq, l_tags)
  
  USE DNS_GLOBAL, ONLY: isize_particle, inb_particle
  
  USE DNS_MPI

  IMPLICIT NONE
  
#include "mpif.h"

  TINTEGER nzone_grid, nzone_south, nzone_north

  TREAL, DIMENSION(isize_particle, inb_particle) :: l_q
  TREAL, DIMENSION(isize_particle, inb_particle) :: l_hq
  TREAL, DIMENSION(isize_particle)               :: l_tags !Attention. Chosen TREAL on purpose. 
  TREAL, DIMENSION(*) :: p_buffer_1, p_buffer_2 !allocation = isize_particle/4*7
  
! -------------------------------------------------------------------
  TINTEGER nzone_send_south, nzone_send_north
  TINTEGER source_south, source_north
  TINTEGER dest_south, dest_north
  TINTEGER l
  TINTEGER mpireq(ims_npro*4)
  TINTEGER status(MPI_STATUS_SIZE, ims_npro*4)
  TINTEGER i, k, k_loc, m

!#######################################################################
  m=(inb_particle*2)+1 !Sending size of the buffer_parts

  ims_size_p(ims_pro+1) = nzone_grid

!###################################################################
!Setting up destination and source + send/recv number of particles which will be send
!###################################################################
  IF (ims_pro_k .EQ. 0) THEN ! Fisrt row
     dest_south= ims_pro_i + ims_npro_i*(ims_npro_k-1)!Destination of the message
     source_south= ims_pro_i + ims_npro_i !Source of the message
     
  ELSE IF (ims_pro_k .EQ. (ims_npro_k-1))THEN !Last row
     dest_south=  ims_pro_i + ims_npro_i*(ims_npro_k-2)
     source_south=  ims_pro_i

  ELSE !Any case
     dest_south= ims_pro_i +  ims_npro_i*(ims_pro_k-1)  !Destination of the message
     source_south= ims_pro_i +  ims_npro_i*(ims_pro_k+1) !Source of the message
  ENDIF

  mpireq(1:ims_npro*4)=MPI_REQUEST_NULL
  l = 4*ims_pro +1
!Amount of particles which will be sent/received
  CALL MPI_ISEND(nzone_south,1,MPI_INTEGER4,dest_south,0,MPI_COMM_WORLD,mpireq(l),ims_err)
  CALL MPI_IRECV(nzone_send_south,1,MPI_INTEGER4,source_south,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err)

!    CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)  ! ALBERTO CONSIDER DELETING THIS ONE

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


!    mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
!    l = 2*ims_pro +1
!Amount of particles which will be sent/received
  CALL MPI_ISEND(nzone_north,1,MPI_INTEGER4,dest_north,0,MPI_COMM_WORLD,mpireq(l+2),ims_err)
  CALL MPI_IRECV(nzone_send_north,1,MPI_INTEGER4,source_north,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+3),ims_err)

  CALL MPI_Waitall(ims_npro*4,mpireq,status,ims_err)  

!#################################################################
! Send to west and get from east
!#################################################################
  ims_size_p(ims_pro+1)=ims_size_p(ims_pro+1)+nzone_send_south

  IF (nzone_south .NE. 0) THEN

!Send particles and co. to south
!Fill p_buffer_2
     k=0
     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_south
           p_buffer_2(i+((k-1)*nzone_south)) = l_q(nzone_grid+i,k_loc)
        END DO

     END DO

     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_south
           p_buffer_2(i+((k-1)*nzone_south)) = l_hq(nzone_grid+i,k_loc)
        END DO

     END DO

     k=k+1
     DO i=1,nzone_south
        p_buffer_2(i+((k-1)*nzone_south)) = l_tags(nzone_grid+i)
     END DO

  END IF

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(p_buffer_2,nzone_south*m,MPI_REAL8,dest_south,0,MPI_COMM_WORLD,mpireq(l),ims_err)
  CALL MPI_IRECV(p_buffer_1,nzone_send_south*m,MPI_REAL8,source_south,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err)

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)  

  IF (nzone_north .NE. 0) THEN      

!Send particles and co. to north
!Fill p_buffer_2 again
     k=0
     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_north
           p_buffer_2(i+((k-1)*nzone_north)) = l_q(nzone_grid+nzone_south+i,k_loc)
        END DO

     END DO

     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_north
!           p_buffer_2(i+((k-1)*nzone_north)) = l_hq(nzone_grid+nzone_east+i,k_loc)
           p_buffer_2(i+((k-1)*nzone_north)) = l_hq(nzone_grid+nzone_south+i,k_loc)
        END DO

     END DO

     k=k+1
     DO i=1,nzone_north
        p_buffer_2(i+((k-1)*nzone_north)) = l_tags(nzone_grid+nzone_south+i)
     END DO

  END IF

  IF (nzone_send_south .NE. 0) THEN

!Received particles and co. from north
     k=0 
     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_send_south
           l_q(nzone_grid+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_south))
        END DO
     END DO

     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_send_south
           l_hq(nzone_grid+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_south))
        END DO

     END DO

     k=k+1
     DO i=1,nzone_send_south
        l_tags(nzone_grid+i) = p_buffer_1(i+((k-1)*nzone_send_south))
     END DO

  END IF

!#################################################################
! Send to east and get from west
!#################################################################
  ims_size_p(ims_pro+1)=ims_size_p(ims_pro+1)+nzone_send_north

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(p_buffer_2,nzone_north*m,MPI_REAL8,dest_north,0,MPI_COMM_WORLD,mpireq(l),ims_err)
  CALL MPI_IRECV(p_buffer_1,nzone_send_north*m,MPI_REAL8,source_north,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err)

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)  

  IF (nzone_send_north .NE. 0) THEN
!Received particles and co. from south
     k=0 
     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_send_north
           l_q(nzone_grid+nzone_send_south+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_north))
        END DO

     END DO

     DO k_loc=1,inb_particle

        k=k+1
        DO i=1,nzone_send_north
           l_hq(nzone_grid+nzone_send_south+i,k_loc) = p_buffer_1(i+((k-1)*nzone_send_north))
        END DO

     END DO

     k=k+1
     DO i=1,nzone_send_north
        l_tags(nzone_grid+nzone_send_south+i) = p_buffer_1(i+((k-1)*nzone_send_north))

     END DO

  END IF

  RETURN
END SUBROUTINE PARTICLE_SEND_RECV_K
