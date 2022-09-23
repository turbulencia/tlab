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
!#           Created
!# 2017/12 - J.P. Mellado
!#           Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Interpolate particle information to a field
!#
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD(l_q, particle_property, field_out, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY: imax,jmax,kmax
  USE PARTICLE_VARS, ONLY: isize_part
#ifdef USE_MPI
  USE MPI
   USE TLAB_MPI_VARS, ONLY: ims_err
#endif

   IMPLICIT NONE

  TREAL, DIMENSION(isize_part,3)    :: l_q
  TREAL, DIMENSION(isize_part)      :: particle_property
  TREAL, DIMENSION(imax,jmax,kmax)      :: field_out
  TREAL, DIMENSION(*)                   :: wrk2d
  TREAL, DIMENSION(imax+1,jmax,kmax+1)  :: wrk3d

!#######################################################################
!Write data to fielf
!Field is an extended field with grid, halo1, halo2 and halo3
!#######################################################################
  wrk3d = C_0_R
  CALL PARTICLE_TO_FIELD_INTERPOLATE(l_q, particle_property, wrk3d)

#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

!#######################################################################
!SEND to the East and RECV from the West
!Sum the most left colum of field
!SEND to the North and RECV from the South
!Sum the lowest row of field
!#######################################################################
 CALL PARTICLE_TO_FIELD_SEND_RECV_EAST (wrk2d(1),wrk2d((jmax*(kmax+1))+1),wrk3d)
 CALL PARTICLE_TO_FIELD_SEND_RECV_NORTH(wrk2d(1),wrk2d(((imax+1)*jmax)+1),wrk3d)

#endif

!#######################################################################
!Put wrk3d into continuous memory field_out
!#######################################################################
  field_out(1:imax,1:jmax,1:kmax)=wrk3d(1:imax,1:jmax,1:kmax)

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD_INTERPOLATE(l_q, particle_property, field)

  USE TLAB_VARS,     ONLY: imax,jmax,kmax
  USE TLAB_VARS,     ONLY: g
  USE PARTICLE_VARS,     ONLY: isize_part
  USE PARTICLE_ARRAYS,ONLY: l_g
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_VARS, ONLY: ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE
#include "integers.h"

  TREAL, DIMENSION(imax+1,jmax,kmax+1) :: field
  TREAL, DIMENSION(isize_part,3)   :: l_q
  TREAL, DIMENSION(isize_part)     :: particle_property

  TREAL length_g_p(6), cube_g_p(4)
  TINTEGER  g_p(6)
  TINTEGER i
  TREAL dx_loc_inv, dz_loc_inv

! ######################################################################
  dx_loc_inv = M_REAL( g(1)%size ) /g(1)%scale
  dz_loc_inv = M_REAL( g(3)%size ) /g(3)%scale

  g_p(5)= 1 ! Default is 2D
  g_p(6)= 1

! ######################################################################
  DO i=1,l_g%np

     length_g_p(1) = l_q(i,1) *dx_loc_inv            ! Local X position
     g_p(1)        = FLOOR( length_g_p(1) )
     length_g_p(1) = length_g_p(1) -M_REAL( g_p(1) )
#ifdef USE_MPI
     g_p(1)        = g_p(1) +1 -ims_offset_i
#else
     g_p(1)        = g_p(1) +1
#endif
     g_p(2)        = g_p(1) +1
     length_g_p(2) = C_1_R -length_g_p(1)

     IF (  g(3)%size .NE. 1 ) THEN
        length_g_p(5) = l_q(i,3) *dz_loc_inv            ! Local Z position
        g_p(5)        = FLOOR( length_g_p(5) )
        length_g_p(5) = length_g_p(5) -M_REAL( g_p(5) )
#ifdef USE_MPI
        g_p(5)        = g_p(5) +1 -ims_offset_k
#else
        g_p(5)        = g_p(5) +1
#endif
        g_p(6)        = g_p(5) +1
        length_g_p(6) = C_1_R -length_g_p(5)
     ENDIF

     g_p(3)        = l_g%nodes(i)                    ! Local Y position
     g_p(4)        = g_p(3) +1
     length_g_p(3) =(l_q(i,2) - g(2)%nodes(g_p(3))) /(g(2)%nodes(g_p(4))-g(2)%nodes(g_p(3)))
     length_g_p(4) = C_1_R -length_g_p(3)

     cube_g_p(1) = length_g_p(1) *length_g_p(3) ! bilear cubes for X and Y
     cube_g_p(2) = length_g_p(1) *length_g_p(4) ! be carefull multiply other side cube of grid for correct interpolation
     cube_g_p(3) = length_g_p(4) *length_g_p(2)
     cube_g_p(4) = length_g_p(2) *length_g_p(3)

!###################################################################
!Two bilinear calculation for each k direction (g_p(5) and g_p(6))
!Then multipled by (1-length) for Trilinear aspect
!###################################################################
     IF ( g(3)%size .EQ. 1 ) THEN
!U
        field(g_p(1),g_p(3),g_p(5))= field(g_p(1),g_p(3),g_p(5)) + particle_property(i)*cube_g_p(3)
        field(g_p(1),g_p(4),g_p(5))= field(g_p(1),g_p(4),g_p(5)) + particle_property(i)*cube_g_p(4)
        field(g_p(2),g_p(4),g_p(5))= field(g_p(2),g_p(4),g_p(5)) + particle_property(i)*cube_g_p(1)
        field(g_p(2),g_p(3),g_p(5))= field(g_p(2),g_p(3),g_p(5)) + particle_property(i)*cube_g_p(2)

     ELSE
!Z 1
        field(g_p(1),g_p(3),g_p(5))= field(g_p(1),g_p(3),g_p(5)) + particle_property(i)*cube_g_p(3)*length_g_p(6)
        field(g_p(1),g_p(4),g_p(5))= field(g_p(1),g_p(4),g_p(5)) + particle_property(i)*cube_g_p(4)*length_g_p(6)
        field(g_p(2),g_p(4),g_p(5))= field(g_p(2),g_p(4),g_p(5)) + particle_property(i)*cube_g_p(1)*length_g_p(6)
        field(g_p(2),g_p(3),g_p(5))= field(g_p(2),g_p(3),g_p(5)) + particle_property(i)*cube_g_p(2)*length_g_p(6)

!Z 2
        field(g_p(1),g_p(3),g_p(6))= field(g_p(1),g_p(3),g_p(6)) + particle_property(i)*cube_g_p(3)*length_g_p(5)
        field(g_p(1),g_p(4),g_p(6))= field(g_p(1),g_p(4),g_p(6)) + particle_property(i)*cube_g_p(4)*length_g_p(5)
        field(g_p(2),g_p(4),g_p(6))= field(g_p(2),g_p(4),g_p(6)) + particle_property(i)*cube_g_p(1)*length_g_p(5)
        field(g_p(2),g_p(3),g_p(6))= field(g_p(2),g_p(3),g_p(6)) + particle_property(i)*cube_g_p(2)*length_g_p(5)

     END IF

  END DO

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD_INTERPOLATE

#ifdef USE_MPI

!########################################################################
!# DESCRIPTION
!#
!# Sends field_halo to neighbouring processors
!# Field has size imax+1 and kmax+1
!# Only halo columns are needed for interpolation to field
!#
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_EAST(f_buffer_1,f_buffer_2, field )

  USE TLAB_VARS, ONLY : imax,jmax,kmax
  USE MPI
  USE TLAB_MPI_VARS

  IMPLICIT NONE

  TINTEGER source_east
  TINTEGER dest_east
  TINTEGER l
  TINTEGER mpireq(ims_npro*2)
  TINTEGER status(MPI_STATUS_SIZE, ims_npro*2)

  TREAL, DIMENSION(imax+1,jmax,kmax+1) :: field
  TREAL, DIMENSION(jmax,kmax+1) :: f_buffer_1, f_buffer_2

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
  f_buffer_1(1:jmax,1:(kmax+1))=field(imax+1,1:jmax,1:(kmax+1))

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(f_buffer_1,jmax*(kmax+1),MPI_REAL8,dest_east,0,MPI_COMM_WORLD,mpireq(l),ims_err) ! Send particles to the west in buffer2
  CALL MPI_IRECV(f_buffer_2,jmax*(kmax+1),MPI_REAL8,source_east,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err) !Get particles from the east in buffer1

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

  field(1,1:jmax,1:(kmax+1))=field(1,1:jmax,1:(kmax+1)) + f_buffer_2(1:jmax,1:(kmax+1))

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_EAST

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_NORTH(f_buffer_1,f_buffer_2, field )

  USE TLAB_VARS, ONLY : imax,jmax,kmax
  USE MPI
  USE TLAB_MPI_VARS

  IMPLICIT NONE

  TINTEGER  source_north
  TINTEGER  dest_north
  TINTEGER l
  TINTEGER mpireq(ims_npro*2)
  TINTEGER status(MPI_STATUS_SIZE, ims_npro*2)

  TREAL, DIMENSION(imax+1,jmax,kmax+1) :: field
  TREAL, DIMENSION(imax+1,jmax)        :: f_buffer_1, f_buffer_2

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
  f_buffer_1(1:(imax+1),1:jmax)=field(1:(imax+1),1:jmax,kmax+1)

!CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

  mpireq(1:ims_npro*2)=MPI_REQUEST_NULL
  l = 2*ims_pro +1
!Actual sending/receiving of particles
  CALL MPI_ISEND(f_buffer_1,(imax+1)*jmax,MPI_REAL8,dest_north,0,MPI_COMM_WORLD,mpireq(l),ims_err) !Send p_buffer_2 to the east
  CALL MPI_IRECV(f_buffer_2,(imax+1)*jmax,MPI_REAL8,source_north,MPI_ANY_TAG,MPI_COMM_WORLD,mpireq(l+1),ims_err) !Get particles from the west into p_buffer_1

  CALL MPI_Waitall(ims_npro*2,mpireq,status,ims_err)

  field(1:(imax+1),1:jmax,1)=field(1:(imax+1),1:jmax,1) + f_buffer_2(1:(imax+1),1:jmax)

  RETURN
END SUBROUTINE PARTICLE_TO_FIELD_SEND_RECV_NORTH

#endif
