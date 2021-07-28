#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_TIME_SUBSTEP(dte, l_q, l_hq, l_comm )    
  
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : isize_particle, inb_part
  USE LAGRANGE_VARS, ONLY : isize_l_comm
  USE LAGRANGE_VARS, ONLY : l_g
#ifdef USE_MPI
  USE LAGRANGE_VARS, ONLY : isize_pbuffer
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL dte
  TREAL, DIMENSION(isize_particle,*)     :: l_q
  TREAL, DIMENSION(isize_particle,*)     :: l_hq
  TREAL, DIMENSION(isize_l_comm), TARGET :: l_comm

! -------------------------------------------------------------------
  TINTEGER is, i

#ifdef USE_MPI
  TREAL, DIMENSION(:), POINTER :: p_buffer_1, p_buffer_2
  TINTEGER nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north
#else
  TREAL x_right, z_right
#endif

!#####################################################################
#ifdef USE_MPI
  p_buffer_1(1:isize_pbuffer)=> l_comm(               1:isize_pbuffer   )
  p_buffer_2(1:isize_pbuffer)=> l_comm(isize_pbuffer +1:isize_pbuffer *2)

#else
  x_right = g(1)%nodes(1) +g(1)%scale
  z_right = g(3)%nodes(1) +g(3)%scale

#endif

!#######################################################################
! Particle new position
!#######################################################################
  DO i = 1,l_g%np
     DO is = 1,inb_part
        l_q(i,is) = l_q(i,is) + dte*l_hq(i,is)

     ENDDO
  END DO

!#####################################################################
! Boundaries
!#####################################################################
#ifdef USE_MPI
! -------------------------------------------------------------------
!Particle sorting for Send/Recv X-Direction
! -------------------------------------------------------------------
  CALL PARTICLE_SORT(1, l_g,l_q, l_hq, nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north)

  IF (ims_pro_i .EQ. 0) THEN !Take care of periodic boundary conditions west
     IF (nzone_west .NE. 0) THEN
        l_q(nzone_grid+1:nzone_grid+nzone_west,1) =&
             l_q(nzone_grid+1:nzone_grid+nzone_west,1) +g(1)%scale
     END IF
  END IF

  IF( ims_pro_i .EQ. (ims_npro_i-1)) THEN !Take care of periodic boundary conditions east
     IF (nzone_east .NE. 0) THEN
        l_q(nzone_grid+nzone_west+1:nzone_grid+nzone_west+nzone_east,1) =&
             l_q(nzone_grid+nzone_west+1:nzone_grid+nzone_west+nzone_east,1) -g(1)%scale
     END IF
  END IF

  CALL PARTICLE_SEND_RECV_I(nzone_grid, nzone_west, nzone_east, & 
       p_buffer_1, p_buffer_2, l_q, l_hq, l_g%tags, l_g%np) 

! -------------------------------------------------------------------
!Particle sorting for Send/Recv Z-Direction
! -------------------------------------------------------------------
  CALL PARTICLE_SORT(3, l_g,l_q, l_hq, nzone_grid, nzone_west, nzone_east,nzone_south,nzone_north)

  IF (ims_pro_k .EQ. 0) THEN !Take care of periodic boundary conditions south
     IF (nzone_south .NE. 0) THEN
        l_q(nzone_grid+1:nzone_grid+nzone_south,3) = &
             l_q(nzone_grid+1:nzone_grid+nzone_south,3) +g(3)%scale
     END IF
  END IF

  IF( ims_pro_k .EQ. (ims_npro_k-1)) THEN !Take care of periodic boundary conditions north
     IF (nzone_north .NE. 0) THEN
        l_q(nzone_grid+nzone_south+1:nzone_grid+nzone_south+nzone_north,3) =&
             l_q(nzone_grid+nzone_south+1:nzone_grid+nzone_south+nzone_north,3) -g(3)%scale
     END IF
  END IF

  CALL PARTICLE_SEND_RECV_K(nzone_grid, nzone_south, nzone_north, & 
       p_buffer_1, p_buffer_2, l_q, l_hq, l_g%tags, l_g%np)

  NULLIFY(p_buffer_1,p_buffer_2)
  
#else
!#######################################################################
! Serial
  DO i = 1,l_g%np
     IF    ( l_q(i,1) .GT. x_right       ) THEN
        l_q(i,1) = l_q(i,1) - g(1)%scale

     ELSEIF( l_q(i,1) .LT. g(1)%nodes(1) ) THEN
        l_q(i,1) = l_q(i,1) + g(1)%scale
        
     END IF
     
     IF    ( l_q(i,3) .GT. z_right       ) THEN
        l_q(i,3) = l_q(i,3) - g(3)%scale
        
     ELSEIF( l_q(i,3) .LT. g(3)%nodes(1) ) THEN
        l_q(i,3) = l_q(i,3) + g(3)%scale
        
     END IF
  END DO
  
#endif
  
!#######################################################################
! Recalculating closest node below in Y direction
!#######################################################################
  CALL PARTICLE_LOCATE_Y( l_g%np, l_q(1,2), l_g%nodes, g(2)%size, g(2)%nodes )

  RETURN
END SUBROUTINE PARTICLE_TIME_SUBSTEP
