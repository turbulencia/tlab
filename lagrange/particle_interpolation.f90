#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_INTERPOLATION &
     (iflag, nvar, data_in, data_out, l_q, grid_start, grid_end)
  
  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers_dt, pointers3d_dt
  USE DNS_GLOBAL,     ONLY : isize_particle
  USE DNS_GLOBAL,     ONLY : g
  USE LAGRANGE_GLOBAL,ONLY : jmin_part
#ifdef USE_MPI
  USE DNS_MPI, ONLY: ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE
#include "integers.h"

  TINTEGER iflag, nvar, grid_start, grid_end
  TYPE(pointers3d_dt), DIMENSION(nvar)             :: data_in     
  TYPE(pointers_dt),   DIMENSION(nvar)             :: data_out
  TREAL,               DIMENSION(isize_particle,3) :: l_q

! -------------------------------------------------------------------
  TREAL length_g_p(6), cube_g_p(4)
  TINTEGER  gridpoint(10), g1loc, g2loc, g5loc, g6loc
  TINTEGER i, j
  TREAL particle_local_grid_posx, particle_local_grid_posy, particle_local_grid_posz
  TREAL dx_loc, dy_loc, dz_loc
  
! ######################################################################
! Setting grid spacings
  dx_loc= g(1)%scale /g(1)%size
!  dy_loc= g(2)%scale/g(2)%size
  dy_loc= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
  dz_loc= g(3)%scale /g(3)%size

! Managing the iflag option outside the loop
  gridpoint(7) =1
  gridpoint(8) =2
  gridpoint(9) =1
  gridpoint(10)=2

  g1loc = 1
  g2loc = 2
  g5loc = 5
  g6loc = 6
  IF ( iflag .EQ. 1 ) THEN
     g1loc = 7
     g2loc = 8
  ELSE IF ( iflag .EQ. 2 ) THEN
     g5loc = 9
     g6loc = 10
  ELSE IF ( iflag .EQ. 3 ) THEN
     g1loc = 7
     g2loc = 8
     g5loc = 9
     g6loc = 10
  ENDIF

! ######################################################################
  IF  ( g(3)%size .NE. 1 ) THEN

     DO i = grid_start,grid_end ! loop over all particles

#ifdef USE_MPI
        particle_local_grid_posx = l_q(i,1)                         /dx_loc +1 - ims_offset_i
        particle_local_grid_posz = l_q(i,3)                         /dz_loc +1 - ims_offset_k
#else
        particle_local_grid_posx = l_q(i,1)                         /dx_loc +1 
        particle_local_grid_posz = l_q(i,3)                         /dz_loc +1
#endif
        particle_local_grid_posy = (l_q(i,2)-g(2)%nodes(jmin_part)) /dy_loc +jmin_part  

        gridpoint(1)= floor(particle_local_grid_posx)                   !position to the left (x1)
        gridpoint(2)= gridpoint(1)+1                                    !to the right (x2)
        gridpoint(3)= (floor((l_q(i,2)-g(2)%nodes(jmin_part))/dy_loc))+jmin_part !to the bottom 
        gridpoint(4)= gridpoint(3)+1                                    !to the top (y2)
        gridpoint(5)= floor(particle_local_grid_posz)                   !front side
        gridpoint(6)= gridpoint(5)+1                                    !back side

        length_g_p(1) = particle_local_grid_posx - gridpoint(1)  !length between x(i) and p
        length_g_p(2) = gridpoint(2) - particle_local_grid_posx
        length_g_p(3) = particle_local_grid_posy - gridpoint(3)
        length_g_p(4) = gridpoint(4) - particle_local_grid_posy
        length_g_p(5) = particle_local_grid_posz - gridpoint(5)  !length between z(i) and p
        length_g_p(6) = gridpoint(6) - particle_local_grid_posz  !length between z(i+1) and p

        cube_g_p(1) = length_g_p(1) *length_g_p(3) ! cubes
        cube_g_p(2) = length_g_p(1) *length_g_p(4) ! be carefull multiply other side cube of grid for correct interpolation
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

        ! IF ( iflag .EQ. 1 ) THEN
        !    gridpoint(1)=1
        !    gridpoint(2)=2
        ! ELSE IF ( iflag .EQ. 2 ) THEN
        !    gridpoint(5)=1
        !    gridpoint(6)=2
        ! ELSE IF ( iflag .EQ. 3 ) THEN
        !    gridpoint(1)=1
        !    gridpoint(2)=2
        !    gridpoint(5)=1
        !    gridpoint(6)=2
        ! ENDIF
        
! -------------------------------------------------------------------
! Trilinear interpolation
! Two bilinear interpolations for each k plane (gridpoint(5) and gridpoint(6)
! Then multipled by (1-length) for trilinear aspect
! -------------------------------------------------------------------
        DO j = 1,nvar
           data_out(j)%field(i) = data_out(j)%field(i) +&
                ((cube_g_p(3) *data_in(j)%field(gridpoint(g1loc),gridpoint(3),gridpoint(g5loc)) &
                 +cube_g_p(4) *data_in(j)%field(gridpoint(g1loc),gridpoint(4),gridpoint(g5loc)) &
                 +cube_g_p(1) *data_in(j)%field(gridpoint(g2loc),gridpoint(4),gridpoint(g5loc)) &
                 +cube_g_p(2) *data_in(j)%field(gridpoint(g2loc),gridpoint(3),gridpoint(g5loc)))*length_g_p(6)) &
               +((cube_g_p(3) *data_in(j)%field(gridpoint(g1loc),gridpoint(3),gridpoint(g6loc)) &
                 +cube_g_p(4) *data_in(j)%field(gridpoint(g1loc),gridpoint(4),gridpoint(g6loc)) &
                 +cube_g_p(1) *data_in(j)%field(gridpoint(g2loc),gridpoint(4),gridpoint(g6loc)) &
                 +cube_g_p(2) *data_in(j)%field(gridpoint(g2loc),gridpoint(3),gridpoint(g6loc)))*length_g_p(5))
        ENDDO
        
     END DO
     
! ######################################################################
  ELSE !2D case

     DO i = grid_start,grid_end

#ifdef USE_MPI
        particle_local_grid_posx = l_q(i,1)/dx_loc + 1 - ims_offset_i
#else
        particle_local_grid_posx = l_q(i,1)/dx_loc + 1 
#endif
        particle_local_grid_posy = ((l_q(i,2)-g(2)%nodes(jmin_part))/dy_loc)+jmin_part  

        gridpoint(1)= floor(particle_local_grid_posx)
        gridpoint(2)= gridpoint(1)+1
        gridpoint(3)= (floor((l_q(i,2)-g(2)%nodes(jmin_part))/dy_loc))+jmin_part
        gridpoint(4)= gridpoint(3)+1
        ! gridpoint(5)=1
        ! gridpoint(6)=1

        length_g_p(1) = particle_local_grid_posx - gridpoint(1)
        length_g_p(2) = gridpoint(2) - particle_local_grid_posx
        length_g_p(3) = particle_local_grid_posy - gridpoint(3)
        length_g_p(4) = gridpoint(4) - particle_local_grid_posy

        cube_g_p(1) = length_g_p(1) *length_g_p(3)
        cube_g_p(2) = length_g_p(1) *length_g_p(4)
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

        ! IF ( iflag .EQ. 1 ) THEN
        !    gridpoint(1)=1
        !    gridpoint(2)=2
        ! ENDIF
        
! -------------------------------------------------------------------
! Bilinear interpolation
! -------------------------------------------------------------------
        DO j = 1,nvar
           data_out(j)%field(i) = data_out(j)%field(i) +&
                (cube_g_p(3) *data_in(j)%field(gridpoint(g1loc),gridpoint(3),1) &
                +cube_g_p(4) *data_in(j)%field(gridpoint(g1loc),gridpoint(4),1) &
                +cube_g_p(1) *data_in(j)%field(gridpoint(g2loc),gridpoint(4),1) &
                +cube_g_p(2) *data_in(j)%field(gridpoint(g2loc),gridpoint(3),1))        
        ENDDO
        
     END DO
     
  ENDIF

  RETURN
END SUBROUTINE PARTICLE_INTERPOLATION
