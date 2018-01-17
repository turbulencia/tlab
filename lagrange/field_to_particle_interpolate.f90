#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!########################################################################
SUBROUTINE FIELD_TO_PARTICLE_INTERPOLATE &
     (iflag, nvar, data_in, data_out, l_q, grid_start, grid_end)
  
  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers_dt, pointers3d_dt
  USE DNS_GLOBAL,     ONLY : isize_particle
  USE DNS_GLOBAL,     ONLY : g
  USE LAGRANGE_GLOBAL,ONLY : l_g !, !jmin_part, l_g
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
  TINTEGER  g_p(10), g1loc, g2loc, g5loc, g6loc
  TINTEGER i, j
!  TREAL particle_local_grid_posx, particle_local_grid_posy, particle_local_grid_posz
!  TREAL dx_loc, dz_loc !, dy_loc
  TREAL dx_loc_inv, dz_loc_inv
  
! ######################################################################
! Setting grid spacings
!  dx_loc= g(1)%scale /g(1)%size
!  dy_loc= g(2)%scale/g(2)%size
!  dy_loc= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
!  dz_loc= g(3)%scale /g(3)%size
  dx_loc_inv = M_REAL( g(1)%size ) /g(1)%scale
  dz_loc_inv = M_REAL( g(3)%size ) /g(3)%scale
  
! Managing the iflag option outside the loop
  g_p(7) =1
  g_p(8) =2
  g_p(9) =1
  g_p(10)=2

  g1loc = 1
  g2loc = 2
  g5loc = 5
  g6loc = 6
  IF      ( iflag .EQ. 1 ) THEN ! halo_zone_x
     g1loc = 7
     g2loc = 8
  ELSE IF ( iflag .EQ. 2 ) THEN ! halo_zone_z
     g5loc = 9
     g6loc = 10
  ELSE IF ( iflag .EQ. 3 ) THEN ! halo_zone_diagonal
     g1loc = 7
     g2loc = 8
     g5loc = 9
     g6loc = 10
  ENDIF

! ######################################################################
  IF  ( g(3)%size .NE. 1 ) THEN

     DO i = grid_start,grid_end ! loop over all particles

! #ifdef USE_MPI
!         particle_local_grid_posx = l_q(i,1)                         /dx_loc +1 - ims_offset_i
!         particle_local_grid_posz = l_q(i,3)                         /dz_loc +1 - ims_offset_k
! #else
!         particle_local_grid_posx = l_q(i,1)                         /dx_loc +1 
!         particle_local_grid_posz = l_q(i,3)                         /dz_loc +1
! #endif
!         particle_local_grid_posy = (l_q(i,2)-g(2)%nodes(jmin_part)) /dy_loc +jmin_part  

!         g_p(1)= floor(particle_local_grid_posx)                   !position to the left (x1)
!         g_p(2)= g_p(1)+1                                    !to the right (x2)
!         g_p(3)= (floor((l_q(i,2)-g(2)%nodes(jmin_part))/dy_loc))+jmin_part !to the bottom 
!         g_p(4)= g_p(3)+1                                    !to the top (y2)
!         g_p(5)= floor(particle_local_grid_posz)                   !front side
!         g_p(6)= g_p(5)+1                                    !back side

!         length_g_p(1) = particle_local_grid_posx - g_p(1)  !length between x(i) and p
!         length_g_p(2) = g_p(2) - particle_local_grid_posx
!         length_g_p(3) = particle_local_grid_posy - g_p(3)
!         length_g_p(4) = g_p(4) - particle_local_grid_posy
!         length_g_p(5) = particle_local_grid_posz - g_p(5)  !length between z(i) and p
!         length_g_p(6) = g_p(6) - particle_local_grid_posz  !length between z(i+1) and p

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

        g_p(3)        = l_g%nodes(i)                    ! Local Y position
        g_p(4)        = g_p(3) +1
        length_g_p(3) =(l_q(i,2) - g(2)%nodes(g_p(3))) /(g(2)%nodes(g_p(4))-g(2)%nodes(g_p(3)))
        length_g_p(4) = C_1_R -length_g_p(3)
        
        cube_g_p(1) = length_g_p(1) *length_g_p(3) ! cubes
        cube_g_p(2) = length_g_p(1) *length_g_p(4) ! be carefull multiply other side cube of grid for correct interpolation
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

! -------------------------------------------------------------------
! Trilinear interpolation
! Two bilinear interpolations for each k plane (g_p(5) and g_p(6)
! Then multipled by (1-length) for trilinear aspect
! -------------------------------------------------------------------
        DO j = 1,nvar
           data_out(j)%field(i) = data_out(j)%field(i) +&
                ((cube_g_p(3) *data_in(j)%field(g_p(g1loc),g_p(3),g_p(g5loc)) &
                 +cube_g_p(4) *data_in(j)%field(g_p(g1loc),g_p(4),g_p(g5loc)) &
                 +cube_g_p(1) *data_in(j)%field(g_p(g2loc),g_p(4),g_p(g5loc)) &
                 +cube_g_p(2) *data_in(j)%field(g_p(g2loc),g_p(3),g_p(g5loc)))*length_g_p(6)) &
               +((cube_g_p(3) *data_in(j)%field(g_p(g1loc),g_p(3),g_p(g6loc)) &
                 +cube_g_p(4) *data_in(j)%field(g_p(g1loc),g_p(4),g_p(g6loc)) &
                 +cube_g_p(1) *data_in(j)%field(g_p(g2loc),g_p(4),g_p(g6loc)) &
                 +cube_g_p(2) *data_in(j)%field(g_p(g2loc),g_p(3),g_p(g6loc)))*length_g_p(5))
        ENDDO
        
     END DO
     
! ######################################################################
  ELSE !2D case

     DO i = grid_start,grid_end

! #ifdef USE_MPI
!         particle_local_grid_posx = l_q(i,1)/dx_loc + 1 - ims_offset_i
! #else
!         particle_local_grid_posx = l_q(i,1)/dx_loc + 1 
! #endif
!         particle_local_grid_posy = ((l_q(i,2)-g(2)%nodes(jmin_part))/dy_loc)+jmin_part  

!         g_p(1)= floor(particle_local_grid_posx)
!         g_p(2)= g_p(1)+1
!         g_p(3)= (floor((l_q(i,2)-g(2)%nodes(jmin_part))/dy_loc))+jmin_part
!         g_p(4)= g_p(3)+1

!         length_g_p(1) = particle_local_grid_posx - g_p(1)
!         length_g_p(2) = g_p(2) - particle_local_grid_posx
!         length_g_p(3) = particle_local_grid_posy - g_p(3)
!         length_g_p(4) = g_p(4) - particle_local_grid_posy

        length_g_p(1) = l_q(i,1) *dx_loc_inv
        g_p(1)        = FLOOR( length_g_p(1) )
        length_g_p(1) = length_g_p(1) -M_REAL( g_p(1) )
#ifdef USE_MPI
        g_p(1)        = g_p(1) +1 -ims_offset_i
#else
        g_p(1)        = g_p(1) +1
#endif
        g_p(2)        = g_p(1) +1
        length_g_p(2) = C_1_R -length_g_p(1) 
        
        g_p(3)        = l_g%nodes(i)
        g_p(4)        = g_p(3) +1
        length_g_p(3) =(l_q(i,2) - g(2)%nodes(g_p(3))) /(g(2)%nodes(g_p(4))-g(2)%nodes(g_p(3)))
        length_g_p(4) = C_1_R -length_g_p(3)

        cube_g_p(1) = length_g_p(1) *length_g_p(3)
        cube_g_p(2) = length_g_p(1) *length_g_p(4)
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

! -------------------------------------------------------------------
! Bilinear interpolation
! -------------------------------------------------------------------
        DO j = 1,nvar
           data_out(j)%field(i) = data_out(j)%field(i) +&
                (cube_g_p(3) *data_in(j)%field(g_p(g1loc),g_p(3),1) &
                +cube_g_p(4) *data_in(j)%field(g_p(g1loc),g_p(4),1) &
                +cube_g_p(1) *data_in(j)%field(g_p(g2loc),g_p(4),1) &
                +cube_g_p(2) *data_in(j)%field(g_p(g2loc),g_p(3),1))        
        ENDDO
        
     END DO
     
  ENDIF

  RETURN
END SUBROUTINE FIELD_TO_PARTICLE_INTERPOLATE

SUBROUTINE PARTICLE_LOCATE_Y( pmax, y_part, j_part, jmax, y_grid )

  IMPLICIT NONE

  TINTEGER, INTENT(IN) :: pmax, jmax
  TREAL,    DIMENSION(pmax), INTENT(IN) :: y_part
  TINTEGER, DIMENSION(pmax), INTENT(OUT) :: j_part
  TREAL,    DIMENSION(jmax), INTENT(IN) :: y_grid
  
  TINTEGER ip, jm, jp, jc

  DO ip = 1,pmax
     jp = jmax
     jm = 1
     jc = ( jm +jp ) /2  
     DO WHILE ( (y_part(ip)-y_grid(jc))*(y_part(ip)-y_grid(jc+1)) .GT. C_0_R .AND. jc .GT. jm )
        IF ( y_part(ip) .LT. y_grid(jc) ) THEN; jp = jc;
        ELSE;                                   jm = jc; END IF
        jc = ( jm +jp ) /2
     END DO
     j_part(ip) = jc
!     WRITE(*,'(i,3f)') ip, y_grid(jc), y_part(ip), y_grid(jc+1) 
  END DO
  
  RETURN
END SUBROUTINE PARTICLE_LOCATE_Y
