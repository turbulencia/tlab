#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!########################################################################
SUBROUTINE PARTICLE_TO_FIELD_INTERPOLATE(l_q, particle_property, field)

  USE DNS_GLOBAL,     ONLY: imax,jmax,kmax
  USE DNS_GLOBAL,     ONLY: g
  USE DNS_GLOBAL,     ONLY: isize_particle
  USE LAGRANGE_GLOBAL,ONLY: l_g, particle_number_local  !jmin_part
#ifdef USE_MPI
  USE DNS_MPI, ONLY: ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(imax+1,jmax,kmax+1) :: field
  TREAL, DIMENSION(isize_particle,3)   :: l_q  
  TREAL, DIMENSION(isize_particle)     :: particle_property 

  TREAL length_g_p(6), cube_g_p(4)
  TINTEGER  g_p(6)
  TINTEGER i 
  ! TREAL particle_local_grid_posx, particle_local_grid_posy, particle_local_grid_posz
  ! TREAL dx_loc, dy_loc, dz_loc
  TREAL dx_loc_inv, dz_loc_inv

! ######################################################################
! Setting grid spacings
!   dx_loc= g(1)%scale /g(1)%size
! !  dy_loc= g(2)%scale/g(2)%size
!   dy_loc= g(2)%nodes(jmin_part+1)-g(2)%nodes(jmin_part)
!   dz_loc= g(3)%scale /g(3)%size
  dx_loc_inv = M_REAL( g(1)%size ) /g(1)%scale
  dz_loc_inv = M_REAL( g(3)%size ) /g(3)%scale

  g_p(5)= 1 ! Default is 2D
  g_p(6)= 1

! ######################################################################
  DO i=1,particle_number_local

! #ifdef USE_MPI
!      particle_local_grid_posx = l_q(i,1)                         /dx_loc +1 - ims_offset_i
!      particle_local_grid_posz = l_q(i,3)                         /dz_loc +1 - ims_offset_k
! #else
!      particle_local_grid_posx = l_q(i,1)                         /dx_loc +1 
!      particle_local_grid_posz = l_q(i,3)                         /dz_loc +1
! #endif
!      particle_local_grid_posy = (l_q(i,2)-g(2)%nodes(jmin_part)) /dy_loc +jmin_part  

! !###################################################################
! !Calculating g_ps AFTER particles are shifted
! !Grid is locally per processor but particle positions is global
! !###################################################################
!      g_p(1)= floor(particle_local_grid_posx)       !tracer position to the left (x1)
!      g_p(2)= g_p(1)+1               !to the right (x2)
!      g_p(3)= (floor((l_q(i,2)-g(2)%nodes(jmin_part))/dy_loc))+jmin_part       !to the bottom 
!      g_p(4)= g_p(3)+1               !to the top (y2)
!      IF     (g(3)%size .NE. 1) THEN
!         g_p(5)= floor(particle_local_grid_posz)       !front side
!         g_p(6)= g_p(5)+1               !back side
!      ELSEIF (g(3)%size .EQ. 1) THEN
!         g_p(5)= 1     !1
!         g_p(6)= 1     !1
!      END IF

! !###################################################################
! !Lenght between g_ps and particles
! !###################################################################
!      length_g_p(1) = particle_local_grid_posx - g_p(1)  !legnth x between x(i) and p
!      length_g_p(2) = g_p(2) - particle_local_grid_posx
!      length_g_p(3) = particle_local_grid_posy-g_p(3)
!      length_g_p(4) = g_p(4)-particle_local_grid_posy
!      IF ( g(3)%size .NE. 1 ) THEN
!         length_g_p(5)=particle_local_grid_posz - g_p(5)  !length between z(i) and p
!         length_g_p(6)=g_p(6) - particle_local_grid_posz  !length between z(i+1) and p
!      END IF

! !###################################################################
! !2 linear cubes for x and y
! !###################################################################
!      cube_g_p(1) = length_g_p(1) *length_g_p(3) ! cubes
!      cube_g_p(2) = length_g_p(1) *length_g_p(4) ! be carefull multiply other side cube of grid for correct interpolation
!      cube_g_p(3) = length_g_p(4) *length_g_p(2)
!      cube_g_p(4) = length_g_p(2) *length_g_p(3)

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
