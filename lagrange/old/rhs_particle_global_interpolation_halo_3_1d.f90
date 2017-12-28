!########################################################################
!# Tool/Library
!#
!########################################################################
!# DESCRIPTION
!#
!# Particle set as tracer.
!# Determing velocity with bilinear interpolation.
!# 1D stands for general interpolation of a scalar variable
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE  RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3_1D &
    (halo_field,l_q,particle_property,y,wrk1d,halo_start, halo_end)

USE DNS_GLOBAL, ONLY: imax,jmax,kmax
USE DNS_GLOBAL, ONLY: isize_particle
USE LAGRANGE_GLOBAL, ONLY: jmin_part
#ifdef USE_MPI
   USE DNS_MPI, ONLY: ims_pro_i, ims_pro_k, ims_pro
#endif



IMPLICIT NONE
#include "integers.h"

  TREAL, DIMENSION(*)                 :: y
  TREAL, DIMENSION(2,jmax,2)          :: halo_field
  TREAL, DIMENSION(isize_particle,3)  :: l_q
  TREAL, DIMENSION(isize_particle)    :: particle_property
  TREAL, DIMENSION(*),intent(in)      :: wrk1d
  TREAL length_g_p(6), cube_g_p(4)
  TINTEGER  gridpoint(6)
  TREAL particle_local_grid_posx, particle_local_grid_posy, particle_local_grid_posz
  TINTEGER i, halo_end, halo_start


  DO i=halo_start,halo_end
 
#ifdef USE_MPI
      particle_local_grid_posx = l_q(i,1)/wrk1d(1) + 1 - ims_pro_i*imax
      particle_local_grid_posy = ((l_q(i,2)-y(jmin_part))/wrk1d(2))+jmin_part  
      particle_local_grid_posz = l_q(i,3)/wrk1d(3) + 1 - ims_pro_k*kmax
#else
      particle_local_grid_posx = l_q(i,1)/wrk1d(1) + 1 
      particle_local_grid_posy = ((l_q(i,2)-y(jmin_part))/wrk1d(2))+jmin_part  
      particle_local_grid_posz = l_q(i,3)/wrk1d(3) + 1 
#endif

     

      !Calculating gridpoints AFTER particles are shifted
      gridpoint(1)= floor(particle_local_grid_posx)      !tracer position to the left (x1)
      gridpoint(2)= gridpoint(1)+1               !to the right (x2)
      gridpoint(3)= (floor((l_q(i,2)-y(jmin_part))/wrk1d(2)))+jmin_part       !to the bottom 
      gridpoint(4)= gridpoint(3)+1               !to the top (y2)
      gridpoint(5)= floor(particle_local_grid_posz)       !front side
      gridpoint(6)= gridpoint(5)+1               !back side

    !#######################################################################
    !Interpolation
    !#######################################################################
    
      length_g_p(1)=particle_local_grid_posx - gridpoint(1) !legnth x between x(i) and p
      length_g_p(2)=gridpoint(2) - particle_local_grid_posx
      length_g_p(3)=particle_local_grid_posy-gridpoint(3)
      length_g_p(4)=gridpoint(4)-particle_local_grid_posy
      length_g_p(5)=particle_local_grid_posz - gridpoint(5)  !length between z(i) and p
      length_g_p(6)=gridpoint(6) - particle_local_grid_posz  !length between z(i+1) and p
   



      cube_g_p(1)=length_g_p(1)*length_g_p(3) ! cubes
      cube_g_p(2)=length_g_p(1)*length_g_p(4)
      cube_g_p(3)=length_g_p(4)*length_g_p(2)
      cube_g_p(4)=length_g_p(2)*length_g_p(3)
  
      !Safety check
      !cube_g_p(5)=cube_g_p(1)+cube_g_p(2)+cube_g_p(3)+cube_g_p(4)
      !IF  (cube_g_p(5) .GT. 1) THEN
      !    print*,'zu grosse wuerfel'
      !END IF

        gridpoint(1)=1
        gridpoint(2)=2
        gridpoint(5)=1
        gridpoint(6)=2
  
     ! ###################################################################
     ! Two bilinear calculation for each k direction (gridpoint(5) and gridpoint(6)
     ! Then multipled by (1-length) for Trilinear aspect
     ! ###################################################################
     
      
     particle_property(i)= particle_property(i) + &
     ((cube_g_p(3)*halo_field(gridpoint(1),gridpoint(3),gridpoint(5)) &
     +cube_g_p(4)*halo_field(gridpoint(1),gridpoint(4),gridpoint(5)) &
     +cube_g_p(1)*halo_field(gridpoint(2),gridpoint(4),gridpoint(5)) &
     +cube_g_p(2)*halo_field(gridpoint(2),gridpoint(3),gridpoint(5)))*length_g_p(6)) &
     +((cube_g_p(3)*halo_field(gridpoint(1),gridpoint(3),gridpoint(6)) &
     +cube_g_p(4)*halo_field(gridpoint(1),gridpoint(4),gridpoint(6)) & 
     +cube_g_p(1)*halo_field(gridpoint(2),gridpoint(4),gridpoint(6)) & 
     +cube_g_p(2)*halo_field(gridpoint(2),gridpoint(3),gridpoint(6)))*length_g_p(5))

 

    
   
  END DO    
    
  RETURN
END SUBROUTINE RHS_PARTICLE_GLOBAL_INTERPOLATION_HALO_3_1D
