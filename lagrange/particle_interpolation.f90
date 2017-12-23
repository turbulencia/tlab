#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Here particles are inside the cpu grid.
!#
!########################################################################
SUBROUTINE PARTICLE_INTERPOLATION &
     (iflag, nx,ny,nz,nvar, data_in, data_out, l_q, y, wrk1d, grid_start, grid_end)
  
  USE DNS_CONSTANTS,  ONLY : efile
  USE DNS_TYPES,      ONLY : pointers_dt
  USE DNS_GLOBAL,     ONLY : isize_particle
  USE DNS_GLOBAL,     ONLY : g
  USE LAGRANGE_GLOBAL,ONLY : jmin_part
#ifdef USE_MPI
  USE DNS_MPI, ONLY: ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE
#include "integers.h"

  TINTEGER iflag, nx,ny,nz,nvar, grid_start, grid_end
  TYPE(pointers_dt), DIMENSION(nvar)             :: data_in     
  TYPE(pointers_dt), DIMENSION(nvar)             :: data_out
  TREAL,             DIMENSION(*)                :: y, wrk1d
  TREAL,             DIMENSION(isize_particle,3) :: l_q

! -------------------------------------------------------------------
  TREAL length_g_p(6), cube_g_p(4)
  TINTEGER  gridpoint(6)
  TINTEGER i, j
  TREAL particle_local_grid_posx, particle_local_grid_posy, particle_local_grid_posz

  TREAL, DIMENSION(:,:,:), POINTER :: tmp

! ######################################################################
  IF  ( g(3)%size .NE. 1) THEN ! 3D case

     DO i = grid_start,grid_end !loop over all particles inside the grid (no halo)

#ifdef USE_MPI
        particle_local_grid_posx = l_q(i,1)/wrk1d(1) + 1 - ims_offset_i !ims_pro_i*imax
        particle_local_grid_posy = ((l_q(i,2)-y(jmin_part))/wrk1d(2))+jmin_part  
        particle_local_grid_posz = l_q(i,3)/wrk1d(3) + 1 - ims_offset_k !ims_pro_k*kmax
#else
        particle_local_grid_posx = l_q(i,1)/wrk1d(1) + 1 
        particle_local_grid_posy = ((l_q(i,2)-y(jmin_part))/wrk1d(2))+jmin_part  
        particle_local_grid_posz = l_q(i,3)/wrk1d(3) + 1
#endif

        gridpoint(1)= floor(particle_local_grid_posx)       !tracer position to the left (x1)
        gridpoint(2)= gridpoint(1)+1                        !to the right (x2)
        gridpoint(3)= (floor((l_q(i,2)-y(jmin_part))/wrk1d(2)))+jmin_part       !to the bottom 
        gridpoint(4)= gridpoint(3)+1                        !to the top (y2)
        gridpoint(5)= floor(particle_local_grid_posz)       !front side
        gridpoint(6)= gridpoint(5)+1                        !back side

! ###################################################################
! Interpolation
! ###################################################################
        length_g_p(1) = particle_local_grid_posx - gridpoint(1)  !legnth x between x(i) and p
        length_g_p(2) = gridpoint(2) - particle_local_grid_posx
        length_g_p(3) = particle_local_grid_posy - gridpoint(3)
        length_g_p(4) = gridpoint(4) - particle_local_grid_posy
        length_g_p(5) = particle_local_grid_posz - gridpoint(5)  !length between z(i) and p
        length_g_p(6) = gridpoint(6) - particle_local_grid_posz  !length between z(i+1) and p

        cube_g_p(1) = length_g_p(1) *length_g_p(3) ! cubes
        cube_g_p(2) = length_g_p(1) *length_g_p(4) !  be carefull multiply other side cube of grid for correct interpolation
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

        IF ( iflag .EQ. 1 ) THEN
           gridpoint(1)=1
           gridpoint(2)=2
        ELSE IF ( iflag .EQ. 2 ) THEN
           gridpoint(5)=1
           gridpoint(6)=2
        ELSE IF ( iflag .EQ. 3 ) THEN
           gridpoint(1)=1
           gridpoint(2)=2
           gridpoint(5)=1
           gridpoint(6)=2
        ENDIF
        
! ###################################################################
! Set the field arrays into the particle
! Two bilinear calculation for each k direction (gridpoint(5) and gridpoint(6)
! Then multipled by (1-length) for Trilinear aspect
! ###################################################################
        DO j = 1,nvar
           tmp(1:nx,1:ny,1:nz) => data_in(j)%field(:)

           data_out(j)%field(i) = data_out(j)%field(i) +&
                ((cube_g_p(3) *tmp(gridpoint(1),gridpoint(3),gridpoint(5)) &
                 +cube_g_p(4) *tmp(gridpoint(1),gridpoint(4),gridpoint(5)) &
                 +cube_g_p(1) *tmp(gridpoint(2),gridpoint(4),gridpoint(5)) &
                 +cube_g_p(2) *tmp(gridpoint(2),gridpoint(3),gridpoint(5)))*length_g_p(6)) &
               +((cube_g_p(3) *tmp(gridpoint(1),gridpoint(3),gridpoint(6)) &
                 +cube_g_p(4) *tmp(gridpoint(1),gridpoint(4),gridpoint(6)) &
                 +cube_g_p(1) *tmp(gridpoint(2),gridpoint(4),gridpoint(6)) &
                 +cube_g_p(2) *tmp(gridpoint(2),gridpoint(3),gridpoint(6)))*length_g_p(5))
        ENDDO
        
     END DO
     
! ######################################################################
  ELSE !2D case

     DO i=grid_start,grid_end !loop over all particles inside the grid (no halo)

#ifdef USE_MPI
        particle_local_grid_posx = l_q(i,1)/wrk1d(1) + 1 - ims_offset_i !ims_pro_i*imax
        particle_local_grid_posy = ((l_q(i,2)-y(jmin_part))/wrk1d(2))+jmin_part  
#else
        particle_local_grid_posx = l_q(i,1)/wrk1d(1) + 1 
        particle_local_grid_posy = ((l_q(i,2)-y(jmin_part))/wrk1d(2))+jmin_part  
#endif

!Calculating gridpoints AFTER particles are shifted
        gridpoint(1)= floor(particle_local_grid_posx)       !tracer position to the left (x1)
        gridpoint(2)= gridpoint(1)+1               !to the right (x2)
        gridpoint(3)= (floor((l_q(i,2)-y(jmin_part))/wrk1d(2)))+jmin_part       !to the bottom 
        gridpoint(4)= gridpoint(3)+1               !to the top (y2)
        gridpoint(5)=1
        gridpoint(6)=1

! ###################################################################
! Own Interpolation
! ###################################################################
        length_g_p(1) = particle_local_grid_posx - gridpoint(1)  !legnth x between x(i) and p
        length_g_p(2) = gridpoint(2) - particle_local_grid_posx
        length_g_p(3) = particle_local_grid_posy - gridpoint(3)
        length_g_p(4) = gridpoint(4) - particle_local_grid_posy

        cube_g_p(1) = length_g_p(1) *length_g_p(3) ! cubes
        cube_g_p(2) = length_g_p(1) *length_g_p(4) !  be carefull multiply other side cube of grid for correct interpolation
        cube_g_p(3) = length_g_p(4) *length_g_p(2)
        cube_g_p(4) = length_g_p(2) *length_g_p(3)

!Safety check
!cube_g_p(5)=cube_g_p(1)+cube_g_p(2)+cube_g_p(3)+cube_g_p(4)
!IF  (cube_g_p(5) .GT. 1) THEN
!    print*,'zu grosse wuerfel'
!END IF

! ###################################################################
! Two bilinear calculation for each k direction (gridpoint(5) and gridpoint(6)
! Then multipled by (1-length) for Trilinear aspect
! ###################################################################
        DO j = 1,nvar
           tmp(1:nx,1:ny,1:nz) => data_in(j)%field(:)

           data_out(j)%field(i) = data_out(j)%field(i) +&
                (cube_g_p(3) *tmp(gridpoint(1),gridpoint(3),gridpoint(5)) &
                +cube_g_p(4) *tmp(gridpoint(1),gridpoint(4),gridpoint(5)) &
                +cube_g_p(1) *tmp(gridpoint(2),gridpoint(4),gridpoint(5)) &
                +cube_g_p(2) *tmp(gridpoint(2),gridpoint(3),gridpoint(5)))        
        ENDDO
        
     END DO
     
  ENDIF

  RETURN
END SUBROUTINE PARTICLE_INTERPOLATION
