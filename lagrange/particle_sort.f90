#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2014/02 - L. Muessle
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!#  Sorting of particles into grid-west-east
!#  Sending to processors in west and east
!#  Sorting of particles into grid-south-north
!#  Sending to processors in south and north
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE PARTICLE_SORT(x,z, nzone_grid, nzone_west, nzone_east,nzone_south,nzone_north,&
                          x_or_z, h_particle, particle_id, particle )    

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_txc_field
  USE DNS_GLOBAL, ONLY: isize_particle, inb_particle
  USE LAGRANGE_GLOBAL, ONLY : particle_number, nzone_max

#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(isize_particle,inb_particle) :: particle
  TREAL, DIMENSION(isize_particle,inb_particle) :: h_particle

  TREAL, DIMENSION(*)     :: x,z
  TREAL                   dx_grid, dz_grid
  

  TREAL dummy, lower_limit, upper_limit
  TINTEGER nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north
  TINTEGER nzone_west_south, nzone_east_north
  TINTEGER counter_swap
  INTEGER(8), DIMENSION(isize_particle)  :: particle_id
  INTEGER(8) dummy_2
  TINTEGER i, j, k, x_or_z

      nzone_west_south=0
      nzone_east_north=0

      nzone_west=0
      nzone_east=0
      nzone_grid=0


      dx_grid= x(2)-x(1)  ! Distance between gridpoints
      dz_grid= z(2)-z(1) 

  IF (x_or_z .EQ. 1) THEN   !Sort in West-East direction
    lower_limit=x(imax*ims_pro_i+1)  ! lower_limit is West
    upper_limit=x(imax*(ims_pro_i+1))+dx_grid  !upper_limit is East
  

  ELSEIF (x_or_z .EQ. 3) THEN !Sort in South-North direction
    lower_limit=z(kmax*ims_pro_k+1) !lower_limit is south
    upper_limit=z(kmax*(ims_pro_k+1))+dz_grid !upper_limit is north
  END IF
  

!#######################################################################
!Sorting structure grid-west-east or grid-south-north
!First Algorythm sorts all grid-particle into first part of particle
!#######################################################################

  i=1 !Starting point of sorting algorythm
  j=ims_size_p(ims_pro+1)  !End point of sorting algorythm

  DO WHILE (i .LT. j )
    IF ( particle(i,x_or_z) .LT. lower_limit) THEN !If particle is out to West
  
      counter_swap=0 !Particle must be swapped
      DO WHILE (counter_swap .eq. 0) !Enter in the swaping region
        IF ( (particle(j,x_or_z) .LT. lower_limit) .OR. (particle(j,x_or_z) .GT. upper_limit)  ) THEN !Partcile here in right area
          j=j-1 !Go to next particle in the swaping region
          IF (i .EQ. j) THEN !You finished your upwards loop
             counter_swap=1
          ENDIF
        ELSE ! Forund a particle which belongs to grid so SWAP
          DO k=1,inb_particle
            dummy=particle(i,k)
            particle(i,k)=particle(j,k)
            particle(j,k)=dummy

            dummy_2=particle_id(i)
            particle_id(i)=particle_id(j)
            particle_id(j)=dummy_2
           
            dummy=h_particle(i,k)
            h_particle(i,k)=h_particle(j,k)
            h_particle(j,k)=dummy
          END DO
  
          j=j-1 
          nzone_grid=nzone_grid+1
          counter_swap=1
        END IF
      END DO
      i=i+1 !Go to next particle

    ELSEIF (particle(i,x_or_z) .GT. upper_limit) THEN ! IF particle is out to the east
      counter_swap=0            
      DO WHILE (counter_swap .eq. 0)  !Same procedure as above for east
      !  IF (i .EQ. j ) THEN
        !  counter_swap=1
        IF ( (particle(j,x_or_z) .LT. lower_limit) .OR. (particle(j,x_or_z) .GT. upper_limit) ) THEN
          j=j-1            
          IF (i .EQ. j) THEN !You finished your upwards loop
             counter_swap=1
          ENDIF
       ELSE                    
          DO k=1,inb_particle
            dummy=particle(i,k)
            particle(i,k)=particle(j,k)
            particle(j,k)=dummy

            dummy_2=particle_id(i)
            particle_id(i)=particle_id(j)
            particle_id(j)=dummy_2
           
            dummy=h_particle(i,k)
            h_particle(i,k)=h_particle(j,k)
            h_particle(j,k)=dummy
          END DO
                                
          j=j-1                 
          nzone_grid=nzone_grid+1
          counter_swap=1
        END IF
      END DO
      i=i+1

    ELSE  ! Particle is in the grid


      i=i+1
      nzone_grid=nzone_grid+1

    END IF

  END DO

!Last particle might not be properly checked. We check it here.
  IF (i .EQ. j) THEN !Probably not needed
     IF ( (particle(i,x_or_z) .LT. lower_limit) .OR. (particle(i,x_or_z) .GT. upper_limit) ) THEN !Particle out the grid
        !Do nothing
     ELSE 
        nzone_grid=nzone_grid+1 !The last particle was not checked but it is in the grid
     END IF

  END IF



  j=ims_size_p(ims_pro+1)
  i=nzone_grid+1 
 
  DO WHILE (i .LT. j )
    IF (particle(i,x_or_z) .LT. lower_limit) THEN
      i=i+1
      nzone_west_south=nzone_west_south+1
    ELSE  !particle is out to the east
      counter_swap=0
      DO WHILE (counter_swap .eq. 0)
        IF (particle(j,x_or_z) .LT. lower_limit) then !if particle is out to west
          DO k=1,inb_particle
            dummy=particle(i,k)
            particle(i,k)=particle(j,k)
            particle(j,k)=dummy

            dummy_2=particle_id(i)
            particle_id(i)=particle_id(j)
            particle_id(j)=dummy_2
           
            dummy=h_particle(i,k)
            h_particle(i,k)=h_particle(j,k)
            h_particle(j,k)=dummy
          END DO
          

          counter_swap=1
          j=j-1
          nzone_west_south=nzone_west_south+1
        
        
        ELSE
          j=j-1
          IF (i .EQ. j) THEN !You finished your upwards loop
             counter_swap=1
          ENDIF
        END IF
      END DO
      i=i+1
    END IF
  END DO

!Last particle might not be properly checked. We check it here.
  IF (i .EQ. j) THEN !Probably not needed
     IF ( particle(i,x_or_z) .GT. upper_limit ) THEN !Particle is out east
        !Do nothing
     ELSE 
          nzone_west_south=nzone_west_south+1
     END IF
 END IF

  !Calculating the number of particles send to east or north
  nzone_east_north=ims_size_p(ims_pro+1) - nzone_grid - nzone_west_south


!  ims_size_p(ims_pro+1) = nzone_grid !For better readability in particle_send_recv

  IF (x_or_z .EQ. 1) THEN
    nzone_west=nzone_west_south
    nzone_east=nzone_east_north 
 
  ELSEIF (x_or_z .EQ. 3) THEN
    nzone_south=nzone_west_south
    nzone_north=nzone_east_north
  END IF

  RETURN
END SUBROUTINE PARTICLE_SORT
