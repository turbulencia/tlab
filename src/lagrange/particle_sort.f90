#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

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
SUBROUTINE PARTICLE_SORT(x_or_z, l_g,l_q, l_hq, &
     nzone_grid,nzone_west,nzone_east,nzone_south,nzone_north)    

  USE DNS_GLOBAL, ONLY : imax,kmax
  USE DNS_GLOBAL, ONLY : isize_particle, inb_part_array, inb_part
  USE DNS_GLOBAL, ONLY : g
  USE LAGRANGE_GLOBAL, ONLY: particle_dt
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER nzone_grid, nzone_west, nzone_east, nzone_south, nzone_north, x_or_z
  TYPE(particle_dt)                  :: l_g
  TREAL, DIMENSION(isize_particle,*) :: l_q
  TREAL, DIMENSION(isize_particle,*) :: l_hq

! -------------------------------------------------------------------
  TREAL dx_grid, dz_grid
  TREAL dummy, lower_limit, upper_limit
  TINTEGER nzone_west_south, nzone_east_north
  TINTEGER counter_swap, idummy
  INTEGER(8) idummy8
  TINTEGER i, j, k
  
!#######################################################################
  nzone_west_south=0
  nzone_east_north=0
  
  nzone_west=0
  nzone_east=0
  nzone_grid=0
  
  dx_grid = g(1)%nodes(2)-g(1)%nodes(1)  ! Distance between gridpoints, to deal with periodicity
  dz_grid = g(3)%nodes(2)-g(3)%nodes(1) 
  
  IF     (x_or_z .EQ. 1) THEN   !Sort in West-East direction
     lower_limit=g(1)%nodes(ims_offset_i +1)             !lower_limit is West
     upper_limit=g(1)%nodes(ims_offset_i +imax) +dx_grid !upper_limit is East
     
  ELSEIF (x_or_z .EQ. 3) THEN !Sort in South-North direction
     lower_limit=g(3)%nodes(ims_offset_k +1)             !lower_limit is south
     upper_limit=g(3)%nodes(ims_offset_k +kmax) +dz_grid !upper_limit is north
  END IF

!#######################################################################
!Sorting structure grid-west-east or grid-south-north
!First Algorythm sorts all grid-particle into first part of particle
!#######################################################################
  i=1 !Starting point of sorting algorythm
  j=l_g%np  !End point of sorting algorythm

  DO WHILE (i .LT. j )
     IF ( l_q(i,x_or_z) .LT. lower_limit) THEN !If particle is out to West

        counter_swap=0 !Particle must be swapped
        DO WHILE (counter_swap .eq. 0) !Enter in the swaping region
           IF ( (l_q(j,x_or_z) .LT. lower_limit) .OR. (l_q(j,x_or_z) .GT. upper_limit)  ) THEN !Partcile here in right area
              j=j-1 !Go to next particle in the swaping region
              IF (i .EQ. j) THEN !You finished your upwards loop
                 counter_swap=1
              ENDIF
           ELSE ! Forund a particle which belongs to grid so SWAP
              idummy      =l_g%nodes(i)
              l_g%nodes(i)=l_g%nodes(j)
              l_g%nodes(j)=idummy
              
              idummy8    =l_g%tags(i)
              l_g%tags(i)=l_g%tags(j)
              l_g%tags(j)=idummy8

              DO k=1,inb_part_array
                 dummy   =l_q(i,k)
                 l_q(i,k)=l_q(j,k)
                 l_q(j,k)=dummy
              ENDDO
              
              DO k=1,inb_part
                 dummy    =l_hq(i,k)
                 l_hq(i,k)=l_hq(j,k)
                 l_hq(j,k)=dummy
              END DO

              j=j-1 
              nzone_grid=nzone_grid+1
              counter_swap=1
           END IF
        END DO
        i=i+1 !Go to next particle

     ELSEIF (l_q(i,x_or_z) .GT. upper_limit) THEN ! IF particle is out to the east
        counter_swap=0            
        DO WHILE (counter_swap .eq. 0)  !Same procedure as above for east
!  IF (i .EQ. j ) THEN
!  counter_swap=1
           IF ( (l_q(j,x_or_z) .LT. lower_limit) .OR. (l_q(j,x_or_z) .GT. upper_limit) ) THEN
              j=j-1            
              IF (i .EQ. j) THEN !You finished your upwards loop
                 counter_swap=1
              ENDIF
           ELSE                    
              idummy      =l_g%nodes(i)
              l_g%nodes(i)=l_g%nodes(j)
              l_g%nodes(j)=idummy
              
              idummy8    =l_g%tags(i)
              l_g%tags(i)=l_g%tags(j)
              l_g%tags(j)=idummy8

              DO k=1,inb_part_array
                 dummy   =l_q(i,k)
                 l_q(i,k)=l_q(j,k)
                 l_q(j,k)=dummy
              ENDDO
              
              DO k=1,inb_part
                 dummy    =l_hq(i,k)
                 l_hq(i,k)=l_hq(j,k)
                 l_hq(j,k)=dummy
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
     IF ( (l_q(i,x_or_z) .LT. lower_limit) .OR. (l_q(i,x_or_z) .GT. upper_limit) ) THEN !Particle out the grid
!Do nothing
     ELSE 
        nzone_grid=nzone_grid+1 !The last particle was not checked but it is in the grid
     END IF

  END IF

! -------------------------------------------------------------------
  j=l_g%np
  i=nzone_grid+1 

  DO WHILE (i .LT. j )
     IF (l_q(i,x_or_z) .LT. lower_limit) THEN
        i=i+1
        nzone_west_south=nzone_west_south+1
     ELSE  !particle is out to the east
        counter_swap=0
        DO WHILE (counter_swap .eq. 0)
           IF (l_q(j,x_or_z) .LT. lower_limit) then !if particle is out to west
              idummy      =l_g%nodes(i)
              l_g%nodes(i)=l_g%nodes(j)
              l_g%nodes(j)=idummy
              
              idummy8    =l_g%tags(i)
              l_g%tags(i)=l_g%tags(j)
              l_g%tags(j)=idummy8

              DO k=1,inb_part_array
                 dummy   =l_q(i,k)
                 l_q(i,k)=l_q(j,k)
                 l_q(j,k)=dummy
              ENDDO
              
              DO k=1,inb_part
                 dummy    =l_hq(i,k)
                 l_hq(i,k)=l_hq(j,k)
                 l_hq(j,k)=dummy
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
     IF ( l_q(i,x_or_z) .GT. upper_limit ) THEN !Particle is out east
!Do nothing
     ELSE 
        nzone_west_south=nzone_west_south+1
     END IF
  END IF

!Calculating the number of particles send to east or north
  nzone_east_north = l_g%np - nzone_grid - nzone_west_south

  IF (x_or_z .EQ. 1) THEN
     nzone_west=nzone_west_south
     nzone_east=nzone_east_north 

  ELSEIF (x_or_z .EQ. 3) THEN
     nzone_south=nzone_west_south
     nzone_north=nzone_east_north
  END IF

  RETURN
END SUBROUTINE PARTICLE_SORT
