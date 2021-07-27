#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Sorting structure grid-halo_x-halo_z-halo_diagonal
!# First algorithm sorts all grid-particle into first part of particle
!########################################################################
SUBROUTINE PARTICLE_SORT_HALO(l_g,l_q, nvar,data, grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal)

  USE DNS_TYPES,      ONLY : pointers_dt
  USE DNS_GLOBAL,     ONLY : isize_particle, inb_part_array
  USE DNS_GLOBAL,     ONLY : g
  USE LAGRANGE_GLOBAL,ONLY : particle_dt

#ifdef USE_MPI
  USE DNS_GLOBAL, ONLY : imax,kmax
  USE TLAB_MPI_VARS,    ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE

  TINTEGER grid_zone, halo_zone_x, halo_zone_z, halo_zone_diagonal, nvar
  TYPE(pointers_dt), DIMENSION(nvar)                        :: data
  TYPE(particle_dt)                                         :: l_g
  TREAL,             DIMENSION(isize_particle,inb_part_array) :: l_q

! -------------------------------------------------------------------
  TREAL dummy, right_limit, upper_limit
  TINTEGER counter_swap, idummy
  INTEGER(8) idummy8
  TINTEGER i, j, k 

!#######################################################################
#ifdef USE_MPI
  right_limit=g(1)%nodes(ims_offset_i+imax)  !right_limit is east
  upper_limit=g(3)%nodes(ims_offset_k+kmax)  !upper_limit is north
#else
  right_limit=g(1)%nodes(g(1)%size)
  upper_limit=g(3)%nodes(g(3)%size)
#endif

! Initialize counters  
  grid_zone=0
  halo_zone_x=0
  halo_zone_z=0
  halo_zone_diagonal=0

!#######################################################################
! Grid zone
  i=1                      !Starting point of sorting algorythm
  j=l_g%np  !End point of sorting algorythm

  DO WHILE (i .LT. j )

     IF ( l_q(i,1) .GT. right_limit) THEN !If particle is out to East
        counter_swap=0 !Particle must be swapped
        DO WHILE (counter_swap .eq. 0) !Enter in the swaping region
           IF ( (l_q(j,1) .GT. right_limit) .OR. (l_q(j,3) .GT. upper_limit)  ) THEN !Partcile here in right area
              j=j-1 !Go to next particle in the swaping region
              IF (i .EQ. j) THEN !You finished your upwards loop
                 counter_swap=1
              ENDIF
           ELSE ! Found a particle which belongs to grid so SWAP
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
              END DO

              DO k=1,nvar
                 dummy           =data(k)%field(i)
                 data(k)%field(i)=data(k)%field(j)
                 data(k)%field(j)=dummy
              END DO

              j=j-1 
              grid_zone=grid_zone+1
              counter_swap=1

           END IF
        END DO
        i=i+1 !Go to next particle

     ELSEIF (l_q(i,3) .GT. upper_limit) THEN ! IF l_q is out to the north
        counter_swap=0            
        DO WHILE (counter_swap .eq. 0)  !Same procedure as above for East
           IF ( (l_q(j,1) .GT. right_limit) .OR. (l_q(j,3) .GT. upper_limit) ) THEN
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
              
              DO k=1,nvar
                 dummy           =data(k)%field(i)
                 data(k)%field(i)=data(k)%field(j)
                 data(k)%field(j)=dummy
              END DO

              j=j-1                 
              grid_zone=grid_zone+1
              counter_swap=1

           END IF
        END DO
        i=i+1

     ELSE  ! Particle is in the grid
        i=i+1
        grid_zone=grid_zone+1

     END IF

  END DO

!Last particle might not be properly checked. We check it here. 
  IF (i .EQ. j) THEN !Probably not needed
     IF ( (l_q(i,1) .GT. right_limit) .OR. (l_q(i,3) .GT. upper_limit) ) THEN !Particle out the grid
!Do nothing
     ELSE 
        grid_zone=grid_zone+1 !The last particle was not checked but it is in the grid
     END IF
  ENDIF
!Possible optimization to avoid if i. EQ. j. Not completed.
! ELSE IF (i-1 .GT. j+1) THEN 
!   j=j+1
!   i=i-1
!   grid_zone=grid_zone-1 
!   IF ( (l_q(j,1) .GT. right_limit) .OR. l_q(j,3) .GT. upper_limit) THEN
!          DO k=1,3
!           dummy=l_q(i,k)
!           l_q(i,k)=l_q(j,k)
!           l_q(j,k)=dummy

!           idummy8=l_g%tags(i)
!           l_g%tags(i)=l_g%tags(j)
!           l_g%tags(j)=idummy8
!         END DO
!   END IF 

! -------------------------------------------------------------------
  i=grid_zone+1 
  j=l_g%np

  DO WHILE (i .LT. j )
     IF (l_q(i,3) .GT. upper_limit) THEN  !If particle is out North
        counter_swap=0
        DO WHILE (counter_swap .eq. 0)
           IF (l_q(j,3) .GT. upper_limit) THEN  !Particle is right here
              j=j-1
              IF (i .EQ. j) THEN !You finished your upwards loop
                 counter_swap=1
              ENDIF

           ELSE  !If particle is only out to East
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
              END DO

              DO k=1,nvar
                 dummy           =data(k)%field(i)
                 data(k)%field(i)=data(k)%field(j)
                 data(k)%field(j)=dummy
              END DO

              counter_swap=1
              j=j-1
              halo_zone_x=halo_zone_x+1

           END IF
        END DO
        i=i+1
     ELSE  
        i=i+1
        halo_zone_x=halo_zone_x+1
     END IF
  END DO

!Last particle might not be properly checked. We check it here.
  IF (i .EQ. j) THEN !Probably not needed
     IF ( l_q(i,3) .GT. upper_limit ) THEN !Particle is out North
!Do nothing
     ELSE 
        halo_zone_x=halo_zone_x+1
     END IF
  END IF

! -------------------------------------------------------------------
  i=grid_zone+halo_zone_x+1 
  j=l_g%np  !End point of sorting algorythm

  DO WHILE (i .LT. j )
     IF (l_q(i,1) .GT. right_limit) THEN  !If particle is out East
        counter_swap=0
        DO WHILE (counter_swap .eq. 0)
           IF (l_q(j,1) .GT. right_limit) then !If particle is out to East
              j=j-1
              IF (i .EQ. j) THEN !You finished your upwards loop
                 counter_swap=1
              ENDIF
           ELSE  !Particle is only out to the North
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
              END DO
              
              DO k=1,nvar
                 dummy           =data(k)%field(i)
                 data(k)%field(i)=data(k)%field(j)
                 data(k)%field(j)=dummy
              END DO

              counter_swap=1
              j=j-1
              halo_zone_z=halo_zone_z+1

           END IF
        END DO
        i=i+1
     ELSE
        halo_zone_z=halo_zone_z+1
        i=i+1
     END IF

  END DO

!Last particle might not be properly checked. We check it here.
  IF (i .EQ. j) THEN !Probably not needed
     IF ( l_q(i,1) .GT. right_limit ) THEN !Particle is out East
!Do nothing
     ELSE 
        halo_zone_z=halo_zone_z+1
     END IF
  END IF

!Calculating the number of particles send to east or north
  halo_zone_diagonal = l_g%np -grid_zone -halo_zone_x -halo_zone_z
  
  RETURN
END SUBROUTINE PARTICLE_SORT_HALO
