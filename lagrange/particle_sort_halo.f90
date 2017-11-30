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
SUBROUTINE PARTICLE_SORT_HALO(x,z, nzone_grid, halo_zone_x, halo_zone_z, halo_zone_diagonal,&
                          l_hq, l_tags, l_q )    

  USE DNS_GLOBAL, ONLY : imax,kmax
  USE DNS_GLOBAL, ONLY: isize_particle, inb_particle
  USE DNS_GLOBAL, ONLY: g
  USE LAGRANGE_GLOBAL, ONLY : particle_number

#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(isize_particle,inb_particle) :: l_q
  TREAL, DIMENSION(isize_particle,inb_particle) :: l_hq

  TREAL, DIMENSION(*)     :: x,z
  

  TREAL dummy, right_limit, upper_limit
  TINTEGER nzone_grid, halo_zone_x, halo_zone_z, halo_zone_diagonal
  TINTEGER counter_swap
  INTEGER(8), DIMENSION(isize_particle)  :: l_tags
  INTEGER(8) dummy_2
  TINTEGER i, j, k 

      nzone_grid=0

      halo_zone_x=0
      halo_zone_z=0
      halo_zone_diagonal=0

#ifdef USE_MPI
    right_limit=x(imax*(ims_pro_i+1))  ! right_limit is east
    upper_limit=z(kmax*(ims_pro_k+1))  !upper_limit is north
#else
    right_limit=x(g(1)%size)  ! right_limit is east
    upper_limit=z(g(3)%size)  
#endif
!    right_limit=x(g(1)%size)  ! right_limit is east

!#######################################################################
!Sorting structure grid-halo_x-halo_z-halo_diagonal
!First Algorythm sorts all grid-particle into first part of particle
!#######################################################################

  i=1 !Starting point of sorting algorythm
#ifdef USE_MPI
  j=particle_vector(ims_pro+1)  !End point of sorting algorythm
#else
  j=particle_number
#endif
  DO WHILE (i .LT. j )
    IF ( l_q(i,1) .GT. right_limit) THEN !If particle is out to East
  
      counter_swap=0 !Particle must be swapped
      DO WHILE (counter_swap .eq. 0) !Enter in the swaping region
        IF ( (l_q(j,1) .GT. right_limit) .OR. (l_q(j,3) .GT. upper_limit)  ) THEN !Partcile here in right area
          j=j-1 !Go to next particle in the swaping region
          IF (i .EQ. j) THEN !You finished your upwards loop
             counter_swap=1
          ENDIF
        ELSE ! Forund a particle which belongs to grid so SWAP
          DO k=1,inb_particle
            dummy=l_q(i,k)
            l_q(i,k)=l_q(j,k)
            l_q(j,k)=dummy

            dummy_2=l_tags(i)
            l_tags(i)=l_tags(j)
            l_tags(j)=dummy_2
           
            dummy=l_hq(i,k)
            l_hq(i,k)=l_hq(j,k)
            l_hq(j,k)=dummy
          END DO
  
          j=j-1 
          nzone_grid=nzone_grid+1
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
          DO k=1,inb_particle
            dummy=l_q(i,k)
            l_q(i,k)=l_q(j,k)
            l_q(j,k)=dummy

            dummy_2=l_tags(i)
            l_tags(i)=l_tags(j)
            l_tags(j)=dummy_2
           
            dummy=l_hq(i,k)
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
     IF ( (l_q(i,1) .GT. right_limit) .OR. (l_q(i,3) .GT. upper_limit) ) THEN !Particle out the grid
        !Do nothing
     ELSE 
        nzone_grid=nzone_grid+1 !The last particle was not checked but it is in the grid
     END IF
  ENDIF
!Possible optimization to avoid if i. EQ. j. Not completed.
  ! ELSE IF (i-1 .GT. j+1) THEN 
  !   j=j+1
  !   i=i-1
  !   nzone_grid=nzone_grid-1 
  !   IF ( (l_q(j,1) .GT. right_limit) .OR. l_q(j,3) .GT. upper_limit) THEN
  !          DO k=1,3
  !           dummy=l_q(i,k)
  !           l_q(i,k)=l_q(j,k)
  !           l_q(j,k)=dummy

  !           dummy_2=l_tags(i)
  !           l_tags(i)=l_tags(j)
  !           l_tags(j)=dummy_2
           
  !           dummy=l_hq(i,k)
  !           l_hq(i,k)=l_hq(j,k)
  !           l_hq(j,k)=dummy
  !         END DO
  !   END IF 
    

  

#ifdef USE_MPI
  j=particle_vector(ims_pro+1)  !End point of sorting algorythm
#else
  j=particle_number
#endif

  i=nzone_grid+1 
 
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
          DO k=1,inb_particle
            dummy=l_q(i,k)
            l_q(i,k)=l_q(j,k)
            l_q(j,k)=dummy

            dummy_2=l_tags(i)
            l_tags(i)=l_tags(j)
            l_tags(j)=dummy_2
           
            dummy=l_hq(i,k)
            l_hq(i,k)=l_hq(j,k)
            l_hq(j,k)=dummy
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


#ifdef USE_MPI
  j=particle_vector(ims_pro+1)  !End point of sorting algorythm
#else
  j=particle_number
#endif

  i=nzone_grid+halo_zone_x+1 
 
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
          DO k=1,inb_particle
            dummy=l_q(i,k)
            l_q(i,k)=l_q(j,k)
            l_q(j,k)=dummy

            dummy_2=l_tags(i)
            l_tags(i)=l_tags(j)
            l_tags(j)=dummy_2
           
            dummy=l_hq(i,k)
            l_hq(i,k)=l_hq(j,k)
            l_hq(j,k)=dummy
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
#ifdef USE_MPI
  halo_zone_diagonal=particle_vector(ims_pro+1) - nzone_grid - halo_zone_x - halo_zone_z
#else
  halo_zone_diagonal=particle_number - nzone_grid - halo_zone_x - halo_zone_z
#endif


  RETURN
END SUBROUTINE PARTICLE_SORT_HALO
