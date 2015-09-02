#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "LAGRANGE_DISPERSION"

!########################################################################
!# Tool/Library STATISTICS
!#
!########################################################################
!# HISTORY
!#
!# 2015/01/06 - L. Muessle
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# Post processing to determine the dispersion of particles in time
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM LAGRANGE_DISPERSION

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------

  TINTEGER  ierr
  TINTEGER  i, j, k 
  TINTEGER  dummy, n_pairs 
  TINTEGER, DIMENSION(1) :: x_seed, wildcard 
  TINTEGER, DIMENSION(:), ALLOCATABLE :: dummy_pos
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: int_rnd_number
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: l_tags
  INTEGER(8), DIMENSION(:,:), ALLOCATABLE ::  tag_neighbor_start
#ifdef USE_MPI
  INTEGER(8), DIMENSION(:,:), ALLOCATABLE :: all_tag_neighbor_start
#endif
  TREAL, DIMENSION(:), ALLOCATABLE ::  rnd_number
  TREAL, DIMENSION(:,:),    ALLOCATABLE :: l_q, l_txc
  TREAL, DIMENSION(:,:), ALLOCATABLE :: l_dispersion
  TREAL, DIMENSION(:,:), ALLOCATABLE :: end_l_dispersion
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: start_particle, end_particle
  TREAL, DIMENSION(:,:,:,:), ALLOCATABLE :: all_particle

  TINTEGER nitera_first,nitera_last

  CHARACTER*32 inifile
  CHARACTER*64 str,fname
  CHARACTER*128 line
  CHARACTER*32 bakfile

  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! Get the local information from the dns.ini
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Start',      '0',  nitera_first)
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'End',        '0',  nitera_last )


#include "dns_alloc_larrays.h"

  n_pairs=num_dispersion  !should be determined in .ini file in FUTURE


  ALLOCATE(dummy_pos(n_pairs))
  ALLOCATE(rnd_number(n_pairs))
  ALLOCATE(int_rnd_number(n_pairs)) 
  ALLOCATE(tag_neighbor_start(2,n_pairs))
#ifdef USE_MPI
  ALLOCATE(all_tag_neighbor_start(ims_npro*n_pairs,2))  !(amount_of_tags,1st_and_2nd)    
  ALLOCATE(all_particle(ims_npro*n_pairs,3,2,2))        !(amount_of_pairs,xyz,start_and_end,1st_and_2nd)
  ALLOCATE(start_particle(ims_npro*n_pairs,3,2))        !(amount_ofpairs,xyz,1st_and_2nd) 
  ALLOCATE(end_particle(ims_npro*n_pairs,3,2))          !(amount_ofpairs,xyz,1st_and_2nd)
  ALLOCATE(end_l_dispersion(ims_npro*n_pairs,2))        !(amount_of_dispersion,1st_and_2nd) 
#else
  ALLOCATE(start_particle(n_pairs,3,2))                 !(amount_of_tags,xyz,1st_and_2nd) 
  ALLOCATE(end_particle(n_pairs,3,2))                   !(amount_of_tags,xyz,1st_and_2nd) 
  ALLOCATE(end_l_dispersion(n_pairs,2))                 !(amount_of_tags,xyz,1st_and_2nd) 
#endif
  !ALLOCATE l_dispersion further down 
  all_particle(:,:,:,:)=C_0_R

  !#######################################################################
  !READ THE (FIRST) FILE
  !#######################################################################
  WRITE(fname,*) nitera_first; fname = "particle_id."//TRIM(ADJUSTL(fname))
  CALL DNS_READ_PARTICLE_TAGS(fname,l_tags)
  
  WRITE(fname,*) nitera_first; fname = "particle."//TRIM(ADJUSTL(fname))
  CALL DNS_READ_PARTICLE(fname,l_q) 

#ifdef USE_MPI
  !ALLOCATE after particle_vector is read
  ALLOCATE(l_dispersion(particle_vector(ims_pro+1),n_pairs)) 
#endif

  !#######################################################################
  !DETERMINE RANDOM PARTICLES FOR DISPERSION
  !#######################################################################
#ifdef USE_MPI
  x_seed=(/ims_pro+1/)
  CALL RANDOM_SEED(PUT=x_seed)
  DO j=1,n_pairs
    CALL RANDOM_NUMBER(rnd_number(j))
    int_rnd_number(j)=ceiling(rnd_number(j)*particle_vector(ims_pro+1))
  ENDDO
#else
  ALLOCATE(l_dispersion(particle_number,n_pairs))
  x_seed=(/1/)
  CALL RANDOM_SEED(PUT=x_seed)
  DO j=1,n_pairs
    CALL RANDOM_NUMBER(rnd_number(j))
    int_rnd_number(j)=ceiling(rnd_number(j)*particle_number)
  ENDDO
#endif
  !#######################################################################
  !CALCULATE NEAREST NEIGHBOR
  !#######################################################################
#ifdef USE_MPI
  DO j=1,n_pairs 
    DO i=1,particle_vector(ims_pro+1)
#else
  DO j=1,n_pairs 
    DO i=1,particle_number
#endif
      IF (i .EQ. int_rnd_number(j)) THEN
        l_dispersion(i,j)=999
      ELSE
        l_dispersion(i,j)=sqrt(((l_q(i,1)-l_q(int_rnd_number(j),1))**2) &
                          + ((l_q(i,2)-l_q(int_rnd_number(j),2))**2) &
                          + ((l_q(i,3)-l_q(int_rnd_number(j),3))**2))
      ENDIF
    ENDDO

#ifdef USE_MPI
  dummy=ims_pro*n_pairs+j
#else
  dummy=j
#endif
  wildcard=MINLOC(l_dispersion(:,j))  !MINLOC needs a 1-dim array
  dummy_pos(j)=wildcard(1)            !therefore use of wildcard
  tag_neighbor_start(1,j)=l_tags(int_rnd_number(j)) 
  tag_neighbor_start(2,j)=l_tags(dummy_pos(j))

  !3-dimension of 1st pair part
  start_particle(dummy,1,1)=l_q(int_rnd_number(j),1)
  start_particle(dummy,2,1)=l_q(int_rnd_number(j),2)
  start_particle(dummy,3,1)=l_q(int_rnd_number(j),3)

  !3-dimension of 2nd pair part
  start_particle(dummy,1,2)=l_q(dummy_pos(j),1)
  start_particle(dummy,2,2)=l_q(dummy_pos(j),2)
  start_particle(dummy,3,2)=l_q(dummy_pos(j),3)

  ENDDO
  !#######################################################################
  !ALL_GATHER TAG INFORMATION
  !#######################################################################
#ifdef USE_MPI
  CALL MPI_ALLGATHER(tag_neighbor_start(1,:),n_pairs,MPI_INTEGER8,all_tag_neighbor_start(:,1),n_pairs,MPI_INTEGER8,MPI_COMM_WORLD,ims_err)
  CALL MPI_ALLGATHER(tag_neighbor_start(2,:),n_pairs,MPI_INTEGER8,all_tag_neighbor_start(:,2),n_pairs,MPI_INTEGER8,MPI_COMM_WORLD,ims_err)

  !#######################################################################
  !COLLECT POSITIONS OF PAIRS (FIRST TIMESTEP)
  !#######################################################################
  CALL MPI_REDUCE(start_particle(1,1,1),all_particle(1,1,1,1),ims_npro*n_pairs*3,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ims_err)
  CALL MPI_REDUCE(start_particle(1,1,2),all_particle(1,1,1,2),ims_npro*n_pairs*3,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ims_err)

#endif
  !#######################################################################
  !READ THE (LAST) FILE
  !#######################################################################
  WRITE(fname,*) nitera_last; fname = "particle_id."//TRIM(ADJUSTL(fname))
  CALL DNS_READ_PARTICLE_TAGS(fname,l_tags)
  
  WRITE(fname,*) nitera_last; fname = "particle."//TRIM(ADJUSTL(fname))
  CALL DNS_READ_PARTICLE(fname,l_q) 

  !#######################################################################
  !SEARCH FOR PARTICLE_TAGS IN LAST TIMESTEP
  !#######################################################################
#ifdef USE_MPI
  DO i=1,particle_vector(ims_pro+1)
    DO k=1,ims_npro*n_pairs
      IF ( l_tags(i) .EQ. all_tag_neighbor_start(k,2)) THEN
          dummy=k
          end_particle(dummy,1,2)=l_q(i,1)
          end_particle(dummy,2,2)=l_q(i,2)
          end_particle(dummy,3,2)=l_q(i,3)
      ENDIF
      IF ( l_tags(i) .EQ. all_tag_neighbor_start(k,1)) THEN
          dummy=k
          end_particle(dummy,1,1)=l_q(i,1)
          end_particle(dummy,2,1)=l_q(i,2)
          end_particle(dummy,3,1)=l_q(i,3)
      ENDIF
    ENDDO
  ENDDO


#else
  DO j=1,n_pairs 
        end_particle(j,1,2)=l_q(dummy_pos(j),1)
        end_particle(j,2,2)=l_q(dummy_pos(j),2)
        end_particle(j,3,2)=l_q(dummy_pos(j),3)

        end_particle(j,1,1)=l_q(int_rnd_number(j),1)
        end_particle(j,2,1)=l_q(int_rnd_number(j),2)
        end_particle(j,3,1)=l_q(int_rnd_number(j),3)
  ENDDO
#endif

  !#######################################################################
  !COLLECT POSITIONS OF PAIRS (LAST TIMESTEP)
  !#######################################################################
#ifdef USE_MPI
  CALL MPI_REDUCE(end_particle(1,1,1),all_particle(1,1,2,1),ims_npro*n_pairs*3,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ims_err)
  CALL MPI_REDUCE(end_particle(1,1,2),all_particle(1,1,2,2),ims_npro*n_pairs*3,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ims_err)
#endif

  !#######################################################################
  !CALCULATE NEW DISPERSION
  !#######################################################################
  
#ifdef USE_MPI
  IF (ims_pro .eq. 0) THEN
    DO j=1,ims_npro*n_pairs
      DO k=1,2
        end_l_dispersion(j,k)=sqrt((all_particle(j,1,k,2)-all_particle(j,1,k,1))**2 &
                              +(all_particle(j,2,k,2)-all_particle(j,2,k,1))**2 &
                              +(all_particle(j,3,k,2)-all_particle(j,3,k,1))**2) 
      ENDDO
    ENDDO
  ENDIF

#else
  DO j=1,n_pairs
    k=1
    end_l_dispersion(j,k)=sqrt((start_particle(j,1,2)-start_particle(j,1,1))**2 &
                          +(start_particle(j,2,2)-start_particle(j,2,1))**2 &
                          +(start_particle(j,3,2)-start_particle(j,3,1))**2) 
    k=2
    end_l_dispersion(j,k)=sqrt((end_particle(j,1,2)-end_particle(j,1,1))**2 &
                          +(end_particle(j,2,2)-end_particle(j,2,1))**2 &
                          +(end_particle(j,3,2)-end_particle(j,3,1))**2) 

  ENDDO

print*, '###'
DO i=1,n_pairs
print*, end_l_dispersion(i,1)
ENDDO
print*, '###'
DO i=1,n_pairs
print*, end_l_dispersion(i,2)
ENDDO
#endif


CALL DNS_END(0)

  STOP
END PROGRAM LAGRANGE_DISPERSION
