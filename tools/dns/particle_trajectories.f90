#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

MODULE PARTICLE_TRAJECTORIES

  IMPLICIT NONE
  SAVE
  
  TREAL,      DIMENSION(:,:,:), ALLOCATABLE :: l_trajectories
  INTEGER(8), DIMENSION(:),     ALLOCATABLE :: l_trajectories_tags
  TINTEGER                                  :: save_point, save_time

CONTAINS
  
SUBROUTINE PARTICLE_TRAJECTORIES_XXX(nitera_last, nitera_save, nitera_first, l_q, l_tags, wrk3d)

  USE DNS_CONSTANTS, ONLY : efile, tag_traj
  USE DNS_GLOBAL,    ONLY : itime, isize_particle
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  
  TINTEGER nitera_last, nitera_save, nitera_first
  TREAL,      DIMENSION(isize_particle,3)                 :: l_q 
  INTEGER(8), DIMENSION(isize_particle)                   :: l_tags
  TREAL,      DIMENSION(inb_trajectory,isize_trajectory,nitera_save) :: wrk3d

! -------------------------------------------------------------------
  CHARACTER(len=32) name
  TINTEGER dummy_ims_npro
  TINTEGER dummy_isize_trajectory

  TINTEGER i,j, ierr

!#######################################################################
! Initialize; track only isize_trajectory particles 
!#######################################################################
  IF ( itime .EQ. nitera_first+1 ) THEN
     
     ALLOCATE(l_trajectories(inb_trajectory,isize_trajectory,nitera_save),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF

     ALLOCATE(l_trajectories_tags(isize_trajectory),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories_tags.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF

     save_point = nitera_first +nitera_save
     IF      ( itrajectory .EQ. LAG_TRAJECTORY_FIRST   ) THEN ! Track first isize_trajectory particles
        DO j=1,isize_trajectory              
           l_trajectories_tags(j) = INT(j, KIND=8)
        ENDDO
        
     ELSE IF ( itrajectory .EQ. LAG_TRAJECTORY_LARGEST ) THEN ! Read file with tags of largest particles, the ones to track
#ifdef USE_MPI
        IF (ims_pro .EQ. 0) THEN
#endif
           WRITE(name,*) nitera_last; name='largest_particle.'//TRIM(ADJUSTL(name))
           OPEN(unit=117, file=name, access='stream', form='unformatted')
           READ(117) dummy_ims_npro                             !is a integer
           READ(117, POS=SIZEOFINT+1) dummy_isize_trajectory  !is an integer
           READ(117, POS=SIZEOFINT*2+1) wrk3d                   !is real(8)
           READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_trajectory) l_trajectories_tags ! attention is integer(8)
           CLOSE(117)
#ifdef USE_MPI
        ENDIF
        CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
        CALL MPI_BCAST(l_trajectories_tags,isize_trajectory,MPI_INTEGER8,0,MPI_COMM_WORLD,ims_err)
#endif
     ENDIF
     l_trajectories = C_0_R
     save_time      = 0 
  ENDIF
  
!#######################################################################
! Accumulate trajectory information
!#######################################################################
  save_time = save_time +1

#ifdef USE_MPI
  DO i=1,ims_size_p(ims_pro+1)
#else
  DO i=1,particle_number
#endif        
     DO j=1,isize_trajectory
        IF ( l_tags(i) .EQ. l_trajectories_tags(j) ) THEN
           l_trajectories(1,j,save_time) = l_q(i,1)
           l_trajectories(2,j,save_time) = l_q(i,2)
           l_trajectories(3,j,save_time) = l_q(i,3)
        ENDIF
     ENDDO
  ENDDO

!#######################################################################
! Reduce and save
!#######################################################################
  IF ( itime .EQ. save_point ) THEN
#ifdef USE_MPI
     wrk3d = 0
     CALL MPI_REDUCE(l_trajectories, wrk3d, inb_trajectory*isize_trajectory*nitera_save, MPI_REAL8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
     CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
     IF(ims_pro .EQ. 0) THEN
#endif
        WRITE(name,*) itime; name =TRIM(ADJUSTL(tag_traj))//TRIM(ADJUSTL(name))
#define LOC_UNIT_ID 115
#define LOC_STATUS 'new'        
#include "dns_open_file.h"
        REWIND(LOC_UNIT_ID)
#ifdef USE_MPI
        WRITE(LOC_UNIT_ID) SNGL(wrk3d)
#else        
        WRITE(LOC_UNIT_ID) SNGL(l_trajectories)
#endif
        CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
     END IF
#endif
     l_trajectories = 0
     save_time      = 0
     save_point     = save_point +nitera_save     
  ENDIF

  RETURN
END SUBROUTINE PARTICLE_TRAJECTORIES_XXX

END MODULE PARTICLE_TRAJECTORIES
