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

  USE DNS_CONSTANTS, ONLY : efile
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
  TREAL,      DIMENSION(3,isize_trajectories,nitera_save) :: wrk3d

! -------------------------------------------------------------------
  CHARACTER(len=32) name
  TINTEGER dummy_ims_npro
  TINTEGER dummy_isize_trajectories

  TINTEGER i,j, ierr

!#######################################################################
! Initialize; track only isize_trajectories particles 
!#######################################################################
  IF ( itime .EQ. nitera_first+1 ) THEN
     
     ALLOCATE(l_trajectories(3,isize_trajectories,nitera_save),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF

     ALLOCATE(l_trajectories_tags(isize_trajectories),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_trajectories_tags.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF

     save_point = nitera_first +nitera_save
     IF      ( itrajectories .EQ. LAG_TRAJECTORY_FIRST   ) THEN ! Track first isize_trajectories particles
        DO j=1,isize_trajectories              
           l_trajectories_tags(j) = j
        ENDDO
        
     ELSE IF ( itrajectories .EQ. LAG_TRAJECTORY_LARGEST ) THEN ! Read file with tags of largest particles, the ones to track
#ifdef USE_MPI
        IF (ims_pro .EQ. 0) THEN
#endif
           WRITE(name,*) nitera_last; name='largest_particle.'//TRIM(ADJUSTL(name))
           OPEN(unit=117, file=name, access='stream', form='unformatted')
           READ(117) dummy_ims_npro                             !is a integer
           READ(117, POS=SIZEOFINT+1) dummy_isize_trajectories  !is an integer
           READ(117, POS=SIZEOFINT*2+1) wrk3d                   !is real(8)
           READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_trajectories) l_trajectories_tags ! attention is integer(8)
           CLOSE(117)
#ifdef USE_MPI
        ENDIF
        CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
        CALL MPI_BCAST(l_trajectories_tags,isize_trajectories,MPI_INTEGER8,0,MPI_COMM_WORLD,ims_err)
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
  DO i=1,particle_vector(ims_pro+1)
#else
  DO i=1,particle_number
#endif        
     DO j=1,isize_trajectories
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
     CALL MPI_REDUCE(l_trajectories, wrk3d, 3*isize_trajectories*nitera_save, MPI_REAL8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
     CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
     IF(ims_pro .EQ. 0) THEN
#endif
!         DO i=1,save_time
!            WRITE(name,*) i +save_point -nitera_save; name ='trajectories.'//TRIM(ADJUSTL(name))//'.vtk'
!            OPEN(unit=115, file=name)
!            DO j=1,isize_trajectories
! #ifdef USE_MPI
!               WRITE (115,*) wrk3d(1,j,i),wrk3d(2,j,i),wrk3d(3,j,i)
! #else
!               WRITE (115,*) l_trajectories(1,j,i),l_trajectories(2,j,i),l_trajectories(3,j,i)
! #endif
!            ENDDO
!            CLOSE(115)
!         ENDDO
        WRITE(name,*) itime; name ='trajectories.'//TRIM(ADJUSTL(name))
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
