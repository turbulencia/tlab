#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif


!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# write particle information for trajectories in a .vtk file
!# 
!########################################################################
!# ARGUMENTS 
!#
!#
!########################################################################
SUBROUTINE DNS_WRITE_TRAJECTORIES(fname, l_q, l_tags, l_trajectories, l_trajectories_tags, wrk3d,txc, itime,nitera_last, nitera_save, nitera_first)

  USE DNS_GLOBAL, ONLY :  isize_field 
  USE DNS_GLOBAL, ONLY: isize_particle
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif


  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  
  CHARACTER(*)  ::  fname
  CHARACTER(len=30)  ::  str

  TREAL, DIMENSION(isize_particle,3) :: l_q 
  INTEGER(8), DIMENSION(isize_particle) :: l_tags
  TREAL, DIMENSION(isize_field) ::wrk3d, txc
  TREAL, DIMENSION(3,num_trajectories,nitera_save) :: l_trajectories
  INTEGER(8), DIMENSION(num_trajectories) :: l_trajectories_tags
  
  TREAL, DIMENSION(num_trajectories) :: dummy_big_overall ! maybe better pointer to wrk3d 
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: big_overall
  TINTEGER  dummy_ims_npro
  TINTEGER  dummy_num_trajectories
  TINTEGER  nitera_last, nitera_save, nitera_first
  
  TINTEGER i,j, itime, version
  TINTEGER,SAVE :: save_point, save_time
  TLONGINTEGER m



  
!  version=1 !0=old 1=new

  !#######################################################################
  !WRITE THE N LARGEST PARTICLES
  !#######################################################################
!  IF (version .EQ. 1) THEN
  
    !#######################################################################
    !Read the file with largest particles
    !#######################################################################
#ifdef USE_MPI
    IF (itime .EQ. nitera_first+1) THEN !READ FILE ONLY FIRST STEP
      save_point=nitera_first+nitera_save
      IF (ims_pro .EQ. 0) THEN
        fname = 'largest_particle'
        WRITE(str,*) nitera_last;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
        OPEN(unit=117, file=str, access='stream', form='unformatted')
        READ(117) dummy_ims_npro   !is a integer
        READ(117, POS=SIZEOFINT+1) dummy_num_trajectories  !is an integer
        READ(117, POS=SIZEOFINT*2+1) dummy_big_overall  !is real(8)
        READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*num_trajectories) l_trajectories_tags ! attention is integer(8)
        CLOSE(117)
      ENDIF
  
      !#######################################################################
      !Broadcast the ID of the largest particles
      !#######################################################################
      CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
      CALL MPI_BCAST(l_trajectories_tags,num_trajectories,MPI_INTEGER8,0,MPI_COMM_WORLD,ims_err)
  
      l_trajectories(:,:,:) = C_0_R
      save_time=1 
    ENDIF 

!      save_time=itime
    !#######################################################################
    !Search for the largest particles in every processor
    !#######################################################################
    DO i=1,particle_vector(ims_pro+1)
      DO j=1,num_trajectories
        IF (l_tags(i) .EQ. l_trajectories_tags(j)) THEN
          l_trajectories(1,j,save_time)=l_q(i,1)
          l_trajectories(2,j,save_time)=l_q(i,2)
          l_trajectories(3,j,save_time)=l_q(i,3)
        ENDIF
      ENDDO
    ENDDO 
    save_time=save_time+1
!DO i=1,16  
!IF (ims_pro .EQ. 6 .OR. ims_pro .EQ. 7) THEN
!print*, 'time step', save_time-1, particle_number
!DO j=1,num_trajectoris
!j=12
!print*, j, i-1, l_trajectories(1,j,save_time-1)
!ENDDO
!ENDIF
!ENDDO
    !#######################################################################
    !Every writing point reduce information to root - only at nitera_save
    !#######################################################################
    IF (itime .EQ. save_point ) THEN
      IF (ims_pro .EQ. 0) THEN
        ALLOCATE(big_overall(3,num_trajectories,nitera_save)) ! DO WRK3D later
        big_overall(:,:,:)=C_0_R
      ENDIF
      CALL MPI_REDUCE(l_trajectories, big_overall, 3*num_trajectories*nitera_save, MPI_REAL8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
    
      !#######################################################################
      !Write particle file with processor root
      !#######################################################################
      IF(ims_pro .EQ. 0) THEN
        DO i=1,save_time-1 !because save_time is last step+1 
          WRITE(fname,*) i+save_point-nitera_save; fname ='trajectories.'//TRIM(ADJUSTL(fname))
          WRITE(str,*) "vtk";  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! adjust the name with the number for direction
          OPEN(unit=115, file=str)
          DO j=1,num_trajectories
            WRITE (115,*) big_overall(1,j,i),big_overall(2,j,i),big_overall(3,j,i)
          ENDDO
          CLOSE(115)
        ENDDO
      END IF
  
      IF (ims_pro .EQ. 0) THEN
        DEALLOCATE(big_overall)
      ENDIF
      l_trajectories(:,:,:) = 0
      save_time=1 
      save_point=save_point + nitera_save

    ENDIF !itime savepoint

#else
    !#######################################################################
    !Read first the file with largest particles
    !#######################################################################
    IF (itime .EQ. nitera_first+1) THEN !READ FILE ONLY FIRST STEP
      save_point=nitera_first+nitera_save
      fname = 'largest_particle'
      WRITE(str,*) nitera_last;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
      OPEN(unit=117, file=str, access='stream', form='unformatted')
      READ(117) dummy_ims_npro   !is a integer
      READ(117, POS=SIZEOFINT+1) dummy_num_trajectories  !is an integer
      READ(117, POS=SIZEOFINT*2+1) dummy_big_overall  !is real(8)
      READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*num_trajectories) l_trajectories_tags ! attention is integer(8)
      CLOSE(117)
      l_trajectories = 0
      save_time=1 
    ENDIF 
    !#######################################################################
    !Search for the largest particles in every processor
    !#######################################################################
    DO i=1,particle_number
      DO j=1,num_trajectories
        IF (l_tags(i) .EQ. l_trajectories_tags(j)) THEN
          l_trajectories(1,j,save_time)=l_q(i,1)
          l_trajectories(2,j,save_time)=l_q(i,2)
          l_trajectories(3,j,save_time)=l_q(i,3)
        ENDIF
      ENDDO
    ENDDO 
    save_time=save_time+1
    !#######################################################################
    !Write particle file with processor root
    !#######################################################################
    IF (itime .EQ. save_point ) THEN
      DO i=1,save_time-1 !because save_time is last step+1 
        WRITE(fname,*) i+save_point-nitera_save; fname = 'trajectories.'//TRIM(ADJUSTL(fname))
        WRITE(str,*) "vtk";  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! adjust the name with the number for direction
        OPEN(unit=115, file=str)
        DO j=1,num_trajectories
          WRITE (115,*) l_trajectories(1,j,i),l_trajectories(2,j,i),l_trajectories(3,j,i)
        ENDDO
        CLOSE(115)
      ENDDO

      l_trajectories(:,:,:) = 0
      save_time=1 
      save_point=save_point + nitera_save

    ENDIF ! itime save point


#endif
  !#######################################################################
  !OLD VERSION - WRITE ALL TRAJECTORIES EVERY TIMESTEP
  !#######################################################################
!  ELSEIF (version .EQ. 0) THEN
!
!#ifdef USE_MPI
!    !Sort particles
!  
!    DO j=1,3 
!      DO i=1,particle_vector(ims_pro+1)
!        wrk3d(l_tags(i),j)= l_q(i,j)
!      END DO
!    END DO
!  
!  
!  
!  
!    CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err) 
!  
!    CALL MPI_REDUCE(wrk3d, txc, particle_number*3, MPI_REAL8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
!  
!   
!    IF(ims_pro .EQ. 0) THEN
!  
!  
!        WRITE(str,*) "vtk";  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! adjust the name with the number for direction
!        OPEN(unit=115, file=str)
!  !      WRITE (115,*)  '# vtk DataFile Version 3.0'
!  !      WRITE (115,*)  'particle'
!  !      WRITE (115,*)  'ASCII'
!  !      WRITE (115,*)  'DATASET UNSTRUCTURED_GRID'
!  !      WRITE (115,*)  '' !maybee delete this line in the header, need to be tested
!  !      WRITE (115,*)  'POINTS 50 float'
!        DO m=1,particle_number
!          WRITE (115,*) txc(m,1), txc(m,2), txc(m,3)
!        END DO
!        CLOSE(115)
!    END IF
!#else
!    DO j=1,3 
!      DO i=1,particle_number
!        wrk3d(l_tags(i),j)= l_q(i,j)
!      END DO
!    END DO
!    WRITE(str,*) "vtk";  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! adjust the name with the number for direction
!    OPEN(unit=115, file=str)
!    DO m=1,particle_number
!      WRITE (115,*) wrk3d(m,1), wrk3d(m,2), wrk3d(m,3)
!    END DO
!    CLOSE(115)
!#endif 
!  ENDIF 
  RETURN
END SUBROUTINE DNS_WRITE_TRAJECTORIES

