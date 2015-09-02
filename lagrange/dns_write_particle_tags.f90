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
!# write particle_id information
!# 
!########################################################################
!# ARGUMENTS 
!#
!#
!########################################################################
SUBROUTINE DNS_WRITE_PARTICLE_TAGS(fname, l_tags)

  USE DNS_GLOBAL, ONLY: isize_particle
  USE LAGRANGE_GLOBAL, ONLY :  particle_number
#ifdef USE_MPI
  USE DNS_MPI, ONLY : particle_vector, ims_pro, ims_npro, ims_err
#endif


  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  
  CHARACTER(*)  ::  fname
#ifdef USE_MPI
  CHARACTER(len=20)  ::  str
  TINTEGER mpio_fh
  INTEGER (KIND=8)  mpio_disp
  TINTEGER status(MPI_STATUS_SIZE)
  TLONGINTEGER, DIMENSION(ims_npro) :: particle_vector_final
  TINTEGER, DIMENSION(:),   ALLOCATABLE :: particle_vector_buffer
#endif



  INTEGER(8), DIMENSION(isize_particle) :: l_tags
  TINTEGER i,  particle_pos
#ifdef USE_MPI

  ALLOCATE(particle_vector_buffer(ims_npro))


!#######################################################################
!GET particle number for each processor
!#######################################################################
  CALL MPI_ALLGATHER(particle_vector(ims_pro+1),1,MPI_INTEGER4,particle_vector_buffer,1,MPI_INTEGER4,MPI_COMM_WORLD,ims_err)
  particle_vector(:)=particle_vector_buffer(:)

!#######################################################################
! HEADER   
!#######################################################################
  IF(ims_pro .EQ. 0) THEN
      OPEN(unit=117, file=fname, access='stream', form='unformatted')
      INQUIRE(UNIT=117, POS=particle_pos)
      WRITE (117)  ims_npro  !header
      INQUIRE(UNIT=117, POS=particle_pos)
      WRITE (117)  particle_number  !header
      INQUIRE(UNIT=117, POS=particle_pos)
      WRITE (117)  particle_vector
      INQUIRE(UNIT=117, POS=particle_pos)
      CLOSE(117)
  END IF

!#######################################################################
! MPI I/O
!#######################################################################
 
  !Calculating positions for processors to write
  !First position is already displacemnet for second processor (ims_pro=1)
  particle_vector_final(1)=INT(particle_vector(1),KIND=8)
  DO i=2,ims_npro
    particle_vector_final(i)=INT(particle_vector(i),KIND=8) + particle_vector_final(i-1)
  END DO

    
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY, MPI_INFO_NULL,mpio_fh,ims_err)
    
    !DISPLACEMENT

    !Header displacement
    mpio_disp = INT(((ims_npro+1)*SIZEOFINT)+SIZEOFLONGINT+1, KIND=8)
    
    !Displacement per processor
    IF(ims_pro .EQ. 0) THEN
      mpio_disp=mpio_disp 
    ELSE
      mpio_disp=mpio_disp + (particle_vector_final(ims_pro))*INT(SIZEOFLONGINT, KIND=8) !due to start point (ims_pro+1 would be normal)
    END IF   

    !Set the view per each processor using displacement
    CALL MPI_FILE_SET_VIEW(mpio_fh,mpio_disp,MPI_INTEGER8,MPI_INTEGER8,'native',MPI_INFO_NULL, ims_err)    

    !WRITE PARTICLE
    CALL MPI_FILE_WRITE_ALL(mpio_fh,l_tags(1:particle_vector(ims_pro+1)),particle_vector(ims_pro+1),MPI_INTEGER8,status,ims_err)

    CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
  DEALLOCATE(particle_vector_buffer)
#else 

!#######################################################################
! SERIAL WRITE
!#######################################################################
    OPEN(unit=117, file=fname, access='stream', form='unformatted')
     
    INQUIRE(UNIT=117, POS=particle_pos)
    WRITE (117)  particle_number  !header
    INQUIRE(UNIT=117, POS=particle_pos)
    DO i=1,particle_number
    WRITE (117) l_tags(i)
    INQUIRE(UNIT=117, POS=particle_pos)
    END DO
    CLOSE(117)


#endif 
 
  RETURN
END SUBROUTINE DNS_WRITE_PARTICLE_TAGS

