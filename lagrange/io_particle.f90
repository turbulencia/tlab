#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!#######################################################################
!#######################################################################
#define LOC_UNIT_ID 117
#define LOC_STATUS 'old'

SUBROUTINE IO_READ_PARTICLE(fname, l_tags, l_q)

  USE DNS_CONSTANTS,   ONLY : lfile, efile
  USE DNS_GLOBAL,      ONLY : isize_particle, inb_particle
  USE DNS_GLOBAL,      ONLY : g
  USE LAGRANGE_GLOBAL, ONLY : particle_number_total, particle_number_local, l_g
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_size_p, ims_pro, ims_npro, ims_err
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif  
 
  CHARACTER*(*)                                                :: fname
  INTEGER(8), DIMENSION(isize_particle)                        :: l_tags
  TREAL,      DIMENSION(isize_particle,inb_particle) :: l_q !, OPTIONAL :: l_q 
  
! -------------------------------------------------------------------
  TINTEGER i
  CHARACTER(len=32) name
#ifdef USE_MPI
  TINTEGER ims_npro_loc
  TINTEGER mpio_fh
  INTEGER (KIND=8) mpio_disp, count
  TINTEGER status(MPI_STATUS_SIZE)
#else
  TINTEGER particle_number_loc
  TINTEGER idummy
#endif

  CALL IO_WRITE_ASCII(lfile, 'Reading field '//TRIM(ADJUSTL(fname))//'...')

#ifdef USE_MPI
!#######################################################################
! Parallel case
!#######################################################################
! -------------------------------------------------------------------
! Let Process 0 handle header
! -------------------------------------------------------------------  
  IF ( ims_pro .EQ. 0 ) THEN
     name = TRIM(ADJUSTL(fname))//".id"
#include "dns_open_file.h"
     READ(LOC_UNIT_ID) ims_npro_loc
     READ(LOC_UNIT_ID) ims_size_p(1:ims_npro_loc)
     CLOSE(LOC_UNIT_ID)
  END IF

! Check
  CALL MPI_BCAST(ims_npro_loc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ims_err)
  IF ( ims_npro .NE. ims_npro_loc) THEN
     CALL IO_WRITE_ASCII(efile, 'IO_PARTICLE. Number-of-processors mismatch.')
     CALL DNS_STOP(DNS_ERROR_PARTICLE)
  ENDIF

! Broadcast number of particles per processor
  CALL MPI_BCAST(ims_size_p,ims_npro,MPI_INTEGER,0,MPI_COMM_WORLD,ims_err)

! Displacement per processor
  mpio_disp = INT( (ims_npro+1)*SIZEOFINT, KIND=8 )

  count = 0
  DO i = 1,ims_pro
     count = count +INT(ims_size_p(i),KIND=8)
  ENDDO
  mpio_disp = mpio_disp +count *INT(SIZEOFLONGINT,KIND=8)

! Check
  DO i = ims_pro+1,ims_npro
     count = count +INT(ims_size_p(i),KIND=8)
  ENDDO
  IF ( particle_number_total .NE. count ) THEN
     CALL IO_WRITE_ASCII(efile, 'IO_PARTICLE. Number-of-particles mismatch.')
     CALL DNS_STOP(DNS_ERROR_PARTICLE)
  ENDIF

! Number of particles in local processor
  particle_number_local = ims_size_p(ims_pro+1)

! -------------------------------------------------------------------
! Use MPI-IO to read particle tags in each processor
! -------------------------------------------------------------------
  name = TRIM(ADJUSTL(fname))//".id"
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
  CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_FILE_READ_ALL(mpio_fh, l_tags, particle_number_local, MPI_INTEGER8, status, ims_err)
  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_particle
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
        CALL MPI_FILE_READ_ALL(mpio_fh, l_q(1,i), particle_number_local, MPI_REAL8, status, ims_err)
        CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
     ENDDO
!  ENDIF
  
#else
! #######################################################################
! Serial case
! #######################################################################
  name = TRIM(ADJUSTL(fname))//".id"
#include "dns_open_file.h"
  READ(LOC_UNIT_ID) idummy                   ! dummy, should be 1 in serial
  READ(LOC_UNIT_ID) particle_number_loc
! Check
  IF ( particle_number_total .NE. INT(particle_number_loc,KIND=8) ) THEN
     CALL IO_WRITE_ASCII(efile, 'IO_PARTICLE. Number-of-particles mismatch.')
     CLOSE(LOC_UNIT_ID)
     CALL DNS_STOP(DNS_ERROR_PARTICLE)
  ENDIF
  READ(LOC_UNIT_ID) l_tags
  CLOSE(LOC_UNIT_ID)

! For homogeneity with MPI version
! If we need more than 4 bytes, we should be using MPI...
  particle_number_local = INT(particle_number_total)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_particle
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
#include "dns_open_file.h"
        READ(LOC_UNIT_ID) idummy             ! dummy, should be 1 in serial
        READ(LOC_UNIT_ID) particle_number_loc
        READ(LOC_UNIT_ID) l_q(:,i)
        CLOSE(LOC_UNIT_ID)
     ENDDO
!  ENDIF
  
#endif

!     CALL PARTICLE_LOCATE_Y( particle_number_local, l_q(:,2), l_g, g(2)%size, g(2)%nodes(:) )
     
  RETURN
END SUBROUTINE IO_READ_PARTICLE

#undef LOC_UNIT_ID
#undef LOC_STATUS

!#######################################################################
!#######################################################################
#define LOC_UNIT_ID 118
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_PARTICLE(fname, l_tags, l_q)

  USE DNS_CONSTANTS,   ONLY : lfile
  USE DNS_GLOBAL,      ONLY : isize_particle, inb_particle
  USE LAGRANGE_GLOBAL, ONLY : particle_number_local
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_size_p, ims_pro, ims_npro, ims_err
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*)                                                :: fname
  INTEGER(8), DIMENSION(isize_particle)                        :: l_tags
  TREAL,      DIMENSION(isize_particle,inb_particle) :: l_q !, OPTIONAL :: l_q 

! -------------------------------------------------------------------  
  TINTEGER i, idummy
  CHARACTER(len=32) name
#ifdef USE_MPI
  TINTEGER mpio_fh
  INTEGER (KIND=8)  mpio_disp, count
  TINTEGER status(MPI_STATUS_SIZE)
#endif

  CALL IO_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(fname))//'...')

#ifdef USE_MPI
!#######################################################################
! Parallel case
!#######################################################################
! -------------------------------------------------------------------
! Let Process 0 handle header
! -------------------------------------------------------------------  
  CALL MPI_ALLGATHER(particle_number_local, 1, MPI_INTEGER4, ims_size_p, 1, MPI_INTEGER4, MPI_COMM_WORLD, ims_err)

  IF ( ims_pro .EQ. 0 ) THEN
     name = TRIM(ADJUSTL(fname))//".id"     
#include "dns_open_file.h"
     WRITE(LOC_UNIT_ID) ims_npro
     WRITE(LOC_UNIT_ID) ims_size_p
     CLOSE(LOC_UNIT_ID)

!     IF ( PRESENT(l_q) ) THEN
        DO i = 1,inb_particle
           WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
#include "dns_open_file.h"
           WRITE(LOC_UNIT_ID) ims_npro
           WRITE(LOC_UNIT_ID) ims_size_p
           CLOSE(LOC_UNIT_ID)
        ENDDO
!     ENDIF
  END IF

! Displacement per processor
  mpio_disp = INT( (ims_npro+1)*SIZEOFINT, KIND=8 )

  count = 0
  DO i = 1,ims_pro
     count = count +INT(ims_size_p(i),KIND=8)
  ENDDO
  mpio_disp = mpio_disp +count *INT(SIZEOFLONGINT,KIND=8)

! -------------------------------------------------------------------
! Use MPI-IO to write particle tags in each processor
! -------------------------------------------------------------------
  name = TRIM(ADJUSTL(fname))//".id"
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
  CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ims_err)    
  CALL MPI_FILE_WRITE_ALL(mpio_fh, l_tags, particle_number_local, MPI_INTEGER8, status, ims_err)
  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_particle
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
        CALL MPI_FILE_WRITE_ALL(mpio_fh, l_q(1,i), particle_number_local, MPI_REAL8, status, ims_err)
        CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
     ENDDO
!  ENDIF
  
#else 
! #######################################################################
! Serial case
! #######################################################################
  idummy = 1
  name = TRIM(ADJUSTL(fname))//".id"
#include "dns_open_file.h"
  WRITE(LOC_UNIT_ID) idummy  
  WRITE(LOC_UNIT_ID) particle_number_local
  WRITE(LOC_UNIT_ID) l_tags
  CLOSE(LOC_UNIT_ID)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_particle
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
#include "dns_open_file.h"
        WRITE(LOC_UNIT_ID) idummy  
        WRITE(LOC_UNIT_ID) particle_number_local
        WRITE(LOC_UNIT_ID) l_q(:,i)
        CLOSE(LOC_UNIT_ID)
     ENDDO
!  ENDIF

#endif 

  RETURN
END SUBROUTINE IO_WRITE_PARTICLE

