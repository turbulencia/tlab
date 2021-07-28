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

SUBROUTINE IO_READ_PARTICLE(fname, l_g, l_q)

  USE TLAB_CONSTANTS,   ONLY : lfile, efile
  USE TLAB_VARS,      ONLY : isize_particle, inb_part_array
  USE TLAB_VARS,      ONLY : g
  USE TLAB_PROCS
  USE LAGRANGE_GLOBAL, ONLY : particle_dt, particle_number_total
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_size_p, ims_pro, ims_npro, ims_err
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*)     fname
  TYPE(particle_dt) l_g
  TREAL, DIMENSION(isize_particle,inb_part_array) :: l_q !, OPTIONAL :: l_q

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

  CALL TLAB_WRITE_ASCII(lfile, 'Reading field '//TRIM(ADJUSTL(fname))//'...')

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
     CALL TLAB_WRITE_ASCII(efile, 'IO_PARTICLE. Number-of-processors mismatch.')
     CALL TLAB_STOP(DNS_ERROR_PARTICLE)
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
     CALL TLAB_WRITE_ASCII(efile, 'IO_PARTICLE. Number-of-particles mismatch.')
     CALL TLAB_STOP(DNS_ERROR_PARTICLE)
  ENDIF

! Number of particles in local processor
  l_g%np = ims_size_p(ims_pro+1)

! -------------------------------------------------------------------
! Use MPI-IO to read particle tags in each processor
! -------------------------------------------------------------------
  name = TRIM(ADJUSTL(fname))//".id"
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
  CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_FILE_READ_ALL(mpio_fh, l_g%tags, l_g%np, MPI_INTEGER8, status, ims_err)
  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_part_array
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
        CALL MPI_FILE_READ_ALL(mpio_fh, l_q(1,i), l_g%np, MPI_REAL8, status, ims_err)
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
     CALL TLAB_WRITE_ASCII(efile, 'IO_PARTICLE. Number-of-particles mismatch.')
     CLOSE(LOC_UNIT_ID)
     CALL TLAB_STOP(DNS_ERROR_PARTICLE)
  ENDIF
  READ(LOC_UNIT_ID) l_g%tags
  CLOSE(LOC_UNIT_ID)

! For homogeneity with MPI version
! If we need more than 4 bytes, we should be using MPI...
  l_g%np = INT(particle_number_total)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_part_array
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
#include "dns_open_file.h"
        READ(LOC_UNIT_ID) idummy             ! dummy, should be 1 in serial
        READ(LOC_UNIT_ID) particle_number_loc
        READ(LOC_UNIT_ID) l_q(:,i)
        CLOSE(LOC_UNIT_ID)
     ENDDO
!  ENDIF

#endif

  CALL PARTICLE_LOCATE_Y( l_g%np, l_q(1,2), l_g%nodes, g(2)%size, g(2)%nodes )

  RETURN
END SUBROUTINE IO_READ_PARTICLE

#undef LOC_UNIT_ID
#undef LOC_STATUS

!#######################################################################
!#######################################################################
#define LOC_UNIT_ID 118
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_PARTICLE(fname, l_g, l_q)

  USE TLAB_CONSTANTS,   ONLY : lfile
  USE TLAB_VARS,      ONLY : isize_particle, inb_part_array
  USE TLAB_PROCS
  USE LAGRANGE_GLOBAL, ONLY : particle_dt
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_size_p, ims_pro, ims_npro, ims_err
#endif

  IMPLICIT NONE
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*)     fname
  TYPE(particle_dt) l_g
  TREAL, DIMENSION(isize_particle,inb_part_array) :: l_q !, OPTIONAL :: l_q

! -------------------------------------------------------------------
  TINTEGER i
  CHARACTER(len=32) name
#ifdef USE_MPI
  TINTEGER mpio_fh
  INTEGER (KIND=8)  mpio_disp, count
  TINTEGER status(MPI_STATUS_SIZE)
#else
  TINTEGER idummy
#endif

  CALL TLAB_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(fname))//'...')

#ifdef USE_MPI
!#######################################################################
! Parallel case
!#######################################################################
! -------------------------------------------------------------------
! Let Process 0 handle header
! -------------------------------------------------------------------
  CALL MPI_ALLGATHER(l_g%np, 1, MPI_INTEGER4, ims_size_p, 1, MPI_INTEGER4, MPI_COMM_WORLD, ims_err)

  IF ( ims_pro .EQ. 0 ) THEN
     name = TRIM(ADJUSTL(fname))//".id"
#include "dns_open_file.h"
     WRITE(LOC_UNIT_ID) ims_npro
     WRITE(LOC_UNIT_ID) ims_size_p
     CLOSE(LOC_UNIT_ID)

!     IF ( PRESENT(l_q) ) THEN
        DO i = 1,inb_part_array
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
  CALL MPI_FILE_WRITE_ALL(mpio_fh, l_g%tags, l_g%np, MPI_INTEGER8, status, ims_err)
  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_part_array
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
        CALL MPI_FILE_WRITE_ALL(mpio_fh, l_q(1,i), l_g%np, MPI_REAL8, status, ims_err)
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
  WRITE(LOC_UNIT_ID) l_g%np
  WRITE(LOC_UNIT_ID) l_g%tags
  CLOSE(LOC_UNIT_ID)

!  IF ( PRESENT(l_q) ) THEN
     DO i = 1,inb_part_array
        WRITE(name,*) i; name = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(name))
#include "dns_open_file.h"
        WRITE(LOC_UNIT_ID) idummy
        WRITE(LOC_UNIT_ID) l_g%np
        WRITE(LOC_UNIT_ID) l_q(:,i)
        CLOSE(LOC_UNIT_ID)
     ENDDO
!  ENDIF

#endif

  RETURN
END SUBROUTINE IO_WRITE_PARTICLE

!#######################################################################
!#######################################################################
SUBROUTINE PARTICLE_LOCATE_Y( pmax, y_part, j_part, jmax, y_grid )

  IMPLICIT NONE

  TINTEGER, INTENT(IN) :: pmax, jmax
  TREAL,    DIMENSION(pmax), INTENT(IN) :: y_part
  TINTEGER, DIMENSION(pmax), INTENT(OUT) :: j_part
  TREAL,    DIMENSION(jmax), INTENT(IN) :: y_grid

  TINTEGER ip, jm, jp, jc

  DO ip = 1,pmax
     jp = jmax
     jm = 1
     jc = ( jm +jp ) /2
     DO WHILE ( (y_part(ip)-y_grid(jc))*(y_part(ip)-y_grid(jc+1)) .GT. C_0_R .AND. jc .GT. jm )
        IF ( y_part(ip) .LT. y_grid(jc) ) THEN; jp = jc;
        ELSE;                                   jm = jc; END IF
        jc = ( jm +jp ) /2
     END DO
     j_part(ip) = jc
!     WRITE(*,'(i,3f)') ip, y_grid(jc), y_part(ip), y_grid(jc+1)
  END DO

  RETURN
END SUBROUTINE PARTICLE_LOCATE_Y
