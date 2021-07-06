#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

SUBROUTINE PARTICLE_ALLOCATE(C_FILE_LOC)
  USE DNS_CONSTANTS, ONLY : lfile, efile
  USE LAGRANGE_GLOBAL
  USE LAGRANGE_ARRAYS
  USE DNS_GLOBAL, ONLY : isize_particle, inb_part_array, inb_part_txc

  IMPLICIT NONE

  CHARACTER(LEN=*) C_FILE_LOC

  ! -------------------------------------------------------------------
  CHARACTER(LEN=128) str, line
  TINTEGER ierr

  ! ###################################################################
  WRITE(str,*) isize_particle; line = 'Allocating array l_g.tags of size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(l_g%tags(isize_particle),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
    CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_tags.')
    CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  WRITE(str,*) isize_particle; line = 'Allocating array l_g.nodes of size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(l_g%nodes(isize_particle),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
    CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_g.')
    CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  WRITE(str,*) isize_particle; line = 'Allocating array l_q of size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) inb_part_array; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(l_q(isize_particle,inb_part_array),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
    CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_q.')
    CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

  IF ( inb_part_txc .GT. 0 ) THEN
    WRITE(str,*) isize_particle; line = 'Allocating array l_txc of size '//TRIM(ADJUSTL(str))//'x'
    WRITE(str,*) inb_part_txc; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    CALL IO_WRITE_ASCII(lfile,line)
    ALLOCATE(l_txc(isize_particle,inb_part_txc),stat=ierr)
    IF ( ierr .NE. 0 ) THEN
      CALL IO_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_txc.')
      CALL DNS_STOP(DNS_ERROR_ALLOC)
    ENDIF
  ENDIF

  RETURN
END SUBROUTINE PARTICLE_ALLOCATE
