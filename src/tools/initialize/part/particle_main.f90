#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "INIPART"

PROGRAM INIPART

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL

  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  ! -------------------------------------------------------------------
  TINTEGER  ierr

  TREAL, DIMENSION(:,:),      ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:),      ALLOCATABLE, SAVE :: q,s,txc
  TREAL, DIMENSION(:),        ALLOCATABLE, SAVE :: wrk1d,wrk2d, wrk3d

  TREAL,      DIMENSION(:,:), ALLOCATABLE, SAVE :: l_q, l_txc
  TREAL,      DIMENSION(:),   ALLOCATABLE, SAVE :: l_comm

  CHARACTER*64 str, line

  !########################################################################
  !########################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(ifile)

  IF ( icalc_part .EQ. 1 ) THEN
    CALL PARTICLE_READ_GLOBAL(ifile)
#ifdef USE_MPI
    CALL DNS_MPI_INITIALIZE
#endif

    ! -------------------------------------------------------------------
    ! Allocating memory space
    ! -------------------------------------------------------------------
    ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
    ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

    inb_flow_array = 0
    inb_scal_array = 0
    isize_wrk3d    = imax*jmax*kmax
    inb_txc        = inb_scal
#include "dns_alloc_arrays.h"
#include "dns_alloc_larrays.h"
    WRITE(str,*) isize_l_comm; line = 'Allocating array l_comm of size '//TRIM(ADJUSTL(str))
    CALL IO_WRITE_ASCII(lfile,line)
    ALLOCATE(l_comm(isize_l_comm), stat=ierr)
    IF ( ierr .NE. 0 ) THEN
      CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_comm.')
      CALL DNS_STOP(DNS_ERROR_ALLOC)
    ENDIF

    ! -------------------------------------------------------------------
    ! Read the grid
    ! -------------------------------------------------------------------
#include "dns_read_grid.h"

    ! -------------------------------------------------------------------
    ! Initialize particle information
    ! -------------------------------------------------------------------
    CALL PARTICLE_RANDOM_POSITION(l_g,l_q,l_txc,l_comm, txc, wrk3d)

    CALL IO_WRITE_PARTICLE(TRIM(ADJUSTL(tag_part))//'ics', l_g, l_q)

  ENDIF

  CALL DNS_STOP(0)
END PROGRAM INIPART
