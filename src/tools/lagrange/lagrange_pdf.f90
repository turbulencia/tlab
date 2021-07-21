#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "LAGRANGE_PDF"

!########################################################################
!########################################################################
PROGRAM LAGRANGE_PDF

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE TLAB_ARRAYS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_PROCS
#endif
  USE LAGRANGE_GLOBAL
  USE LAGRANGE_ARRAYS

  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
! Additional local arrays
  TREAL, DIMENSION(:),       ALLOCATABLE, SAVE :: l_comm

  TINTEGER nitera_first, nitera_last, nitera_save
  TINTEGER ierr, i

  CHARACTER*64 fname, str
  CHARACTER*128 line
  CHARACTER*32 bakfile

  bakfile = TRIM(ADJUSTL(ifile))//'.bak'

  CALL TLAB_START

  CALL DNS_READ_GLOBAL(ifile)
  IF ( icalc_part .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif
!  CALL DNS_READ_LOCAL(ifile) !for nitera stuff

! Get the local information from the dns.ini
  CALL SCANINIINT(bakfile, ifile, 'Iteration', 'Start',      '0',  nitera_first)
  CALL SCANINIINT(bakfile, ifile, 'Iteration', 'End',        '0',  nitera_last )
  CALL SCANINIINT(bakfile, ifile, 'Iteration', 'Restart',    '50', nitera_save )

  inb_part_txc = 1

  CALL PARTICLE_ALLOCATE(C_FILE_LOC)

  isize_wrk3d = imax*jmax*kmax
  isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
  isize_wrk3d = MAX(isize_wrk3d,(jmax*(kmax+1)*inb_particle_interp*2))
  isize_wrk3d = MAX(isize_wrk3d,(jmax*(imax+1)*inb_particle_interp*2))

  isize_wrk2d = MAX(isize_wrk2d, jmax*inb_particle_interp)

  ! IF (jmax_part .EQ. 1) THEN
  !    jmax_part   = jmax ! 1 by default
  ! ENDIF
! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d,1))
  ALLOCATE(wrk3d(isize_wrk3d))

  ALLOCATE(s(isize_field, inb_scal_array))

  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !Allocte memory to read fields
     ALLOCATE(txc(isize_field,3))
  ENDIF

   WRITE(str,*) isize_l_comm; line = 'Allocating array l_comm of size '//TRIM(ADJUSTL(str))
  CALL TLAB_WRITE_ASCII(lfile,line)
  ALLOCATE(l_comm(isize_l_comm), stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL TLAB_WRITE_ASCII(efile,'DNS. Not enough memory for l_comm.')
     CALL TLAB_STOP(DNS_ERROR_ALLOC)
  ENDIF

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
CALL FDM_INITIALIZE(x, g(1), wrk1d)
CALL FDM_INITIALIZE(y, g(2), wrk1d)
CALL FDM_INITIALIZE(z, g(3), wrk1d)

!#######################################################################
!Setup the file properties
!#######################################################################
  DO i = nitera_first, nitera_last, nitera_save

!#######################################################################
!READ ALL FILES
!#######################################################################
     WRITE(fname,*) i;  fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
     CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
     CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))

     WRITE(fname,*) i; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
     CALL IO_READ_PARTICLE(fname, l_g, l_q)

! ######################################################################
! Save particle pathlines for particle_pdf
! ######################################################################
     WRITE(fname,*) i; fname = "particle_pdf."//TRIM(ADJUSTL(fname))
     CALL PARTICLE_PDF(fname,s, l_g,l_q,l_txc,l_comm, wrk3d)

  ENDDO

  CALL TLAB_STOP(0)
END PROGRAM
