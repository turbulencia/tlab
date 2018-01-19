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
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif
    
  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
  TINTEGER  ierr,isize_wrk3d, i
  TREAL, DIMENSION(:,:),     ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:),       ALLOCATABLE       :: wrk1d,wrk2d, wrk3d
  TREAL, DIMENSION(:,:),     ALLOCATABLE       :: txc
  TREAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE :: s  

  TREAL, DIMENSION(:,:),     ALLOCATABLE       :: l_q, l_txc
  TREAL, DIMENSION(:),       ALLOCATABLE, SAVE :: l_comm
  TLONGINTEGER, DIMENSION(:,:),ALLOCATABLE     :: particle_bins
  TREAL, DIMENSION(:),       ALLOCATABLE       :: counter_interval
  TREAL y_pdf_max, y_pdf_min

  TINTEGER nitera_first,nitera_last,nitera_save

  CHARACTER*32 inifile
  CHARACTER*64 fname, str
  CHARACTER*128 line
  CHARACTER*32 bakfile

  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_part .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif
!  CALL DNS_READ_LOCAL(inifile) !for nitera stuff

! Get the local information from the dns.ini
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Start',      '0',  nitera_first)
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'End',        '0',  nitera_last )
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Restart',    '50', nitera_save )

  number_of_bins = INT(particle_pdf_max/particle_pdf_interval)
  inb_particle_txc = 1 

#include "dns_alloc_larrays.h"
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
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d))
  ALLOCATE(wrk3d(isize_wrk3d))

  ALLOCATE(s(imax,jmax,kmax, inb_scal_array))

  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !Allocte memory to read fields
     ALLOCATE(txc(isize_field,3))
  ENDIF
   WRITE(str,*) isize_l_comm; line = 'Allocating array l_comm of size '//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(l_comm(isize_l_comm), stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS. Not enough memory for l_comm.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
  
  ALLOCATE(particle_bins(number_of_bins,3)) 
  ALLOCATE(counter_interval(number_of_bins)) 

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

!#######################################################################
!Setup the file properties
!#######################################################################
  y_pdf_max=y_particle_pdf_pos+0.5*y_particle_pdf_width
  y_pdf_min=y_particle_pdf_pos-0.5*y_particle_pdf_width
  particle_bins=int(0,KIND=8)

  DO i = nitera_first, nitera_last, nitera_save

!#######################################################################
!READ ALL FILES
!#######################################################################
     WRITE(fname,*) i;  fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
     CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
     CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,1,1,inb_scal_array))

     WRITE(fname,*) i; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
     CALL IO_READ_PARTICLE(fname, l_g, l_q)

! ######################################################################
! Save particle pathlines for particle_pdf
! ######################################################################
     number_of_bins = INT(particle_pdf_max/particle_pdf_interval)

     WRITE(fname,*) i; fname = "particle_pdf."//TRIM(ADJUSTL(fname))
     CALL PARTICLE_PDF(fname,s, l_g,l_q,l_txc,l_comm, wrk2d,wrk3d)

  ENDDO

  CALL DNS_END(0)

  STOP
END PROGRAM
