#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

PROGRAM TRANSFORMAT

  USE DNS_GLOBAL, ONLY : imax, jmax, kmax, itime
  USE DNS_GLOBAL, ONLY : icalc_flow, icalc_scal, inb_flow, inb_scal
  USE DNS_GLOBAL, ONLY : imode_files
  USE DNS_CONSTANTS
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: s, q
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk3d

  CHARACTER*32 fname
  TINTEGER isize_wrk3d, ierr, iopt, iheader

  TINTEGER itime_size_max, itime_size, i
  PARAMETER(itime_size_max=128)
  TINTEGER itime_vec(itime_size_max)

  TINTEGER iopt_size_max, iopt_size
  PARAMETER(iopt_size_max=20)
  TREAL opt_vec(iopt_size_max)

  CHARACTER*512 sRes

#ifdef USE_MPI
  INTEGER icount
#endif

! ###################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL('dns.ini')

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(q(imax*jmax*kmax,inb_flow),STAT=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'FORMAT_MAIN. Not enough memory, mark1.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
  
  ALLOCATE(s(imax*jmax*kmax,inb_scal),STAT=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'FORMAT_MAIN. Not enough memory, mark2.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
     
  isize_wrk3d = imax*jmax*kmax
  ALLOCATE(wrk3d(isize_wrk3d))
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'FORMAT_MAIN. Not enough memory, mark4.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CALL SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'Files', '-1', sRes)
     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'Integral Iterations ?'
        READ(*,'(A512)') sRes
     ENDIF
     itime_size = itime_size_max
     CALL LIST_INTEGER(sRes, itime_size, itime_vec)
#ifdef USE_MPI
  ENDIF
  CALL MPI_BCAST(itime_size, 1,      MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  icount = itime_size
  CALL MPI_BCAST(itime_vec,  icount, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Additional options
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CALL SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'ParamFFormat', '-1', sRes)
     iopt_size = iopt_size_max
     CALL LIST_REAL(sRes, iopt_size, opt_vec)

     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'FFormat tool. Option ?'
        WRITE(*,*) '1. From RawArray to RawSplit'
        WRITE(*,*) '2. From RawSplit to RawArray'
        WRITE(*,*) '3. Drop header'
        READ(*,*) iopt
     ELSE
        iopt = DINT(opt_vec(1))
     ENDIF

#ifdef USE_MPI
  ENDIF
  CALL MPI_BCAST(iopt, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif


! ###################################################################
! Processing files
! ###################################################################
  DO i = 1,itime_size
     itime = itime_vec(i)

     WRITE(sRes,*) itime; sRes = 'Processing iteration It'//TRIM(ADJUSTL(sRes))
     CALL IO_WRITE_ASCII(lfile, sRes)

! -------------------------------------------------------------------
! Flow fields
! -------------------------------------------------------------------
     IF ( icalc_flow .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        iheader = 2

! Read original format
        IF      ( iopt .EQ. 1 ) THEN; imode_files = DNS_FILE_RAWARRAY
        ELSE IF ( iopt .EQ. 2 ) THEN; imode_files = DNS_FILE_RAWSPLIT; ENDIF
        CALL DNS_READ_FIELDS(fname, iheader, imax,jmax,kmax, inb_flow,i0, isize_wrk3d, q, wrk3d)

! Write new format
        IF      ( iopt .EQ. 1 ) THEN; imode_files = DNS_FILE_RAWSPLIT
        ELSE IF ( iopt .EQ. 2 ) THEN; imode_files = DNS_FILE_RAWARRAY
        ELSE IF ( iopt .EQ. 3 ) THEN; iheader = 0; fname = TRIM(ADJUSTL(fname))//'.trn'; ENDIF
        CALL DNS_WRITE_FIELDS(fname, iheader, imax,jmax,kmax, inb_flow, isize_wrk3d, q, wrk3d)

     ENDIF

! -------------------------------------------------------------------
! Scalar fields
! -------------------------------------------------------------------
     IF ( icalc_scal .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        iheader = 1

! Read original format
        IF      ( iopt .EQ. 1 ) THEN; imode_files = DNS_FILE_RAWARRAY
        ELSE IF ( iopt .EQ. 2 ) THEN; imode_files = DNS_FILE_RAWSPLIT; ENDIF
        CALL DNS_READ_FIELDS(fname, iheader, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)

! Write new format
        IF      ( iopt .EQ. 1 ) THEN; imode_files = DNS_FILE_RAWSPLIT
        ELSE IF ( iopt .EQ. 2 ) THEN; imode_files = DNS_FILE_RAWARRAY
        ELSE IF ( iopt .EQ. 3 ) THEN; iheader = 0; fname = TRIM(ADJUSTL(fname))//'.trn'; ENDIF
        CALL DNS_WRITE_FIELDS(fname, iheader, imax,jmax,kmax, inb_scal, isize_wrk3d, s, wrk3d)

     ENDIF

  ENDDO

END PROGRAM TRANSFORMAT
