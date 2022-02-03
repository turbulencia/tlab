#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# iheader    In      Header control flag:
!#                    0 No header
!#                    1 Scalar fields
!#                    2 Flow fields
!# nfield     In      Number of fields in the file
!# iread      In      Control flag, if set to 0 read all fields (array a
!#                    w/ space for nfield), if not read only iread (array
!#                    a only 1 field)
!#
!########################################################################
SUBROUTINE DNS_READ_FIELDS(fname, iheader, nx,ny,nz, nfield, iread, itxc, a, txc)

  USE TLAB_VARS, ONLY : imode_files
  USE TLAB_VARS, ONLY : itime, rtime, visc
  USE TLAB_CONSTANTS, ONLY : efile, lfile
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS,    ONLY : ims_npro_i, ims_npro_k, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname
  TINTEGER itxc, iheader, nfield, iread, nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*) :: a
  TREAL, DIMENSION(itxc)       :: txc

  ! -------------------------------------------------------------------
  CHARACTER*32 :: str
  CHARACTER*128 :: line
  TINTEGER nx_total, ny_total, nz_total
  TINTEGER ifield, iz

  TINTEGER isize_max, isize
  PARAMETER(isize_max=20)
  TREAL params(isize_max)

  ! ###################################################################
  IF ( iread .GT. nfield ) THEN
    CALL TLAB_WRITE_ASCII(efile, 'DNS_READ_FIELDS. Array size error')
    CALL TLAB_STOP(DNS_ERROR_SCALFORMAT)
  ENDIF

#ifdef USE_MPI
  nx_total = nx*ims_npro_i
  ny_total = ny
  nz_total = nz*ims_npro_k
#else
  nx_total = nx
  ny_total = ny
  nz_total = nz
#endif

  line = 'Reading field '//TRIM(ADJUSTL(fname))//' of size'
  WRITE(str,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(str))
  WRITE(str,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))
  WRITE(str,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))//'...'

  CALL TLAB_WRITE_ASCII(lfile, line)

  ! ###################################################################
  SELECT CASE( imode_files )

  CASE( DNS_NOFILE )        ! Do nothing
    a(:,1:nfield) = C_0_R

  CASE( DNS_FILE_NETCDF )   ! To be implemented

  CASE DEFAULT              ! One file with header per field
#ifdef USE_MPI
#ifdef USE_MPI_IO
    IF ( itxc .LT. nx*ny*nz ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'DNS_READ_FIELDS. Work array size error')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF
#endif
#endif

    ! -------------------------------------------------------------------
    ! read data
    DO ifield = 1,nfield
      IF ( iread .EQ. 0 .OR. iread .EQ. ifield ) THEN
        WRITE(str,'(I2)') ifield
        str=TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(str))
        isize = isize_max
        IF ( iread .EQ. 0 ) THEN; iz = ifield
        ELSE;                     iz = 1
        ENDIF
        CALL IO_READ_ONE_FIELD(str, iheader, nx,ny,nz,itime, isize,params, a(1,iz),txc)

      ENDIF
    ENDDO

    ! -------------------------------------------------------------------
    ! process header info
    IF ( isize .GT. isize_max ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'DNS_READ_FIELDS. Parameters array size error')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF

    IF      ( iheader .EQ. 1 ) THEN
      rtime = params(1)
    ELSE IF ( iheader .EQ. 2 ) THEN
      rtime = params(1)
      visc  = params(2)
    ENDIF
#ifdef USE_MPI
    CALL MPI_BCAST(itime, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
    CALL MPI_BCAST(rtime, 1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
    IF ( iheader .EQ. 2 ) THEN ! Flow type
      CALL MPI_BCAST(visc,  1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
    ENDIF
#endif

  END SELECT

  RETURN
END SUBROUTINE DNS_READ_FIELDS
