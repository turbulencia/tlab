#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define SIZEOFBYTE 1

#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'

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
!# To be renamed io_read_fields
!########################################################################
SUBROUTINE DNS_READ_FIELDS(fname, iheader, nx,ny,nz, nfield, iread, itxc, a, txc)

  USE TLAB_VARS, ONLY : imode_files
  USE TLAB_VARS, ONLY : itime, rtime, visc
  USE TLAB_CONSTANTS, ONLY : efile, lfile
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k
#endif
  USE IO_FIELDS

  IMPLICIT NONE

  CHARACTER(LEN=*) fname
  TINTEGER,  INTENT(IN   ) :: itxc, iheader, nfield, iread, nx,ny,nz
  TREAL,     INTENT(  OUT) :: a(nx*ny*nz,*)
  TREAL,     INTENT(INOUT) :: txc(itxc)

  ! -------------------------------------------------------------------
  CHARACTER*32 :: name
  CHARACTER*128 :: line
  TINTEGER nx_total, ny_total, nz_total, header_offset
  TINTEGER ifield, iz

  TINTEGER isize_max, isize
  PARAMETER(isize_max=20)
  TREAL params(isize_max)

#ifdef USE_MPI
#include "mpif.h"
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp
  TINTEGER                :: subarray, ndims
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset
#endif

  ! ###################################################################
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
  WRITE(name,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(name))
  WRITE(name,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(name))
  WRITE(name,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(name))//'...'
  CALL TLAB_WRITE_ASCII(lfile, line)

  ! ###################################################################
  SELECT CASE( imode_files )

  CASE( DNS_NOFILE )        ! Do nothing
    a(:,1:nfield) = C_0_R

  CASE( DNS_FILE_NETCDF )   ! To be implemented

  CASE DEFAULT              ! One file with header per field
#ifdef USE_MPI
    ndims = 3
    sizes(1)   = nx*ims_npro_i; sizes(2)   = ny;           sizes(3)   = nz*ims_npro_k
    locsize(1) = nx;            locsize(2) = ny;           locsize(3) = nz
    offset(1)  = ims_offset_i;  offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

    CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
        MPI_ORDER_FORTRAN, MPI_REAL8, subarray, ims_err)
    CALL MPI_Type_commit(subarray, ims_err)

    IF ( itxc .LT. nx*ny*nz ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'DNS_READ_FIELDS. Work array size error')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF
#endif

    ! -------------------------------------------------------------------
    ! read data
    iz = 0
    DO ifield = 1,nfield
      IF ( iread .EQ. 0 .OR. iread .EQ. ifield ) THEN
        WRITE(name,'(I2)') ifield
        name=TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(name))
        isize = isize_max

        ! -------------------------------------------------------------------
        ! header
#ifdef USE_MPI
        IF ( ims_pro .EQ. 0 ) THEN
#endif

          header_offset = 0
#include "dns_open_file.h"
          REWIND(LOC_UNIT_ID)
          IF ( iheader .GT. 0 ) THEN
            CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,itime, params)
            isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
          ENDIF
          CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
        ENDIF
        CALL MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
        CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

        ! -------------------------------------------------------------------
        ! field
        iz = iz +1
        ! CALL IO_READ_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a(1,iz),txc)
#ifdef USE_MPI
        mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field
        mpio_locsize = nx*ny*nz
        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
        CALL MPI_File_read_all(mpio_fh, a(1,iz), mpio_locsize, MPI_REAL8, status, ims_err)
        CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
        READ(LOC_UNIT_ID,POS=header_offset+1) a(:,iz)
        CLOSE(LOC_UNIT_ID)

#endif

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

#undef LOC_UNIT_ID
#undef LOC_STATUS
