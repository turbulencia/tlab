#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define SIZEOFBYTE 1

!########################################################################
!#
!# Read & write a file of size (nx*ims_npro_i)x(ny*ims_npro_y)x(nz*ims_npro_z)
!#
!# Unformatted records
!# No embedded record information
!#
!########################################################################

!########################################################################
!########################################################################
#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'

SUBROUTINE IO_READ_SUBARRAY_FIELD(name, iheader, nx,ny,nz,nt, isize,params, a)

#ifdef USE_MPI
  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_PROCS
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER(LEN=*) name
  TINTEGER, INTENT(IN   ) :: iheader, nx,ny,nz,nt
  TINTEGER, INTENT(INOUT) :: isize
  TREAL,    INTENT(INOUT) :: params(isize)
  TREAL,    INTENT(  OUT) :: a(nx*ny*nz)

  ! -------------------------------------------------------------------
  TINTEGER header_offset
  TINTEGER nx_total, ny_total, nz_total

#ifdef USE_MPI
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

  ! -------------------------------------------------------------------
  ! header
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

    header_offset = 0
#include "dns_open_file.h"
    REWIND(LOC_UNIT_ID)
    IF ( iheader .GT. 0 ) THEN
      CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,nt, params)
      isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
    ENDIF
    CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
  CALL MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)

  ! -------------------------------------------------------------------
  ! field
  CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field

  ndims = 3
  sizes(1)   = nx_total;     sizes(2)   = ny_total;     sizes(3)   = nz_total
  locsize(1) = nx;           locsize(2) = ny;           locsize(3) = nz
  offset(1)  = ims_offset_i; offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
      MPI_ORDER_FORTRAN, MPI_REAL8, subarray, ims_err)
  CALL MPI_Type_commit(subarray, ims_err)

  mpio_locsize = nx*ny*nz
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)

  CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_File_read_all(mpio_fh, a, mpio_locsize, MPI_REAL8, status, ims_err)
  CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
  header_offset = header_offset +1
  READ(LOC_UNIT_ID,POS=header_offset) a
  CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_READ_SUBARRAY_FIELD

#undef LOC_UNIT_ID
#undef LOC_STATUS

!########################################################################
!########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_SUBARRAY_FIELD(name, iheader, nx,ny,nz,nt, isize,params, a)

#ifdef USE_MPI
  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_PROCS
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER(LEN=*) name
  TINTEGER,  INTENT(IN   ) :: iheader, nx,ny,nz,nt, isize
  TREAL,     INTENT(IN   ) :: params(isize)
  TREAL,     INTENT(IN   ) :: a(nx*ny*nz)

  ! -------------------------------------------------------------------
  TINTEGER header_offset
  TINTEGER nx_total, ny_total, nz_total

#ifdef USE_MPI
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp, mpio_locoff
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

  ! -------------------------------------------------------------------
  ! header
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

    header_offset = 0
#include "dns_open_file.h"
    IF ( iheader .GT. 0 ) THEN
      CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,nt, params)
      header_offset = 5*SIZEOFINT + isize*SIZEOFREAL
    ENDIF
    CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF

  ! -------------------------------------------------------------------
  ! field
  CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  mpio_disp = header_offset*SIZEOFBYTE

  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

  ndims = 3
  sizes(1)   = nx_total;     sizes(2)   = ny_total;     sizes(3)   = nz_total
  locsize(1) = nx;           locsize(2) = ny;           locsize(3) = nz
  offset(1)  = ims_offset_i; offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
      MPI_ORDER_FORTRAN, MPI_REAL8, subarray, ims_err)
  CALL MPI_Type_commit(subarray, ims_err)

  mpio_locsize = nx*ny*nz
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
  CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_File_write_all(mpio_fh, a, mpio_locsize, MPI_REAL8, status, ims_err)
  CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
  header_offset = header_offset +1
  WRITE(LOC_UNIT_ID,POS=header_offset) a
  CLOSE(LOC_UNIT_ID)
#endif

  RETURN
END SUBROUTINE IO_WRITE_SUBARRAY_FIELD

!########################################################################
!########################################################################
SUBROUTINE IO_READ_HEADER(unit, offset, nx,ny,nz,nt, params)

  USE TLAB_CONSTANTS, ONLY : efile, wfile
  USE TLAB_PROCS

  IMPLICIT NONE

  TINTEGER unit, offset, nx,ny,nz,nt
  TREAL, DIMENSION(*) :: params

  ! -------------------------------------------------------------------
  TINTEGER isize, nx_loc, ny_loc, nz_loc, nt_loc

  !########################################################################
  READ(unit) offset, nx_loc, ny_loc, nz_loc, nt_loc

  isize = offset - 5*SIZEOFINT
  IF ( isize .GT. 0 .AND. MOD(isize,SIZEOFREAL) .EQ. 0 ) THEN
    isize = isize/SIZEOFREAL
    READ(unit) params(1:isize)

  ELSE
    CALL TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Header format incorrect.')
    CALL TLAB_STOP(DNS_ERROR_RECLEN)

  ENDIF

  ! Check
  IF ( nx .NE. nx_loc .OR. ny .NE. ny_loc .OR. nz .NE. nz_loc ) THEN
    CLOSE(unit)
    CALL TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Grid size mismatch.')
    CALL TLAB_STOP(DNS_ERROR_DIMGRID)
  ENDIF

  IF ( nt .NE. nt_loc ) THEN
    CALL TLAB_WRITE_ASCII(wfile, 'IO_READ_HEADER. ItNumber mismatch. Filename value ignored.')
    !     nt = nt_loc
  ENDIF

  RETURN
END SUBROUTINE IO_READ_HEADER

!########################################################################
!########################################################################
SUBROUTINE IO_WRITE_HEADER(unit, isize, nx,ny,nz,nt, params)

  IMPLICIT NONE

  TINTEGER unit, isize, nx,ny,nz,nt
  TREAL, DIMENSION(isize) :: params

  ! -------------------------------------------------------------------
  TINTEGER offset

  !########################################################################
  offset = 5*SIZEOFINT + isize*SIZEOFREAL

  WRITE(unit) offset, nx, ny, nz, nt
  WRITE(unit) params(1:isize)

  RETURN
END SUBROUTINE IO_WRITE_HEADER
