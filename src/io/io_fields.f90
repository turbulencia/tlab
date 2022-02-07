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
!# with header
!#
!# Unformatted records
!# No embedded record information
!#
!########################################################################

MODULE IO_FIELDS

  USE TLAB_CONSTANTS, ONLY : lfile, wfile, efile
  USE TLAB_PROCS, ONLY : TLAB_STOP, TLAB_WRITE_ASCII
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k
  !
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

  PUBLIC :: IO_READ_FIELD_INT1, IO_WRITE_FIELD_INT1
  PUBLIC :: IO_READ_FIELD_XPENCIL, IO_WRITE_FIELD_XPENCIL
  PUBLIC :: IO_READ_HEADER, IO_WRITE_HEADER

  PRIVATE

  TINTEGER nx_total, ny_total, nz_total
  CHARACTER*32 str
  CHARACTER*128 line

#ifdef USE_MPI
#include "mpif.h"
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp
  TINTEGER                :: subarray, ndims
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset
  !
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_locoff
  TREAL, DIMENSION(:), POINTER :: p_read, p_write
  TINTEGER id, npage
#endif

CONTAINS

#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'

  !########################################################################
  !########################################################################
    SUBROUTINE IO_READ_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a, wrk)
    IMPLICIT NONE

#include "integers.h"

    CHARACTER(LEN=*) name
    TINTEGER, INTENT(IN   ) :: header_offset, nx,ny,nz
    TREAL,    INTENT(  OUT) :: a(nx*ny*nz)
    TREAL,    INTENT(INOUT) :: wrk(nx*ny*nz)

    TARGET a, wrk

    ! ###################################################################
#ifdef USE_MPI
    mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field

    IF ( ims_npro_i .GT. 1 ) THEN
      ! We always initialize types here. For the general field files, we could
      ! use TLAB_MPI_I_PARTIAL, but we use this routine for other files.
      CALL TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for reading in IO_READ_FIELDS_SPLIT.')
      id = TLAB_MPI_I_AUX1
      npage = nz*ny
      CALL TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

      p_read => wrk

    ELSE
      p_read => a

    ENDIF

    mpio_locsize = nx*ny*nz
    mpio_locoff  = mpio_locsize         ! mpio_locoff might be of type larger than INT4
    mpio_locoff  = ims_pro*mpio_locoff  ! mpio_locoff might be of type larger than INT4
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
    CALL MPI_FILE_READ_AT_ALL(mpio_fh, mpio_locoff, p_read, mpio_locsize, MPI_REAL8, status, ims_err)
    CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

    IF ( ims_npro_i .GT. 1 ) THEN
      CALL TLAB_MPI_TRPB_I(p_read, a, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
    ENDIF

    NULLIFY(p_read)

#else
#include "dns_open_file.h"
    READ(LOC_UNIT_ID,POS=header_offset +1) a
    CLOSE(LOC_UNIT_ID)

#endif

    RETURN
  END SUBROUTINE IO_READ_FIELD_XPENCIL

#undef LOC_UNIT_ID
#undef LOC_STATUS

#define LOC_UNIT_ID 34
#define LOC_STATUS 'old'
  !########################################################################
  !########################################################################
  SUBROUTINE IO_READ_FIELD_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)
    IMPLICIT NONE

    CHARACTER(LEN=*) name
    TINTEGER,   INTENT(IN   ) :: iheader, nx,ny,nz,nt
    TINTEGER,   INTENT(INOUT) :: isize
    TREAL,      INTENT(INOUT) :: params(isize)
    INTEGER(1), INTENT(  OUT) :: a(nx*ny*nz)

    TINTEGER header_offset

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

    line = 'Reading field '//TRIM(ADJUSTL(name))//' of size'
    WRITE(str,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(str))
    WRITE(str,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))
    WRITE(str,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))//'...'
    CALL TLAB_WRITE_ASCII(lfile, line)

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
        MPI_ORDER_FORTRAN, MPI_INTEGER1, subarray, ims_err)
    CALL MPI_Type_commit(subarray, ims_err)

    mpio_locsize = nx*ny*nz
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_INTEGER1, subarray, 'native', MPI_INFO_NULL, ims_err)
    CALL MPI_File_read_all(mpio_fh, a, mpio_locsize, MPI_INTEGER1, status, ims_err)
    CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
    header_offset = header_offset +1
    READ(LOC_UNIT_ID,POS=header_offset) a
    CLOSE(LOC_UNIT_ID)

#endif

    RETURN
  END SUBROUTINE IO_READ_FIELD_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

  !########################################################################
  !########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

  SUBROUTINE IO_WRITE_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a, wrk)
    IMPLICIT NONE

#include "integers.h"

    CHARACTER(LEN=*) name
    TINTEGER,  INTENT(IN   ) :: header_offset, nx,ny,nz
    TREAL,     INTENT(IN   ) :: a(nx*ny*nz)
    TREAL,     INTENT(INOUT) :: wrk(nx*ny*nz)

    TARGET a, wrk

    ! ###################################################################
#ifdef USE_MPI
    mpio_disp = header_offset*SIZEOFBYTE

    CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    IF ( ims_npro_i .GT. 1 ) THEN
      ! We always initialize types here. For the general field files, we could
      ! use TLAB_MPI_I_PARTIAL, but we use this routine for other files.
      CALL TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for writing in IO_WRITE_FIELDS_SPLIT.')
      id = TLAB_MPI_I_AUX1
      npage = nz*ny
      CALL TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

      CALL TLAB_MPI_TRPF_I(a, wrk, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
      p_write => wrk

    ELSE
      p_write => a

    ENDIF

    mpio_locsize = nx*ny*nz
    mpio_locoff  = mpio_locsize         ! reclen might be of type larger than INT4
    mpio_locoff  = ims_pro*mpio_locoff  ! reclen might be of type larger than INT4
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
    CALL MPI_FILE_WRITE_AT_ALL(mpio_fh, mpio_locoff, p_write, mpio_locsize, MPI_REAL8, status, ims_err)
    CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
    NULLIFY(p_write)

#else
#include "dns_open_file.h"
    WRITE(LOC_UNIT_ID,POS=header_offset+1) a
    CLOSE(LOC_UNIT_ID)
#endif

    RETURN
  END SUBROUTINE IO_WRITE_FIELD_XPENCIL

#undef LOC_UNIT_ID
#undef LOC_STATUS

#define LOC_UNIT_ID 35
#define LOC_STATUS 'unknown'
  !########################################################################
  !########################################################################
  SUBROUTINE IO_WRITE_FIELD_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)
    IMPLICIT NONE

    CHARACTER(LEN=*) name
    TINTEGER,   INTENT(IN   ) :: iheader, nx,ny,nz,nt, isize
    TREAL,      INTENT(IN   ) :: params(isize)
    INTEGER(1), INTENT(IN   ) :: a(nx*ny*nz)

    TINTEGER header_offset

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

    line = 'Writing field '//TRIM(ADJUSTL(name))//' of size'
    WRITE(str,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(str))
    WRITE(str,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))
    WRITE(str,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))//'...'
    CALL TLAB_WRITE_ASCII(lfile, line)

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
        MPI_ORDER_FORTRAN, MPI_INTEGER1, subarray, ims_err)
    CALL MPI_Type_commit(subarray, ims_err)

    mpio_locsize = nx*ny*nz
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
    CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_INTEGER1, subarray, 'native', MPI_INFO_NULL, ims_err)
    CALL MPI_File_write_all(mpio_fh, a, mpio_locsize, MPI_INTEGER1, status, ims_err)
    CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
    header_offset = header_offset +1
    WRITE(LOC_UNIT_ID,POS=header_offset) a
    CLOSE(LOC_UNIT_ID)
#endif

    RETURN
  END SUBROUTINE IO_WRITE_FIELD_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

  !########################################################################
  !########################################################################
  SUBROUTINE IO_READ_HEADER(unit, offset, nx,ny,nz,nt, params)
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

END MODULE IO_FIELDS
