#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define SIZEOFBYTE 1

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2010/04/04 - J.P. Mellado
!#              Created
!# 2012/12/01 - J.P. Mellado
!#              Not using global variables {imax,jmax,kmax}_total any more
!#
!########################################################################
!# DESCRIPTION
!#
!# Read & write a file of size (nx*ims_npro_i)x(ny*ims_npro_y)x(nz*ims_npro_z)
!#
!# There are PARALLEL (only MPI_IO) and SERIAL modes
!#
!# IEEE double precision floating point representation.
!# Unformatted records
!# Byte ordering is big-endian
!# No embedded record information (different to io_fields_array)
!#
!########################################################################

!########################################################################
!########################################################################
#define LOC_UNIT_ID i54
#define LOC_STATUS 'old'

SUBROUTINE IO_READ_FIELDS_SPLIT(name, iheader, nx,ny,nz,nt, isize,params, a, wrk)

#ifdef USE_MPI
  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_PROCS
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) name
  TINTEGER iheader, nx,ny,nz,nt, isize
  TREAL, DIMENSION(isize)            :: params
  TREAL, DIMENSION(nx*ny*nz), TARGET :: a, wrk

! -------------------------------------------------------------------
  TINTEGER header_offset
  TINTEGER nx_total, ny_total, nz_total

#ifdef USE_MPI
#ifdef USE_MPI_IO
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp, mpio_locoff
  TREAL, DIMENSION(:), POINTER :: p_read
  TINTEGER id, npage
#endif
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

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
  IF ( ims_npro_i .GT. 1 ) THEN
! We always initialize types here. For the general field files, we could
! use TLAB_MPI_I_PARTIAL, but we use this routine for other files like
! buffer regions of transformed fields.
     CALL TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for reading in IO_READ_FIELDS_SPLIT.')
     id = TLAB_MPI_I_AUX1
     npage = nz*ny
     CALL TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( iheader .GT. 0 ) THEN
     IF ( ims_pro .EQ. 0 ) THEN
#include "dns_open_file.h"
        REWIND(LOC_UNIT_ID)
        CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,nt, params)
        CLOSE(LOC_UNIT_ID)
     ENDIF
     CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)

! Displacement to start of field
     mpio_disp = header_offset*SIZEOFBYTE

! Size of array params
     isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL
     CALL MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)

  ELSE
     mpio_disp = 0

  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)

  IF ( ims_npro_i .GT. 1 ) THEN; p_read => wrk
  ELSE;                          p_read => a; ENDIF

  mpio_locsize = nx*ny*nz
  mpio_locoff  = mpio_locsize         ! mpio_locoff might be of type larger than INT4
  mpio_locoff  = ims_pro*mpio_locoff  ! mpio_locoff might be of type larger than INT4
  CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_FILE_READ_AT_ALL(mpio_fh, mpio_locoff, p_read, mpio_locsize, MPI_REAL8, status, ims_err)

  IF ( ims_npro_i .GT. 1 ) THEN
     CALL TLAB_MPI_TRPB_I(p_read, a, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF

  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
  NULLIFY(p_read)

#endif

#else
! ###################################################################
! Serial case
! ###################################################################
#include "dns_open_file.h"
  REWIND(LOC_UNIT_ID)

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( iheader .GT. 0 ) THEN
     CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,nt, params)
! Size of array params
     isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  READ(LOC_UNIT_ID) a

  CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_READ_FIELDS_SPLIT

#undef LOC_UNIT_ID
#undef LOC_STATUS

!########################################################################
!########################################################################
#define LOC_UNIT_ID i55
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_FIELDS_SPLIT(name, iheader, nx,ny,nz,nt, isize,params, a, wrk)

#ifdef USE_MPI
  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_PROCS
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) name
  TINTEGER iheader, isize, nx,ny,nz,nt
  TREAL, DIMENSION(isize)            :: params
  TREAL, DIMENSION(nx*ny*nz), TARGET :: a, wrk

! -------------------------------------------------------------------
  TINTEGER nx_total, ny_total, nz_total

#ifdef USE_MPI
#ifdef USE_MPI_IO
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE), header_offset
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp, mpio_locoff
  TREAL, DIMENSION(:), POINTER :: p_write
  TINTEGER id, npage
#endif

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

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
  IF ( ims_npro_i .GT. 1 ) THEN
! We always initialize types here. For the general field files, we could
! use TLAB_MPI_I_PARTIAL, but we use this routine for other files like
! buffer regions of transformed fields.
     CALL TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for writing in IO_WRITE_FIELDS_SPLIT.')
     id = TLAB_MPI_I_AUX1
     npage = nz*ny
     CALL TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
          ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( ims_pro .EQ. 0 ) THEN
#include "dns_open_file.h"
     IF ( iheader .GT. 0 ) THEN
        CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,nt, params)

! Displacement to start of field
        header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

     ELSE
        header_offset = 0

     ENDIF
     CLOSE(LOC_UNIT_ID)
  ENDIF

  CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  mpio_disp = header_offset*SIZEOFBYTE

! This is an MPI barrier ! do not remove
  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL TLAB_MPI_TRPF_I(a, wrk, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     p_write => wrk
  ELSE
     p_write => a
  ENDIF

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)

  mpio_locsize = nx*ny*nz
  mpio_locoff  = mpio_locsize         ! reclen might be of type larger than INT4
  mpio_locoff  = ims_pro*mpio_locoff  ! reclen might be of type larger than INT4

  CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_FILE_WRITE_AT_ALL(mpio_fh, mpio_locoff, p_write, mpio_locsize, MPI_REAL8, status, ims_err)
  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
  NULLIFY(p_write)

#endif

#else
! ###################################################################
! Serial case
! ###################################################################
#include "dns_open_file.h"

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( iheader .GT. 0 ) THEN
     CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,nt, params)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  WRITE(LOC_UNIT_ID) a

  CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_WRITE_FIELDS_SPLIT

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
