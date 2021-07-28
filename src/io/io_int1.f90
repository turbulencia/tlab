#include "types.h"
#include "dns_error.h"

#define SIZEOFBYTE 1

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2008/04/02 - J.P. Mellado
!#              Created
!# 2015/03/05 - J.P. Mellado
!#              Updated to new general structure
!#
!########################################################################
!# DESCRIPTION
!#
!# Writing gate field. Generalized to writing fields of type integer(1).
!#
!########################################################################

!########################################################################
!########################################################################
#define LOC_UNIT_ID i34
#define LOC_STATUS 'old'

SUBROUTINE IO_READ_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)

  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k, ims_err
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*),                   INTENT(IN)    :: name
  TINTEGER,                        INTENT(IN)    :: iheader, nx,ny,nz,nt
  TINTEGER,                        INTENT(INOUT) :: isize
  TREAL,      DIMENSION(isize),    INTENT(OUT)   :: params
  INTEGER(1), DIMENSION(nx*ny*nz), INTENT(OUT)   :: a

! -------------------------------------------------------------------
  CHARACTER*32 str
  CHARACTER*128 line
  TINTEGER header_offset
  TINTEGER nx_total,ny_total,nz_total

#ifdef USE_MPI
#ifdef USE_MPI_IO
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp
  TINTEGER                :: subarray, ndims
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset
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

  line = 'Reading field '//TRIM(ADJUSTL(name))//' of size'
  WRITE(str,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(str))
  WRITE(str,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))
  WRITE(str,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))//'...'

  CALL TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
#ifdef USE_MPI_IO
  ndims = 3
  sizes(1)   = nx_total;     sizes(2)   = ny_total;     sizes(3)   = nz_total
  locsize(1) = nx;           locsize(2) = ny;           locsize(3) = nz
  offset(1)  = ims_offset_i; offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
       MPI_ORDER_FORTRAN, MPI_INTEGER1, subarray, ims_err)
  CALL MPI_Type_commit(subarray, ims_err)
#endif
#endif

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
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
  mpio_locsize = nx*ny*nz

  CALL MPI_File_open(MPI_COMM_WORLD,TRIM(ADJUSTL(name)),&
       MPI_MODE_RDONLY,MPI_INFO_NULL,mpio_fh, ims_err)

  CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_INTEGER1, subarray, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_File_read_all(mpio_fh, a, mpio_locsize, MPI_INTEGER1, status, ims_err)
  CALL MPI_File_close(mpio_fh, ims_err)
#endif

#else
! ###################################################################
! Serial case
! ###################################################################
#include "dns_open_file.h"
  REWIND(LOC_UNIT_ID)
  IF ( iheader .GT. 0 ) THEN
     CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,nt, params)
     isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
  ENDIF
  READ(LOC_UNIT_ID) a
  CLOSE(LOC_UNIT_ID)
#endif

  RETURN
END SUBROUTINE IO_READ_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

!########################################################################
!########################################################################
#define LOC_UNIT_ID i35
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)

  USE TLAB_CONSTANTS, ONLY : lfile
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k, ims_err
#endif
  USE TLAB_PROCS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*),                   INTENT(IN) :: name
  TINTEGER,                        INTENT(IN) :: iheader, nx,ny,nz,nt, isize
  TREAL,      DIMENSION(isize),    INTENT(IN) :: params
  INTEGER(1), DIMENSION(nx*ny*nz), INTENT(IN) :: a

! -------------------------------------------------------------------
  CHARACTER*32  str
  CHARACTER*128 line
  TINTEGER nx_total,ny_total,nz_total

#ifdef USE_MPI
#ifdef USE_MPI_IO
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE), header_offset
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp
  TINTEGER                :: subarray, ndims
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset
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

  line = 'Writing field '//TRIM(ADJUSTL(name))//' of size'
  WRITE(str,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(str))
  WRITE(str,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))
  WRITE(str,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))//'...'

  CALL TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
#ifdef USE_MPI_IO
  ndims = 3
  sizes(1)   = nx_total;     sizes(2)   = ny_total;     sizes(3)   = nz_total
  locsize(1) = nx;           locsize(2) = ny;           locsize(3) = nz
  offset(1)  = ims_offset_i; offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
       MPI_ORDER_FORTRAN, MPI_INTEGER1, subarray, ims_err)
  CALL MPI_Type_commit(subarray, ims_err)
#endif
#endif

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( ims_pro .EQ. 0 ) THEN
#include "dns_open_file.h"
     IF ( iheader .GT. 0 ) THEN
        CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,nt, params)
        header_offset = 5*SIZEOFINT + isize*SIZEOFREAL ! Displacement to start of field
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
  mpio_locsize = nx*ny*nz

  CALL MPI_File_open(MPI_COMM_WORLD,TRIM(ADJUSTL(name)),&
       IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,mpio_fh, ims_err)

  CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_INTEGER1, subarray, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_File_write_all(mpio_fh, a, mpio_locsize, MPI_INTEGER1, status, ims_err)
  CALL MPI_File_close(mpio_fh, ims_err)

#endif

#else
! ###################################################################
! Serial case
! ###################################################################
#include "dns_open_file.h"
  IF ( iheader .GT. 0 ) THEN
     CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,nt, params)
  ENDIF
  WRITE(LOC_UNIT_ID) a
  CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_WRITE_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS
