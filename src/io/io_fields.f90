#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define SIZEOFBYTE 1

!########################################################################
!#
!# Read/write files of size (nx*ims_npro_i)x(ny*ims_npro_y)x(nz*ims_npro_z)
!# With header/metadata
!# Unformatted records
!# No embedded record information
!#
!########################################################################

MODULE IO_FIELDS

  USE TLAB_CONSTANTS, ONLY : lfile, wfile, efile
  USE TLAB_CONSTANTS, ONLY : sp,dp
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

  PUBLIC :: IO_READ_FIELDS, IO_WRITE_FIELDS
  PUBLIC :: IO_READ_FIELD_INT1, IO_WRITE_FIELD_INT1
  PUBLIC :: IO_FLOW, IO_SCAL

  PRIVATE

  TINTEGER, PARAMETER :: IO_SCAL = 1 ! Header of scalar field
  TINTEGER, PARAMETER :: IO_FLOW = 2 ! Header of flow field

  TINTEGER nx_total, ny_total, nz_total
  CHARACTER(LEN=32 ) str, name
  CHARACTER(LEN=128) line

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
  SUBROUTINE IO_READ_FIELDS(fname, iheader, nx,ny,nz, nfield,iread, a, txc)
    USE TLAB_VARS, ONLY : imode_files, imode_precision_files
    USE TLAB_VARS, ONLY : itime, rtime, visc
    USE, INTRINSIC :: ISO_C_binding, ONLY : c_f_pointer, c_loc

    CHARACTER(LEN=*) fname
    TINTEGER,  INTENT(IN   ) :: iheader         ! 1 for scalar header, 2 flow header
    TINTEGER,  INTENT(IN   ) :: nfield, iread   ! iread=0 reads all nfields, otherwise iread field
    TINTEGER,  INTENT(IN   ) :: nx,ny,nz
    TREAL,     INTENT(  OUT) :: a(nx*ny*nz,*)
    TREAL,     INTENT(INOUT) :: txc(nx*ny*nz)

    TARGET txc

    ! -------------------------------------------------------------------
    TINTEGER header_offset
    TINTEGER ifield, iz
    REAL(sp), POINTER :: s_wrk(:) => NULL()

    TINTEGER isize_max, isize
    PARAMETER(isize_max=20)
    TREAL params(isize_max)

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

    ! Pass memory address from double precision array to single precision array
    CALL c_f_POINTER(c_LOC(txc), s_wrk, shape=[nx*ny*nz])

    ! ###################################################################
    SELECT CASE( imode_files )

    CASE( IO_NOFILE )        ! Do nothing
      a(:,1:nfield) = C_0_R

    CASE( IO_NETCDF )   ! To be implemented

    CASE DEFAULT              ! One file with header per field
#ifdef USE_MPI
      ndims = 3
      sizes(1)   = nx*ims_npro_i; sizes(2)   = ny;           sizes(3)   = nz*ims_npro_k
      locsize(1) = nx;            locsize(2) = ny;           locsize(3) = nz
      offset(1)  = ims_offset_i;  offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

      IF ( imode_precision_files == IO_TYPE_SINGLE ) THEN ! to be finished; here just as an idea
        CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
            MPI_ORDER_FORTRAN, MPI_REAL4, subarray, ims_err)
      ELSE
        CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
            MPI_ORDER_FORTRAN, MPI_REAL8, subarray, ims_err)
      END IF
      CALL MPI_Type_commit(subarray, ims_err)
#endif

      ! -------------------------------------------------------------------
      ! read data
      iz = 0
      DO ifield = 1,nfield
        IF ( iread == 0 .OR. iread == ifield ) THEN
          WRITE(name,'(I2)') ifield
          name=TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(name))
          isize = isize_max

          ! -------------------------------------------------------------------
          ! header
#ifdef USE_MPI
          IF ( ims_pro == 0 ) THEN
#endif

            header_offset = 0
#include "dns_open_file.h"
            REWIND(LOC_UNIT_ID)
            CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,itime, params)
            isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
            CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
          END IF
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
          IF ( imode_precision_files == IO_TYPE_SINGLE ) THEN ! to be finished; here just as an idea
            CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
            CALL MPI_File_read_all(mpio_fh, s_wrk, mpio_locsize, MPI_REAL4, status, ims_err)
            a(:,iz) = REAL(s_wrk(:),dp)
          ELSE
            CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
            CALL MPI_File_read_all(mpio_fh, a(1,iz), mpio_locsize, MPI_REAL8, status, ims_err)
          END IF
          CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
          IF ( imode_precision_files == IO_TYPE_SINGLE ) THEN ! to be finished; here just as an idea
            READ(LOC_UNIT_ID,POS=header_offset+1) s_wrk(:)
            a(:,iz) = REAL(s_wrk(:),dp)
    ELSE
            READ(LOC_UNIT_ID,POS=header_offset+1) a(:,iz)
          END IF
          CLOSE(LOC_UNIT_ID)

#endif

        END IF
      END DO

      ! -------------------------------------------------------------------
      ! process header info
      IF ( isize > isize_max ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'IO_READ_FIELDS. Parameters array size error')
        CALL TLAB_STOP(DNS_ERROR_ALLOC)
      END IF

      rtime = params(1)
      IF ( iheader == IO_FLOW ) THEN
        visc  = params(2)
      END IF
#ifdef USE_MPI
      CALL MPI_BCAST(itime, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
      CALL MPI_BCAST(rtime, 1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
      IF ( iheader == IO_FLOW ) THEN
        CALL MPI_BCAST(visc,  1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
      END IF
#endif

    END SELECT

    RETURN
  END SUBROUTINE IO_READ_FIELDS

  !########################################################################
  !########################################################################
  SUBROUTINE IO_READ_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a, wrk)

#include "integers.h"

    CHARACTER(LEN=*) name
    TINTEGER, INTENT(IN   ) :: header_offset, nx,ny,nz
    TREAL,    INTENT(  OUT) :: a(nx*ny*nz)
    TREAL,    INTENT(INOUT) :: wrk(nx*ny*nz)

    TARGET a, wrk

    ! ###################################################################
#ifdef USE_MPI
    mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field

    IF ( ims_npro_i > 1 ) THEN
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

    END IF

    mpio_locsize = nx*ny*nz
    mpio_locoff  = mpio_locsize         ! mpio_locoff might be of type larger than INT4
    mpio_locoff  = ims_pro*mpio_locoff  ! mpio_locoff might be of type larger than INT4
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
    CALL MPI_FILE_READ_AT_ALL(mpio_fh, mpio_locoff, p_read, mpio_locsize, MPI_REAL8, status, ims_err)
    CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

    IF ( ims_npro_i > 1 ) THEN
      CALL TLAB_MPI_TRPB_I(p_read, a, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
    END IF

    NULLIFY(p_read)

#else
#include "dns_open_file.h"
    READ(LOC_UNIT_ID,POS=header_offset +1) a
    CLOSE(LOC_UNIT_ID)

#endif

    RETURN
  END SUBROUTINE IO_READ_FIELD_XPENCIL

  !########################################################################
  !########################################################################
  SUBROUTINE IO_READ_FIELD_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)

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

#ifdef USE_MPI
    ndims = 3
    sizes(1)   = nx_total;     sizes(2)   = ny_total;     sizes(3)   = nz_total
    locsize(1) = nx;           locsize(2) = ny;           locsize(3) = nz
    offset(1)  = ims_offset_i; offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

    CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
        MPI_ORDER_FORTRAN, MPI_INTEGER1, subarray, ims_err)
    CALL MPI_Type_commit(subarray, ims_err)
#endif

    ! -------------------------------------------------------------------
    ! header
#ifdef USE_MPI
    IF ( ims_pro == 0 ) THEN
#endif

      header_offset = 0
#include "dns_open_file.h"
      REWIND(LOC_UNIT_ID)
      IF ( iheader > 0 ) THEN
        CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,nt, params)
        isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
      END IF
      CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
    END IF
    CALL MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
    CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)

    ! -------------------------------------------------------------------
    ! field
    mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field
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

  SUBROUTINE IO_WRITE_FIELDS(fname, iheader, nx,ny,nz, nfield, a, txc)

    USE TLAB_VARS, ONLY : imode_files, imode_precision_files, imode_eqns
    USE TLAB_VARS, ONLY : itime, rtime
    USE TLAB_VARS, ONLY : visc, froude, rossby, damkohler, prandtl, mach
    USE TLAB_VARS, ONLY : schmidt
    USE THERMO_VARS, ONLY : gama0
    USE, INTRINSIC :: ISO_C_binding, ONLY : c_f_pointer, c_loc

    CHARACTER(LEN=*) fname
    TINTEGER,  INTENT(IN   ) :: iheader         ! Scalar or Flow headers
    TINTEGER,  INTENT(IN   ) :: nfield
    TINTEGER,  INTENT(IN   ) :: nx,ny,nz
    TREAL,     INTENT(IN   ) :: a(nx*ny*nz,nfield)
    TREAL,     INTENT(INOUT) :: txc(nx*ny*nz)

    TARGET txc

    ! -------------------------------------------------------------------
    TINTEGER header_offset
    TINTEGER ifield
    REAL(sp), POINTER :: s_wrk(:) => NULL()

    TINTEGER isize_max, isize
    PARAMETER(isize_max=20)
    TREAL params(isize_max)

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

    line = 'Writing field '//TRIM(ADJUSTL(fname))//' of size'
    WRITE(name,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(name))
    WRITE(name,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(name))
    WRITE(name,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(name))//'...'
    CALL TLAB_WRITE_ASCII(lfile, line)

    ! Pass memory address from double precision array to single precision array
    CALL c_f_POINTER(c_LOC(txc), s_wrk, shape=[nx*ny*nz])

    ! ###################################################################
    SELECT CASE( imode_files )

    CASE( IO_NOFILE )        ! Do nothing

    CASE( IO_NETCDF )   ! To be implemented

    CASE DEFAULT              ! One file with header per field
#ifdef USE_MPI
      ndims = 3
      sizes(1)   = nx*ims_npro_i; sizes(2)   = ny;           sizes(3)   = nz*ims_npro_k
      locsize(1) = nx;            locsize(2) = ny;           locsize(3) = nz
      offset(1)  = ims_offset_i;  offset(2)  = ims_offset_j; offset(3)  = ims_offset_k

      IF ( imode_precision_files == IO_TYPE_SINGLE ) THEN ! to be finished; here just as an idea
        CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
            MPI_ORDER_FORTRAN, MPI_REAL4, subarray, ims_err)
      ELSE
        CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
            MPI_ORDER_FORTRAN, MPI_REAL8, subarray, ims_err)
      END IF
      CALL MPI_Type_commit(subarray, ims_err)
#endif

      ! -------------------------------------------------------------------
      ! process header info
      isize = 0
      isize = isize+1; params(isize) = rtime
      isize = isize+1; params(isize) = visc ! inverse of reynolds
      IF      ( iheader == IO_SCAL ) THEN
        isize = isize+1+1                     ! prepare space for schmidt and damkohler

      ELSE IF ( iheader == IO_FLOW ) THEN
        isize = isize+1; params(isize) = froude
        isize = isize+1; params(isize) = rossby
        IF ( imode_eqns == DNS_EQNS_INTERNAL .OR. imode_eqns == DNS_EQNS_TOTAL ) THEN
          isize = isize+1; params(isize) = gama0
          isize = isize+1; params(isize) = prandtl
          isize = isize+1; params(isize) = mach
        END IF
      END IF

      IF ( isize > isize_max ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'IO_WRITE_FIELDS. Parameters array size error.')
        CALL TLAB_STOP(DNS_ERROR_ALLOC)
      END IF

      ! -------------------------------------------------------------------
      ! write data
      DO ifield = 1,nfield
        IF ( iheader == IO_SCAL ) params(isize-1) = schmidt(ifield)   ! Scalar header
        IF ( iheader == IO_SCAL ) params(isize  ) = damkohler(ifield) ! Scalar header
        WRITE(name,'(I2)') ifield
        name=TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(name))

        ! -------------------------------------------------------------------
        ! header
#ifdef USE_MPI
        IF ( ims_pro == 0 ) THEN
#endif
          header_offset = 0
#include "dns_open_file.h"
          CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,itime, params)
          header_offset = 5*SIZEOFINT + isize*SIZEOFREAL
          CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
        END IF
        CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

        ! -------------------------------------------------------------------
        ! field
        ! CALL IO_WRITE_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a(1,ifield),txc)
#ifdef USE_MPI
        CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

        mpio_disp = header_offset*SIZEOFBYTE
        mpio_locsize = nx*ny*nz
        CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        IF ( imode_precision_files == IO_TYPE_SINGLE ) THEN ! to be finished; here just as an idea
          CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
          s_wrk(:) = REAL(a(:,ifield),sp)
          CALL MPI_File_write_all(mpio_fh, s_wrk, mpio_locsize, MPI_REAL4, status, ims_err)
        ELSE
          CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
          CALL MPI_File_write_all(mpio_fh, a(1,ifield), mpio_locsize, MPI_REAL8, status, ims_err)
        END IF
        CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
        IF ( imode_precision_files == IO_TYPE_SINGLE ) THEN ! to be finished; here just as an idea
          s_wrk(:) = REAL(a(:,ifield),sp)
          WRITE(LOC_UNIT_ID,POS=header_offset+1) s_wrk(:)
        ELSE
          WRITE(LOC_UNIT_ID,POS=header_offset+1) a(:,ifield)
        END IF
        CLOSE(LOC_UNIT_ID)
#endif

      END DO

    END SELECT

    RETURN
  END SUBROUTINE IO_WRITE_FIELDS

  !########################################################################
  !########################################################################
  SUBROUTINE IO_WRITE_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a, wrk)

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

    IF ( ims_npro_i > 1 ) THEN
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

    END IF

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

  !########################################################################
  !########################################################################
  SUBROUTINE IO_WRITE_FIELD_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)

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
    IF ( ims_pro == 0 ) THEN
#endif

      header_offset = 0
#include "dns_open_file.h"
      IF ( iheader > 0 ) THEN
        CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,nt, params)
        header_offset = 5*SIZEOFINT + isize*SIZEOFREAL
      END IF
      CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
    END IF

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
    TINTEGER unit, offset, nx,ny,nz,nt
    TREAL, DIMENSION(*) :: params

    ! -------------------------------------------------------------------
    TINTEGER isize, nx_loc, ny_loc, nz_loc, nt_loc

    !########################################################################
    READ(unit) offset, nx_loc, ny_loc, nz_loc, nt_loc

    isize = offset - 5*SIZEOFINT
    IF ( isize > 0 .AND. MOD(isize,SIZEOFREAL) == 0 ) THEN
      isize = isize/SIZEOFREAL
      READ(unit) params(1:isize)

    ELSE
      CALL TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Header format incorrect.')
      CALL TLAB_STOP(DNS_ERROR_RECLEN)

    END IF

    ! Check
    IF ( nx .NE. nx_loc .OR. ny .NE. ny_loc .OR. nz .NE. nz_loc ) THEN
      CLOSE(unit)
      CALL TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Grid size mismatch.')
      CALL TLAB_STOP(DNS_ERROR_DIMGRID)
    END IF

    IF ( nt .NE. nt_loc ) THEN
      CALL TLAB_WRITE_ASCII(wfile, 'IO_READ_HEADER. ItNumber mismatch. Filename value ignored.')
      !     nt = nt_loc
    END IF

    RETURN
  END SUBROUTINE IO_READ_HEADER

  !########################################################################
  !########################################################################
  SUBROUTINE IO_WRITE_HEADER(unit, isize, nx,ny,nz,nt, params)
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
