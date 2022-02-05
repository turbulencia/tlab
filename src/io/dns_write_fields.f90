#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define SIZEOFBYTE 1

#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

!########################################################################
!#
!# nfield     In      Number of fields in the file
!# iheader    In      Header control flag:
!#                    0 No header
!#                    1 Scalar fields
!#                    2 Flow fields
!#
!# To be renamed io_write_fields
!########################################################################
SUBROUTINE DNS_WRITE_FIELDS(fname, iheader, nx,ny,nz, nfield, itxc, a, txc)

  USE TLAB_VARS, ONLY : imode_files, imode_eqns
  USE TLAB_CONSTANTS, ONLY : efile, lfile
  USE TLAB_VARS, ONLY : itime, rtime
  USE TLAB_VARS, ONLY : visc, froude, rossby, damkohler, prandtl, mach
  USE TLAB_VARS, ONLY : schmidt
  USE THERMO_VARS, ONLY : gama0
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_j, ims_offset_k
#endif
  USE IO_FIELDS

  IMPLICIT NONE

  CHARACTER(LEN=*) fname
  TINTEGER,  INTENT(IN   ) :: itxc, iheader, nfield, nx,ny,nz
  TREAL,     INTENT(IN   ) :: a(nx*ny*nz,nfield)
  TREAL,     INTENT(INOUT) :: txc(itxc)

  ! -------------------------------------------------------------------
  CHARACTER*32 :: name
  CHARACTER*128 :: line
  TINTEGER nx_total, ny_total, nz_total, header_offset
  TINTEGER ifield

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

  line = 'Writing field '//TRIM(ADJUSTL(fname))//' of size'
  WRITE(name,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(name))
  WRITE(name,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(name))
  WRITE(name,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(name))//'...'
  CALL TLAB_WRITE_ASCII(lfile, line)

  ! ###################################################################
  SELECT CASE( imode_files )

  CASE( DNS_NOFILE )        ! Do nothing

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
      CALL TLAB_WRITE_ASCII(efile, 'DNS_WRITE_FIELDS. Work array size error')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF
#endif

    ! -------------------------------------------------------------------
    ! process header info
    isize = 0
    IF      ( iheader .EQ. 1 ) THEN ! Scalar
      isize = isize+1; params(isize) = rtime
      isize = isize+1; params(isize) = visc ! inverse of reynolds
      isize = isize+1+1                     ! prepare space for schmidt and damkohler

    ELSE IF ( iheader .EQ. 2 ) THEN ! Flow
      isize = isize+1; params(isize) = rtime
      isize = isize+1; params(isize) = visc ! inverse of reynolds
      isize = isize+1; params(isize) = froude
      isize = isize+1; params(isize) = rossby
      IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
        isize = isize+1; params(isize) = gama0
        isize = isize+1; params(isize) = prandtl
        isize = isize+1; params(isize) = mach
      ENDIF
    ENDIF

    IF ( isize .GT. isize_max ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'DNS_WRITE_FIELDS. Parameters array size error.')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF

    ! -------------------------------------------------------------------
    ! write data
    DO ifield = 1,nfield
      IF ( iheader .EQ. 1 ) params(isize-1) = schmidt(ifield)   ! Scalar header
      IF ( iheader .EQ. 1 ) params(isize  ) = damkohler(ifield) ! Scalar header
      WRITE(name,'(I2)') ifield
      name=TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(name))

      ! -------------------------------------------------------------------
      ! header
#ifdef USE_MPI
      IF ( ims_pro .EQ. 0 ) THEN
#endif
        header_offset = 0
#include "dns_open_file.h"
        IF ( iheader .GT. 0 ) THEN
          CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,itime, params)
          header_offset = 5*SIZEOFINT + isize*SIZEOFREAL
        ENDIF
        CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
      ENDIF
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
      CALL MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
      CALL MPI_File_write_all(mpio_fh, a(1,ifield), mpio_locsize, MPI_REAL8, status, ims_err)
      CALL MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
      WRITE(LOC_UNIT_ID,POS=header_offset+1) a(:,ifield)
      CLOSE(LOC_UNIT_ID)
#endif

    ENDDO

  END SELECT

  RETURN
END SUBROUTINE DNS_WRITE_FIELDS

#undef LOC_UNIT_ID
#undef LOC_STATUS
