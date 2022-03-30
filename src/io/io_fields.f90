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

module IO_FIELDS

  use TLAB_CONSTANTS, only : lfile, wfile, efile
  use TLAB_CONSTANTS, only : sp,dp
  use TLAB_PROCS, only : TLAB_STOP, TLAB_WRITE_ASCII
#ifdef USE_MPI
  use MPI
  use TLAB_MPI_VARS, only : ims_err
  use TLAB_MPI_VARS, only : ims_pro, ims_npro_i, ims_npro_k
  use TLAB_MPI_VARS, only : ims_offset_i, ims_offset_j, ims_offset_k
#endif

  implicit none

  public :: IO_READ_FIELDS, IO_WRITE_FIELDS
  public :: IO_READ_FIELD_INT1, IO_WRITE_FIELD_INT1
  public :: IO_FLOW, IO_SCAL
#ifdef USE_MPI
  public :: IO_CREATE_SUBARRAY_XOY, IO_CREATE_SUBARRAY_XOZ, IO_CREATE_SUBARRAY_ZOY
#endif

  private

  TINTEGER, parameter :: IO_SCAL = 1 ! Header of scalar field
  TINTEGER, parameter :: IO_FLOW = 2 ! Header of flow field

  TINTEGER nx_total, ny_total, nz_total
  character(LEN=32 ) str, name
  character(LEN=128) line

#ifdef USE_MPI
  integer mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
  integer(KIND=MPI_OFFSET_KIND) mpio_disp
  TINTEGER subarray
#endif

contains

#ifdef USE_MPI
  !########################################################################
  !########################################################################
  function IO_CREATE_SUBARRAY_XOY( nx,ny, mpi_type ) result( subarray )
    TINTEGER, intent(in   ) :: nx,ny
    integer,  intent(in   ) :: mpi_type

    integer  :: subarray
    TINTEGER, parameter :: ndims = 2
    TINTEGER :: sizes(ndims), locsize(ndims), offset(ndims)

    sizes   = [ nx*ims_npro_i, ny           ]
    locsize = [ nx,            ny           ]
    offset  = [ ims_offset_i,  ims_offset_j ]

    call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
         MPI_ORDER_FORTRAN, mpi_type, subarray, ims_err)
    call MPI_Type_commit(subarray, ims_err)

  end function IO_CREATE_SUBARRAY_XOY

  function IO_CREATE_SUBARRAY_XOZ(nx,ny,nz,mpi_type) result(subarray)
    TINTEGER, intent(in   ) :: nx,ny,nz
    integer,  intent(in   ) :: mpi_type

    integer  :: subarray
    TINTEGER, parameter :: ndims = 3
    TINTEGER :: sizes(ndims), locsize(ndims), offset(ndims)

    sizes   = [ nx*ims_npro_i, ny,           nz*ims_npro_k ]
    locsize = [ nx,            ny,           nz            ]
    offset  = [ ims_offset_i,  ims_offset_j, ims_offset_k  ]

    call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
         MPI_ORDER_FORTRAN, mpi_type, subarray, ims_err)
    call MPI_Type_commit(subarray, ims_err)

  end function IO_CREATE_SUBARRAY_XOZ

  function IO_CREATE_SUBARRAY_ZOY(ny,nz, mpi_type) result(subarray)
    TINTEGER, intent(in   ) :: ny,nz
    integer,  intent(in   ) :: mpi_type

    integer  :: subarray
    TINTEGER, parameter :: ndims = 2
    TINTEGER :: sizes(ndims), locsize(ndims), offset(ndims)

    sizes   = [ ny,           nz*ims_npro_k ]
    locsize = [ ny,           nz            ]
    offset  = [ ims_offset_j, ims_offset_k  ]

    call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
         MPI_ORDER_FORTRAN, mpi_type, subarray, ims_err)
    call MPI_Type_commit(subarray, ims_err)

  end function IO_CREATE_SUBARRAY_ZOY
#endif

  !########################################################################
  !########################################################################
#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'

  subroutine IO_READ_FIELDS(fname, iheader, nx,ny,nz, nfield,iread, a, txc)
    use TLAB_VARS, only : imode_files, imode_precision_files
    use TLAB_VARS, only : itime, rtime, visc
    use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

    character(LEN=*) fname
    TINTEGER,  intent(IN   ) :: iheader         ! 1 for scalar header, 2 flow header
    TINTEGER,  intent(IN   ) :: nfield, iread   ! iread=0 reads all nfields, otherwise iread field
    TINTEGER,  intent(IN   ) :: nx,ny,nz
    TREAL,     intent(  OUT) :: a(nx*ny*nz,*)
    TREAL,     intent(INOUT) :: txc(nx*ny*nz)

    target txc

    ! -------------------------------------------------------------------
    TINTEGER header_offset
    TINTEGER ifield, iz
    real(sp), pointer :: s_wrk(:) => null()

    TINTEGER, parameter :: isize_max = 20
    TREAL params(isize_max)
    TINTEGER isize

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

    line = 'Reading field '//trim(adjustl(fname))//' of size'
    write(name,*) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(name))
    write(name,*) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))
    write(name,*) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))//'...'
    call TLAB_WRITE_ASCII(lfile, line)

    ! Pass memory address from double precision array to single precision array
    call c_f_pointer(c_loc(txc), s_wrk, shape=[nx*ny*nz])

    ! ###################################################################
    select case( imode_files )

    case( IO_NOFILE )         ! Do nothing
      a(:,1:nfield) = C_0_R

    case( IO_NETCDF )         ! To be implemented

    case DEFAULT              ! One file with header per field
#ifdef USE_MPI
      if ( imode_precision_files == IO_TYPE_SINGLE ) then ! to be finished; here just as an idea
        subarray = IO_CREATE_SUBARRAY_XOZ( nx,ny,nz, MPI_REAL4 )
      else
        subarray = IO_CREATE_SUBARRAY_XOZ( nx,ny,nz, MPI_REAL8 )
      end if
#endif

      ! -------------------------------------------------------------------
      ! read data
      iz = 0
      do ifield = 1,nfield
        if ( iread == 0 .or. iread == ifield ) then
          iz = iz +1
          write(name,'(I2)') ifield
          name=trim(adjustl(fname))//'.'//trim(adjustl(name))

          ! -------------------------------------------------------------------
          ! header
#ifdef USE_MPI
          if ( ims_pro == 0 ) then
#endif
#include "dns_open_file.h"
            rewind(LOC_UNIT_ID)
            call IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,itime, params)
            close(LOC_UNIT_ID)

#ifdef USE_MPI
          end if
          call MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
          ! isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
          ! if ( isize > 0 ) call MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)

          ! -------------------------------------------------------------------
          ! field
          ! CALL IO_READ_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a(1,iz),txc)
          mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field
          mpio_locsize = nx*ny*nz
          call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
          if ( imode_precision_files == IO_TYPE_SINGLE ) then ! to be finished; here just as an idea
            call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
            call MPI_File_read_all(mpio_fh, s_wrk, mpio_locsize, MPI_REAL4, status, ims_err)
            a(:,iz) = real(s_wrk(:),dp)
          else
            call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
            call MPI_File_read_all(mpio_fh, a(1,iz), mpio_locsize, MPI_REAL8, status, ims_err)
          end if
          call MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
          if ( imode_precision_files == IO_TYPE_SINGLE ) then ! to be finished; here just as an idea
            read(LOC_UNIT_ID,POS=header_offset+1) s_wrk(:)
            a(:,iz) = real(s_wrk(:),dp)
          else
            read(LOC_UNIT_ID,POS=header_offset+1) a(:,iz)
          end if
          close(LOC_UNIT_ID)

#endif

        end if
      end do

      ! -------------------------------------------------------------------
      ! process header info
      isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
      if ( isize > isize_max ) then
        call TLAB_WRITE_ASCII(efile, 'IO_READ_FIELDS. Parameters array size error')
        call TLAB_STOP(DNS_ERROR_ALLOC)
      end if

      rtime = params(1)
      if ( iheader == IO_FLOW ) then
        visc  = params(2)
      end if
#ifdef USE_MPI
      call MPI_BCAST(itime, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
      call MPI_BCAST(rtime, 1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
      if ( iheader == IO_FLOW ) then
        call MPI_BCAST(visc,  1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
      end if
#endif

    end select

    return
  end subroutine IO_READ_FIELDS

  !########################################################################
  !########################################################################
  subroutine IO_READ_FIELD_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)

    character(LEN=*) name
    TINTEGER,   intent(IN   ) :: iheader, nx,ny,nz,nt
    TINTEGER,   intent(INOUT) :: isize
    TREAL,      intent(INOUT) :: params(isize)
    integer(1), intent(  OUT) :: a(nx*ny*nz)

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

    line = 'Reading field '//trim(adjustl(name))//' of size'
    write(str,*) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(str))
    write(str,*) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))
    write(str,*) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))//'...'
    call TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
    subarray = IO_CREATE_SUBARRAY_XOZ( nx,ny,nz, MPI_INTEGER1 )
#endif

    ! -------------------------------------------------------------------
    ! header
#ifdef USE_MPI
    if ( ims_pro == 0 ) then
#endif
      header_offset = 0
#include "dns_open_file.h"
      rewind(LOC_UNIT_ID)
      if ( iheader > 0 ) then
        call IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total,ny_total,nz_total,nt, params)
      end if
      close(LOC_UNIT_ID)

#ifdef USE_MPI
    end if
    call MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)

    ! -------------------------------------------------------------------
    ! field
    mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field
    mpio_locsize = nx*ny*nz
    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    call MPI_File_set_view(mpio_fh, mpio_disp, MPI_INTEGER1, subarray, 'native', MPI_INFO_NULL, ims_err)
    call MPI_File_read_all(mpio_fh, a, mpio_locsize, MPI_INTEGER1, status, ims_err)
    call MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
    header_offset = header_offset +1
    read(LOC_UNIT_ID,POS=header_offset) a
    close(LOC_UNIT_ID)

#endif

    ! -------------------------------------------------------------------
    ! process header info
    isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
#ifdef USE_MPI
    if ( isize > 0 ) call MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

    return
  end subroutine IO_READ_FIELD_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

  !########################################################################
  !########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

  subroutine IO_WRITE_FIELDS(fname, iheader, nx,ny,nz, nfield, a, txc)

    use TLAB_VARS, only : imode_files, imode_precision_files, imode_eqns
    use TLAB_VARS, only : itime, rtime
    use TLAB_VARS, only : visc, froude, rossby, damkohler, prandtl, mach
    use TLAB_VARS, only : schmidt
    use THERMO_VARS, only : gama0
    use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

    character(LEN=*) fname
    TINTEGER,  intent(IN   ) :: iheader         ! Scalar or Flow headers
    TINTEGER,  intent(IN   ) :: nfield
    TINTEGER,  intent(IN   ) :: nx,ny,nz
    TREAL,     intent(IN   ) :: a(nx*ny*nz,nfield)
    TREAL,     intent(INOUT) :: txc(nx*ny*nz)

    target txc

    ! -------------------------------------------------------------------
    TINTEGER header_offset
    TINTEGER ifield
    real(sp), pointer :: s_wrk(:) => null()

    TINTEGER, parameter :: isize_max = 20
    TREAL params(isize_max)
    TINTEGER isize

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

    line = 'Writing field '//trim(adjustl(fname))//' of size'
    write(name,*) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(name))
    write(name,*) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))
    write(name,*) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))//'...'
    call TLAB_WRITE_ASCII(lfile, line)

    ! Pass memory address from double precision array to single precision array
    call c_f_pointer(c_loc(txc), s_wrk, shape=[nx*ny*nz])

    ! ###################################################################
    select case( imode_files )

    case( IO_NOFILE )         ! Do nothing

    case( IO_NETCDF )         ! To be implemented

    case DEFAULT              ! One file with header per field
#ifdef USE_MPI
      if ( imode_precision_files == IO_TYPE_SINGLE ) then ! to be finished; here just as an idea
        subarray = IO_CREATE_SUBARRAY_XOZ( nx,ny,nz, MPI_REAL4 )
      else
        subarray = IO_CREATE_SUBARRAY_XOZ( nx,ny,nz, MPI_REAL8 )
      end if
#endif

      ! -------------------------------------------------------------------
      ! process header info
      isize = 0
      isize = isize+1; params(isize) = rtime
      isize = isize+1; params(isize) = visc ! inverse of reynolds
      if      ( iheader == IO_SCAL ) then
        isize = isize+1+1                     ! prepare space for schmidt and damkohler

      else if ( iheader == IO_FLOW ) then
        isize = isize+1; params(isize) = froude
        isize = isize+1; params(isize) = rossby
        if ( imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL ) then
          isize = isize+1; params(isize) = gama0
          isize = isize+1; params(isize) = prandtl
          isize = isize+1; params(isize) = mach
        end if
      end if

      if ( isize > isize_max ) then
        call TLAB_WRITE_ASCII(efile, 'IO_WRITE_FIELDS. Parameters array size error.')
        call TLAB_STOP(DNS_ERROR_ALLOC)
      end if

      header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

      ! -------------------------------------------------------------------
      ! write data
      do ifield = 1,nfield
        if ( iheader == IO_SCAL ) params(isize-1) = schmidt(ifield)   ! Scalar header
        if ( iheader == IO_SCAL ) params(isize  ) = damkohler(ifield) ! Scalar header
        write(name,'(I2)') ifield
        name=trim(adjustl(fname))//'.'//trim(adjustl(name))

        ! -------------------------------------------------------------------
        ! header
#ifdef USE_MPI
        if ( ims_pro == 0 ) then
#endif
#include "dns_open_file.h"
          call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,itime, params)
          close(LOC_UNIT_ID)
#ifdef USE_MPI
        end if
#endif

        ! -------------------------------------------------------------------
        ! field
        ! CALL IO_WRITE_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a(1,ifield),txc)
#ifdef USE_MPI
        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

        mpio_disp = header_offset*SIZEOFBYTE
        mpio_locsize = nx*ny*nz
        call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        if ( imode_precision_files == IO_TYPE_SINGLE ) then ! to be finished; here just as an idea
          call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
          s_wrk(:) = real(a(:,ifield),sp)
          call MPI_File_write_all(mpio_fh, s_wrk, mpio_locsize, MPI_REAL4, status, ims_err)
        else
          call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
          call MPI_File_write_all(mpio_fh, a(1,ifield), mpio_locsize, MPI_REAL8, status, ims_err)
        end if
        call MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
        if ( imode_precision_files == IO_TYPE_SINGLE ) then ! to be finished; here just as an idea
          s_wrk(:) = real(a(:,ifield),sp)
          write(LOC_UNIT_ID,POS=header_offset+1) s_wrk(:)
        else
          write(LOC_UNIT_ID,POS=header_offset+1) a(:,ifield)
        end if
        close(LOC_UNIT_ID)
#endif

      end do

    end select

    return
  end subroutine IO_WRITE_FIELDS

  !########################################################################
  !########################################################################
  subroutine IO_WRITE_FIELD_INT1(name, iheader, nx,ny,nz,nt, isize,params, a)

    character(LEN=*) name
    TINTEGER,   intent(IN   ) :: iheader, nx,ny,nz,nt, isize
    TREAL,      intent(IN   ) :: params(isize)
    integer(1), intent(IN   ) :: a(nx*ny*nz)

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

    line = 'Writing field '//trim(adjustl(name))//' of size'
    write(str,*) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(str))
    write(str,*) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))
    write(str,*) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))//'...'
    call TLAB_WRITE_ASCII(lfile, line)

    ! ###################################################################
#ifdef USE_MPI
    subarray = IO_CREATE_SUBARRAY_XOZ( nx,ny,nz, MPI_INTEGER1 )
#endif

    ! -------------------------------------------------------------------
    ! header
    if ( iheader > 0 ) then
      header_offset = 5*SIZEOFINT + isize*SIZEOFREAL
#ifdef USE_MPI
      if ( ims_pro == 0 ) then
#endif
#include "dns_open_file.h"
        call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total,ny_total,nz_total,nt, params)
        close(LOC_UNIT_ID)
#ifdef USE_MPI
      end if
#endif

    else
      header_offset = 0

    end if

    ! -------------------------------------------------------------------
    ! field
#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    mpio_disp = header_offset*SIZEOFBYTE
    mpio_locsize = nx*ny*nz
    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, ior(MPI_MODE_WRONLY,MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
    call MPI_File_set_view(mpio_fh, mpio_disp, MPI_INTEGER1, subarray, 'native', MPI_INFO_NULL, ims_err)
    call MPI_File_write_all(mpio_fh, a, mpio_locsize, MPI_INTEGER1, status, ims_err)
    call MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
    header_offset = header_offset +1
    write(LOC_UNIT_ID,POS=header_offset) a
    close(LOC_UNIT_ID)
#endif

    return
  end subroutine IO_WRITE_FIELD_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

  !########################################################################
  !########################################################################
  subroutine IO_READ_HEADER(unit, offset, nx,ny,nz,nt, params)
    TINTEGER unit, offset, nx,ny,nz,nt
    TREAL, dimension(*) :: params

    ! -------------------------------------------------------------------
    TINTEGER isize, nx_loc, ny_loc, nz_loc, nt_loc

    !########################################################################
    read(unit) offset, nx_loc, ny_loc, nz_loc, nt_loc

    isize = offset - 5*SIZEOFINT
    if ( isize > 0 .and. mod(isize,SIZEOFREAL) == 0 ) then
      isize = isize/SIZEOFREAL
      read(unit) params(1:isize)

    else
      call TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Header format incorrect.')
      call TLAB_STOP(DNS_ERROR_RECLEN)

    end if

    ! Check
    if ( nx /= nx_loc .or. ny /= ny_loc .or. nz /= nz_loc ) then
      close(unit)
      call TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Grid size mismatch.')
      call TLAB_STOP(DNS_ERROR_DIMGRID)
    end if

    if ( nt /= nt_loc ) then
      call TLAB_WRITE_ASCII(wfile, 'IO_READ_HEADER. ItNumber mismatch. Filename value ignored.')
      !     nt = nt_loc
    end if

    return
  end subroutine IO_READ_HEADER

  !########################################################################
  !########################################################################
  subroutine IO_WRITE_HEADER(unit, isize, nx,ny,nz,nt, params)
    TINTEGER unit, isize, nx,ny,nz,nt
    TREAL, dimension(isize) :: params

    ! -------------------------------------------------------------------
    TINTEGER offset

    !########################################################################
    offset = 5*SIZEOFINT + isize*SIZEOFREAL

    write(unit) offset, nx, ny, nz, nt
    write(unit) params(1:isize)

    return
  end subroutine IO_WRITE_HEADER

  !########################################################################
  !########################################################################
#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'

  subroutine IO_READ_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a, wrk)
#ifdef USE_MPI
    use TLAB_MPI_VARS, only : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLAB_MPI_PROCS
#endif
#include "integers.h"

    character(LEN=*) name
    TINTEGER, intent(IN   ) :: header_offset, nx,ny,nz
    TREAL,    intent(  OUT) :: a(nx*ny*nz)
    TREAL,    intent(INOUT) :: wrk(nx*ny*nz)

    target a, wrk

#ifdef USE_MPI
    integer(KIND=MPI_OFFSET_KIND) mpio_locoff
    TREAL, dimension(:), pointer :: p_read, p_write
    TINTEGER id, npage
#endif

    ! ###################################################################
#ifdef USE_MPI
    mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field

    if ( ims_npro_i > 1 ) then
      ! We always initialize types here. For the general field files, we could
      ! use TLAB_MPI_I_PARTIAL, but we use this routine for other files.
      call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for reading in IO_READ_FIELDS_SPLIT.')
      id = TLAB_MPI_I_AUX1
      npage = nz*ny
      call TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
           ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

      p_read => wrk

    else
      p_read => a

    end if

    mpio_locsize = nx*ny*nz
    mpio_locoff  = mpio_locsize         ! mpio_locoff might be of type larger than INT4
    mpio_locoff  = ims_pro*mpio_locoff  ! mpio_locoff might be of type larger than INT4
    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
    call MPI_FILE_READ_AT_ALL(mpio_fh, mpio_locoff, p_read, mpio_locsize, MPI_REAL8, status, ims_err)
    call MPI_FILE_CLOSE(mpio_fh, ims_err)

    if ( ims_npro_i > 1 ) then
      call TLAB_MPI_TRPB_I(p_read, a, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
    end if

    nullify(p_read)

#else
#include "dns_open_file.h"
    read(LOC_UNIT_ID,POS=header_offset +1) a
    close(LOC_UNIT_ID)

#endif

    return
  end subroutine IO_READ_FIELD_XPENCIL

#undef LOC_UNIT_ID
#undef LOC_STATUS

  !########################################################################
  !########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

  subroutine IO_WRITE_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a, wrk)
#ifdef USE_MPI
    use TLAB_MPI_VARS, only : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLAB_MPI_PROCS
#endif
#include "integers.h"

    character(LEN=*) name
    TINTEGER,  intent(IN   ) :: header_offset, nx,ny,nz
    TREAL,     intent(IN   ) :: a(nx*ny*nz)
    TREAL,     intent(INOUT) :: wrk(nx*ny*nz)

    target a, wrk

#ifdef USE_MPI
    integer(KIND=MPI_OFFSET_KIND) mpio_locoff
    TREAL, dimension(:), pointer :: p_read, p_write
    TINTEGER id, npage
#endif

    ! ###################################################################
#ifdef USE_MPI
    mpio_disp = header_offset*SIZEOFBYTE

    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

    if ( ims_npro_i > 1 ) then
      ! We always initialize types here. For the general field files, we could
      ! use TLAB_MPI_I_PARTIAL, but we use this routine for other files.
      call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for writing in IO_WRITE_FIELDS_SPLIT.')
      id = TLAB_MPI_I_AUX1
      npage = nz*ny
      call TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
           ims_size_i(id), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

      call TLAB_MPI_TRPF_I(a, wrk, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
      p_write => wrk

    else
      p_write => a

    end if

    mpio_locsize = nx*ny*nz
    mpio_locoff  = mpio_locsize         ! reclen might be of type larger than INT4
    mpio_locoff  = ims_pro*mpio_locoff  ! reclen might be of type larger than INT4
    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
    call MPI_FILE_WRITE_AT_ALL(mpio_fh, mpio_locoff, p_write, mpio_locsize, MPI_REAL8, status, ims_err)
    call MPI_FILE_CLOSE(mpio_fh, ims_err)
    nullify(p_write)

#else
#include "dns_open_file.h"
    write(LOC_UNIT_ID,POS=header_offset+1) a
    close(LOC_UNIT_ID)
#endif

    return
  end subroutine IO_WRITE_FIELD_XPENCIL

#undef LOC_UNIT_ID
#undef LOC_STATUS

end module IO_FIELDS
