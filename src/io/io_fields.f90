#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define USE_ACCESS_STREAM

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
    use TLAB_CONSTANTS, only: lfile, wfile, efile, wp, wi, sp, dp
    use TLAB_PROCS, only: TLAB_STOP, TLAB_WRITE_ASCII
    use TLAB_ARRAYS, only: wrk3d
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only: ims_err
    use TLAB_MPI_VARS, only: ims_pro, ims_npro_i, ims_npro_k, ims_pro_i, ims_pro_k
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_j, ims_offset_k
    use TLAB_MPI_PROCS, only: TLAB_MPI_PANIC
#endif
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    public :: IO_READ_FIELDS, IO_WRITE_FIELDS
    public :: IO_READ_FIELD_INT1, IO_WRITE_FIELD_INT1
    integer, parameter, public :: IO_SCAL = 1 ! Header of scalar field
    integer, parameter, public :: IO_FLOW = 2 ! Header of flow field

    type subarray_dt
        sequence
        integer :: precision = IO_TYPE_DOUBLE
#ifdef USE_MPI
        logical active, lpadding(3)
        integer communicator
        integer subarray
        integer(KIND=MPI_OFFSET_KIND) offset
#else
        integer offset
#endif
    end type subarray_dt

    public :: subarray_dt
#ifdef USE_MPI
    public :: IO_CREATE_SUBARRAY_XOY, IO_CREATE_SUBARRAY_XOZ, IO_CREATE_SUBARRAY_ZOY
#endif
    public :: IO_WRITE_SUBARRAY, IO_READ_SUBARRAY

    ! Global IO subarrays; to be transformed into local variables in corrresponding apps
    integer, parameter, public :: IO_SUBARRAY_VISUALS_XOY = 1
    integer, parameter, public :: IO_SUBARRAY_VISUALS_ZOY = 2
    integer, parameter, public :: IO_SUBARRAY_VISUALS_XOZ = 3
    integer, parameter, public :: IO_SUBARRAY_PLANES_XOY  = 4
    integer, parameter, public :: IO_SUBARRAY_PLANES_ZOY  = 5
    integer, parameter, public :: IO_SUBARRAY_PLANES_XOZ  = 6
    integer, parameter, public :: IO_SUBARRAY_BUFFER_ZOY  = 7
    integer, parameter, public :: IO_SUBARRAY_BUFFER_XOZ  = 8
    integer, parameter, public :: IO_SUBARRAY_SPECTRA_X   = 9
    integer, parameter, public :: IO_SUBARRAY_SPECTRA_Z   = 10
    integer, parameter, public :: IO_SUBARRAY_SPECTRA_XZ  = 11
    integer, parameter, public :: IO_SUBARRAY_ENVELOPES   = 12
    integer, parameter, public :: IO_SUBARRAY_AUX         = 13
    integer, parameter, public :: IO_SUBARRAY_SIZE        = 13
    type(subarray_dt), public :: io_aux(IO_SUBARRAY_SIZE)

    integer(wi) nx_total, ny_total, nz_total
    character(len=64) str, name
    character(len=128) line

    real(sp), pointer :: s_wrk(:) => null()

    integer, parameter :: sizeofreal = sizeof(1.0_wp)
    integer, parameter :: sizeofint = sizeof(1_wi)

#ifdef USE_MPI
    integer mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
    integer(KIND=MPI_OFFSET_KIND) mpio_disp
    integer subarray
#endif

contains

#ifdef USE_MPI
    !########################################################################
    !########################################################################
    function IO_CREATE_SUBARRAY_XOY(nx, ny, mpi_type) result(subarray)
        integer(wi), intent(in) :: nx, ny
        integer, intent(in) :: mpi_type

        integer :: subarray
        integer, parameter :: ndims = 2
        integer(wi) :: sizes(ndims), locsize(ndims), offset(ndims)

        sizes = [nx*ims_npro_i, ny]
        locsize = [nx, ny]
        offset = [nx*ims_pro_i, 0]

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, mpi_type, subarray, ims_err)
        call MPI_Type_commit(subarray, ims_err)

    end function IO_CREATE_SUBARRAY_XOY

    function IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, mpi_type) result(subarray)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: mpi_type

        integer :: subarray
        integer, parameter :: ndims = 3
        integer(wi) :: sizes(ndims), locsize(ndims), offset(ndims)

        sizes = [nx*ims_npro_i, ny, nz*ims_npro_k]
        locsize = [nx, ny, nz]
        offset = [nx*ims_pro_i, 0, nz*ims_pro_k]

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, mpi_type, subarray, ims_err)
        call MPI_Type_commit(subarray, ims_err)

    end function IO_CREATE_SUBARRAY_XOZ

    function IO_CREATE_SUBARRAY_ZOY(ny, nz, mpi_type) result(subarray)
        integer(wi), intent(in) :: ny, nz
        integer, intent(in) :: mpi_type

        integer :: subarray
        integer, parameter :: ndims = 2
        integer(wi) :: sizes(ndims), locsize(ndims), offset(ndims)

        sizes = [ny, nz*ims_npro_k]
        locsize = [ny, nz]
        offset = [0, nz*ims_pro_k]

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, mpi_type, subarray, ims_err)
        call MPI_Type_commit(subarray, ims_err)

    end function IO_CREATE_SUBARRAY_ZOY
#endif

    !########################################################################
    !########################################################################
#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'

    subroutine IO_READ_FIELDS(fname, iheader, nx, ny, nz, nfield, iread, a)
        use TLAB_VARS, only: imode_files, imode_precision_files
        use TLAB_VARS, only: itime, rtime, visc

        character(LEN=*) fname
        integer, intent(in) :: iheader         ! 1 for scalar header, 2 flow header
        integer, intent(in) :: nfield, iread   ! iread=0 reads all nfields, otherwise iread field
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(out) :: a(nx*ny*nz, *)

        ! -------------------------------------------------------------------
        integer(wi) header_offset
        integer ifield, iz

        integer, parameter :: isize_max = 20
        real(wp) params(isize_max)
        integer isize

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

        if (imode_precision_files == IO_TYPE_SINGLE) then
            line = 'Reading single precision field '//trim(adjustl(fname))//' of size'
            ! Pass memory address from double precision array to single precision array
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[nx*ny*nz])
        else
            line = 'Reading double precision field '//trim(adjustl(fname))//' of size'
        end if
        write (name, *) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(name))
        write (name, *) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))
        write (name, *) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))//'...'
        call TLAB_WRITE_ASCII(lfile, line)

        ! ###################################################################
        select case (imode_files)

        case (IO_NOFILE)         ! Do nothing
            a(:, 1:nfield) = 0.0_wp

        case (IO_NETCDF)         ! To be implemented

        case DEFAULT              ! One file with header per field
#ifdef USE_MPI
            if (imode_precision_files == IO_TYPE_SINGLE) then
                subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL4)
            else
                subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL8)
            end if
#endif

            ! -------------------------------------------------------------------
            ! read data
            iz = 0
            do ifield = 1, nfield
                if (iread == 0 .or. iread == ifield) then
                    iz = iz + 1
                    write (name, '(I2)') ifield
                    name = trim(adjustl(fname))//'.'//trim(adjustl(name))

                    ! -------------------------------------------------------------------
                    ! header
#ifdef USE_MPI
                    if (ims_pro == 0) then
#endif
#include "dns_open_file.h"
                        rewind (LOC_UNIT_ID)
                        call IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total, ny_total, nz_total, itime, params)
                        close (LOC_UNIT_ID)

#ifdef USE_MPI
                    end if
                    call MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)

                    ! -------------------------------------------------------------------
                    ! field
                    ! CALL IO_READ_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a(1,iz),wrk3d)
                    mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field
                    mpio_locsize = nx*ny*nz
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
                    if (imode_precision_files == IO_TYPE_SINGLE) then
                        call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
                        call MPI_File_read_all(mpio_fh, s_wrk, mpio_locsize, MPI_REAL4, status, ims_err)
                        a(:, iz) = real(s_wrk(:), dp)
                    else
                        call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
                        call MPI_File_read_all(mpio_fh, a(1, iz), mpio_locsize, MPI_REAL8, status, ims_err)
                    end if
                    call MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
                    if (imode_precision_files == IO_TYPE_SINGLE) then
                        read (LOC_UNIT_ID, POS=header_offset + 1) s_wrk(:)
                        a(:, iz) = real(s_wrk(:), dp)
                    else
                        read (LOC_UNIT_ID, POS=header_offset + 1) a(:, iz)
                    end if
                    close (LOC_UNIT_ID)

#endif

                end if
            end do

            ! -------------------------------------------------------------------
            ! process header info
            isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
            if (isize > isize_max) then
                call TLAB_WRITE_ASCII(efile, 'IO_READ_FIELDS. Parameters array size error')
                call TLAB_STOP(DNS_ERROR_ALLOC)
            end if

            rtime = params(1)
            if (iheader == IO_FLOW) then
                visc = params(2)
            end if
#ifdef USE_MPI
            call MPI_BCAST(itime, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
            call MPI_BCAST(rtime, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
            if (iheader == IO_FLOW) then
                call MPI_BCAST(visc, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
            end if
#endif

        end select

        return
    end subroutine IO_READ_FIELDS

    !########################################################################
    !########################################################################
    subroutine IO_READ_FIELD_INT1(name, iheader, nx, ny, nz, nt, isize, params, a)
        character(len=*) name
        integer, intent(in) :: iheader
        integer(wi), intent(in) :: nx, ny, nz, nt
        integer, intent(inout) :: isize
        real(wp), intent(inout) :: params(isize)
        integer(1), intent(out) :: a(nx*ny*nz)

        integer(wi) header_offset

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
        write (str, *) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(str))
        write (str, *) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))
        write (str, *) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))//'...'
        call TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
        subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_INTEGER1)
#endif

        ! -------------------------------------------------------------------
        ! header
#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
            header_offset = 0
#include "dns_open_file.h"
            rewind (LOC_UNIT_ID)
            if (iheader > 0) then
                call IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total, ny_total, nz_total, nt, params)
            end if
            close (LOC_UNIT_ID)

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
        header_offset = header_offset + 1
        read (LOC_UNIT_ID, POS=header_offset) a
        close (LOC_UNIT_ID)

#endif

        ! -------------------------------------------------------------------
        ! process header info
        isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
#ifdef USE_MPI
        if (isize > 0) call MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

        return
    end subroutine IO_READ_FIELD_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

    !########################################################################
    !########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

    subroutine IO_WRITE_FIELDS(fname, iheader, nx, ny, nz, nfield, a)
        use TLAB_VARS, only: imode_files, imode_precision_files, imode_eqns
        use TLAB_VARS, only: itime, rtime
        use TLAB_VARS, only: visc, froude, rossby, damkohler, prandtl, mach
        use TLAB_VARS, only: schmidt
        use THERMO_VARS, only: gama0

        character(len=*) fname
        integer, intent(in) :: iheader         ! Scalar or Flow headers
        integer, intent(in) :: nfield
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: a(nx*ny*nz, nfield)

        ! -------------------------------------------------------------------
        integer(wi) header_offset
        integer ifield

        integer, parameter :: isize_max = 20
        real(wp) params(isize_max)
        integer isize

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

        if (imode_precision_files == IO_TYPE_SINGLE) then
            line = 'Writing single precision field '//trim(adjustl(fname))//' of size'
            ! Pass memory address from double precision array to single precision array
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[nx*ny*nz])
        else
            line = 'Writing double precision field '//trim(adjustl(fname))//' of size'
        end if
        write (name, *) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(name))
        write (name, *) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))
        write (name, *) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))//'...'
        call TLAB_WRITE_ASCII(lfile, line)

        ! ###################################################################
        select case (imode_files)

        case (IO_NOFILE)         ! Do nothing

        case (IO_NETCDF)         ! To be implemented

        case DEFAULT              ! One file with header per field
#ifdef USE_MPI
            if (imode_precision_files == IO_TYPE_SINGLE) then
                subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL4)
            else
                subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL8)
            end if
#endif

            ! -------------------------------------------------------------------
            ! process header info
            isize = 0
            isize = isize + 1; params(isize) = rtime
            isize = isize + 1; params(isize) = visc ! inverse of reynolds
            if (iheader == IO_SCAL) then
                isize = isize + 1 + 1                     ! prepare space for schmidt and damkohler

            else if (iheader == IO_FLOW) then
                isize = isize + 1; params(isize) = froude
                isize = isize + 1; params(isize) = rossby
                if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) then
                    isize = isize + 1; params(isize) = gama0
                    isize = isize + 1; params(isize) = prandtl
                    isize = isize + 1; params(isize) = mach
                end if
            end if

            if (isize > isize_max) then
                call TLAB_WRITE_ASCII(efile, 'IO_WRITE_FIELDS. Parameters array size error.')
                call TLAB_STOP(DNS_ERROR_ALLOC)
            end if

            header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

            ! -------------------------------------------------------------------
            ! write data
            do ifield = 1, nfield
                if (iheader == IO_SCAL) params(isize - 1) = schmidt(ifield)   ! Scalar header
                if (iheader == IO_SCAL) params(isize) = damkohler(ifield) ! Scalar header
                write (name, '(I2)') ifield
                name = trim(adjustl(fname))//'.'//trim(adjustl(name))

                ! -------------------------------------------------------------------
                ! header
#ifdef USE_MPI
                if (ims_pro == 0) then
#endif
#include "dns_open_file.h"
                    call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total, ny_total, nz_total, itime, params)
                    close (LOC_UNIT_ID)
#ifdef USE_MPI
                end if
#endif

                ! -------------------------------------------------------------------
                ! field
                ! CALL IO_WRITE_FIELD_XPENCIL(name, header_offset, nx,ny,nz, a(1,ifield),wrk3d)
#ifdef USE_MPI
                call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

                mpio_disp = header_offset*SIZEOFBYTE
                mpio_locsize = nx*ny*nz
                call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
                if (imode_precision_files == IO_TYPE_SINGLE) then
                    call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
                    s_wrk(:) = real(a(:, ifield), sp)
                    call MPI_File_write_all(mpio_fh, s_wrk, mpio_locsize, MPI_REAL4, status, ims_err)
                else
                    call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
                    call MPI_File_write_all(mpio_fh, a(1, ifield), mpio_locsize, MPI_REAL8, status, ims_err)
                end if
                call MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
                if (imode_precision_files == IO_TYPE_SINGLE) then
                    s_wrk(:) = real(a(:, ifield), sp)
                    write (LOC_UNIT_ID, POS=header_offset + 1) s_wrk(:)
                else
                    write (LOC_UNIT_ID, POS=header_offset + 1) a(:, ifield)
                end if
                close (LOC_UNIT_ID)
#endif

            end do

        end select

        return
    end subroutine IO_WRITE_FIELDS

    !########################################################################
    !########################################################################
    subroutine IO_WRITE_FIELD_INT1(name, iheader, nx, ny, nz, nt, isize, params, a)
        character(len=*) name
        integer, intent(in) :: iheader, isize
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in) :: params(isize)
        integer(1), intent(in) :: a(nx*ny*nz)

        integer(wi) header_offset

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
        write (str, *) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(str))
        write (str, *) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))
        write (str, *) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))//'...'
        call TLAB_WRITE_ASCII(lfile, line)

        ! ###################################################################
#ifdef USE_MPI
        subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_INTEGER1)
#endif

        ! -------------------------------------------------------------------
        ! header
        if (iheader > 0) then
            header_offset = 5*SIZEOFINT + isize*SIZEOFREAL
#ifdef USE_MPI
            if (ims_pro == 0) then
#endif
#include "dns_open_file.h"
                call IO_WRITE_HEADER(LOC_UNIT_ID, isize, nx_total, ny_total, nz_total, nt, params)
                close (LOC_UNIT_ID)
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
        call MPI_FILE_OPEN(MPI_COMM_WORLD, name, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
        call MPI_File_set_view(mpio_fh, mpio_disp, MPI_INTEGER1, subarray, 'native', MPI_INFO_NULL, ims_err)
        call MPI_File_write_all(mpio_fh, a, mpio_locsize, MPI_INTEGER1, status, ims_err)
        call MPI_File_close(mpio_fh, ims_err)

#else
#include "dns_open_file.h"
        header_offset = header_offset + 1
        write (LOC_UNIT_ID, POS=header_offset) a
        close (LOC_UNIT_ID)
#endif

        return
    end subroutine IO_WRITE_FIELD_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

    !########################################################################
    !########################################################################
    subroutine IO_READ_HEADER(unit, offset, nx, ny, nz, nt, params)
        integer, intent(in) :: unit
        integer(wi), intent(out) :: offset
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(inout) :: params(:)

        ! -------------------------------------------------------------------
        integer isize
        integer(wi) nx_loc, ny_loc, nz_loc, nt_loc

        !########################################################################
        read (unit) offset, nx_loc, ny_loc, nz_loc, nt_loc

        ! Check
        if (nx /= nx_loc .or. ny /= ny_loc .or. nz /= nz_loc) then
            call TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Grid size mismatch.')
            call TLAB_STOP(DNS_ERROR_DIMGRID)
        end if

        if (nt /= nt_loc) then
            call TLAB_WRITE_ASCII(wfile, 'IO_READ_HEADER. ItNumber mismatch. Filename value ignored.')
            !     nt = nt_loc
        end if

        isize = offset - 5*SIZEOFINT
        if (isize > 0 .and. mod(isize, SIZEOFREAL) == 0) then
            isize = isize/SIZEOFREAL
            read (unit) params(1:isize)
        elseif (isize == 0) then
            continue ! no params to read; header format is correct
        else
            call TLAB_WRITE_ASCII(efile, 'IO_READ_HEADER. Header format incorrect.')
            call TLAB_STOP(DNS_ERROR_RECLEN)
        end if

        return

    end subroutine IO_READ_HEADER

    !########################################################################
    !########################################################################
    subroutine IO_WRITE_HEADER(unit, isize, nx, ny, nz, nt, params)
        integer, intent(in) :: unit, isize
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in) :: params(isize)

        ! -------------------------------------------------------------------
        integer(wi) offset

        !########################################################################
        offset = 5*SIZEOFINT + isize*SIZEOFREAL

        write (unit) offset, nx, ny, nz, nt

        if (isize > 0) then   ! do not write params to file if there are none
            write (unit) params(1:isize)
        end if

        return
    end subroutine IO_WRITE_HEADER

!########################################################################
!########################################################################
    subroutine IO_WRITE_SUBARRAY(aux, fname, varname, data, sizes)
        type(subarray_dt), intent(in) :: aux
        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
        character*32, intent(in) :: varname(sizes(5))
        real(wp), intent(in) :: data(sizes(1), sizes(5))

        ! -----------------------------------------------------------------------
        integer(wi) iv, isize

#ifdef USE_MPI
#else
        integer(wi) :: ioffset_local
#endif

        ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

        isize = (sizes(3) - sizes(2))/sizes(4) + 1

        if (aux%precision == IO_TYPE_SINGLE) then
            line = 'Writing single precision field'
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[isize])
        else
            line = 'Writing double precision field'
        end if

#ifdef USE_MPI
        if (aux%active) then
#endif

            do iv = 1, sizes(5)
                name = trim(adjustl(fname))
                if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))
                call TLAB_WRITE_ASCII(lfile, trim(adjustl(line))//' '//trim(adjustl(name))//'...')

#ifdef USE_MPI
                call MPI_File_open(aux%communicator, trim(adjustl(name)), &
                                   ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                if (aux%precision == IO_TYPE_SINGLE) then
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL4, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                    s_wrk(1:isize) = real(data(sizes(2):sizes(3):sizes(4), iv), sp)
                    call MPI_File_write_all(mpio_fh, s_wrk, isize, MPI_REAL4, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                else
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL8, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                    wrk3d(1:isize) = data(sizes(2):sizes(3):sizes(4), iv)
                    call MPI_File_write_all(mpio_fh, wrk3d, isize, MPI_REAL8, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                end if
                call MPI_File_close(mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)

#else
#include "dns_open_file.h"
                ioffset_local = aux%offset + 1
                if (aux%precision == IO_TYPE_SINGLE) then
                    write (LOC_UNIT_ID, POS=ioffset_local) real(data(sizes(2):sizes(3):sizes(4), iv), sp)
                else
                    write (LOC_UNIT_ID, POS=ioffset_local) data(sizes(2):sizes(3):sizes(4), iv)
                end if
                close (LOC_UNIT_ID)
#endif

            end do

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine IO_WRITE_SUBARRAY

!########################################################################
!########################################################################
    subroutine IO_READ_SUBARRAY(aux, fname, varname, data, sizes)
        type(subarray_dt), intent(in) :: aux
        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
        character*32, intent(in) :: varname(sizes(5))
        real(wp), intent(out) :: data(sizes(1), sizes(5))

        ! -----------------------------------------------------------------------
        integer(wi) iv, isize

#ifdef USE_MPI
#else
        integer(wi) :: ioffset_local
#endif

        ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

        isize = (sizes(3) - sizes(2))/sizes(4) + 1

        if (aux%precision == IO_TYPE_SINGLE) then
            line = 'Reading single precision field'
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[isize])
        else
            line = 'Reading double precision field'
        end if

#ifdef USE_MPI
        if (aux%active) then
#endif

            do iv = 1, sizes(5)
                name = trim(adjustl(fname))
                if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))
                call TLAB_WRITE_ASCII(lfile, trim(adjustl(line))//' '//trim(adjustl(name))//'...')

#ifdef USE_MPI
                call MPI_File_open(aux%communicator, trim(adjustl(name)), &
                                   MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                if (aux%precision == IO_TYPE_SINGLE) then
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL4, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                    call MPI_File_read_all(mpio_fh, s_wrk, isize, MPI_REAL4, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                else
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL8, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                    call MPI_File_read_all(mpio_fh, wrk3d, isize, MPI_REAL8, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
                end if
                call MPI_File_close(mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
#else
#include "dns_open_file.h"
                ioffset_local = aux%offset + 1
                if (aux%precision == IO_TYPE_SINGLE) then
                    read(LOC_UNIT_ID, POS=ioffset_local) s_wrk(1:isize)
                else
                    read(LOC_UNIT_ID, POS=ioffset_local) wrk3d(1:isize)
                end if
                close (LOC_UNIT_ID)
#endif
                if (aux%precision == IO_TYPE_SINGLE) then
                    data(sizes(2):sizes(3):sizes(4), iv) = real(s_wrk(1:isize), wp)
                else
                    data(sizes(2):sizes(3):sizes(4), iv) = wrk3d(1:isize)
                end if

            end do

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine IO_READ_SUBARRAY

!     !########################################################################
!     !########################################################################
! #define LOC_UNIT_ID 54
! #define LOC_STATUS 'old'

!     subroutine IO_READ_FIELD_XPENCIL(name, header_offset, nx, ny, nz, a, wrk)
! #ifdef USE_MPI
!         use TLAB_MPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
!         use TLAB_MPI_PROCS
! #endif

!         character(LEN=*) name
!         integer(wi), intent(in) :: header_offset, nx, ny, nz
!         real(wp), intent(out) :: a(nx*ny*nz)
!         real(wp), intent(inout) :: wrk(nx*ny*nz)

!         target a, wrk

! #ifdef USE_MPI
!         integer(KIND=MPI_OFFSET_KIND) mpio_locoff
!         real(wp), dimension(:), pointer :: p_read, p_write
!         integer(wi) id, npage
! #endif

!         ! ###################################################################
! #ifdef USE_MPI
!         mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field

!         if (ims_npro_i > 1) then
!             ! We always initialize types here. For the general field files, we could
!             ! use TLAB_MPI_I_PARTIAL, but we use this routine for other files.
!             call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for reading in IO_READ_FIELDS_SPLIT.')
!             id = TLAB_MPI_I_AUX1
!             npage = nz*ny
!             call TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
!                                  ims_size_i(id), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

!             p_read => wrk

!         else
!             p_read => a

!         end if

!         mpio_locsize = nx*ny*nz
!         mpio_locoff = mpio_locsize         ! mpio_locoff might be of type larger than INT4
!         mpio_locoff = ims_pro*mpio_locoff  ! mpio_locoff might be of type larger than INT4
!         call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
!         call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
!         call MPI_FILE_READ_AT_ALL(mpio_fh, mpio_locoff, p_read, mpio_locsize, MPI_REAL8, status, ims_err)
!         call MPI_FILE_CLOSE(mpio_fh, ims_err)

!         if (ims_npro_i > 1) then
!             call TLAB_MPI_TRPB_I(p_read, a, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
!         end if

!         nullify (p_read)

! #else
! #include "dns_open_file.h"
!         read (LOC_UNIT_ID, POS=header_offset + 1) a
!         close (LOC_UNIT_ID)

! #endif

!         return
!     end subroutine IO_READ_FIELD_XPENCIL

! #undef LOC_UNIT_ID
! #undef LOC_STATUS

!     !########################################################################
!     !########################################################################
! #define LOC_UNIT_ID 55
! #define LOC_STATUS 'unknown'

!     subroutine IO_WRITE_FIELD_XPENCIL(name, header_offset, nx, ny, nz, a, wrk)
! #ifdef USE_MPI
!         use TLAB_MPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
!         use TLAB_MPI_PROCS
! #endif

!         character(LEN=*) name
!         integer(wi), intent(in) :: header_offset, nx, ny, nz
!         real(wp), intent(in) :: a(nx*ny*nz)
!         real(wp), intent(inout) :: wrk(nx*ny*nz)

!         target a, wrk

! #ifdef USE_MPI
!         integer(KIND=MPI_OFFSET_KIND) mpio_locoff
!         real(wp), dimension(:), pointer :: p_read, p_write
!         integer(wi) id, npage
! #endif

!         ! ###################################################################
! #ifdef USE_MPI
!         mpio_disp = header_offset*SIZEOFBYTE

!         call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

!         if (ims_npro_i > 1) then
!             ! We always initialize types here. For the general field files, we could
!             ! use TLAB_MPI_I_PARTIAL, but we use this routine for other files.
!             call TLAB_WRITE_ASCII(lfile, 'Initializing MPI types for writing in IO_WRITE_FIELDS_SPLIT.')
!             id = TLAB_MPI_I_AUX1
!             npage = nz*ny
!             call TLAB_MPI_TYPE_I(ims_npro_i, nx, npage, i1, i1, i1, i1, &
!                                  ims_size_i(id), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))

!             call TLAB_MPI_TRPF_I(a, wrk, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
!             p_write => wrk

!         else
!             p_write => a

!         end if

!         mpio_locsize = nx*ny*nz
!         mpio_locoff = mpio_locsize         ! reclen might be of type larger than INT4
!         mpio_locoff = ims_pro*mpio_locoff  ! reclen might be of type larger than INT4
!         call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
!         call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
!         call MPI_FILE_WRITE_AT_ALL(mpio_fh, mpio_locoff, p_write, mpio_locsize, MPI_REAL8, status, ims_err)
!         call MPI_FILE_CLOSE(mpio_fh, ims_err)
!         nullify (p_write)

! #else
! #include "dns_open_file.h"
!         write (LOC_UNIT_ID, POS=header_offset + 1) a
!         close (LOC_UNIT_ID)
! #endif

!         return
!     end subroutine IO_WRITE_FIELD_XPENCIL

! #undef LOC_UNIT_ID
! #undef LOC_STATUS

end module IO_FIELDS
