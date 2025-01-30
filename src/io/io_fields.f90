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
    use TLab_Constants, only: lfile, wfile, efile, wp, wi, sp, dp, sizeofint, sizeofreal, MAX_PARS, MAX_VARS
    use TLab_WorkFlow, only: TLab_Stop, TLab_Write_ASCII
    use TLab_Arrays, only: wrk3d
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS
    use TLabMPI_PROCS, only: TLabMPI_Panic
#endif
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    integer, public :: io_fileformat                ! files format
    integer, public :: io_datatype                  ! single or double precision

    type :: io_header                               ! header information
        sequence
        integer size
        real(wp) params(MAX_PARS)
    end type
    type(io_header), public :: io_header_q(1), io_header_s(MAX_VARS)

    public :: IO_READ_FIELDS, IO_WRITE_FIELDS
    public :: IO_READ_FIELD_INT1, IO_WRITE_FIELD_INT1

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
    public :: TLabMPI_WRITE_PE0_SINGLE
#endif
    public :: IO_WRITE_SUBARRAY, IO_READ_SUBARRAY

    ! Global IO subarrays; to be transformed into local variables in corrresponding apps
    integer, parameter, public :: IO_SUBARRAY_VISUALS_XOY = 1
    integer, parameter, public :: IO_SUBARRAY_VISUALS_ZOY = 2
    integer, parameter, public :: IO_SUBARRAY_VISUALS_XOZ = 3
    integer, parameter, public :: IO_SUBARRAY_PLANES_XOY = 4
    integer, parameter, public :: IO_SUBARRAY_PLANES_ZOY = 5
    integer, parameter, public :: IO_SUBARRAY_PLANES_XOZ = 6
    integer, parameter, public :: IO_SUBARRAY_BUFFER_ZOY = 7
    integer, parameter, public :: IO_SUBARRAY_BUFFER_XOZ = 8
    integer, parameter, public :: IO_SUBARRAY_SPECTRA_X = 9
    integer, parameter, public :: IO_SUBARRAY_SPECTRA_Z = 10
    integer, parameter, public :: IO_SUBARRAY_SPECTRA_XZ = 11
    integer, parameter, public :: IO_SUBARRAY_ENVELOPES = 12
    integer, parameter, public :: IO_SUBARRAY_AUX = 13
    integer, parameter, public :: IO_SUBARRAY_SIZE = 13
    type(subarray_dt), public :: io_aux(IO_SUBARRAY_SIZE)

    integer(wi) nx_total, ny_total, nz_total
    character(len=64) str, name
    character(len=128) line

    real(sp), pointer :: s_wrk(:) => null()

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

    subroutine IO_READ_FIELDS(fname, nx, ny, nz, nt, nfield, iread, a, params)
        character(LEN=*) fname
        integer, intent(in) :: nfield, iread   ! iread=0 reads all nfields, otherwise iread field
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(out) :: a(nx*ny*nz, *)
        real(wp), intent(inout) :: params(:)

        ! -------------------------------------------------------------------
        integer(wi) header_offset, isize
        integer ifield, iz

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

        if (io_datatype == IO_TYPE_SINGLE) then
            line = 'Reading single precision field '//trim(adjustl(fname))//' of size'
            ! Pass memory address from double precision array to single precision array
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[nx*ny*nz])
        else
            line = 'Reading double precision field '//trim(adjustl(fname))//' of size'
        end if
        write (name, *) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(name))
        write (name, *) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))
        write (name, *) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))//'...'
        call TLab_Write_ASCII(lfile, line)

        ! ###################################################################
        select case (io_fileformat)

        case (IO_NOFILE)         ! Do nothing
            a(:, 1:nfield) = 0.0_wp

        case (IO_NETCDF)         ! To be implemented

        case DEFAULT              ! One file with header per field
#ifdef USE_MPI
            if (io_datatype == IO_TYPE_SINGLE) then
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
                        call IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total, ny_total, nz_total, nt, params)
                        close (LOC_UNIT_ID)

#ifdef USE_MPI
                    end if
                    call MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)

                    ! -------------------------------------------------------------------
                    ! field
                    mpio_disp = header_offset*SIZEOFBYTE ! Displacement to start of field
                    mpio_locsize = nx*ny*nz
                    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
                    if (io_datatype == IO_TYPE_SINGLE) then
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
                    if (io_datatype == IO_TYPE_SINGLE) then
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
#ifdef USE_MPI
            if (isize > 0) call MPI_BCAST(params, size(params), MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

        end select

        return
    end subroutine IO_READ_FIELDS

    !########################################################################
    !########################################################################
    subroutine IO_READ_FIELD_INT1(name, nx, ny, nz, nt, a, params)
        character(len=*) name
        integer(wi), intent(in) :: nx, ny, nz, nt
        integer(1), intent(out) :: a(nx*ny*nz)
        real(wp), intent(inout) :: params(:)

        ! -------------------------------------------------------------------
        integer(wi) header_offset, isize

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
        call TLab_Write_ASCII(lfile, line)

#ifdef USE_MPI
        subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_INTEGER1)
#endif

        ! -------------------------------------------------------------------
        ! header
#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
#include "dns_open_file.h"
            rewind (LOC_UNIT_ID)
            call IO_READ_HEADER(LOC_UNIT_ID, header_offset, nx_total, ny_total, nz_total, nt, params)
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
        if (isize > 0) call MPI_BCAST(params, size(params), MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

        return
    end subroutine IO_READ_FIELD_INT1

#undef LOC_UNIT_ID
#undef LOC_STATUS

    !########################################################################
    !########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

    subroutine IO_WRITE_FIELDS(fname, nx, ny, nz, nt, nfield, a, locHeader)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: nfield
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in) :: a(nx*ny*nz, nfield)
        type(io_header), intent(in), optional :: locHeader(:)

        ! -------------------------------------------------------------------
        integer(wi) header_offset
        integer ifield, ih

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

        if (io_datatype == IO_TYPE_SINGLE) then
            line = 'Writing single precision field '//trim(adjustl(fname))//' of size'
            ! Pass memory address from double precision array to single precision array
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[nx*ny*nz])
        else
            line = 'Writing double precision field '//trim(adjustl(fname))//' of size'
        end if
        write (name, *) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(name))
        write (name, *) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))
        write (name, *) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(name))//'...'
        call TLab_Write_ASCII(lfile, line)

        ! ###################################################################
        select case (io_fileformat)

        case (IO_NOFILE)         ! Do nothing

        case (IO_NETCDF)         ! To be implemented

        case DEFAULT              ! One file with header per field
#ifdef USE_MPI
            if (io_datatype == IO_TYPE_SINGLE) then
                subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL4)
            else
                subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL8)
            end if
#endif

            do ifield = 1, nfield
                write (name, '(I2)') ifield
                name = trim(adjustl(fname))//'.'//trim(adjustl(name))

                ! -------------------------------------------------------------------
                ! header
                header_offset = 5*SIZEOFINT
                if (present(locHeader)) then
                    ih = min(size(locHeader), ifield)      ! use always the 1 if you enter with one, or the ifield
                    header_offset = header_offset + locHeader(ih)%size*SIZEOFREAL
                end if

#ifdef USE_MPI
                if (ims_pro == 0) then
#endif
#include "dns_open_file.h"
                    if (present(locHeader)) then
                        call IO_WRITE_HEADER(LOC_UNIT_ID, nx_total, ny_total, nz_total, nt, locHeader(ih)%params(1:locHeader(ih)%size))
                    else
                        call IO_WRITE_HEADER(LOC_UNIT_ID, nx_total, ny_total, nz_total, nt)
                    end if
                    close (LOC_UNIT_ID)
#ifdef USE_MPI
                end if
#endif

                ! -------------------------------------------------------------------
                ! field
#ifdef USE_MPI
                call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

                mpio_disp = header_offset*SIZEOFBYTE
                mpio_locsize = nx*ny*nz
                call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
                if (io_datatype == IO_TYPE_SINGLE) then
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
                if (io_datatype == IO_TYPE_SINGLE) then
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
    subroutine IO_WRITE_FIELD_INT1(name, nx, ny, nz, nt, a, params)
        character(len=*) name
        integer(wi), intent(in) :: nx, ny, nz, nt
        integer(1), intent(in) :: a(nx*ny*nz)
        real(wp), intent(in), optional :: params(:)

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
        call TLab_Write_ASCII(lfile, line)

        ! ###################################################################
#ifdef USE_MPI
        subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_INTEGER1)
#endif

        ! -------------------------------------------------------------------
        ! header
        header_offset = 5*SIZEOFINT
        if (present(params)) then
            header_offset = header_offset + size(params)*SIZEOFREAL
        end if

#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
#include "dns_open_file.h"
            call IO_WRITE_HEADER(LOC_UNIT_ID, nx_total, ny_total, nz_total, nt, params(:))
            close (LOC_UNIT_ID)
#ifdef USE_MPI
        end if
#endif

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
            call TLab_Write_ASCII(efile, 'IO_READ_HEADER. Grid size mismatch.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        if (nt /= nt_loc) then
            call TLab_Write_ASCII(wfile, 'IO_READ_HEADER. ItNumber mismatch. Filename value ignored.')
            ! nt = nt_loc
        end if

        isize = offset - 5*SIZEOFINT
        if (isize >= 0 .and. mod(isize, SIZEOFREAL) == 0) then
            ! isize = isize/SIZEOFREAL
            read (unit) params(:)
            ! elseif (isize == 0) then
            !     continue ! no params to read; header format is correct
        else
            call TLab_Write_ASCII(efile, 'IO_READ_HEADER. Header format incorrect.')
            call TLab_Stop(DNS_ERROR_RECLEN)
        end if

        return

    end subroutine IO_READ_HEADER

    !########################################################################
    !########################################################################
    subroutine IO_WRITE_HEADER(unit, nx, ny, nz, nt, params)
        integer, intent(in) :: unit
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in), optional :: params(:)

        ! -------------------------------------------------------------------
        integer(wi) offset

        !########################################################################
        offset = 5*SIZEOFINT
        if (present(params)) then
            offset = offset + size(params)*SIZEOFREAL
        end if

        write (unit) offset, nx, ny, nz, nt

        if (present(params)) then
            write (unit) params(:)
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
                call TLab_Write_ASCII(lfile, trim(adjustl(line))//' '//trim(adjustl(name))//'...')

#ifdef USE_MPI
                call MPI_File_open(aux%communicator, trim(adjustl(name)), &
                                   ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                if (aux%precision == IO_TYPE_SINGLE) then
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL4, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                    s_wrk(1:isize) = real(data(sizes(2):sizes(3):sizes(4), iv), sp)
                    call MPI_File_write_all(mpio_fh, s_wrk, isize, MPI_REAL4, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                else
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL8, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                    wrk3d(1:isize) = data(sizes(2):sizes(3):sizes(4), iv)
                    call MPI_File_write_all(mpio_fh, wrk3d, isize, MPI_REAL8, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                end if
                call MPI_File_close(mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)

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
                call TLab_Write_ASCII(lfile, trim(adjustl(line))//' '//trim(adjustl(name))//'...')

#ifdef USE_MPI
                call MPI_File_open(aux%communicator, trim(adjustl(name)), &
                                   MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                if (aux%precision == IO_TYPE_SINGLE) then
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL4, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                    call MPI_File_read_all(mpio_fh, s_wrk, isize, MPI_REAL4, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                else
                    call MPI_File_set_view(mpio_fh, aux%offset, MPI_REAL8, aux%subarray, 'native', MPI_INFO_NULL, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                    call MPI_File_read_all(mpio_fh, wrk3d, isize, MPI_REAL8, status, ims_err)
                    if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
                end if
                call MPI_File_close(mpio_fh, ims_err)
                if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
#else
#include "dns_open_file.h"
                ioffset_local = aux%offset + 1
                if (aux%precision == IO_TYPE_SINGLE) then
                    read (LOC_UNIT_ID, POS=ioffset_local) s_wrk(1:isize)
                else
                    read (LOC_UNIT_ID, POS=ioffset_local) wrk3d(1:isize)
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

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine TLabMPI_WRITE_PE0_SINGLE(iunit, nx, ny, nz, subdomain, u, tmp1, tmp2)
        use TLabMPI_Transpose

        integer ims_tag
        integer(wi) iunit, nx, ny, nz, subdomain(6)
        real(wp), dimension(nx*ny*nz), target :: u, tmp1
        real(wp), dimension(nx*ims_npro_i, *), target :: tmp2

        ! -------------------------------------------------------------------
        integer(wi) nx_total, ny_total, nz_total
        integer(wi) nx_min, nx_max, ny_min, ny_max, nz_min, nz_max
        integer(wi) nyz

        integer(wi) ip_i, ip_k, joffset_loc, koffset_loc, id
        integer(wi) i, jk, j_loc, k_loc
        integer mpio_size, mpio_ip
        integer status(MPI_STATUS_SIZE)

        real(wp), dimension(:), pointer :: p_org

        ! ###################################################################
        nx_total = nx*ims_npro_i
        ny_total = ny
        nz_total = nz*ims_npro_k

        nx_min = subdomain(1); nx_max = subdomain(2)
        ny_min = subdomain(3); ny_max = subdomain(4)
        nz_min = subdomain(5); nz_max = subdomain(6)

        koffset_loc = 0
        joffset_loc = 0

        ! id = TLAB_MPI_TRP_I_PARTIAL

        ! -------------------------------------------------------------------
        ! Transposing along Ox
        ! -------------------------------------------------------------------
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Forward(u, tmp1, ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL))
            p_org => tmp1
            ! nyz = ims_size_i(id)
            nyz = ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL)%nlines
        else
            p_org => u
            nyz = ny*nz
        end if
        mpio_size = nyz*nx_total

        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

        ! -------------------------------------------------------------------
        ! Passing all data through PE#0
        ! -------------------------------------------------------------------
        if (ims_pro == 0) then

            do ip_k = 1, ims_npro_k
                koffset_loc = nz*(ip_k - 1)

                do ip_i = 1, ims_npro_i
                    joffset_loc = nyz*(ip_i - 1) ! Remember that data is Ox-transposed

                    mpio_ip = ims_npro_i*(ip_k - 1) + ip_i - 1
                    if (mpio_ip == 0) then
                        tmp2(1:mpio_size, 1) = p_org(1:mpio_size)
                    else
                        call MPI_RECV(tmp2, mpio_size, MPI_REAL8, mpio_ip, ims_tag, MPI_COMM_WORLD, status, ims_err)
                    end if

                    do jk = 1, nyz
                        j_loc = mod((jk - 1 + joffset_loc), ny_total) + 1
                        k_loc = ((jk - 1 + joffset_loc)/ny_total) + 1 + koffset_loc

                        if ((j_loc >= ny_min) .and. (j_loc <= ny_max) .and. &
                            (k_loc >= nz_min) .and. (k_loc <= nz_max)) then
                            write (iunit) (SNGL(tmp2(i, jk)), i=nx_min, nx_max)
                        end if

                    end do

                end do
            end do

        else
            call MPI_SEND(p_org, mpio_size, MPI_REAL8, 0, ims_tag, MPI_COMM_WORLD, ims_err)
        end if

        return
    end subroutine TLabMPI_WRITE_PE0_SINGLE

#endif
end module IO_FIELDS
