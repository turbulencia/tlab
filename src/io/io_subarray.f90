#include "types.h"
#include "dns_error.h"
#include "dns_const_mpi.h"

#define SIZEOFBYTE 1

! The offset should be converted into logical and,
! when .TRUE., read the first integer in file as offset

subroutine IO_WRITE_SUBARRAY4(iflag_mode, fname, varname, data, sizes, work)

    use TLAB_TYPES, only: subarray_dt
    use TLAB_CONSTANTS, only: lfile
    use TLAB_VARS, only: io_aux
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS, only: TLAB_MPI_PANIC
#endif
    implicit none

    TINTEGER, intent(IN) :: iflag_mode
    character*(*), intent(IN) :: fname
    TINTEGER, intent(IN) :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
    character*32, dimension(sizes(5)), intent(IN) :: varname
    TREAL, dimension(sizes(1), sizes(5)), intent(IN) :: data
    real(4), dimension(sizes(1)), intent(INOUT) :: work

    ! -----------------------------------------------------------------------
    TINTEGER iv, isize
    character*64 name

#ifdef USE_MPI
    TINTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh, ims_err
#else
    TINTEGER :: ioffset_local
#endif

    ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

    isize = (sizes(3) - sizes(2))/sizes(4) + 1

#ifdef USE_MPI
    if (io_aux(iflag_mode)%active) then
#endif

        do iv = 1, sizes(5)
            name = trim(adjustl(fname))
            if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))

            call TLAB_WRITE_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')

            work(1:isize) = SNGL(data(sizes(2):sizes(3):sizes(4), iv))

#ifdef USE_MPI
            call MPI_File_open(io_aux(iflag_mode)%communicator, trim(adjustl(name)), &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_set_view(mpio_fh, io_aux(iflag_mode)%offset, MPI_REAL4, io_aux(iflag_mode)%subarray, &
                                   'native', MPI_INFO_NULL, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_write_all(mpio_fh, work, isize, MPI_REAL4, mpio_status, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_close(mpio_fh, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)

#else
#include "dns_open_file.h"
            ioffset_local = io_aux(iflag_mode)%offset + 1
            write (LOC_UNIT_ID, POS=ioffset_local) work(1:isize)
            close (LOC_UNIT_ID)

#endif

        end do

#ifdef USE_MPI
    end if
#endif

    return
end subroutine IO_WRITE_SUBARRAY4

!########################################################################
!########################################################################
subroutine IO_READ_SUBARRAY8(iflag_mode, fname, varname, data, sizes, work)

    use TLAB_TYPES, only: subarray_dt
    use TLAB_CONSTANTS, only: lfile, efile
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
#endif
    use TLAB_VARS, only: io_aux
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only: ims_err
    use TLAB_MPI_PROCS, only: TLAB_MPI_PANIC
#endif

    implicit none

    TINTEGER, intent(IN) :: iflag_mode
    character*(*), intent(IN) :: fname
    TINTEGER, intent(IN) :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
    character*32, dimension(sizes(5)), intent(IN) :: varname
    TREAL, dimension(sizes(1), sizes(5)), intent(OUT) :: data
    TREAL, dimension(sizes(1)), intent(INOUT) :: work

    ! -----------------------------------------------------------------------
    TINTEGER iv, isize
    character*64 name

#ifdef USE_MPI
    integer :: mpio_status(MPI_STATUS_SIZE), mpio_fh
#else
    TINTEGER :: ioffset_local
#endif

    ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING IO_READ_SUBARRAY8')
#endif
    isize = (sizes(3) - sizes(2))/sizes(4) + 1

    do iv = 1, sizes(5)
        name = trim(adjustl(fname))
        if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))
        call TLAB_WRITE_ASCII(lfile, 'Reading field '//trim(adjustl(name))//'...')
    end do

#ifdef USE_MPI
    if (io_aux(iflag_mode)%active) then
#endif

        do iv = 1, sizes(5)
            name = trim(adjustl(fname))
            if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))
#ifdef USE_MPI
            call MPI_File_open(io_aux(iflag_mode)%communicator, trim(adjustl(name)), &
                               MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_set_view(mpio_fh, io_aux(iflag_mode)%offset, MPI_REAL8, io_aux(iflag_mode)%subarray, &
                                   'native', MPI_INFO_NULL, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_read_all(mpio_fh, work, isize, MPI_REAL8, mpio_status, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_close(mpio_fh, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
#else
#include "dns_open_file.h"
            ioffset_local = io_aux(iflag_mode)%offset + 1
            read (LOC_UNIT_ID, POS=ioffset_local) work(1:isize)
            close (LOC_UNIT_ID)
#endif
            data(sizes(2):sizes(3):sizes(4), iv) = work(1:isize)
        end do

#ifdef USE_MPI
    end if
#endif
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING IO_READ_SUBARRAY8')
#endif

    return
end subroutine IO_READ_SUBARRAY8

!########################################################################
!########################################################################
subroutine IO_WRITE_SUBARRAY8(iflag_mode, fname, varname, data, sizes, work)

    use TLAB_TYPES, only: subarray_dt
    use TLAB_CONSTANTS, only: lfile
    use TLAB_VARS, only: io_aux
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_PROCS, only: TLAB_MPI_PANIC
#endif

    implicit none

    TINTEGER, intent(IN) :: iflag_mode
    character*(*), intent(IN) :: fname
    TINTEGER, intent(IN) :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
    character*32, dimension(sizes(5)), intent(IN) :: varname
    TREAL, dimension(sizes(1), sizes(5)), intent(IN) :: data
    TREAL, dimension(sizes(1)), intent(INOUT) :: work

    ! -----------------------------------------------------------------------
    TINTEGER iv, isize
    character*64 name

#ifdef USE_MPI
    TINTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh, ims_err
#else
    TINTEGER :: ioffset_local
#endif

    ! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

    isize = (sizes(3) - sizes(2))/sizes(4) + 1

#ifdef USE_MPI
    if (io_aux(iflag_mode)%active) then
#endif

        do iv = 1, sizes(5)
            name = trim(adjustl(fname))
            if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))

            call TLAB_WRITE_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')

            work(1:isize) = data(sizes(2):sizes(3):sizes(4), iv)

#ifdef USE_MPI
            call MPI_File_open(io_aux(iflag_mode)%communicator, trim(adjustl(name)), &
                               ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_set_view(mpio_fh, io_aux(iflag_mode)%offset, MPI_REAL8, io_aux(iflag_mode)%subarray, &
                                   'native', MPI_INFO_NULL, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_write_all(mpio_fh, work, isize, MPI_REAL8, mpio_status, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)
            call MPI_File_close(mpio_fh, ims_err)
            if (ims_err /= MPI_SUCCESS) call TLAB_MPI_PANIC(__FILE__, ims_err)

#else
#include "dns_open_file.h"
            ioffset_local = io_aux(iflag_mode)%offset + 1
            write (LOC_UNIT_ID, POS=ioffset_local) work(1:isize)
            close (LOC_UNIT_ID)

#endif

        end do

#ifdef USE_MPI
    end if
#endif

    return
end subroutine IO_WRITE_SUBARRAY8
