#include "types.h"
#include "dns_const.h"

#define LOC_UNIT_ID i55
#define LOC_STATUS 'unknown'

subroutine IO_WRITE_VISUALS(fname, iformat, nx, ny, nz, nfield, subdomain, field, txc)

    use TLAB_TYPES, only: subarray_dt
    use TLAB_VARS, only: g, isize_txc_field, io_aux
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only: ims_pro
    use TLAB_MPI_PROCS
#endif
    use IO_FIELDS

    implicit none

#include "integers.h"

    TINTEGER iformat, nx, ny, nz, nfield, subdomain(6)
    TREAL, dimension(isize_txc_field, nfield) :: field
    TREAL, dimension(nx*ny*nz, 2) :: txc
    character*(*) fname

    ! -------------------------------------------------------------------
    TINTEGER sizes(5), nx_aux, ny_aux, nz_aux, ifield, i
    character*32 varname(16), name
    TINTEGER iflag_mode

    ! ###################################################################
    sizes(5) = nfield

    nx_aux = subdomain(2) - subdomain(1) + 1
    ny_aux = subdomain(4) - subdomain(3) + 1
    nz_aux = subdomain(6) - subdomain(5) + 1

    iflag_mode = 0 ! default
    sizes(1) = isize_txc_field     ! array size
    sizes(2) = 1                   ! lower bound
    if (subdomain(2) - subdomain(1) + 1 == g(1)%size .and. &
        subdomain(6) - subdomain(5) + 1 == 1) then! xOy plane
        iflag_mode = IO_SUBARRAY_VISUALS_XOY
        sizes(3) = ny_aux*nx     ! upper bound
        sizes(4) = 1              ! stride

    else if (subdomain(6) - subdomain(5) + 1 == g(3)%size .and. &
             subdomain(2) - subdomain(1) + 1 == 1) then! zOy plane
        iflag_mode = IO_SUBARRAY_VISUALS_ZOY
        sizes(3) = ny_aux*nx*nz ! upper bound
        sizes(4) = nx             ! stride

    else if (subdomain(2) - subdomain(1) + 1 == g(1)%size .and. &
             subdomain(6) - subdomain(5) + 1 == g(3)%size) then! xOz blocks
        iflag_mode = IO_SUBARRAY_VISUALS_XOZ
        sizes(3) = ny_aux*nx*nz ! upper bound
        sizes(4) = 1              ! stride

    end if

    ! ###################################################################
    ! We need to rearrange the arrays here, because
    ! IO_WRITE_FIELDS and ENSIGHT_FIELD expects field to be aligned by nx*ny*nz
    ! (instead of isize_txc_field)
    if (iformat < 2. .and. &
        nfield > 1 .and. isize_txc_field > nx*ny*nz) then
        do ifield = 2, nfield
            do i = 1, nx*ny*nz
                field((ifield - 1)*nx*ny*nz + i, 1) = field(i, ifield)
            end do
        end do
    end if

    ! ###################################################################

    if (iformat == 0) then ! standard scalar format

        call IO_WRITE_FIELDS(fname, IO_SCAL, nx, ny, nz, nfield, field)

        ! -------------------------------------------------------------------
    else if (iformat == 1) then  ! ensight; to be removed
        call ENSIGHT_FIELD(fname, i1, nx, ny, nz, nfield, subdomain, field, txc)

        ! -------------------------------------------------------------------
    else if (iformat == 2 .and. iflag_mode > 0) then  ! single precision, using MPI_IO
        if (ny_aux /= ny) then
            do ifield = 1, nfield
                call REDUCE_BLOCK_INPLACE(nx, ny, nz, i1, subdomain(3), i1, nx, ny_aux, nz, field(1, ifield), txc)
            end do
        end if

        varname = ''
        if (nfield > 1) then
            do ifield = 1, nfield; write (varname(ifield), *) ifield; varname(ifield) = trim(adjustl(varname(ifield)))
            end do
        end if
!        call IO_WRITE_SUBARRAY4(iflag_mode, fname, varname, field, sizes, txc)
        call IO_WRITE_SUBARRAY(io_aux(iflag_mode), fname, varname, field, sizes)

        ! -------------------------------------------------------------------
    else                                                     ! single precision, through PE0
        do ifield = 1, nfield
            if (nfield > 1) then
                write (name, *) ifield; name = trim(adjustl(fname))//'.'//trim(adjustl(name))
            else
                name = fname
            end if

#ifdef USE_MPI
            if (ims_pro == 0) then
#endif
#include "dns_open_file.h"

#ifdef USE_MPI
            end if
            call TLAB_MPI_WRITE_PE0_SINGLE(LOC_UNIT_ID, nx, ny, nz, subdomain, field(1, ifield), txc(1, 1), txc(1, 2))
            if (ims_pro == 0) then
#else
      call REDUCE_BLOCK_INPLACE(nx, ny, nz, subdomain(1), subdomain(3), subdomain(5), nx_aux, ny_aux, nz_aux, field(1, ifield), txc)
                write (LOC_UNIT_ID) SNGL(field(1:nx_aux*ny_aux*nz_aux, ifield))
#endif
                close (LOC_UNIT_ID)
#ifdef USE_MPI
            end if
#endif
        end do

    end if

    return
end subroutine IO_WRITE_VISUALS

!########################################################################
! Writing data in Ensight Gold Variable File Format
!########################################################################
subroutine ENSIGHT_FIELD(name, iheader, nx, ny, nz, nfield, subdomain, field, tmp_mpi)

#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_pro
    use TLAB_MPI_PROCS
#endif

    implicit none

#include "integers.h"

    character*(*) name
    TINTEGER, intent(IN) :: iheader ! 0 no header; 1 header

    TINTEGER, intent(IN) :: nx, ny, nz, nfield, subdomain(6)
    TREAL, dimension(nx, ny, nz, nfield) :: field
    TREAL, dimension(nx*ny*nz, 2) :: tmp_mpi

    ! -------------------------------------------------------------------
    character*80 line
#ifdef USE_MPI
    TINTEGER ifield
#else
    TINTEGER i, j, k, ifield
#endif

    ! ###################################################################
    ! Header in Ensight Gold Variable File Format
    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
#include "dns_open_file.h"

        if (iheader == 1) then
            line = 'description line              '
            write (LOC_UNIT_ID) line
            line = 'part                          '
            write (LOC_UNIT_ID) line
            write (LOC_UNIT_ID) i1
            line = 'block                         '
            write (LOC_UNIT_ID) line
        end if

#ifdef USE_MPI
    end if
#endif

    ! ###################################################################
    ! Body
    ! ###################################################################
    ! -------------------------------------------------------------------
    ! parallel
    ! -------------------------------------------------------------------
#ifdef USE_MPI
    do ifield = 1, nfield
        call TLAB_MPI_WRITE_PE0_SINGLE(LOC_UNIT_ID, nx, ny, nz, subdomain, field, tmp_mpi(1, 1), tmp_mpi(1, 2))
    end do

    ! -------------------------------------------------------------------
    ! serial
    ! -------------------------------------------------------------------
#else
    do ifield = 1, nfield
        do k = subdomain(5), subdomain(6)
            do j = subdomain(3), subdomain(4)
                write (LOC_UNIT_ID) (SNGL(field(i, j, k, ifield)), i=subdomain(1), subdomain(2))
            end do
        end do
    end do

#endif

    ! ###################################################################
    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif
        close (LOC_UNIT_ID)
#ifdef USE_MPI
    end if
#endif

    return
end subroutine ENSIGHT_FIELD

!########################################################################
! Writing data in Ensight Gold Geometry File Format
! Note that record-length information has to be avoided
!########################################################################
subroutine ENSIGHT_GRID(name, nx, ny, nz, subdomain, x, y, z)

    implicit none

#include "integers.h"

    TINTEGER nx, ny, nz, subdomain(6)
    TREAL x(nx), y(ny), z(nz)
    character*(*) name

    ! -------------------------------------------------------------------
    character*80 line
    TINTEGER ij

    ! ###################################################################
#include "dns_open_file.h"

    line = 'Fortran Binary                '
    write (LOC_UNIT_ID) line

    line = 'description line 1            '
    write (LOC_UNIT_ID) line
    line = 'description line 2            '
    write (LOC_UNIT_ID) line
    line = 'node id off                   '
    write (LOC_UNIT_ID) line
    line = 'element id off                '
    write (LOC_UNIT_ID) line

    line = 'part                          '
    write (LOC_UNIT_ID) line
    write (LOC_UNIT_ID) i1
    line = 'description line              '
    write (LOC_UNIT_ID) line
    line = 'block rectilinear             '
    write (LOC_UNIT_ID) line
    write (LOC_UNIT_ID) subdomain(2) - subdomain(1) + 1, subdomain(4) - subdomain(3) + 1, subdomain(6) - subdomain(5) + 1
    write (LOC_UNIT_ID) (SNGL(x(ij)), ij=subdomain(1), subdomain(2))
    write (LOC_UNIT_ID) (SNGL(y(ij)), ij=subdomain(3), subdomain(4))
    write (LOC_UNIT_ID) (SNGL(z(ij)), ij=subdomain(5), subdomain(6))

    close (LOC_UNIT_ID)

    return
end subroutine ENSIGHT_GRID

! ###################################################################
! ###################################################################
#ifdef USE_MPI

subroutine VISUALS_MPIO_AUX(opt_format, subdomain)

    use TLAB_VARS, only: imax, kmax, io_aux
    use TLAB_MPI_VARS
    use MPI
    use IO_FIELDS
    implicit none

    TINTEGER, intent(IN) :: opt_format, subdomain(6)

    ! -----------------------------------------------------------------------
    TINTEGER id, ny_loc

    ! #######################################################################
    io_aux(:)%active = .false.
    io_aux(:)%offset = 0
    io_aux(:)%precision = IO_TYPE_SINGLE
    if (opt_format == 1) io_aux(:)%offset = 244 ! # bytes of ensight header

    ny_loc = subdomain(4) - subdomain(3) + 1

    ! ###################################################################
    ! Saving full vertical xOy planes; using subdomain(5) to define the plane
    id = IO_SUBARRAY_VISUALS_XOY
    if (ims_pro_k == ((subdomain(5) - 1)/kmax)) io_aux(id)%active = .true.
    io_aux(id)%communicator = ims_comm_x
    io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOY(imax, ny_loc, MPI_REAL4)

    ! ndims = 2
    ! sizes(1)   = imax *ims_npro_i; sizes(2)   = subdomain(4)-subdomain(3)+1
    ! locsize(1) = imax;             locsize(2) = subdomain(4)-subdomain(3)+1
    ! offset(1)  = ims_offset_i;     offset(2)  = 0
    !
    ! CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
    !     MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
    ! CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ! Saving full vertical zOy planes; using subiddomain(1) to define the plane
    id = IO_SUBARRAY_VISUALS_ZOY
    if (ims_pro_i == ((subdomain(1) - 1)/imax)) io_aux(id)%active = .true.
    io_aux(id)%communicator = ims_comm_z
    io_aux(id)%subarray = IO_CREATE_SUBARRAY_ZOY(ny_loc, kmax, MPI_REAL4)

    ! ndims = 2
    ! sizes(1)   = subdomain(4)-subdomain(3)+1; sizes(2)   = kmax *ims_npro_k
    ! locsize(1) = subdomain(4)-subdomain(3)+1; locsize(2) = kmax
    ! offset(1)  = 0;                           offset(2)  = ims_offset_k
    !
    ! CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
    !     MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
    ! CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    ! Saving full blocks xOz planes
    id = IO_SUBARRAY_VISUALS_XOZ
    io_aux(id)%active = .true.
    io_aux(id)%communicator = MPI_COMM_WORLD
    io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOZ(imax, ny_loc, kmax, MPI_REAL4)

    ! ndims = 3
    ! sizes(1)   = imax *ims_npro_i; sizes(2)   = subdomain(4)-subdomain(3)+1; sizes(3)   = kmax *ims_npro_k
    ! locsize(1) = imax;             locsize(2) = subdomain(4)-subdomain(3)+1; locsize(3) = kmax
    ! offset(1)  = ims_offset_i;     offset(2)  = 0;                           offset(3)  = ims_offset_k
    !
    ! CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
    !     MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(id)%subarray, ims_err)
    ! CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

    return
end subroutine VISUALS_MPIO_AUX

#endif
