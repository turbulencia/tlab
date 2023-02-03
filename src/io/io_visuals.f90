#include "dns_const.h"

#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

subroutine IO_WRITE_VISUALS(fname, iformat, nx, ny, nz, nfield, subdomain, field, txc)
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_TYPES, only: subarray_dt
    use TLAB_VARS, only: g, isize_txc_field, io_aux
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only: ims_pro
    use TLAB_MPI_PROCS
#endif
    use IO_FIELDS

    implicit none

    integer(wi) iformat, nx, ny, nz, nfield, subdomain(6)
    real(wp), dimension(isize_txc_field, nfield) :: field
    real(wp), dimension(nx*ny*nz, 2) :: txc
    character*(*) fname

    ! -------------------------------------------------------------------
    integer(wi) sizes(5), nx_aux, ny_aux, nz_aux, ifield, i
    character*32 varname(16), name
    integer(wi) iflag_mode
    integer, parameter :: i1 = 1

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
    use TLAB_CONSTANTS, only: wp, wi
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_pro
    use TLAB_MPI_PROCS
#endif

    implicit none

    character*(*) name
    integer(wi), intent(IN) :: iheader ! 0 no header; 1 header

    integer(wi), intent(IN) :: nx, ny, nz, nfield, subdomain(6)
    real(wp), dimension(nx, ny, nz, nfield) :: field
    real(wp), dimension(nx*ny*nz, 2) :: tmp_mpi

    ! -------------------------------------------------------------------
    character*80 line
#ifdef USE_MPI
    integer(wi) ifield
#else
    integer(wi) i, j, k, ifield
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
            write (LOC_UNIT_ID) 1_wi
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
    use TLAB_CONSTANTS, only: wp, wi

    implicit none

    integer(wi) nx, ny, nz, subdomain(6)
    real(wp) x(nx), y(ny), z(nz)
    character*(*) name

    ! -------------------------------------------------------------------
    character*80 line
    integer(wi) ij

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
    write (LOC_UNIT_ID) 1_wi
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
    use TLAB_CONSTANTS, only: wp, wi

    use TLAB_VARS, only: imax, kmax, io_aux
    use TLAB_MPI_VARS
    use MPI
    use IO_FIELDS
    implicit none

    integer(wi), intent(IN) :: opt_format, subdomain(6)

    ! -----------------------------------------------------------------------
    integer(wi) id, ny_loc

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

    ! Saving full vertical zOy planes; using subiddomain(1) to define the plane
    id = IO_SUBARRAY_VISUALS_ZOY
    if (ims_pro_i == ((subdomain(1) - 1)/imax)) io_aux(id)%active = .true.
    io_aux(id)%communicator = ims_comm_z
    io_aux(id)%subarray = IO_CREATE_SUBARRAY_ZOY(ny_loc, kmax, MPI_REAL4)

    ! Saving full blocks xOz planes
    id = IO_SUBARRAY_VISUALS_XOZ
    io_aux(id)%active = .true.
    io_aux(id)%communicator = MPI_COMM_WORLD
    io_aux(id)%subarray = IO_CREATE_SUBARRAY_XOZ(imax, ny_loc, kmax, MPI_REAL4)

    return
end subroutine VISUALS_MPIO_AUX

#endif
