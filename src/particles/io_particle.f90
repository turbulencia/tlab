#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI

#endif

#define USE_ACCESS_STREAM

!#######################################################################
!#######################################################################
#define LOC_UNIT_ID 117
#define LOC_STATUS 'old'

subroutine IO_READ_PARTICLE(fname, l_g, l_q)
    use TLab_Constants, only: wp, wi, longi, lfile, efile, sizeofint, sizeofreal, sizeoflongint
    use FDM, only: g
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use PARTICLE_VARS, only: isize_part, inb_part_array, isize_part_total
    use PARTICLE_TYPES, only: particle_dt
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_pro, ims_npro, ims_err
    use PARTICLE_ARRAYS, only: ims_np_all
#endif

    implicit none

    character*(*) fname
    type(particle_dt) l_g
    real(wp), dimension(isize_part, inb_part_array) :: l_q !, OPTIONAL :: l_q

! -------------------------------------------------------------------
    integer(wi) i
    character(len=32) name
#ifdef USE_MPI
    integer(wi) ims_npro_loc
    type(MPI_File) mpio_fh
    integer(longi) mpio_disp, count
    type(MPI_Status) status
#else
    integer(wi) idummy
#endif

    call TLab_Write_ASCII(lfile, 'Reading field '//trim(adjustl(fname))//'...')

#ifdef USE_MPI
!#######################################################################
! Parallel case
!#######################################################################
! -------------------------------------------------------------------
! Let Process 0 handle header
! -------------------------------------------------------------------
    if (ims_pro == 0) then
        name = trim(adjustl(fname))//".id"
#include "dns_open_file.h"
        read (LOC_UNIT_ID) ims_npro_loc
        read (LOC_UNIT_ID) ims_np_all(1:ims_npro_loc)
        close (LOC_UNIT_ID)
    end if

! Check
    call MPI_BCAST(ims_npro_loc, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ims_err)
    if (ims_npro /= ims_npro_loc) then
        call TLab_Write_ASCII(efile, 'IO_PARTICLE. Number-of-processors mismatch.')
        call TLab_Stop(DNS_ERROR_PARTICLE)
    end if

! Broadcast number of particles per processor
    call MPI_BCAST(ims_np_all, ims_npro, MPI_INTEGER, 0, MPI_COMM_WORLD, ims_err)

! Number of particles in local processor
    l_g%np = ims_np_all(ims_pro + 1)

! Displacement per processor
    mpio_disp = int((ims_npro + 1)*SIZEOFINT, longi)

    count = 0
    do i = 1, ims_pro
        count = count + int(ims_np_all(i), longi)
    end do
    mpio_disp = mpio_disp + count*int(SIZEOFLONGINT, longi)

! Check
    do i = ims_pro + 1, ims_npro
        count = count + int(ims_np_all(i), longi)
    end do
    if (isize_part_total /= count) then
        call TLab_Write_ASCII(efile, 'IO_PARTICLE. Number-of-particles mismatch.')
        call TLab_Stop(DNS_ERROR_PARTICLE)
    end if

! -------------------------------------------------------------------
! Use MPI-IO to read particle tags in each processor
! -------------------------------------------------------------------
    name = trim(adjustl(fname))//".id"
    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ims_err)
    call MPI_FILE_READ_ALL(mpio_fh, l_g%tags, l_g%np, MPI_INTEGER8, status, ims_err)
    call MPI_FILE_CLOSE(mpio_fh, ims_err)

!  IF ( PRESENT(l_q) ) THEN
    do i = 1, inb_part_array
        write (name, *) i; name = trim(adjustl(fname))//"."//trim(adjustl(name))
        call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
        call MPI_FILE_READ_ALL(mpio_fh, l_q(1, i), l_g%np, MPI_REAL8, status, ims_err)
        call MPI_FILE_CLOSE(mpio_fh, ims_err)
    end do
!  ENDIF

#else
! #######################################################################
! Serial case
! #######################################################################
    name = trim(adjustl(fname))//".id"
#include "dns_open_file.h"
    read (LOC_UNIT_ID) idummy                   ! dummy, should be 1 in serial
    read (LOC_UNIT_ID) l_g%np
! Check
    if (isize_part_total /= int(l_g%np, longi)) then
        call TLab_Write_ASCII(efile, 'IO_PARTICLE. Number-of-particles mismatch.')
        close (LOC_UNIT_ID)
        call TLab_Stop(DNS_ERROR_PARTICLE)
    end if
    read (LOC_UNIT_ID) l_g%tags
    close (LOC_UNIT_ID)

!  IF ( PRESENT(l_q) ) THEN
    do i = 1, inb_part_array
        write (name, *) i; name = trim(adjustl(fname))//"."//trim(adjustl(name))
#include "dns_open_file.h"
        read (LOC_UNIT_ID) idummy             ! dummy, should be 1 in serial
        read (LOC_UNIT_ID) l_g%np
        read (LOC_UNIT_ID) l_q(:, i)
        close (LOC_UNIT_ID)
    end do
!  ENDIF

#endif

    call LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, g(2)%size, g(2)%nodes)

    return
end subroutine IO_READ_PARTICLE

#undef LOC_UNIT_ID
#undef LOC_STATUS

!#######################################################################
!#######################################################################
#define LOC_UNIT_ID 118
#define LOC_STATUS 'unknown'

subroutine IO_WRITE_PARTICLE(fname, l_g, l_q)

    use TLab_Constants, only: wp, wi, longi, lfile, sizeofint, sizeoflongint
    use PARTICLE_VARS, only: isize_part, inb_part_array
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use PARTICLE_TYPES, only: particle_dt
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS, only: ims_pro, ims_npro, ims_err
    use PARTICLE_ARRAYS, only: ims_np_all
#endif

    implicit none

    character(len=*) fname
    type(particle_dt) l_g
    real(wp), dimension(isize_part, inb_part_array) :: l_q !, OPTIONAL :: l_q

! -------------------------------------------------------------------
    integer(wi) i
    character(len=32) name
#ifdef USE_MPI
    type(MPI_File) mpio_fh
    integer(KIND=8) mpio_disp, count
    type(MPI_Status) status
#else
    integer(wi) idummy
#endif

    call TLab_Write_ASCII(lfile, 'Writing field '//trim(adjustl(fname))//'...')

#ifdef USE_MPI
!#######################################################################
! Parallel case
!#######################################################################
! -------------------------------------------------------------------
! Let Process 0 handle header
! -------------------------------------------------------------------
    call MPI_ALLGATHER(l_g%np, 1, MPI_INTEGER4, ims_np_all, 1, MPI_INTEGER4, MPI_COMM_WORLD, ims_err)

    if (ims_pro == 0) then
        name = trim(adjustl(fname))//".id"
#include "dns_open_file.h"
        write (LOC_UNIT_ID) ims_npro
        write (LOC_UNIT_ID) ims_np_all
        close (LOC_UNIT_ID)

!     IF ( PRESENT(l_q) ) THEN
        do i = 1, inb_part_array
            write (name, *) i; name = trim(adjustl(fname))//"."//trim(adjustl(name))
#include "dns_open_file.h"
            write (LOC_UNIT_ID) ims_npro
            write (LOC_UNIT_ID) ims_np_all
            close (LOC_UNIT_ID)
        end do
!     ENDIF
    end if

! Displacement per processor
    mpio_disp = int((ims_npro + 1)*SIZEOFINT, longi)

    count = 0
    do i = 1, ims_pro
        count = count + int(ims_np_all(i), longi)
    end do
    mpio_disp = mpio_disp + count*int(SIZEOFLONGINT, longi)

! -------------------------------------------------------------------
! Use MPI-IO to write particle tags in each processor
! -------------------------------------------------------------------
    name = trim(adjustl(fname))//".id"
    call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
    call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_INTEGER8, MPI_INTEGER8, 'native', MPI_INFO_NULL, ims_err)
    call MPI_FILE_WRITE_ALL(mpio_fh, l_g%tags, l_g%np, MPI_INTEGER8, status, ims_err)
    call MPI_FILE_CLOSE(mpio_fh, ims_err)

!  IF ( PRESENT(l_q) ) THEN
    do i = 1, inb_part_array
        write (name, *) i; name = trim(adjustl(fname))//"."//trim(adjustl(name))
        call MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        call MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
        call MPI_FILE_WRITE_ALL(mpio_fh, l_q(1, i), l_g%np, MPI_REAL8, status, ims_err)
        call MPI_FILE_CLOSE(mpio_fh, ims_err)
    end do
!  ENDIF

#else
! #######################################################################
! Serial case
! #######################################################################
    idummy = 1
    name = trim(adjustl(fname))//".id"
#include "dns_open_file.h"
    write (LOC_UNIT_ID) idummy
    write (LOC_UNIT_ID) l_g%np
    write (LOC_UNIT_ID) l_g%tags
    close (LOC_UNIT_ID)

!  IF ( PRESENT(l_q) ) THEN
    do i = 1, inb_part_array
        write (name, *) i; name = trim(adjustl(fname))//"."//trim(adjustl(name))
#include "dns_open_file.h"
        write (LOC_UNIT_ID) idummy
        write (LOC_UNIT_ID) l_g%np
        write (LOC_UNIT_ID) l_q(:, i)
        close (LOC_UNIT_ID)
    end do
!  ENDIF

#endif

    return
end subroutine IO_WRITE_PARTICLE
