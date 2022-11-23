#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "PARTICLE_PDF"

!########################################################################
!########################################################################
program PARTICLE_BUILD_PDF

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
    use IO_FIELDS
#ifdef USE_MPI
    use TLAB_MPI_PROCS
#endif
    use PARTICLE_VARS
    use PARTICLE_ARRAYS

    implicit none
#include "integers.h"

! -------------------------------------------------------------------

    TINTEGER nitera_first, nitera_last, nitera_save
    TINTEGER ierr, i

    character*64 fname, str
    character*128 line
    character*32 bakfile

    bakfile = TRIM(ADJUSTL(ifile))//'.bak'

    call TLAB_START

    call IO_READ_GLOBAL(ifile)
    call PARTICLE_READ_GLOBAL('dns.ini')

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif
!  CALL DNS_READ_LOCAL(ifile) !for nitera stuff

! Get the local information from the dns.ini
    call SCANINIINT(bakfile, ifile, 'Iteration', 'Start', '0', nitera_first)
    call SCANINIINT(bakfile, ifile, 'Iteration', 'End', '0', nitera_last)
    call SCANINIINT(bakfile, ifile, 'Iteration', 'Restart', '50', nitera_save)

    inb_part_txc = 1

    call PARTICLE_ALLOCATE(C_FILE_LOC)

    isize_wrk3d = imax*jmax*kmax
    isize_wrk3d = MAX(isize_wrk3d, (imax + 1)*jmax*(kmax + 1))
    isize_wrk3d = MAX(isize_wrk3d, (jmax*(kmax + 1)*inb_part_interp*2))
    isize_wrk3d = MAX(isize_wrk3d, (jmax*(imax + 1)*inb_part_interp*2))

    isize_wrk2d = MAX(isize_wrk2d, jmax*inb_part_interp)

    ! IF (jmax_part .EQ. 1) THEN
    !    jmax_part   = jmax ! 1 by default
    ! ENDIF
! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (wrk1d(isize_wrk1d, inb_wrk1d))
    allocate (wrk2d(isize_wrk2d, 1))
    allocate (wrk3d(isize_wrk3d))

    allocate (s(isize_field, inb_scal_array))

    if (imode_part == PART_TYPE_BIL_CLOUD_3 .or. imode_part == PART_TYPE_BIL_CLOUD_4) then !Allocte memory to read fields
        allocate (txc(isize_field, 3))
    end if

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

!#######################################################################
!Setup the file properties
!#######################################################################
    do i = nitera_first, nitera_last, nitera_save

!#######################################################################
!READ ALL FILES
!#######################################################################
        write (fname, *) i; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, i0, s, wrk3d)
        call THERMO_AIRWATER_LINEAR(imax, jmax, kmax, s, s(1, inb_scal_array))

        write (fname, *) i; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
        call IO_READ_PARTICLE(fname, l_g, l_q)

! ######################################################################
! Save particle pathlines for particle_pdf
! ######################################################################
        write (fname, *) i; fname = "particle_pdf."//TRIM(ADJUSTL(fname))
        call PARTICLE_PDF(fname, s, l_g, l_q, l_txc, l_comm, wrk3d)

    end do

    call TLAB_STOP(0)
end program
