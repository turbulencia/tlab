#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "PARTICLE_PDF"

!########################################################################
!########################################################################
program PARTICLE_BUILD_PDF

    use TLab_Constants
    use TLAB_VARS
    use TLab_Arrays
    use TLab_WorkFlow
    use IO_FIELDS
#ifdef USE_MPI
    use TLabMPI_PROCS
#endif
    use Thermodynamics
    use THERMO_AIRWATER
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS

    implicit none

! -------------------------------------------------------------------

    integer(wi) nitera_first, nitera_last, nitera_save
    integer(wi) i

    character*64 fname
    character*32 bakfile

    bakfile = trim(adjustl(ifile))//'.bak'

    call TLAB_START

    call IO_READ_GLOBAL(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif
    call Thermodynamics_Initialize_Parameters(ifile)
    call Particle_Initialize_Parameters('tlab.ini')

!  CALL DNS_READ_LOCAL(ifile) !for nitera stuff

! Get the local information from the tlab.ini
    call SCANINIINT(bakfile, ifile, 'Iteration', 'Start', '0', nitera_first)
    call SCANINIINT(bakfile, ifile, 'Iteration', 'End', '0', nitera_last)
    call SCANINIINT(bakfile, ifile, 'Iteration', 'Restart', '50', nitera_save)

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (wrk1d(isize_wrk1d, inb_wrk1d))
    allocate (wrk2d(isize_wrk2d, 1))
    allocate (wrk3d(isize_wrk3d))

    allocate (s(isize_field, inb_scal_array))

    if (part%type == PART_TYPE_BIL_CLOUD_3 .or. part%type == PART_TYPE_BIL_CLOUD_4) then !Allocte memory to read fields
        allocate (txc(isize_field, 3))
    end if

    inb_part_txc = 1

    call Particle_Initialize_Memory(C_FILE_LOC)

    isize_wrk2d = max(isize_wrk2d, jmax*inb_part_interp)

    ! IF (jmax_part .EQ. 1) THEN
    !    jmax_part   = jmax ! 1 by default
    ! ENDIF

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z)
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
        write (fname, *) i; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, 0, s)
        call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, inb_scal_array))

        write (fname, *) i; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
        call IO_READ_PARTICLE(fname, l_g, l_q)

! ######################################################################
! Save particle pathlines for particle_pdf
! ######################################################################
        write (fname, *) i; fname = "particle_pdf."//trim(adjustl(fname))
        call PARTICLE_PDF(fname, s, l_g, l_q, l_txc)

    end do

    call TLAB_STOP(0)
end program
