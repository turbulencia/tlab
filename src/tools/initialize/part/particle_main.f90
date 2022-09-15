#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "INIPART"

program INIPART
    use TLAB_TYPES, only: cp, ci
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_PROCS
#endif
    use PARTICLE_VARS
    use PARTICLE_ARRAYS

    implicit none

    ! -------------------------------------------------------------------
    integer(ci) ierr
    real(cp), dimension(:), allocatable, save :: l_comm
    character*64 str, line

    !########################################################################
    !########################################################################
    call TLAB_START

    call IO_READ_GLOBAL(ifile)

    call PARTICLE_READ_GLOBAL(ifile)

    if (imode_part /= PART_TYPE_NONE) then
#ifdef USE_MPI
        call TLAB_MPI_INITIALIZE
#endif

        ! -------------------------------------------------------------------
        ! Allocating memory space
        ! -------------------------------------------------------------------
        inb_flow_array = 0
        inb_scal_array = 0
        isize_wrk3d = imax*jmax*kmax
        inb_txc = inb_scal

        call TLAB_ALLOCATE(C_FILE_LOC)

        call PARTICLE_ALLOCATE(C_FILE_LOC)

        write (str, *) isize_l_comm; line = 'Allocating array l_comm of size '//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (l_comm(isize_l_comm), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'Not enough memory for l_comm.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

        ! -------------------------------------------------------------------
        ! Read the grid
        ! -------------------------------------------------------------------
        call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
        call FDM_INITIALIZE(x, g(1), wrk1d)
        call FDM_INITIALIZE(y, g(2), wrk1d)
        call FDM_INITIALIZE(z, g(3), wrk1d)

        call FI_BACKGROUND_INITIALIZE(wrk1d)

        ! -------------------------------------------------------------------
        ! Initialize particle information
        ! -------------------------------------------------------------------
        call PARTICLE_RANDOM_POSITION(l_q, l_txc, l_comm, txc, wrk3d)

        call IO_WRITE_PARTICLE(TRIM(ADJUSTL(tag_part))//'ics', l_g, l_q)

    end if

    call TLAB_STOP(0)
end program INIPART
