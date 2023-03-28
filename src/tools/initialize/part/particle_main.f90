#define C_FILE_LOC "INIPART"

program INIPART
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_PROCS
#endif
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS

    implicit none

    !########################################################################
    !########################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call THERMO_INITIALIZE()
    call PARTICLE_READ_GLOBAL(ifile)

    if (part%type /= PART_TYPE_NONE) then
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

        ! -------------------------------------------------------------------
        ! Read the grid
        ! -------------------------------------------------------------------
        call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
        call FDM_INITIALIZE(x, g(1), wrk1d)
        call FDM_INITIALIZE(y, g(2), wrk1d)
        call FDM_INITIALIZE(z, g(3), wrk1d)

        call FI_BACKGROUND_INITIALIZE()

        call PARTICLE_INITIALIZE()

        ! -------------------------------------------------------------------
        ! Initialize particle information
        ! -------------------------------------------------------------------
        call PARTICLE_RANDOM_POSITION(l_q, l_txc, txc)

        call IO_WRITE_PARTICLE(TRIM(ADJUSTL(tag_part))//'ics', l_g, l_q)

    end if

    call TLAB_STOP(0)
end program INIPART
