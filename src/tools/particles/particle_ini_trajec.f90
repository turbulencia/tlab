#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "INI_TRAJEC"

!########################################################################
!# Tool/Library INIT/PARTICLE
!#
!########################################################################
!# HISTORY
!#
!# 2015/01/16 - L. Muessle
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# Reads position data of the largest files for the start of simulation.
!#
!# Position data from the start of simulation is provided by program
!# l_pos_trajec.x. The file is called pos_largest_file_start.
!#
!# Less particles are created by setting smaller number for
!# particle_number in tlab.ini. Afterwards particles are exchanged with
!# position of the largest particles determined before with l_trajec.x
!# and l_pos_trajec.x
!#
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
program PARTICLE_INI_TRAJEC

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only: ims_err
    use TLAB_MPI_VARS, only: ims_pro, ims_npro
    use TLAB_MPI_PROCS
#endif
    use Thermodynamics
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS

    implicit none

! -------------------------------------------------------------------
! Additional local arrays

    TINTEGER, dimension(:), allocatable :: dummy_proc
    integer(8), dimension(:), allocatable :: l_traj_tags, fake_l_traj_tags
#ifdef USE_MPI
    TINTEGER dummy_ims_npro
    TINTEGER dummy_isize_traj
    integer(8), dimension(:), allocatable :: all_fake_l_traj_tags
#endif
    TREAL, dimension(:, :), allocatable :: l_traj
    TREAL, dimension(:), allocatable :: fake_liquid, all_fake_liquid

    TINTEGER nitera_first

    character*32 bakfile

#ifdef USE_MPI
    character*64 fname, str
    TINTEGER particle_pos, i
    TLONGINTEGER dummy
#endif

    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call Thermodynamics_Initialize(ifile)
    call PARTICLE_READ_GLOBAL(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    bakfile = TRIM(ADJUSTL(ifile))//'.bak'

! Get the local information from the tlab.ini
    call SCANINIINT(bakfile, ifile, 'Iteration', 'Start', '0', nitera_first)

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (wrk1d(isize_wrk1d, inb_wrk1d))
    allocate (wrk2d(isize_wrk2d, 1))
    allocate (wrk3d(isize_wrk3d))

    if (part%type == PART_TYPE_BIL_CLOUD_3 .or. part%type == PART_TYPE_BIL_CLOUD_4) then !Allocte memory to read fields
        allocate (txc(isize_field, 3))
    end if

    call PARTICLE_ALLOCATE(C_FILE_LOC)

    allocate (dummy_proc(isize_traj))
    allocate (l_traj_tags(isize_traj))
    allocate (fake_l_traj_tags(isize_traj))
    allocate (l_traj(3, isize_traj))
    allocate (fake_liquid(isize_traj))
    allocate (all_fake_liquid(isize_traj))
#ifdef USE_MPI
    allocate (all_fake_l_traj_tags(isize_traj))
    all_fake_l_traj_tags(:) = C_0_R
#endif
    fake_l_traj_tags(:) = C_0_R
    fake_liquid(:) = C_0_R
    all_fake_liquid(:) = C_0_R

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    !#######################################################################
    !CREATE THE RANDOM PARTICLE FIELD
    !#######################################################################
    call PARTICLE_RANDOM_POSITION(l_q, l_txc, txc)

#ifdef USE_MPI
    !#######################################################################
    !READ THE LARGEST PARTICLE FILE
    !#######################################################################

    if (ims_pro == 0) then
        fname = 'pos_largest_particle_start'
        write (str, *) nitera_first; str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
        open (unit=117, file=str, access='stream', form='unformatted')
        read (117) dummy_ims_npro   !is a integer
        read (117, POS=SIZEOFINT + 1) dummy_isize_traj  !is an integer
        read (117, POS=SIZEOFINT*2 + 1) l_traj_tags  !is INTEGER(8)
        read (117, POS=(SIZEOFINT*2 + 1) + SIZEOFREAL*isize_traj) l_traj ! attention is integer(8)
        read (117, POS=(SIZEOFINT*2 + 1) + SIZEOFREAL*isize_traj + 3*isize_traj*SIZEOFREAL) dummy_proc ! attention is integer(8)
        close (117)
    end if

    !#######################################################################
    !BROADCAST INFORMATION OF LARGEST PARTICLES
    !#######################################################################
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(l_traj_tags, isize_traj, MPI_INTEGER8, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(l_traj, isize_traj*3, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(dummy_proc, isize_traj, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)

    !#######################################################################
    !REPLACE THE FIRST PARTICLES WITH THE LARGEST
    !CORRESPONDING TO THE PROCESSORS
    !#######################################################################
    dummy = 1
    do i = 1, isize_traj
        if (ims_pro == dummy_proc(i)) then
            l_q(dummy, 1) = l_traj(1, i)
            l_q(dummy, 2) = l_traj(2, i)
            l_q(dummy, 3) = l_traj(3, i)
            fake_l_traj_tags(i) = l_g%tags(dummy)
            fake_liquid(i) = l_q(dummy, 5) ! just to be consistent with other output file
            dummy = dummy + 1
        end if
    end do

    !#######################################################################
    !WRITE FILE WITH NEW FAKE LARGEST IDs
    !#######################################################################
  CALL MPI_REDUCE(fake_l_traj_tags, all_fake_l_traj_tags, isize_traj, MPI_INTEGER8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
    call MPI_REDUCE(fake_liquid, all_fake_liquid, isize_traj, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)

    !#######################################################################
    !WRITE FILE WITH NEW FAKE LARGEST IDs
    !#######################################################################
    if (ims_pro == 0) then
        fname = 'fake_largest_particle'
        open (unit=15, file=fname, access='stream', form='unformatted')
        inquire (UNIT=15, POS=particle_pos) !would be 1
        write (15) ims_npro  !header
        inquire (UNIT=15, POS=particle_pos) !would be 5
        write (15) isize_traj  !header
        inquire (UNIT=15, POS=particle_pos)  !would be 9
        write (15) fake_liquid ! just to be consistent with other output file
        inquire (UNIT=15, POS=particle_pos)  !would be 409 with 50 numbers
        write (15) all_fake_l_traj_tags
        close (15)
    end if

    !#######################################################################
    !WRITE THE PARTICLE FILE
    !#######################################################################
    call IO_WRITE_PARTICLE('part2traj.ics', l_g, l_q)

#endif
    call TLAB_STOP(0)
end program PARTICLE_INI_TRAJEC
