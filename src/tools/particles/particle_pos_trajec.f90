#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "PARTICLE_POS_TRAJEC"

!########################################################################
!# Tool/Library PLOT
!#
!########################################################################
!# HISTORY
!#
!# 2014/01/16 - L. Muessle
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# This program reads the ID(TAG) of the largest particle file from the
!# the program l_trajec.x and then determines the position of these
!# largest particles at the start point of the simulation.
!#
!#
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
program PARTICLE_POS_TRAJEC

    use TLab_Constants
    use TLAB_VARS
    use TLab_WorkFlow
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_err
    use TLabMPI_VARS, only: ims_pro, ims_npro
    use TLabMPI_PROCS
#endif
    use Thermodynamics
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS

    implicit none

! -------------------------------------------------------------------

#ifdef USE_MPI
    integer(wi) ierr, i, j, particle_pos

    integer(wi) dummy_ims_npro
    integer(wi) dummy_isize_traj
    integer(wi), dimension(:), allocatable :: dummy_proc, all_dummy_proc
    integer(8), dimension(:), allocatable :: l_traj_tags
    real(wp), dimension(:), allocatable :: dummy_big_overall
    real(wp), dimension(:, :), allocatable :: l_traj, all_l_traj

    integer(wi) nitera_first, nitera_last
    integer(wi) :: isize_traj                       ! # of saved trajectories

    character*64 str, fname
    character*128 line
    character*32 bakfile

    bakfile = trim(adjustl(ifile))//'.bak'

    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Particle_Initialize_Parameters('tlab.ini')

#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif

! Get the local information from the tlab.ini
    call SCANINIINT(bakfile, ifile, 'Particle', 'TrajNumber', '0', isize_traj)
    call SCANINIINT(bakfile, ifile, 'Iteration', 'Start', '0', nitera_first)

    call Particle_Initialize_Memory(C_FILE_LOC)

    allocate (dummy_proc(isize_traj))
    allocate (all_dummy_proc(isize_traj))
    allocate (dummy_big_overall(isize_traj))
    allocate (l_traj_tags(isize_traj))
    allocate (l_traj(3, isize_traj))
    allocate (all_l_traj(3, isize_traj))

    l_traj(:, :) = C_0_R
    all_l_traj(:, :) = C_0_R
    dummy_proc(:) = C_0_R
    all_dummy_proc(:) = C_0_R

!#######################################################################
!READ THE (FIRST) FILE
!#######################################################################
    write (fname, *) nitera_first; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
    call IO_READ_PARTICLE(fname, l_g, l_q)

    if (ims_pro == 0) then
        fname = 'largest_particle'
        write (str, *) nitera_last; str = trim(adjustl(fname))//"."//trim(adjustl(str)) ! name with the number for direction and scalar
        open (unit=117, file=str, access='stream', form='unformatted')
        read (117) dummy_ims_npro   !is a integer
        read (117, POS=SIZEOFINT + 1) dummy_isize_traj  !is an integer
        read (117, POS=SIZEOFINT*2 + 1) dummy_big_overall  !is real(8)
        read (117, POS=(SIZEOFINT*2 + 1) + SIZEOFREAL*isize_traj) l_traj_tags ! attention is integer(8)
        close (117)
    end if

!#######################################################################
!BROADCAST THE ID OF THE LARGEST PARTICLES
!#######################################################################
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    call MPI_BCAST(l_traj_tags, isize_traj, MPI_INTEGER8, 0, MPI_COMM_WORLD, ims_err)

!#######################################################################
!SEARCH FOR LARGEST PARTICLES
!#######################################################################
    do i = 1, l_g%np
        do j = 1, isize_traj
            if (l_g%tags(i) == l_traj_tags(j)) then
                l_traj(1, j) = l_q(i, 1)
                l_traj(2, j) = l_q(i, 2)
                l_traj(3, j) = l_q(i, 3)
                dummy_proc(j) = ims_pro
            end if
        end do
    end do

!#######################################################################
!REDUCE ALL INFORMATION TO ROOT
!#######################################################################
    call MPI_REDUCE(l_traj, all_l_traj, 3*isize_traj, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
    call MPI_REDUCE(dummy_proc, all_dummy_proc, isize_traj, MPI_INTEGER4, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

!#######################################################################
!WRITE DATA WITH ROOT
!#######################################################################

    if (ims_pro == 0) then
        fname = 'pos_largest_particle_start'
        write (str, *) nitera_first; str = trim(adjustl(fname))//"."//trim(adjustl(str)) ! name with the number for direction and scalar
        open (unit=15, file=str, access='stream', form='unformatted')
        inquire (UNIT=15, POS=particle_pos) !would be 1
        write (15) ims_npro  !header
        inquire (UNIT=15, POS=particle_pos) !would be 5
        write (15) isize_traj  !header
        inquire (UNIT=15, POS=particle_pos)  !would be 9
        write (15) l_traj_tags
        inquire (UNIT=15, POS=particle_pos)  !409
        write (15) all_l_traj
        inquire (UNIT=15, POS=particle_pos)  !would be 1609 with 50 numbers
        write (15) all_dummy_proc
        close (15)
    end if
!    !Just for testing and as template
!      IF (ims_pro .EQ. 0) THEN
!    OPEN(unit=117, file=str, access='stream', form='unformatted')
!    READ(117) test1   !is a integer
!    READ(117, POS=SIZEOFINT+1) test2  !is an integer
!    READ(117, POS=SIZEOFINT*2+1) test3  !is real(8)
!    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_traj) test4 ! attention is integer(8)
!    CLOSE(117)
!    ENDIF

#endif

    call TLAB_STOP(0)
end program
