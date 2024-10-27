#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "PARTICLE_TRAJEC"

!########################################################################
!# Tool/Library PLOT
!#
!########################################################################
!# HISTORY
!#
!# 2014/12/15 - L. Muessle
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# Post processing to find the largest particles
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
program PARTICLE_TRAJEC

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

    integer(wi) i, j, k, particle_pos
    real(wp) temp
    integer(8) itemp
    logical :: swapped
    real(wp), dimension(:), allocatable :: big_part
#ifdef USE_MPI
    real(wp), dimension(:), allocatable :: big_all
    real(wp), dimension(:), allocatable :: big_overall
    integer(8), dimension(:), allocatable :: tag_big_all, tag_big_overall
#else
    integer(wi) dummy
#endif
    integer(8), dimension(:), allocatable :: tag_big_part

    integer(wi) nitera_last
    integer(wi) :: isize_traj                       ! # of saved trajectories

    character*64 str, fname
    character*32 bakfile

!  integer(wi) test1, test2
!  real(wp) test3(50)
!  INTEGER(8) test4(50)
    bakfile = trim(adjustl(ifile))//'.bak'

    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize()
#endif
    call Particle_Initialize_Parameters(ifile)

    call Thermodynamics_Initialize_Parameters(ifile)
! Get the local information from the tlab.ini
    call SCANINIINT(bakfile, ifile, 'Particle', 'TrajNumber', '0', isize_traj)
    call SCANINIINT(bakfile, ifile, 'Iteration', 'End', '0', nitera_last)

    call Particle_Initialize_Memory(C_FILE_LOC)

    allocate (big_part(isize_traj))
    allocate (tag_big_part(isize_traj))
#ifdef USE_MPI
    if (ims_pro == 0) then
        allocate (big_all(isize_traj*ims_npro))
        allocate (big_overall(isize_traj))
        allocate (tag_big_all(isize_traj*ims_npro))
        allocate (tag_big_overall(isize_traj))
    end if

#endif

!#######################################################################
!READ THE (LAST) FILE
!#######################################################################
    write (fname, *) nitera_last; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
    call IO_READ_PARTICLE(fname, l_g, l_q)

!#######################################################################
!Every processor searches for the largest particles
!#######################################################################

    big_part(1:isize_traj) = l_q(1:isize_traj, 5)
    tag_big_part(1:isize_traj) = l_g%tags(1:isize_traj)

!  IF (ims_pro .EQ. 0) THEN
!  print*, tag_big_part
!  print*, big_part
!  ENDIF

    do j = isize_traj - 1, 1, -1
        swapped = .false.
        do i = 1, j
            if (big_part(i) > big_part(i + 1)) then !if 49... bigger than 50...
!START SWAPPING
                temp = big_part(i)
                big_part(i) = big_part(i + 1)
                big_part(i + 1) = temp

                itemp = tag_big_part(i)
                tag_big_part(i) = tag_big_part(i + 1)
                tag_big_part(i + 1) = itemp
                swapped = .true.
            end if
        end do
        if (.not. swapped) exit
    end do

!DO k=isize_traj+1,l_g%np,1
    do k = 1, l_g%np
        if (l_q(k, 5) > big_part(1)) then
            big_part(1) = l_q(k, 5)
            tag_big_part(1) = l_g%tags(k)

            swapped = .false.
            do i = 1, isize_traj - 1
                if (big_part(i) > big_part(i + 1)) then !if 49... bigger than 50...
!START SWAPPING
                    temp = big_part(i)
                    big_part(i) = big_part(i + 1)
                    big_part(i + 1) = temp

                    itemp = tag_big_part(i)
                    tag_big_part(i) = tag_big_part(i + 1)
                    tag_big_part(i + 1) = itemp
                    swapped = .true.
                end if
                if (.not. swapped) exit
            end do

        end if
    end do

!  IF (ims_pro .EQ. 0) THEN
!  print*, 'delimiter'
!  print*, tag_big_part
!  print*, big_part
!  END IF

!#######################################################################
!Send all particles to one processor
!#######################################################################
#ifdef USE_MPI
    call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
    call MPI_GATHER(big_part, isize_traj, MPI_REAL8, big_all, isize_traj, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
    call MPI_GATHER(tag_big_part, isize_traj, MPI_INTEGER8, tag_big_all, isize_traj, MPI_INTEGER8, 0, MPI_COMM_WORLD, ims_err)

    if (ims_pro == 0) then

!#######################################################################
!Determine the oerall biggest 50 of the biggest 50 of each processor
!#######################################################################

        big_overall = big_part
        tag_big_overall = tag_big_part
        do k = isize_traj + 1, isize_traj*ims_npro, 1
            if (big_all(k) > big_overall(1)) then
                big_overall(1) = big_all(k)
                tag_big_overall(1) = tag_big_all(k)

                swapped = .false.
                do i = 1, isize_traj - 1
                    if (big_overall(i) > big_overall(i + 1)) then !if 49... bigger than 50...
!START SWAPPING
                        temp = big_overall(i)
                        big_overall(i) = big_overall(i + 1)
                        big_overall(i + 1) = temp

                        itemp = tag_big_overall(i)
                        tag_big_overall(i) = tag_big_overall(i + 1)
                        tag_big_overall(i + 1) = itemp
                        swapped = .true.
                    end if
                    if (.not. swapped) exit
                end do

            end if
        end do

!#######################################################################
!Write data file with root processor
!#######################################################################

        fname = 'largest_particle'
        write (str, *) nitera_last; str = trim(adjustl(fname))//"."//trim(adjustl(str)) ! name with the number for direction and scalar
        open (unit=15, file=str, access='stream', form='unformatted')
        inquire (UNIT=15, POS=particle_pos) !would be 1
        write (15) ims_npro  !header
        inquire (UNIT=15, POS=particle_pos) !would be 5
        write (15) isize_traj  !header
        inquire (UNIT=15, POS=particle_pos)  !would be 9
        write (15) big_overall
        inquire (UNIT=15, POS=particle_pos)  !would be 409 with 50 numbers
        write (15) tag_big_overall
        close (15)

!    !Just for testing and as template
!    OPEN(unit=117, file=str, access='stream', form='unformatted')
!    READ(117) test1   !is a integer
!    READ(117, POS=SIZEOFINT+1) test2  !is an integer
!    READ(117, POS=SIZEOFINT*2+1) test3  !is real(8)
!    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_traj) test4 ! attention is integer(8)
!    CLOSE(117)

    end if

#else
    dummy = 1 !amount of processors
    fname = 'largest_particle'
    write (str, *) nitera_last; str = trim(adjustl(fname))//"."//trim(adjustl(str)) ! name with the number for direction and scalar
    open (unit=15, file=str, access='stream', form='unformatted')
    inquire (UNIT=15, POS=particle_pos) !would be 1
    write (15) dummy  !only 1
    inquire (UNIT=15, POS=particle_pos) !would be 5
    write (15) isize_traj  !header
    inquire (UNIT=15, POS=particle_pos)  !would be 9
    write (15) big_part
    inquire (UNIT=15, POS=particle_pos)  !would be 409 with 50 numbers
    write (15) tag_big_part
    close (15)

!    !Just for testing and as template
!    OPEN(unit=117, file=str, access='stream', form='unformatted')
!    READ(117) test1   !is a integer
!    READ(117, POS=SIZEOFINT+1) test2  !is an integer
!    READ(117, POS=SIZEOFINT*2+1) test3  !is real(8)
!    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_traj) test4 ! attention is integer(8)
!    CLOSE(117)
!    print*, test1
!    print*, test2
!    print*, test3
!    print*, test4

#endif

!CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
!print*, 'here2', ims_pro

    call TLAB_STOP(0)
end program
