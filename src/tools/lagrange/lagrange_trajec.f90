#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "LAGRANGE_TRAJEC"

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
PROGRAM LAGRANGE_TRAJEC

  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_pro, ims_npro
  USE TLAB_MPI_PROCS
#endif
  USE LAGRANGE_VARS
  USE LAGRANGE_ARRAYS

  IMPLICIT NONE
#include "integers.h"

! -------------------------------------------------------------------

  TINTEGER  i, j, k, particle_pos
  TREAL temp
  INTEGER(8) itemp
  LOGICAL :: swapped
  TREAL, DIMENSION(:), ALLOCATABLE :: big_part
#ifdef USE_MPI
  TREAL, DIMENSION(:), ALLOCATABLE :: big_all
  TREAL, DIMENSION(:), ALLOCATABLE :: big_overall
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: tag_big_all, tag_big_overall
#else
  TINTEGER dummy
#endif
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: tag_big_part

  TINTEGER nitera_last

  CHARACTER*64 str,fname
  CHARACTER*32 bakfile

!  TINTEGER test1, test2
!  TREAL test3(50)
!  INTEGER(8) test4(50)
  bakfile = TRIM(ADJUSTL(ifile))//'.bak'

  CALL TLAB_START()

  CALL DNS_READ_GLOBAL(ifile)
  IF ( icalc_part .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL(ifile)
  ENDIF
#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

! Get the local information from the dns.ini
  CALL SCANINIINT(bakfile, ifile, 'Iteration', 'End',        '0',  nitera_last )

  CALL PARTICLE_ALLOCATE(C_FILE_LOC)

  ALLOCATE(big_part(isize_trajectory))
  ALLOCATE(tag_big_part(isize_trajectory))
#ifdef USE_MPI
  IF (ims_pro .EQ. 0) THEN
     ALLOCATE(big_all(isize_trajectory*ims_npro))
     ALLOCATE(big_overall(isize_trajectory))
     ALLOCATE(tag_big_all(isize_trajectory*ims_npro))
     ALLOCATE(tag_big_overall(isize_trajectory))
  ENDIF

#endif


!#######################################################################
!READ THE (LAST) FILE
!#######################################################################
  WRITE(fname,*) nitera_last; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
  CALL IO_READ_PARTICLE(fname, l_g, l_q)

!#######################################################################
!Every processor searches for the largest particles
!#######################################################################

  big_part(1:isize_trajectory) = l_q(1:isize_trajectory,5)
  tag_big_part(1:isize_trajectory) = l_g%tags(1:isize_trajectory)

!  IF (ims_pro .EQ. 0) THEN
!  print*, tag_big_part
!  print*, big_part
!  ENDIF

  DO j= isize_trajectory-1,1,-1
     swapped = .FALSE.
     DO i = 1,j
        IF ( big_part(i) .GT. big_part(i+1)) THEN !if 49... bigger than 50...
!START SWAPPING
           temp = big_part(i)
           big_part(i) = big_part(i+1)
           big_part(i+1) = temp

           itemp = tag_big_part(i)
           tag_big_part(i) = tag_big_part(i+1)
           tag_big_part(i+1) = itemp
           swapped = .TRUE.
        ENDIF
     ENDDO
     IF (.NOT. swapped) EXIT
  ENDDO

!DO k=isize_trajectory+1,l_g%np,1
  DO k=1,l_g%np
     IF (l_q(k,5) .GT. big_part(1))THEN
        big_part(1) = l_q(k,5)
        tag_big_part(1) = l_g%tags(k)

        swapped = .FALSE.
        DO i = 1,isize_trajectory-1
           IF ( big_part(i) .GT. big_part(i+1)) THEN !if 49... bigger than 50...
!START SWAPPING
              temp = big_part(i)
              big_part(i) = big_part(i+1)
              big_part(i+1) = temp

              itemp = tag_big_part(i)
              tag_big_part(i) = tag_big_part(i+1)
              tag_big_part(i+1) = itemp
              swapped = .TRUE.
           ENDIF
           IF (.NOT. swapped) EXIT
        ENDDO

     ENDIF
  ENDDO

!  IF (ims_pro .EQ. 0) THEN
!  print*, 'delimiter'
!  print*, tag_big_part
!  print*, big_part
!  END IF

!#######################################################################
!Send all particles to one processor
!#######################################################################
#ifdef USE_MPI
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  CALL MPI_GATHER(big_part, isize_trajectory, MPI_REAL8, big_all, isize_trajectory, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_GATHER(tag_big_part, isize_trajectory, MPI_INTEGER8, tag_big_all, isize_trajectory, MPI_INTEGER8, 0, MPI_COMM_WORLD, ims_err)


  IF (ims_pro .EQ. 0) THEN

!#######################################################################
!Determine the oerall biggest 50 of the biggest 50 of each processor
!#######################################################################

     big_overall=big_part
     tag_big_overall=tag_big_part
     DO k=isize_trajectory+1,isize_trajectory*ims_npro,1
        IF (big_all(k) .GT. big_overall(1))THEN
           big_overall(1) = big_all(k)
           tag_big_overall(1) = tag_big_all(k)

           swapped = .FALSE.
           DO i = 1,isize_trajectory-1
              IF ( big_overall(i) .GT. big_overall(i+1)) THEN !if 49... bigger than 50...
!START SWAPPING
                 temp = big_overall(i)
                 big_overall(i) = big_overall(i+1)
                 big_overall(i+1) = temp

                 itemp = tag_big_overall(i)
                 tag_big_overall(i) = tag_big_overall(i+1)
                 tag_big_overall(i+1) = itemp
                 swapped = .TRUE.
              ENDIF
              IF (.NOT. swapped) EXIT
           ENDDO

        ENDIF
     ENDDO



!#######################################################################
!Write data file with root processor
!#######################################################################

     fname = 'largest_particle'
     WRITE(str,*) nitera_last;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
     OPEN(unit=15, file=str, access='stream', form='unformatted')
     INQUIRE(UNIT=15, POS=particle_pos) !would be 1
     WRITE (15)  ims_npro  !header
     INQUIRE(UNIT=15, POS=particle_pos) !would be 5
     WRITE (15)  isize_trajectory  !header
     INQUIRE(UNIT=15, POS=particle_pos)  !would be 9
     WRITE (15)  big_overall
     INQUIRE(UNIT=15, POS=particle_pos)  !would be 409 with 50 numbers
     WRITE (15)  tag_big_overall
     CLOSE(15)

!    !Just for testing and as template
!    OPEN(unit=117, file=str, access='stream', form='unformatted')
!    READ(117) test1   !is a integer
!    READ(117, POS=SIZEOFINT+1) test2  !is an integer
!    READ(117, POS=SIZEOFINT*2+1) test3  !is real(8)
!    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_trajectory) test4 ! attention is integer(8)
!    CLOSE(117)


  ENDIF

#else
  dummy=1 !amount of processors
  fname = 'largest_particle'
  WRITE(str,*) nitera_last;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
  OPEN(unit=15, file=str, access='stream', form='unformatted')
  INQUIRE(UNIT=15, POS=particle_pos) !would be 1
  WRITE (15)  dummy  !only 1
  INQUIRE(UNIT=15, POS=particle_pos) !would be 5
  WRITE (15)  isize_trajectory  !header
  INQUIRE(UNIT=15, POS=particle_pos)  !would be 9
  WRITE (15)  big_part
  INQUIRE(UNIT=15, POS=particle_pos)  !would be 409 with 50 numbers
  WRITE (15)  tag_big_part
  CLOSE(15)


!    !Just for testing and as template
!    OPEN(unit=117, file=str, access='stream', form='unformatted')
!    READ(117) test1   !is a integer
!    READ(117, POS=SIZEOFINT+1) test2  !is an integer
!    READ(117, POS=SIZEOFINT*2+1) test3  !is real(8)
!    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_trajectory) test4 ! attention is integer(8)
!    CLOSE(117)
!    print*, test1
!    print*, test2
!    print*, test3
!    print*, test4

#endif



!CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
!print*, 'here2', ims_pro


  CALL TLAB_STOP(0)
END PROGRAM
