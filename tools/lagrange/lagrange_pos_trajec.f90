#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "LAGRANGE_POS_TRAJEC"

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
PROGRAM LAGRANGE_POS_TRAJEC

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------

#ifdef USE_MPI
  TINTEGER  ierr, i, j, k, particle_pos

  TINTEGER  dummy_ims_npro
  TINTEGER  dummy_isize_trajectory
  TINTEGER, DIMENSION(:), ALLOCATABLE :: dummy_proc, all_dummy_proc
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: l_trajectories_tags
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: l_tags
  TREAL, DIMENSION(:), ALLOCATABLE :: dummy_big_overall
  TREAL, DIMENSION(:,:),    ALLOCATABLE :: l_q, l_txc
  TREAL, DIMENSION(:,:), ALLOCATABLE :: l_trajectories, all_l_trajectories

  TINTEGER nitera_first, nitera_last

  CHARACTER*32 inifile
  CHARACTER*64 str,fname
  CHARACTER*128 line
  CHARACTER*32 bakfile

  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_part .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! Get the local information from the dns.ini
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Start','0',  nitera_first)
    

#include "dns_alloc_larrays.h"
  

  ALLOCATE(dummy_proc(isize_trajectory))
  ALLOCATE(all_dummy_proc(isize_trajectory)) 
  ALLOCATE(dummy_big_overall(isize_trajectory))
  ALLOCATE(l_trajectories_tags(isize_trajectory))
  ALLOCATE(l_trajectories(3,isize_trajectory))
  ALLOCATE(all_l_trajectories(3,isize_trajectory)) 

  l_trajectories(:,:) = C_0_R
  all_l_trajectories(:,:)=C_0_R
  dummy_proc(:) = C_0_R
  all_dummy_proc(:) = C_0_R
  

  !#######################################################################
  !READ THE (FIRST) FILE
  !#######################################################################
  WRITE(fname,*) nitera_first; fname = TRIM(ADJUSTL(tag_part))//TRIM(ADJUSTL(fname))
  CALL IO_READ_PARTICLE(fname, l_tags, l_q)

  IF (ims_pro .EQ. 0) THEN
    fname = 'largest_particle'
    WRITE(str,*) nitera_last;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
    OPEN(unit=117, file=str, access='stream', form='unformatted')
    READ(117) dummy_ims_npro   !is a integer
    READ(117, POS=SIZEOFINT+1) dummy_isize_trajectory  !is an integer
    READ(117, POS=SIZEOFINT*2+1) dummy_big_overall  !is real(8)
    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_trajectory) l_trajectories_tags ! attention is integer(8)
    CLOSE(117)
  ENDIF

  !#######################################################################
  !BROADCAST THE ID OF THE LARGEST PARTICLES
  !#######################################################################
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  CALL MPI_BCAST(l_trajectories_tags,isize_trajectory,MPI_INTEGER8,0,MPI_COMM_WORLD,ims_err)

  !#######################################################################
  !SEARCH FOR LARGEST PARTICLES
  !#######################################################################
  DO i=1,ims_size_p(ims_pro+1)
    DO j=1,isize_trajectory
      IF (l_tags(i) .EQ. l_trajectories_tags(j)) THEN
        l_trajectories(1,j)=l_q(i,1)
        l_trajectories(2,j)=l_q(i,2)
        l_trajectories(3,j)=l_q(i,3)
        dummy_proc(j)=ims_pro
      ENDIF
    ENDDO
  ENDDO 

  !#######################################################################
  !REDUCE ALL INFORMATION TO ROOT
  !#######################################################################
  CALL MPI_REDUCE(l_trajectories, all_l_trajectories, 3*isize_trajectory, MPI_REAL8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
  CALL MPI_REDUCE(dummy_proc, all_dummy_proc, isize_trajectory, MPI_INTEGER4, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)

  !#######################################################################
  !WRITE DATA WITH ROOT
  !#######################################################################

  IF (ims_pro .EQ. 0) THEN
    fname = 'pos_largest_particle_start'
    WRITE(str,*) nitera_first;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
    OPEN(unit=15, file=str, access='stream', form='unformatted')
    INQUIRE(UNIT=15, POS=particle_pos) !would be 1
    WRITE (15)  ims_npro  !header
    INQUIRE(UNIT=15, POS=particle_pos) !would be 5 
    WRITE (15)  isize_trajectory  !header
    INQUIRE(UNIT=15, POS=particle_pos)  !would be 9
    WRITE (15)  l_trajectories_tags
    INQUIRE(UNIT=15, POS=particle_pos)  !409
    WRITE (15)  all_l_trajectories
    INQUIRE(UNIT=15, POS=particle_pos)  !would be 1609 with 50 numbers
    WRITE (15)  all_dummy_proc
    CLOSE(15)
  ENDIF
!    !Just for testing and as template
!      IF (ims_pro .EQ. 0) THEN
!    OPEN(unit=117, file=str, access='stream', form='unformatted')
!    READ(117) test1   !is a integer
!    READ(117, POS=SIZEOFINT+1) test2  !is an integer
!    READ(117, POS=SIZEOFINT*2+1) test3  !is real(8)
!    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*isize_trajectory) test4 ! attention is integer(8)
!    CLOSE(117)
!    ENDIF

#endif

CALL DNS_END(0)

  STOP
END PROGRAM LAGRANGE_POS_TRAJEC
