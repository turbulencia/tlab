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
!# particle_number in dns.ini. Afterwards particles are exchanged with
!# position of the largest particles determined before with l_trajec.x
!# and l_pos_trajec.x
!#
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM INI_TRAJEC

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

  TINTEGER  ierr,isize_wrk3d, i, dummy, particle_pos
  TREAL, DIMENSION(:,:),    ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:),      ALLOCATABLE :: wrk1d,wrk2d, wrk3d
  TREAL, DIMENSION(:,:),    ALLOCATABLE :: txc
  
  TREAL, DIMENSION(:,:),    ALLOCATABLE :: l_q, l_txc, l_hq
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: l_tags

  TINTEGER  dummy_ims_npro
  TINTEGER  dummy_num_trajectories
  TINTEGER, DIMENSION(:), ALLOCATABLE :: dummy_proc
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: l_trajectories_tags, fake_l_trajectories_tags
#ifdef USE_MPI
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: all_fake_l_trajectories_tags
#endif
  TREAL, DIMENSION(:,:), ALLOCATABLE :: l_trajectories
  TREAL, DIMENSION(:),   ALLOCATABLE :: fake_liquid, all_fake_liquid

  TREAL, DIMENSION(:,:), POINTER           :: dx, dy, dz

  TINTEGER nitera_first

  CHARACTER*32 inifile
  CHARACTER*64 str, line, fname
  CHARACTER*32 bakfile

#ifdef USE_MPI
  TINTEGER  particle_number_each
  TLONGINTEGER dummy_int1, dummy_int2 ,dummy_int3
#endif
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! Get the local information from the dns.ini
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Start',      '0',  nitera_first)

  
  inb_particle_txc = 0 ! so far, not needed

#include "dns_alloc_larrays.h"
  isize_wrk3d = imax*jmax*kmax
  IF (jmax_part .EQ. 1) THEN
     jmax_part   = jmax ! 1 by default
  ENDIF
  
! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------      
  ALLOCATE(x(imax_total,inb_grid))
  ALLOCATE(y(jmax_total,inb_grid))
  ALLOCATE(z(kmax_total,inb_grid))
  ! ALLOCATE(dx(imax_total,inb_grid))
  ! ALLOCATE(dy(jmax_total,inb_grid))
  ! ALLOCATE(dz(kmax_total,inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d))
  ALLOCATE(wrk3d(isize_wrk3d))

  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !Allocte memory to read fields
     ALLOCATE(txc(isize_field,3))
!     ALLOCATE(txc(isize_txc_field,6))
     ALLOCATE(l_hq(isize_particle,inb_particle)) !Rubish information. Just to run FIELD_TO_PARTICLE properly
  ENDIF
  ALLOCATE(dummy_proc(num_trajectories))
  ALLOCATE(l_trajectories_tags(num_trajectories))
  ALLOCATE(fake_l_trajectories_tags(num_trajectories))
  ALLOCATE(l_trajectories(3,num_trajectories))
  ALLOCATE(fake_liquid(num_trajectories))
  ALLOCATE(all_fake_liquid(num_trajectories))
#ifdef USE_MPI
  ALLOCATE(all_fake_l_trajectories_tags(num_trajectories))
  all_fake_l_trajectories_tags(:) = C_0_R
#endif
fake_l_trajectories_tags(:) = C_0_R
fake_liquid(:) = C_0_R
all_fake_liquid(:) = C_0_R


! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"


  !#######################################################################
  !CREATE THE RANDOM PARTICLE FIELD
  !#######################################################################
  CALL PARTICLE_RANDOM_POSITION(l_q,l_hq,l_tags,x,y,z,isize_wrk3d,wrk1d,wrk2d,wrk3d,txc)

  !#######################################################################
  !CREATE THE CORRESPONDING TAGS
  !#######################################################################
#ifdef USE_MPI
  particle_number_each=int(particle_number/INT(ims_npro, KIND=8)) 
  
  dummy_int1=INT(ims_pro, KIND=8)
  dummy_int2=particle_number_each

  dummy_int3=dummy_int1*dummy_int2

  DO i=1,particle_number_each
    l_tags(i)=INT(i, KIND=8)+dummy_int3
  END DO
#else

  DO i=1,particle_number
    l_tags(i)=INT(i, KIND=8)
  END DO
#endif

#ifdef USE_MPI
  !#######################################################################
  !READ THE LARGEST PARTICLE FILE
  !#######################################################################

  IF (ims_pro .EQ. 0) THEN
    fname = 'pos_largest_particle_start'
    WRITE(str,*) nitera_first;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str)) ! name with the number for direction and scalar
    OPEN(unit=117, file=str, access='stream', form='unformatted')
    READ(117) dummy_ims_npro   !is a integer
    READ(117, POS=SIZEOFINT+1) dummy_num_trajectories  !is an integer
    READ(117, POS=SIZEOFINT*2+1) l_trajectories_tags  !is INTEGER(8)
    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*num_trajectories) l_trajectories ! attention is integer(8)
    READ(117, POS=(SIZEOFINT*2+1)+SIZEOFREAL*num_trajectories+3*num_trajectories*SIZEOFREAL) dummy_proc ! attention is integer(8)
    CLOSE(117)
  ENDIF
  
  !#######################################################################
  !BROADCAST INFORMATION OF LARGEST PARTICLES
  !#######################################################################
  CALL MPI_BARRIER(MPI_COMM_WORLD,ims_err)
  CALL MPI_BCAST(l_trajectories_tags,num_trajectories,MPI_INTEGER8,0,MPI_COMM_WORLD,ims_err)
  CALL MPI_BCAST(l_trajectories,num_trajectories*3,MPI_REAL8,0,MPI_COMM_WORLD,ims_err)
  CALL MPI_BCAST(dummy_proc,num_trajectories,MPI_INTEGER4,0,MPI_COMM_WORLD,ims_err)

  !#######################################################################
  !REPLACE THE FIRST PARTICLES WITH THE LARGEST 
  !CORRESPONDING TO THE PROCESSORS
  !#######################################################################
  dummy=1
  DO i=1,num_trajectories
    IF (ims_pro .EQ. dummy_proc(i))THEN
      l_q(dummy,1)=l_trajectories(1,i)
      l_q(dummy,2)=l_trajectories(2,i)
      l_q(dummy,3)=l_trajectories(3,i)
      fake_l_trajectories_tags(i)=l_tags(dummy)
      fake_liquid(i)=l_q(dummy,5) ! just to be consistent with other output file
      dummy=dummy+1
    ENDIF
  ENDDO

  !#######################################################################
  !WRITE FILE WITH NEW FAKE LARGEST IDs
  !#######################################################################
  CALL MPI_REDUCE(fake_l_trajectories_tags, all_fake_l_trajectories_tags, num_trajectories, MPI_INTEGER8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)
  CALL MPI_REDUCE(fake_liquid, all_fake_liquid, num_trajectories, MPI_REAL8, MPI_SUM,0, MPI_COMM_WORLD, ims_err)

  !#######################################################################
  !WRITE FILE WITH NEW FAKE LARGEST IDs
  !#######################################################################
  IF (ims_pro .EQ. 0) THEN
    fname = 'fake_largest_particle'
    OPEN(unit=15, file=fname, access='stream', form='unformatted')
    INQUIRE(UNIT=15, POS=particle_pos) !would be 1
    WRITE (15)  ims_npro  !header
    INQUIRE(UNIT=15, POS=particle_pos) !would be 5 
    WRITE (15)  num_trajectories  !header
    INQUIRE(UNIT=15, POS=particle_pos)  !would be 9
    WRITE (15)  fake_liquid ! just to be consistent with other output file
    INQUIRE(UNIT=15, POS=particle_pos)  !would be 409 with 50 numbers
    WRITE (15)  all_fake_l_trajectories_tags
    CLOSE(15)
  ENDIF

  !#######################################################################
  !WRITE THE PARTICLE FILE
  !#######################################################################
  CALL DNS_WRITE_PARTICLE('fake_particle.ics',l_q)
  CALL DNS_WRITE_PARTICLE_TAGS('fake_particle_id.ics',l_tags)


#endif
CALL DNS_END(0)

  STOP
END PROGRAM INI_TRAJEC
