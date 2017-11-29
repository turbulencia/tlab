#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "INIPART"

!########################################################################
!# HISTORY
!#
!# 2013/10/28 - L. Muessle
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# Create an initial condition for the particle field
!#
!########################################################################
PROGRAM INIPART

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro, ims_npro
#endif
  
  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
  TINTEGER  ierr,isize_wrk3d, i

  TREAL, DIMENSION(:,:),    ALLOCATABLE,SAVE,TARGET :: x,y,z
  TREAL, DIMENSION(:),      ALLOCATABLE             :: wrk1d,wrk2d, wrk3d
  TREAL, DIMENSION(:,:),    ALLOCATABLE             :: txc

  TREAL, DIMENSION(:,:),    ALLOCATABLE             :: l_q, l_txc, l_hq
  INTEGER(8), DIMENSION(:), ALLOCATABLE             :: l_tags

  CHARACTER*32 inifile
  CHARACTER*64 str, line

#ifdef USE_MPI
  TINTEGER  particle_number_each
  TLONGINTEGER dummy_int1, dummy_int2, partcile_offset
#endif
  
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL(inifile)
  ELSE
     CALL DNS_END(0)
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  inb_particle_txc = 0 ! so far, not needed

#include "dns_alloc_larrays.h"
  isize_wrk3d = imax*jmax*kmax
  IF (jmax_part .EQ. 1) THEN
     jmax_part   = jmax ! 1 by default
  ENDIF

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------      
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d))
  ALLOCATE(wrk3d(isize_wrk3d))

  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !Allocte memory to read fields
     ALLOCATE(txc(isize_field,3))
     ALLOCATE(l_hq(isize_particle,inb_particle)) !Rubish information. Just to run FIELD_TO_PARTICLE properly
  ENDIF

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

  CALL PARTICLE_RANDOM_POSITION(l_q,l_hq,l_tags,isize_wrk3d,wrk1d,wrk2d,wrk3d,txc)

  CALL DNS_WRITE_PARTICLE('particle.ics',l_q)

#ifdef USE_MPI
  particle_number_each = INT( particle_number /INT(ims_npro, KIND=8) ) 

  dummy_int1 = INT(ims_pro, KIND=8)
  dummy_int2 = INT(particle_number_each, KIND=8)

  partcile_offset = dummy_int1 *dummy_int2

  DO i = 1,particle_number_each
     l_tags(i) = INT(i, KIND=8) +partcile_offset
  END DO
  CALL DNS_WRITE_PARTICLE_TAGS('particle.ics.id',l_tags)
#else


  DO i=1,particle_number
     l_tags(i) = INT(i, KIND=8)
  END DO
  CALL DNS_WRITE_PARTICLE_TAGS('particle.ics.id',l_tags)

#endif

  CALL DNS_END(0)

  STOP
END PROGRAM INIPART
