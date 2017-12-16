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
  TREAL, DIMENSION(:,:),    ALLOCATABLE             :: q,s,txc
  TREAL, DIMENSION(:),      ALLOCATABLE             :: wrk1d,wrk2d, wrk3d

  TREAL,      DIMENSION(:,:),   ALLOCATABLE, SAVE :: l_q, l_hq, l_txc
!  TREAL,      DIMENSION(:),     ALLOCATABLE, SAVE :: l_comm
  INTEGER(8), DIMENSION(:),     ALLOCATABLE, SAVE :: l_tags

  CHARACTER*32 inifile
  CHARACTER*64 str, line

#ifdef USE_MPI
  TINTEGER  particle_number_each
  TLONGINTEGER dummy_int1, dummy_int2, partcile_offset
#endif
  
!########################################################################
!########################################################################
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_part .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL(inifile)
  ELSE
     CALL DNS_END(0)
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Definitions
! -------------------------------------------------------------------
  inb_particle_txc = 0 ! so far, not needed

! #include "dns_alloc_larrays.h"
  IF ( jmax_part .EQ. 1 ) THEN
     jmax_part   = jmax ! 1 by default
  ENDIF

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------      
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

  inb_flow_array = 0
  inb_scal_array = 0
  isize_wrk3d    = imax*jmax*kmax
  inb_txc        = inb_scal
#include "dns_alloc_arrays.h"
#include "dns_alloc_larrays.h"
  IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4 ) THEN
     ALLOCATE(l_hq(isize_particle,inb_particle)) ! Aux space to run FIELD_TO_PARTICLE properly
  ENDIF
  
! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! -------------------------------------------------------------------
! Initialize particle information
! -------------------------------------------------------------------
  CALL PARTICLE_RANDOM_POSITION(l_q,l_hq,l_tags, txc, wrk1d,wrk2d,wrk3d)

#ifdef USE_MPI
  particle_number_each = INT( particle_number /INT(ims_npro, KIND=8) ) 

  dummy_int1 = INT(ims_pro, KIND=8)
  dummy_int2 = INT(particle_number_each, KIND=8)

  partcile_offset = dummy_int1 *dummy_int2

  DO i = 1,particle_number_each
     l_tags(i) = INT(i, KIND=8) +partcile_offset
  END DO

#else
  DO i=1,particle_number
     l_tags(i) = INT(i, KIND=8)
  END DO

#endif

  CALL IO_WRITE_PARTICLE(TRIM(ADJUSTL(tag_part))//'ics', l_tags, l_q)

  CALL DNS_END(0)

  STOP
END PROGRAM INIPART
