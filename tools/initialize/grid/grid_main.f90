!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - S. Stanley
!#              Created
!# 2000/06/01 - J.P. Mellado
!#              Modified
!#
!########################################################################
!# DESCRIPTION
!#
!# Grid generation tool
!# Origin is set always to (0,0,0)
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

PROGRAM INIGRID

  USE GRID_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER                         :: imax,jmax,kmax, nmax
  TREAL, DIMENSION(:), ALLOCATABLE :: x,y,z, work1,work2
  TINTEGER                         :: idir, iseg, isize_wrk1d

! #######################################################################
! Reading the input file
! #######################################################################
  ifile = 'dns.ini'
  ofile = 'dns.log'
  sfile = 'grid.sts'
  ffile = 'grid'

  CALL GRID_READ_LOCAL(ifile,i1) ! read data of Ox direction
  CALL GRID_READ_LOCAL(ifile,i2) ! read data of Oy direction
  CALL GRID_READ_LOCAL(ifile,i3) ! read data of Oz direction

! #######################################################################
! Calculating the dimensions
! #######################################################################
  DO idir = 1,3

! Add points of all segments
     nmax = isegdim(1,idir)
     DO iseg = 2,idir_opts(1,idir)
        nmax = nmax + isegdim(iseg,idir) - 1
     ENDDO

! Mirrored
     IF (idir_opts(3,idir).eq.1) nmax = 2*nmax - 2

! Store data
     IF ( idir .EQ. 1 ) imax = nmax
     IF ( idir .EQ. 2 ) jmax = nmax
     IF ( idir .EQ. 3 ) kmax = nmax

  ENDDO

! #######################################################################
! Allocation of arrays
! #######################################################################
  isize_wrk1d = MAX(imax,MAX(jmax,kmax))
  ALLOCATE(x(imax))
  ALLOCATE(y(jmax))
  ALLOCATE(z(kmax))
  ALLOCATE(work1(isize_wrk1d))
  ALLOCATE(work2(isize_wrk1d))

! #######################################################################
! Direction Ox
! #######################################################################
  IF      ( iseg_opts(1,1,1) .LE. 4 ) THEN; CALL BLD_GEN( i1,x,imax,scale(1))
  ELSE IF ( iseg_opts(1,1,1) .EQ. 5 ) THEN; CALL BLD_TANH(i1,x,imax,scale(1))
  ELSE IF ( iseg_opts(1,1,1) .EQ. 6 ) THEN; CALL BLD_EXP( i1,x,imax,scale(1)); ENDIF

  IF (idir_opts(2,1).eq.1) imax = imax - 1 ! periodic case

  IF ( idir_opts(3,1) .EQ. 1 ) THEN
     IF ( iseg_opts(1,1,1) .EQ. 5 .OR. iseg_opts(1,1,1) .EQ. 6 ) THEN; CALL GRID_MIRROR(i2,imax,x,scale(1))
     ELSE;                                                             CALL GRID_MIRROR(i1,imax,x,scale(1))
     ENDIF
  ENDIF

! #######################################################################
! Direction Oy
! #######################################################################
  IF      ( iseg_opts(1,1,2) .LE. 4 ) THEN; CALL BLD_GEN( i2,y,jmax,scale(2))
  ELSE IF ( iseg_opts(1,1,2) .EQ. 5 ) THEN; CALL BLD_TANH(i2,y,jmax,scale(2))
  ELSE IF ( iseg_opts(1,1,2) .EQ. 6 ) THEN; CALL BLD_EXP( i2,y,jmax,scale(2)); ENDIF

  IF (idir_opts(2,2).eq.1) jmax = jmax - 1 ! periodic case

  IF ( idir_opts(3,2) .EQ. 1 ) THEN
     IF ( iseg_opts(1,1,2) .EQ. 5 .OR. iseg_opts(1,1,2) .EQ. 6 ) THEN; CALL GRID_MIRROR(i2,jmax,y,scale(2))
     ELSE;                                                             CALL GRID_MIRROR(i1,jmax,y,scale(2))
     ENDIF
  ENDIF

! #######################################################################
! Direction Oz
! #######################################################################
  IF      ( iseg_opts(1,1,3) .LE. 4 ) THEN; CALL BLD_GEN( i3,z,kmax,scale(3))
  ELSE IF ( iseg_opts(1,1,3) .EQ. 5 ) THEN; CALL BLD_TANH(i3,z,kmax,scale(3))
  ELSE IF ( iseg_opts(1,1,3) .EQ. 6 ) THEN; CALL BLD_EXP( i3,z,kmax,scale(3)); ENDIF

  IF (idir_opts(2,3).eq.1) kmax = kmax - 1 ! periodic case

  IF ( idir_opts(3,3) .EQ. 1 ) THEN
     IF ( iseg_opts(1,1,3) .EQ. 5 .OR. iseg_opts(1,1,3) .EQ. 6 ) THEN; CALL GRID_MIRROR(i2,kmax,z,scale(3))
     ELSE;                                                             CALL GRID_MIRROR(i1,kmax,z,scale(3))
     ENDIF
  ENDIF

! #######################################################################
! Final step
! #######################################################################
! Statistics
  CALL GRID_STAT(sfile, imax,jmax,kmax, x,y,z, scale(1),scale(2),scale(3), work1,work2)

! Writing data
  CALL IO_WRITE_GRID(ffile, imax,jmax,kmax, scale(1),scale(2),scale(3), x,y,z)

  STOP

END PROGRAM INIGRID
