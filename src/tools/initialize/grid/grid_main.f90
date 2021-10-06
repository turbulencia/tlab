#include "types.h"

!# Grid generation tool. Origin is set always to (0,0,0)
PROGRAM INIGRID

  USE TLAB_TYPES, ONLY     : grid_dt
  USE TLAB_CONSTANTS, ONLY : gfile, ifile
  USE TLAB_PROCS
  USE GRID_LOCAL
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif
  IMPLICIT NONE

#include "integers.h"

  CHARACTER*32  sfile
  TYPE(grid_dt), DIMENSION(3)              :: g
  TINTEGER                                 :: nmax
  TREAL, DIMENSION(:), ALLOCATABLE, TARGET :: x,y,z
  TREAL, DIMENSION(:), ALLOCATABLE         :: work1,work2
  TINTEGER idir, iseg, isize_wrk1d, n
  TREAL dxmx, dxmn, axmx, axmn
  TREAL scale_old, scale_new

#ifdef USE_MPI
#include "mpif.h"
#endif

! #######################################################################
! Initialize
! #######################################################################
  sfile = TRIM(ADJUSTL(gfile))//'.sts'

  g(1)%name = 'x'
  g(2)%name = 'y'
  g(3)%name = 'z'

#ifdef USE_MPI
  CALL TLAB_START()
  IF ( ims_pro .EQ. 0 ) THEN
#endif

  DO idir = 1,3

! Read data
     CALL GRID_READ_LOCAL(ifile, idir, g(idir)%scale, g(idir)%periodic, g(idir)%fixed_scale)

! Add points of all segments
     nmax = isegdim(1,idir)
     DO iseg = 2,idir_opts(1,idir)
        nmax = nmax + isegdim(iseg,idir) - 1
     ENDDO

! Mirrored
     IF ( idir_opts(3,idir) .EQ. 1 ) nmax = 2*nmax - 2

! Size
     g(idir)%size = nmax

  ENDDO

! #######################################################################
! Allocation of arrays
! #######################################################################
  ALLOCATE(x(g(1)%size)); g(1)%nodes => x
  ALLOCATE(y(g(2)%size)); g(2)%nodes => y
  ALLOCATE(z(g(3)%size)); g(3)%nodes => z

  isize_wrk1d = MAX(g(1)%size,MAX(g(2)%size,g(3)%size))
  ALLOCATE(work1(isize_wrk1d))
  ALLOCATE(work2(isize_wrk1d))

! #######################################################################
! Construct grid
! #######################################################################
  DO idir = 1,3

     SELECT CASE(iseg_opts(1,1,idir))

     CASE(:4)
        CALL BLD_GEN( idir, g(idir)%nodes, g(idir)%size, g(idir)%scale)
        IF ( idir_opts(3,idir) .EQ. 1 ) CALL GRID_MIRROR(i1, g(idir)%size, g(idir)%nodes, g(idir)%scale)

     CASE(5)
        CALL BLD_TANH(idir, g(idir)%nodes, g(idir)%size, g(idir)%scale)
        IF ( idir_opts(3,idir) .EQ. 1 ) CALL GRID_MIRROR(i2, g(idir)%size, g(idir)%nodes, g(idir)%scale)

     CASE(6)
        CALL BLD_EXP( idir, g(idir)%nodes, g(idir)%size, g(idir)%scale)
        IF ( idir_opts(3,idir) .EQ. 1 ) CALL GRID_MIRROR(i2, g(idir)%size, g(idir)%nodes, g(idir)%scale)

     END SELECT

     IF ( g(idir)%periodic ) g(idir)%size = g(idir)%size - 1

     IF (g(idir)%fixed_scale .GT. C_0_R) THEN                    ! rescale on exact fixed value
        scale_new     =  g(idir)%fixed_scale; scale_old = g(idir)%scale
        g(idir)%nodes = (g(idir)%nodes / scale_old) * scale_new  ! rescale nodes
        g(idir)%nodes(g(idir)%size) =                 scale_new  ! avoid rounding error
        g(idir)%scale =                               scale_new  ! update scale
     ENDIF

  ENDDO

! #######################################################################
! Statistics
! #######################################################################
  OPEN(20,file=sfile)

  DO idir = 1,3
     WRITE(20,3000) '['//TRIM(ADJUSTL(g(idir)%name))//'-direction]'

     IF ( g(idir)%size .GT. 1 ) THEN
        work1(2) = g(idir)%nodes(2) - g(idir)%nodes(1)
        DO n = 3,g(idir)%size
           work1(n) = g(idir)%nodes(n)-g(idir)%nodes(n-1)
           work2(n) = work1(n)/work1(n-1)
        ENDDO
        dxmx = MAXVAL(work1(2:g(idir)%size)); dxmn = MINVAL(work1(2:g(idir)%size))
        axmx = MAXVAL(work2(3:g(idir)%size)); axmn = MINVAL(work2(3:g(idir)%size))

        WRITE(20,2000) 'number of points .......: ',g(idir)%size
        WRITE(20,1000) 'origin .................: ',g(idir)%nodes(1)
        WRITE(20,1000) 'end point ..............: ',g(idir)%nodes(g(idir)%size)
        WRITE(20,1000) 'scale ..................: ',g(idir)%scale
        WRITE(20,1000) 'minimum step ...........: ',dxmn
        WRITE(20,1000) 'maximum step ...........: ',dxmx
        WRITE(20,1000) 'minimum stretching .....: ',axmn
        WRITE(20,1000) 'maximum stretching .....: ',axmx

     ELSE
        WRITE(20,'(a7)') '2D case'

     ENDIF

  ENDDO

  CLOSE(20)

! #######################################################################
! Writing data
! #######################################################################
  CALL IO_WRITE_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, g(1)%nodes,g(2)%nodes,g(3)%nodes)

#ifdef USE_MPI
ENDIF
CALL TLAB_STOP(0)
#endif

STOP

1000 FORMAT(a25,e12.5)
2000 FORMAT(a25,i5)
3000 FORMAT(a13)

END PROGRAM INIGRID