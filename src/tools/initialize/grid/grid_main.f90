#include "types.h"

!# Grid generation tool. Origin is set always to (0,0,0)
PROGRAM INIGRID

  USE TLAB_TYPES, ONLY     : grid_dt
  USE TLAB_CONSTANTS, ONLY : gfile, ifile, lfile
  USE TLAB_PROCS
  USE GRID_LOCAL
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_pro
#endif
  IMPLICIT NONE

  CHARACTER*32  sfile
  TYPE(grid_dt)              :: g(3)
  TREAL, ALLOCATABLE, TARGET :: x(:),y(:),z(:)
  TREAL, ALLOCATABLE         :: work1(:),work2(:)
  TINTEGER idir, iseg, isize_wrk1d, n
  TREAL scale_old, scale_new

  ! #######################################################################
  ! Initialize
  ! #######################################################################
  sfile = TRIM(ADJUSTL(gfile))//'.sts'

  g(1)%name = 'x'
  g(2)%name = 'y'
  g(3)%name = 'z'

  CALL TLAB_START()

  DO idir = 1,3

    CALL GRID_READ_LOCAL(ifile, idir, g(idir)%scale, g(idir)%periodic)

    ! Calculate total number of points
    g(idir)%size = g_build(idir)%SIZE(1)
    DO iseg = 2,g_build(idir)%nseg
      g(idir)%size = g(idir)%size + g_build(idir)%SIZE(iseg) - 1
    ENDDO
    IF ( g_build(idir)%mirrored ) g(idir)%size = 2*g(idir)%size - 2

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

    SELECT CASE(g_build(idir)%opts(1,1))

    CASE(:4)
      CALL BLD_GEN( idir, g(idir)%nodes, g(idir)%size, g(idir)%scale)

    CASE( GTYPE_TANH )
      CALL BLD_TANH(idir, g(idir)%nodes, g(idir)%size, g(idir)%scale, work1)

    CASE( GTYPE_EXP )
      CALL BLD_EXP( idir, g(idir)%nodes, g(idir)%size, g(idir)%scale)

    END SELECT

    IF ( g_build(idir)%mirrored ) CALL GRID_MIRROR(g(idir)%size, g(idir)%nodes, g(idir)%scale)

    IF ( g(idir)%periodic ) g(idir)%size = g(idir)%size - 1

    IF ( g_build(idir)%fixed_scale > C_0_R ) THEN               ! rescale on exact fixed value
      scale_new     =  g_build(idir)%fixed_scale; scale_old = g(idir)%scale
      g(idir)%nodes = (g(idir)%nodes / scale_old) * scale_new   ! rescale nodes
      g(idir)%nodes(g(idir)%size) =                 scale_new   ! avoid rounding error
      g(idir)%scale =                               scale_new   ! update scale
    ENDIF

  ENDDO

  ! #######################################################################
  ! Statistics
  ! #######################################################################
#ifdef USE_MPI
  IF ( ims_pro == 0 ) THEN
#endif
    OPEN(20,file=sfile)

    DO idir = 1,3
      WRITE(20,3000) '['//TRIM(ADJUSTL(g(idir)%name))//'-direction]'

      IF ( g(idir)%size > 1 ) THEN
        work1(2) = g(idir)%nodes(2) - g(idir)%nodes(1)
        DO n = 3,g(idir)%size
          work1(n) = g(idir)%nodes(n)-g(idir)%nodes(n-1)
          work2(n) = work1(n)/work1(n-1)
        ENDDO

        WRITE(20,2000) 'number of points .......: ',g(idir)%size
        WRITE(20,1000) 'origin .................: ',g(idir)%nodes(1)
        WRITE(20,1000) 'end point ..............: ',g(idir)%nodes(g(idir)%size)
        WRITE(20,1000) 'scale ..................: ',g(idir)%scale
        WRITE(20,1000) 'minimum step ...........: ',MINVAL(work1(2:g(idir)%size))
        WRITE(20,1000) 'maximum step ...........: ',MAXVAL(work1(2:g(idir)%size))
        WRITE(20,1000) 'minimum stretching .....: ',MINVAL(work2(3:g(idir)%size))
        WRITE(20,1000) 'maximum stretching .....: ',MAXVAL(work2(3:g(idir)%size))

      ELSE
        WRITE(20,'(a7)') '2D case'

      ENDIF

    ENDDO

    CLOSE(20)

    ! #######################################################################
    ! Writing data
    ! #######################################################################
    CALL TLAB_WRITE_ASCII(lfile, 'Writing grid.')
    CALL IO_WRITE_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, g(1)%nodes,g(2)%nodes,g(3)%nodes)

#ifdef USE_MPI
  ENDIF
#endif

  CALL TLAB_STOP(0)

1000 FORMAT(a25,e12.5)
2000 FORMAT(a25,i5)
3000 FORMAT(a13)

CONTAINS
  ! #######################################################################
  ! #######################################################################
  SUBROUTINE GRID_MIRROR(imax, x, scale)
    IMPLICIT NONE

    TINTEGER, INTENT(IN   ) :: imax
    TREAL,    INTENT(INOUT) :: x(imax), scale

    ! -----------------------------------------------------------------------
    TINTEGER i
    TREAL offset

    ! #######################################################################
    ! Offset for even number of points
    offset = (x(imax/2+1)-x(imax/2)) / C_2_R
    DO i = imax/2,imax
      x(i) = x(i) - offset
    ENDDO

    ! Mirroring
    DO i = 1,imax/2-1
      x(i) = - x(imax+1-i)
    ENDDO

    ! Global translation
    offset = x(1)
    x = x -offset

    scale = x(imax)-x(1)

    RETURN
  END SUBROUTINE GRID_MIRROR

END PROGRAM INIGRID
