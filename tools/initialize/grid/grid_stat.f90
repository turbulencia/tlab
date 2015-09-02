#include "types.h"

SUBROUTINE GRID_STAT(name, imax,jmax,kmax, x,y,z, scalex,scaley,scalez, work1,work2 )
      
  USE GRID_LOCAL

  IMPLICIT NONE

  CHARACTER*(*) name
  TINTEGER imax, jmax, kmax
  TREAL x(*), y(*), z(*), work1(*), work2(*)
  TREAL scalex, scaley, scalez

  TREAL dxmx, dxmn, axmx, axmn
  TREAL dymx, dymn, aymx, aymn
  TREAL dzmx, dzmn, azmx, azmn
  TINTEGER n !, idir, iseg
!  TREAL dummy1, dummy2, h_ini

! ##################
! # Initialization #
! ##################

  OPEN(20,file=name)

! ###############
! # x-direction #
! ###############

  WRITE(20,3000) '[x-direction]'

  IF (imax.gt.1) THEN

     work1(2) = x(2)-x(1)
     DO n = 3,imax
        work1(n) = x(n)-x(n-1)
        work2(n) = work1(n)/work1(n-1)
     ENDDO
     dxmx = MAXVAL(work1(2:imax)); dxmn = MINVAL(work1(2:imax))
     axmx = MAXVAL(work2(3:imax)); axmn = MINVAL(work2(3:imax))

     WRITE(20,2000) 'number of points .......: ',imax
     WRITE(20,1000) 'origin .................: ',x(1)
     WRITE(20,1000) 'end point ..............: ',x(imax)
     WRITE(20,1000) 'scale ..................: ',scalex
     WRITE(20,1000) 'minimum step ...........: ',dxmn
     WRITE(20,1000) 'maximum step ...........: ',dxmx
     WRITE(20,1000) 'minimum stretching .....: ',axmn
     WRITE(20,1000) 'maximum stretching .....: ',axmx

  ELSE

     WRITE(20,'(a7)') '2D case'

  ENDIF

! ###############
! # y-direction #
! ###############


  WRITE(20,3000) '[y-direction]'

  IF (jmax.gt.1) THEN

     work1(2) = y(2)-y(1)
     DO n = 3,jmax
        work1(n) = y(n)-y(n-1)
        work2(n) = work1(n)/work1(n-1)
     ENDDO
     dymx = MAXVAL(work1(2:jmax)); dymn = MINVAL(work1(2:jmax))
     aymx = MAXVAL(work2(3:jmax)); aymn = MINVAL(work2(3:jmax))

     WRITE(20,2000) 'number of points .......: ',jmax
     WRITE(20,1000) 'origin .................: ',y(1)
     WRITE(20,1000) 'end point ..............: ',y(jmax)
     WRITE(20,1000) 'scale ..................: ',scaley
     WRITE(20,1000) 'minimum step ...........: ',dymn
     WRITE(20,1000) 'maximum step ...........: ',dymx
     WRITE(20,1000) 'minimum stretching .....: ',aymn
     WRITE(20,1000) 'maximum stretching .....: ',aymx

  ELSE

     WRITE(20,'(a7)') '2D case'

  ENDIF

! ###############
! # z-direction #
! ###############

  WRITE(20,3000) '[z-direction]'

  IF (kmax.gt.1) THEN

     work1(2) = z(2)-z(1)
     DO n = 3,kmax
        work1(n) = z(n)-z(n-1)
        work2(n) = work1(n)/work1(n-1)
     ENDDO
     dzmx = MAXVAL(work1(2:kmax)); dzmn = MINVAL(work1(2:kmax))
     azmx = MAXVAL(work2(3:kmax)); azmn = MINVAL(work2(3:kmax))

     WRITE(20,2000) 'number of points .......: ',kmax
     WRITE(20,1000) 'origin .................: ',z(1)
     WRITE(20,1000) 'end point ..............: ',z(kmax)
     WRITE(20,1000) 'scale ..................: ',scalez
     WRITE(20,1000) 'minimum step ...........: ',dzmn
     WRITE(20,1000) 'maximum step ...........: ',dzmx
     WRITE(20,1000) 'minimum stretching .....: ',azmn
     WRITE(20,1000) 'maximum stretching .....: ',azmx

  ELSE

     WRITE(20,'(a7)') '2D case'

  ENDIF

! ###########################
! # Geometrical progression #
! ###########################

!  WRITE(20,*) ' '
!  WRITE(20,4000) '[Geometric progression information]'
!
!  DO idir = 1,3
!     WRITE(20,'(a12,i1)') 'Direction # ',idir
!     DO iseg = 1,idir_opts(1,idir)
!        IF ( val1(iseg,idir) .NE. C_1_R ) THEN
!           val0(iseg,idir) = (C_1_R-val1(iseg,idir)**FLOAT(isegdim(iseg,idir)-1))/&
!                (C_1_R-val1(iseg,idir))
!        ELSE
!           val0(iseg,idir) = FLOAT(isegdim(iseg,idir)-1)
!        ENDIF
!     ENDDO
!     
!     dummy1 = C_0_R
!     DO iseg = 1,idir_opts(1,idir)
!        dummy2 = val0(iseg,idir)
!        DO n = 1,iseg-1
!           dummy2 = dummy2*val1(n,idir)**FLOAT(isegdim(n,idir)-1)
!        ENDDO
!        dummy1 = dummy1 + dummy2
!     ENDDO
!     
!     IF ( dummy1 .GT. C_0_R ) THEN
!        h_ini = isegend(idir_opts(1,idir),idir) / dummy1
!     ENDIF
!     
!! dummy1 for the step size at the begining of the segment
!! dummy2 foR the length
!! -------------------------------------------------------
!
!     dummy1 = h_ini
!     dummy2 = dummy1*val0(1,idir)
!     WRITE(20,6000) 1, dummy2
!     DO iseg = 2,idir_opts(1,idir)
!        dummy1 = dummy1*val1(iseg-1,idir)**FLOAT(isegdim(iseg-1,idir)-1)
!        dummy2 = dummy2 + dummy1*val0(iseg,idir) 
!        WRITE(20,6000) iseg, dummy2
!     ENDDO
!     
!  ENDDO

  CLOSE(20)

  RETURN

! ###########
! # Formats #
! ###########

1000 FORMAT(a25,e12.5)
2000 FORMAT(a25,i5)
3000 FORMAT(a13)
!4000 FORMAT(a35)
!6000 FORMAT('Segment ',i2,' .............:',e12.5)

END SUBROUTINE GRID_STAT
