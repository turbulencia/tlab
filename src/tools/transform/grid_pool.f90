#include "types.h"

! ###################################################################
! ###################################################################
SUBROUTINE RD_GRIDSIZE(name, imax,jmax,kmax, scalex,scaley,scalez)

  IMPLICIT NONE

  CHARACTER*(*) name
  TINTEGER imax, jmax, kmax
  TREAL scalex, scaley, scalez

  TREAL scale(3)

  OPEN(50,file=name, status='old',form='unformatted')
  REWIND(50)
  READ(50) imax,jmax,kmax
  READ(50) scale
  scalex = scale(1)
  scaley = scale(2)
  scalez = scale(3)
  CLOSE(50)

  RETURN
END SUBROUTINE RD_GRIDSIZE

! ###################################################################
! ###################################################################
SUBROUTINE TRANS_DROP_PLANES(nmax, a,  scale)

  IMPLICIT NONE

  TINTEGER nmax
  TREAL a(nmax), scale

  TINTEGER nplanes, option, n
  TREAL correction, tolerance

  tolerance = C_1EM10_R

  WRITE(*,'(20(/),"Drop planes",/)')
  WRITE(*,'(a)')   '1. Drop planes symmetrically'
  WRITE(*,'(a)')   '2. Drop planes at the beginning'
  WRITE(*,'(a)')   '3. Drop planes at the end'
  WRITE(*,'(a)')   '4. Drop intermediate planes'
  WRITE(*,'(a,/)') '9. Exit'
  READ(*,*) option 

  IF ( option .LT. 4 ) THEN
     WRITE(*,'(/,"Total number of planes to drop ? ", $)')
     READ(*,*) nplanes
     IF ( nplanes .GE. nmax ) THEN
        WRITE(*,*) 'Error. Trying to drop equal/more planes than exist.'
        STOP
     ENDIF
  ENDIF
  
  SELECT CASE(option)

  CASE(1) ! Drop planes symmetrically
     nplanes = nplanes /2
     correction = scale - (a(nmax)-a(1)) ! The variable correction takes care of the periodic case
     scale = a(nmax-nplanes) - a(nplanes+1) + correction
     nmax = nmax - 2*nplanes
     DO n = 1,nmax
        a(n) = a(n+nplanes)
     ENDDO

  CASE(2) ! Drop planes at the beginning
     correction = scale - (a(nmax)-a(1))
     scale = a(nmax) - a(nplanes+1) + correction
     DO n = nplanes+1, nmax
        a(n-nplanes) = a(n)
     ENDDO
     nmax = nmax - nplanes
     IF ( nmax .EQ. 1 ) scale = C_1_R

  CASE(3) ! Drop planes at the end
     correction = scale - (a(nmax)-a(1))
     scale = a(nmax-nplanes) - a(1) + correction
     nmax = nmax - nplanes
     IF ( nmax .EQ. 1 ) scale = C_1_R

  CASE(4) ! Drop intermediate planes
     correction = scale - (a(nmax)-a(1))
     DO n = 1,nmax/2
        a(n) = a(2*n-1)
     ENDDO
     nmax = nmax/2
     IF (correction .LT. tolerance) THEN
        DO n = 1,nmax
           a(n) = a(1) + (a(n)-a(1))/(a(nmax)-a(1)) * scale
        ENDDO
     ENDIF

  END SELECT

  RETURN
END SUBROUTINE TRANS_DROP_PLANES

! ###################################################################
! ###################################################################
SUBROUTINE TRANS_ADD_PLANES(nmax, a, scale)

  IMPLICIT NONE

  TINTEGER nmax
  TREAL a(2*nmax), scale

  TINTEGER nplanes, option, n, inipos
  TREAL correction, tolerance, deltaup, deltadown, posnew

  tolerance = C_1EM10_R

  WRITE(*,'(20(/),"Introduce planes and exit",/)')
  WRITE(*,'(a)')   '1. Introduce planes symmetrically'
  WRITE(*,'(a)')   '2. Introduce planes at the beginning'
  WRITE(*,'(a)')   '3. Introduce planes at the end'
  WRITE(*,'(a)')   '4. Introduce planes at midpositions'
  WRITE(*,'(a)')   '5. Introduce zone of planes'
  WRITE(*,'(a)')   '6. Introduce particular planes by stdin'
  WRITE(*,'(a,/)') '9. Exit'
  READ(*,*) option 

  IF ( option .LT. 4 .OR. option .EQ. 5 ) THEN
     WRITE(*,'(/,"Total number of planes to introduce? ", $)')
     READ(*,*) nplanes
     IF ( nplanes .GT. nmax ) THEN
        WRITE(*,*) 'Error: Not enough space in grid arrays.'
        STOP
     ENDIF
  ENDIF
  
  SELECT CASE(option)

  CASE(1) ! Introduce planes symmetrically
     nplanes = nplanes /2
     WRITE(*,'(/,"Initial position ? ", $)')
     READ(*,*) deltadown
     WRITE(*,'(/,"Final position ? ", $)')
     READ(*,*) deltaup

     correction = scale - a(nmax) ! The variable correction takes care of the periodic case
     IF ( deltaup .GT. scale ) THEN; deltaup = (deltaup-scale) / M_REAL(nplanes) 
     ELSE;                           deltaup = a(nmax) - a(nmax-1)
     ENDIF
     IF ( deltadown .GT. a(1) ) THEN; deltadown = a(2) - a(1)
     ELSE;                            deltadown = (a(1)-deltadown) / M_REAL(nplanes)
     ENDIF
     DO n = nmax,1,-1
        a(n+nplanes) = a(n)
     ENDDO
     DO n = 1,nplanes
        a(nmax+nplanes+n) = a(nmax+nplanes+n-1) + deltaup
        a(nplanes+1-n)    = a(nplanes+1-n+1)    - deltadown
     ENDDO
     nmax = nmax + 2*nplanes
     scale = a(nmax) - a(1) + correction 
     
  CASE(2) ! Introduce planes at the beginning
     WRITE(*,'(/,"Initial position ? ", $)')
     READ(*,*) deltadown

     correction = scale - a(nmax)
     IF ( deltadown .GT. a(1) ) THEN; deltadown = a(2) - a(1)
     ELSE;                            deltadown = (a(1)-deltadown) / M_REAL(nplanes)
     ENDIF
     DO n = nmax,1,-1
        a(n+nplanes) = a(n)
     ENDDO
     DO n = 1,nplanes
        a(nplanes+1-n)    = a(nplanes+1-n+1)    - deltadown
     ENDDO
     nmax = nmax + nplanes
     scale = a(nmax) - a(1) + correction 

  CASE(3) ! Introduce planes at the end
     WRITE(*,'(/,"Final position ? ", $)')
     READ(*,*) deltaup

     correction = scale - a(nmax)
     IF ( deltaup .GT. scale ) THEN; deltaup = (deltaup-scale) / M_REAL(nplanes) 
     ELSE;                           deltaup = a(nmax) - a(nmax-1)
     ENDIF
     DO n = 1,nplanes
        a(nmax+n) = a(nmax+n-1) + deltaup
     ENDDO
     nmax = nmax + nplanes
     scale = a(nmax) - a(1) + correction 

  CASE(4) ! Introduce planes at midpositions
     correction = scale - a(nmax)
     DO n = nmax,1,-1
        a(2*n-1) = a(n)
     ENDDO
     DO n = 2,2*nmax-2,2
        a(n) = (a(n+1)+a(n-1))*C_05_R
     ENDDO
     a(2*nmax) = a(2*nmax-1) + (a(2*nmax-1)-a(2*nmax-2))
     nmax = 2*nmax
     IF ( correction .LT. tolerance ) THEN
        DO n = 1,nmax
           a(n) = a(1) + (a(n)-a(1))/(a(nmax)-a(1)) * scale
        ENDDO
     ENDIF

  CASE(5) ! Introduce zone of planes
     WRITE(*,'(/,"Initial position ? ", $)')
     READ(*,*) inipos

     correction = scale - a(nmax)
     deltadown = a(inipos) - a(inipos-1)
     DO n = nmax,inipos,-1
        a(n+nplanes) = a(n) + deltadown*M_REAL(nplanes)
     ENDDO
     DO n = 1,nplanes
        a(inipos+n-1) = a(inipos+n-2) + deltadown
     ENDDO
     nmax = nmax + nplanes
     scale = a(nmax) - a(1) + correction 

  CASE(6) ! Introduce particular planes by stdin
     WRITE(*,'(/,"Plane position ? ", $)') ! to be written in terms of a list_real
     READ(*,*) posnew
     
     correction = scale - a(nmax)
     IF ( correction .GT. tolerance ) THEN
        WRITE(*,*) "ERROR 2: Periodic direction. Use other option."
        STOP
     ENDIF
     DO n = nmax,1,-1
        IF ( a(n) .GT. posnew ) THEN
           a(n+1) = a(n)
        ELSE
           a(n+1) = posnew
           EXIT
        ENDIF
     ENDDO
     IF ( a(1) .GT. posnew ) THEN
        a(1) = posnew
     ENDIF
     nmax = nmax + 1
     scale = a(nmax) - a(1) 

  END SELECT

  RETURN

END SUBROUTINE TRANS_ADD_PLANES

! ###################################################################
! ###################################################################
SUBROUTINE TRANS_DATA(name, grid, work1,work2)
  
  USE DNS_TYPES, ONLY : grid_dt
  IMPLICIT NONE

  TYPE(grid_dt) grid
  
  CHARACTER*(*) name
  TREAL work1(*), work2(*)

  TREAL dxmx, dxmn, axmx, axmn
  TINTEGER n

  OPEN(20,file=name,status='unknown',position='append')

  WRITE(20,'(a13)') '['//TRIM(ADJUSTL(grid%name))//'-direction]'

  IF ( grid%size .GT. 1 ) THEN

     IF ( MOD(grid%size,2) .NE. 0 ) THEN
        WRITE(20,*) 'Warning. Not an even number of points.'
     ENDIF
     
     work1(2) = grid%nodes(2)-grid%nodes(1)
     DO n = 3,grid%size
        work1(n) = grid%nodes(n)-grid%nodes(n-1)
        work2(n) = work1(n)/work1(n-1)
     ENDDO
     dxmx = MAXVAL(work1(2:grid%size)); dxmn = MINVAL(work1(2:grid%size))
     axmx = MAXVAL(work2(3:grid%size)); axmn = MINVAL(work2(3:grid%size))
     
     WRITE(20,'(a25,i5)')    'number of points .......: ',grid%size
     WRITE(20,'(a25,e12.5)') 'origin .................: ',grid%nodes(1)
     WRITE(20,'(a25,e12.5)') 'end point ..............: ',grid%nodes(grid%size)
     WRITE(20,'(a25,e12.5)') 'scale ..................: ',grid%scale
     WRITE(20,'(a25,e12.5)') 'minimum step ...........: ',dxmn
     WRITE(20,'(a25,e12.5)') 'maximum step ...........: ',dxmx
     WRITE(20,'(a25,e12.5)') 'minimum stretching .....: ',axmn
     WRITE(20,'(a25,e12.5)') 'maximum stretching .....: ',axmx

  ELSE
     WRITE(20,'(a7)') '2D case'

  ENDIF

  CLOSE(20)

  RETURN
END SUBROUTINE TRANS_DATA

