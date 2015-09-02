PROGRAM VEFILTER2

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "types.h"  
#include "integers.h"

  TREAL, DIMENSION(:),   ALLOCATABLE :: x, y, z, dx, dy, dz
  TREAL, DIMENSION(:),   ALLOCATABLE :: a, cx, cy, cz
  TREAL, DIMENSION(:,:), ALLOCATABLE :: wrk3d
  TINTEGER :: i

! ###################################################################
  CALL DNS_INITIALIZE
  
  CALL DNS_READ_GLOBAL('dns.ini')

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
  ALLOCATE(x(imax_total))
  ALLOCATE(y(jmax_total))
  ALLOCATE(z(kmax_total))
  ALLOCATE(dx(imax_total*inb_grid))
  ALLOCATE(dy(jmax_total*inb_grid))
  ALLOCATE(dz(kmax_total*inb_grid))

  ALLOCATE(wrk3d(imax*jmax*kmax,2),a(imax*jmax*kmax))
  ALLOCATE(cx(imax*5),cy(jmax*5),cz(kmax_total*5))
  
! ###################################################################
  CALL IO_READ_GRID(gfile, imax,jmax,kmax_total, scalex,scaley,scalez, x,y,z)

  CALL FILT4E_INI(imax,       i1bc, scalex, x, cx)
  CALL FILT4E_INI(jmax,       j1bc, scaley, y, cy)
  CALL FILT4E_INI(kmax_total, k1bc, scalez, z, cz)

  CALL DNS_READ_FIELDS('field.inp', i1, imax,jmax,kmax, i1,i0, i1, a, wrk3d)

!  CALL OPR_FILTER(i4, imax, jmax, kmax,  i1bc, j1bc, k1bc, i1, i1, i1, i1, a, &
!       cx, cy, cz, wrk3d) 

  CALL OPR_FILTER(i3, imax, jmax, kmax, i1bc, j1bc, k1bc, i1, i1, i1, i1, a, &
       cx, cy, cz, wrk3d(1,1)) 
  CALL OPR_FILTER(i3, imax, jmax, kmax,  i1bc, j1bc, k1bc, i1, i1, i1, i1, wrk3d(1,1), &
       cx, cy, cz, wrk3d(1,2)) 
  DO i = 1,imax*jmax*kmax
     wrk3d(i,2) = wrk3d(i,2) + C_3_R*( a(i) - wrk3d(i,1) )
  ENDDO
  CALL OPR_FILTER(i3, imax, jmax, kmax,  i1bc, j1bc, k1bc, i1, i1, i1, i1, wrk3d(1,2), &
       cx, cy, cz, a) 

  CALL DNS_WRITE_FIELDS('field.out', i1, imax,jmax,kmax, i1, i1, a, wrk3d)

  STOP
END PROGRAM VEFILTER2
