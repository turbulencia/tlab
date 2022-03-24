PROGRAM VEFILTER2

  USE TLAB_VARS
  USE IO_FIELDS

  IMPLICIT NONE

#include "types.h"
#include "integers.h"

  TREAL, DIMENSION(:,:), ALLOCATABLE :: x, y, z
  TREAL, DIMENSION(:),   ALLOCATABLE :: a, cx, cy, cz
  TREAL, DIMENSION(:,:), ALLOCATABLE :: wrk3d
  TINTEGER :: i

! ###################################################################
  CALL DNS_START

  CALL DNS_READ_GLOBAL('dns.ini')

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk3d(imax*jmax*kmax,2),a(imax*jmax*kmax))
  ALLOCATE(cx(imax*5),cy(jmax*5),cz(kmax_total*5))

! ###################################################################
  CALL IO_READ_GRID(gfile, imax,jmax,kmax_total, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z)

  ! CALL FLT4E_INI(g(1)%scale, x, cx)
  ! CALL FLT4E_INI(g(2)%scale, y, cy)
  ! CALL FLT4E_INI(g(3)%scale, z, cz)

  CALL IO_READ_FIELDS('field.inp', IO_SCAL, imax,jmax,kmax, i1,i0, i1, a, wrk3d)

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

  CALL IO_WRITE_FIELDS('field.out', IO_SCAL, imax,jmax,kmax, i1, i1, a, wrk3d)

  STOP
END PROGRAM VEFILTER2
