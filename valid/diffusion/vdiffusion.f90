#include "types.h"  
#include "dns_const.h"

PROGRAM VDIFFUSION
  
  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:), ALLOCATABLE :: q, s, s_r
  TREAL, DIMENSION(:),   ALLOCATABLE :: wrk1d, wrk2d, wrk3d
  
  TINTEGER i, j, ij, iopt
  TINTEGER isize_wrk3d
  TREAL dummy, error, pi_loc, factor, wavenumber, x_loc
  CHARACTER*(32) fname

  TREAL, DIMENSION(:,:), POINTER :: dx, dy, dz

! ###################################################################
  CALL DNS_INITIALIZE
  
  CALL DNS_READ_GLOBAL('dns.ini')

  isize_wrk3d = isize_field

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))
  ALLOCATE(wrk3d(isize_wrk3d))
  ALLOCATE(  q(isize_field,3))
  ALLOCATE(  s(isize_field,1))
  ALLOCATE(s_r(isize_field,1))

#include "dns_read_grid.h"

! ###################################################################
  wavenumber = C_1_R
  pi_loc     = ACOS(-C_1_R)

  WRITE(*,*) '1-ICs / 2-Error ?'; READ(*,*) iopt

  IF      ( iopt .EQ. 1 ) THEN
  DO j = 1,jmax; DO i = 1,imax
     ij = i + imax*(j-1)
     s(ij,1) = SIN(C_2_R*pi_loc/scalex*wavenumber*x(i))
  ENDDO; ENDDO
  fname = TRIM(ADJUSTL(tag_scal))//'0'
  CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, i1,i0, s, wrk3d)

  q(:,1) = mean_u; q(:,2) = C_0_R; q(:,3) = C_0_R
  fname = TRIM(ADJUSTL(tag_flow))//'0'
  CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, q, wrk3d)

! ###################################################################
  ELSE IF ( iopt .EQ. 2 ) THEN
  WRITE(*,*) 'Iteration?'; READ(*,*) itime

  WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
  CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, i1,i0, isize_wrk3d, s, wrk3d)

! Theoretical
  factor = EXP(-visc*rtime*(C_2_R*pi_loc/scalex*wavenumber)**2)
  DO j = 1,jmax; DO i = 1,imax
     ij = i + imax*(j-1)
     x_loc = x(i) - mean_u*rtime; 
     s_r(ij,1) = factor*SIN(C_2_R*pi_loc/scalex*wavenumber*x_loc)
  ENDDO; ENDDO

! Error
  error = C_0_R
  dummy = C_0_R
  DO ij = 1,isize_field
     wrk3d(ij) = s(ij,1)-s_r(ij,1)
     error = error + wrk3d(ij)*wrk3d(ij)
     dummy = dummy + s_r(ij,1)*s_r(ij,1)
  ENDDO

  IF ( dummy .GT. C_0_R ) THEN
     WRITE(*,*) 'L-infinity .................: ', MAXVAL(ABS(wrk3d))
     WRITE(*,*) 'Absolute error .............: ', sqrt(error/M_REAL(imax*jmax))
     WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
     fname = 'error'
     CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, i1, i1, wrk3d, wrk3d)
  ENDIF
  
  ENDIF

  STOP
END PROGRAM VDIFFUSION

