#include "types.h"  
#include "dns_const.h"

PROGRAM VPOISSON
  
  USE DNS_CONSTANTS
  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: a, b, c, d, e, f
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: txc
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: wrk1d, wrk2d, bcs_hb, bcs_ht 
  TREAL, DIMENSION(:),     ALLOCATABLE :: cx, cy, cz, wrk3d
  
  TINTEGER i, j, k,  ibc_x(4), ibc_y(4), ibc_z(4)
  TINTEGER isize_wrk3d, itype
  TREAL dummy, error, falpha, lambda

  TREAL, DIMENSION(:,:), POINTER :: dx, dy, dz

! ###################################################################
  CALL DNS_INITIALIZE
  
  CALL DNS_READ_GLOBAL('dns.ini')

  isize_wrk3d = isize_txc_field

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d+1))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(bcs_ht(imax,kmax),bcs_hb(imax,kmax))
  ALLOCATE(a(imax,jmax,kmax),b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(d(imax,jmax,kmax),e(imax,jmax,kmax),f(imax,jmax,kmax))
  ALLOCATE(txc(isize_txc_field,2),wrk3d(isize_wrk3d))
  ALLOCATE(cx(6*imax),cy(6*jmax),cz(6*kmax_total))

#include "dns_read_grid.h"

  dx => x(:,2:) ! to be removed
  dy => y(:,2:)
  dz => z(:,2:)

  falpha = 0.49d0
  CALL FILT4C_INI(imax_total, i1bc, falpha, dx, cx)
  CALL FILT4C_INI(jmax_total, j1bc, falpha, dy, cy)
  CALL FILT4C_INI(kmax_total, k1bc, falpha, dz, cz)
! BCs for the filters (see routine FILTER)
  ibc_x(1) = 1; ibc_x(2) = i1bc; ibc_x(3) = 0; ibc_x(4) = 0
  ibc_y(1) = 0; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0 
  ibc_z(1) = 1; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0 

  CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  
! ###################################################################
! Define forcing term
! ###################################################################
  CALL DNS_READ_FIELDS('field.inp', i1, imax,jmax,kmax, i1,i0, i1, a, wrk3d)

! remove 2\Delta x wave
!  CALL OPR_FILTER(i1,imax,jmax,kmax, ibc_x,ibc_y,ibc_z, i1, a, cx,cy,cz,wrk1d,wrk2d,wrk3d)
! -------------------------------------------------------------------
  ! DO j = 1,jmax
  !    mean = AVG_IK(imax,jmax,kmax, j, a, dx,dz, area)
  !    a(:,j,:) = a(:,j,:) - mean
  ! ENDDO
! -------------------------------------------------------------------
! DC level at lower boundary set to zero
  ! mean = AVG_IK(imax,jmax,kmax, i1, a, dx,dz, area)
  ! a = a - mean

! ###################################################################
  f = a; bcs_hb = C_0_R; bcs_ht = C_0_R
  itype = 2

  IF      ( itype .EQ. 1 ) THEN
     CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, i3, &
          a,c, txc(1,1),txc(1,2), bcs_hb,bcs_ht, wrk1d,wrk1d(1,5),wrk3d)
     e = c ! save dp/dy
  ELSE IF ( itype .EQ. 2 ) THEN
     WRITE(*,*) 'Eigenvalue ?'
     READ(*,*) lambda
     ! CALL OPR_HELMHOLTZ_FXZ(imax,jmax,kmax, i0, lambda,&
     !      a, txc(1,1),txc(1,2), bcs_hb,bcs_ht, wrk1d,wrk1d(1,5),wrk3d)
     CALL OPR_HELMHOLTZ_FXZ_2(imax,jmax,kmax, g, i0, lambda,&
          a, txc(1,1),txc(1,2), bcs_hb,bcs_ht, wrk1d,wrk1d(1,5),wrk3d)
  ENDIF

! -------------------------------------------------------------------
  ! CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, a, c, i0,i0, wrk1d,wrk2d,wrk3d)
  ! CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, c, b, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1,iunifx,imode_fdm, imax,jmax,kmax, i1bc, dx, a,b, i0,i0, i0,i0, c, wrk1d,wrk2d,wrk3d)

  ! CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, a, c, i0,i0, wrk1d,wrk2d,wrk3d)
  ! CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, c, d, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_ZZ(i1,iunifz,imode_fdm, imax,jmax,kmax, k1bc, dz, a,d, i0,i0, i0,i0, c, wrk1d,wrk2d,wrk3d)
  b = b + d
     
  ! CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, a, c, i0,i0, wrk1d,wrk2d,wrk3d)
  ! CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, c, d, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1,iunify,imode_fdm, imax,jmax,kmax, j1bc, dy, a,d, i0,i0, i0,i0, c, wrk1d,wrk2d,wrk3d)
  b = b + d
!  bcs_hb(:,:) = c(:,1,   :); bcs_ht(:,:) = c(:,jmax,:) ! Neumann BCs 

! ###################################################################
  IF ( itype .EQ. 2 ) THEN
     b = b + lambda*a
  ENDIF

! ###################################################################
! solve poisson eqn
! ###################################################################
!  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, i3, &
!       b,c,d,e, txc(1,1),txc(1,2), bcs_hb,bcs_ht, wrk1d,wrk1d(1,5),wrk3d)
  CALL DNS_WRITE_FIELDS('field.out', i1, imax,jmax,kmax, i1, i1, b, wrk3d)

  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, a, c, i0,i0, wrk1d,wrk2d,wrk3d)
! -------------------------------------------------------------------
  a = f ! rhs
  d = e ! dp/dy

! ###################################################################
! Error
! ###################################################################
  error = C_0_R
  dummy = C_0_R
  DO k = 1,kmax
!     DO j = 1,jmax
     DO j = 2,jmax-1
        DO i = 1,imax
           e(i,j,k) = b(i,j,k)-a(i,j,k)
           error = error + e(i,j,k)*e(i,j,k)
           dummy = dummy + a(i,j,k)*a(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, i1, e, wrk3d)
  
! first derivative
  error = C_0_R
  dummy = C_0_R
  DO k = 1,kmax
!     DO j = 1,jmax
     DO j = 2,jmax-1
        DO i = 1,imax
           e(i,j,k) = d(i,j,k)-c(i,j,k)
           error = error + e(i,j,k)*e(i,j,k)
           dummy = dummy + c(i,j,k)*c(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  WRITE(*,*) 'Relative error in df/dy ....: ', sqrt(error)/sqrt(dummy)
!  CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, i1, e, wrk3d)
  
  STOP
END PROGRAM VPOISSON
