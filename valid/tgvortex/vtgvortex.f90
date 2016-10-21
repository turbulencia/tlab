#include "types.h"  
#include "dns_const.h"

PROGRAM VTGVORTEX
  
  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:), ALLOCATABLE :: txc, q
  TREAL, DIMENSION(:),   ALLOCATABLE :: wrk1d, wrk2d, wrk3d
  
  TINTEGER ij, iv, iopt
  TINTEGER isize_wrk3d
  TREAL dummy, error
  CHARACTER*(32) fname

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

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d* inb_wrk2d))
  ALLOCATE(wrk3d(isize_wrk3d))
  ALLOCATE(  q(isize_field,    4))
  ALLOCATE(txc(isize_txc_field,4))

#include "dns_read_grid.h"

! ###################################################################
  CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  
  WRITE(*,*) '1-ICs / 2-Error ?'; READ(*,*) iopt

  IF      ( iopt .EQ. 1 ) THEN
  rtime = C_0_R
  CALL FLOW_TAYLORGREEN(imax,jmax,kmax, rtime,visc, x,y,z, q(1,1),q(1,2),q(1,3),q(1,4))
  q(:,1) = q(:,1) + mean_u
  fname = 'pt0'
  CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, q,      wrk3d)

  q(:,4) = C_0_R
  fname = 'sc0'
  CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, i1,i0, q(1,4), wrk3d)

! ###################################################################
  ELSE IF ( iopt .EQ. 2 ) THEN
  WRITE(*,*) 'Iteration?'
  READ(*,*) itime

  WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
  CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, isize_wrk3d, q, wrk3d)
  txc(:,1) = C_0_R; txc(:,4) = C_0_R
!  CALL FI_FORCING_1(iunifx,iunify, imode_fdm, imax,jmax,kmax, i1bc,j1bc, & 
!       rtime,visc, x,y,dx,dy, txc(1,1),txc(1,4), q(1,1),q(1,2),q(1,3),q(1,4), wrk1d,wrk2d,wrk3d)
!  CALL FI_FORCING_0(imax,jmax,kmax, rtime,visc, x,y, q(1,1),q(1,2), txc(1,1),txc(1,4))
!  CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, isize_wrk3d, q, wrk3d)

  CALL FI_PRESSURE_BOUSSINESQ(q(1,1),q(1,2),q(1,3),txc(1,4), q(1,4), &
       txc(1,1),txc(1,2),txc(1,3), wrk1d,wrk2d,wrk3d)

! Theoretical Taylor-Green in array txc
  x = x - mean_u*rtime
  CALL FLOW_TAYLORGREEN(imax,jmax,kmax, rtime,visc, x,y,z, txc(1,1),txc(1,2),txc(1,3),txc(1,4))
  txc(:,1) = txc(:,1) + mean_u
  
! Error
  DO iv = 1,4
     error = C_0_R
     dummy = C_0_R
     DO ij = 1,isize_field !-imax
        dummy = dummy + txc(ij,iv)*txc(ij,iv)
        txc(ij,iv) = q(ij,iv)-txc(ij,iv)
        error = error + txc(ij,iv)*txc(ij,iv)
     ENDDO
!     print*,dummy*dx(1)*dy(1)

     IF ( dummy .GT. C_0_R ) THEN
     WRITE(*,*) 'L-infinity error ...: ', MAXVAL(ABS(txc(1:isize_field-imax,iv)))
     WRITE(*,*) 'L-2 error ..........: ', sqrt(error*g(1)%jac(1,1)*g(2)%jac(1,1))
     WRITE(*,*) 'Relative error .....: ', sqrt(error)/sqrt(dummy)
     WRITE(fname,*) iv; fname = 'error'//TRIM(ADJUSTL(fname))
     CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, i1, i0, txc(1,iv), wrk3d)
     ENDIF

  ENDDO

  ENDIF

  STOP
END PROGRAM VTGVORTEX

!########################################################################
!########################################################################

SUBROUTINE FLOW_TAYLORGREEN(nx,ny,nz, rtime,visc, x,y,z, u,v,w,p)

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL rtime, visc
  TREAL, DIMENSION(*)        :: x,y,z
  TREAL, DIMENSION(nx,ny,nz) :: u,v,w,p

! -----------------------------------------------------------------------
  TINTEGER i,j,k
  TREAL pi_loc, omega, factor, sigma

! #######################################################################
  pi_loc= ACOS(-C_1_R); omega = C_2_R*pi_loc

  DO k = 1,nz
     DO j = 1,ny; DO i = 1,nx
!        u(i,j,k) = SIN(x(i)*omega)*COS(y(j)*omega)
!        v(i,j,k) =-COS(x(i)*omega)*SIN(y(j)*omega)
        u(i,j,k) = C_05_R*(SIN((x(i)+y(j))*omega) + SIN((x(i)-y(j))*omega))
        v(i,j,k) =-C_05_R*(SIN((x(i)+y(j))*omega) - SIN((x(i)-y(j))*omega))
        p(i,j,k) =(COS(C_2_R*omega*x(i))+COS(C_2_R*omega*y(j))-C_1_R)*C_025_R

!        u(i,j,k) = sin(omega*x(i))*       sin(C_2_R*omega*y(j))
!        v(i,j,k) =-cos(omega*x(i))*(C_1_R-cos(C_2_R*omega*y(j)))*C_05_R
!        p(i,j,k) = cos(C_2_R*omega*x(i))*(C_2_R-cos(C_2_R*omega*y(j)))/C_8_R &
!             - C_05_R*(sin(omega*y(j)))**4

!        u(i,j,k) = sin(pi_loc*x(i))*sin(pi_loc*x(i))*sin(C_2_R*pi_loc*y(j))
!        v(i,j,k) =-sin(C_2_R*pi_loc*x(i))*sin(pi_loc*y(j))*sin(pi_loc*y(j))
!        p(i,j,k) = sin(C_2_R*pi_loc*x(i))*sin(C_2_R*pi_loc*y(j))

     ENDDO; ENDDO
  ENDDO
  w = C_0_R

! -----------------------------------------------------------------------
  sigma = omega*omega*C_2_R*visc
  factor = EXP(-sigma*rtime)
!  factor = cos(omega*rtime) - sigma/omega*sin(omega*rtime)
!  factor = EXP( -sigma*rtime + C_1_R/omega*(C_1_R-cos(omega*rtime)) )
!  factor = cos(omega*rtime)

  u = u*factor; v = v*factor; w = w*factor
  p = p*factor*factor

  RETURN
END SUBROUTINE FLOW_TAYLORGREEN
