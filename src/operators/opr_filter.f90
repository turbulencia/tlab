#include "types.h"
#include "dns_const.h"

SUBROUTINE OPR_FILTER(nx,ny,nz, f, u, wrk1d,wrk2d,txc)
        
  USE DNS_TYPES,  ONLY : filter_dt
  USE TLAB_VARS, ONLY : isize_txc_field, g
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz
  TYPE(filter_dt), DIMENSION(3),                 INTENT(IN)    :: f
  TREAL,           DIMENSION(nx,ny,nz),          INTENT(INOUT) :: u   ! Inplace operation
  TREAL,           DIMENSION(ny,*),              INTENT(INOUT) :: wrk1d
  TREAL,           DIMENSION(nx,nz,2),           INTENT(INOUT) :: wrk2d
  TREAL,           DIMENSION(isize_txc_field,*), INTENT(INOUT) :: txc ! size 2 if ADM
                                                                      ! size 4 if SPECTRAL, HELMHOLTZ
! -------------------------------------------------------------------
  TREAL dummy
  TINTEGER k, flag_bcs, n, bcs(2,2), nxy, ip_b, ip_t

! ###################################################################
  nxy = nx*ny

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

! Global filters
  SELECT CASE( f(1)%type )

  CASE( DNS_FILTER_HELMHOLTZ )
     IF      ( f(2)%BcsMin .EQ. DNS_FILTER_BCS_DIRICHLET ) THEN
        DO k = 1,nz
           wrk2d(:,k,1) = u(:,1, k)
           wrk2d(:,k,2) = u(:,ny,k)
        ENDDO
        flag_bcs = 0
     ELSE IF ( f(2)%BcsMin .EQ. DNS_FILTER_BCS_NEUMANN   ) THEN
        CALL OPR_PARTIAL_Y(OPR_P1, nx,ny,nz, bcs, g(2), u,txc(1,1), txc(1,2), wrk2d,txc(1,2))
        ip_b =             1
        ip_t = nx*(ny-1) + 1
        DO k = 1,nz
           wrk2d(1:nx,k,1) = txc(ip_b:ip_b+nx-1,1); ip_b = ip_b + nxy ! bottom
           wrk2d(1:nx,k,2) = txc(ip_t:ip_t+nx-1,1); ip_t = ip_t + nxy ! top
        ENDDO
        flag_bcs = 3
     ELSE IF ( f(2)%BcsMin .EQ. DNS_FILTER_BCS_SOLID    ) THEN
        wrk2d(:,:,1) = C_0_R
        wrk2d(:,:,2) = C_0_R
        flag_bcs = 3
     ENDIF
     
     txc(1:nx*ny*nz,1) = u(1:nx*ny*nz,1,1) *f(1)%parameters(2) ! I need extended arrays
     CALL OPR_HELMHOLTZ_FXZ(nx,ny,nz, g, flag_bcs, f(1)%parameters(2),&
          txc(1,1), txc(1,2),txc(1,3), &
          wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),txc(1,4))
     u(1:nx*ny*nz,1,1) = txc(1:nx*ny*nz,1)
     
  CASE( DNS_FILTER_BAND )
     dummy = C_1_R /M_REAL( f(1)%size *f(3)%size )
     txc(1:nx*ny*nz,1) = u(1:nx*ny*nz,1,1) ! I need extended arrays
     CALL OPR_FOURIER_F(i2, nx,ny,nz, txc(1,1),txc(1,2), txc(1,3),wrk2d,txc(1,4)) 
     CALL OPR_FILTER_BAND_2D(nx,ny,nz, f(1)%parameters, txc(1,2)) 
     CALL OPR_FOURIER_B(i2, nx,ny,nz, txc(1,2), txc(1,3), txc(1,4))
     u(1:nx*ny*nz,1,1) = txc(1:nx*ny*nz,3) *dummy

  CASE( DNS_FILTER_ERF )
     dummy = C_1_R /M_REAL( f(1)%size *f(3)%size )
     txc(1:nx*ny*nz,1) = u(1:nx*ny*nz,1,1) ! I need extended arrays
     CALL OPR_FOURIER_F(i2, nx,ny,nz, txc(1,1),txc(1,2), txc(1,3),wrk2d,txc(1,4)) 
     CALL OPR_FILTER_ERF_2D(nx,ny,nz,f(1)%parameters,txc(1,2))  
     CALL OPR_FOURIER_B(i2, nx,ny,nz, txc(1,2), txc(1,3), txc(1,4))
     u(1:nx*ny*nz,1,1) = txc(1:nx*ny*nz,3) *dummy
     
  CASE DEFAULT
     
! ###################################################################
! Directional filters
     IF ( f(1)%type .NE. DNS_FILTER_NONE ) THEN
        DO n = 1,f(1)%repeat
           CALL OPR_FILTER_X(nx,ny,nz, f(1), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
        END DO
     ENDIF
     
     IF ( f(2)%type .NE. DNS_FILTER_NONE ) THEN
        DO n = 1,f(2)%repeat
           CALL OPR_FILTER_Y(nx,ny,nz, f(2), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
        END DO
     ENDIF
     
     IF ( f(3)%type .NE. DNS_FILTER_NONE ) THEN
        DO n = 1,f(3)%repeat
           CALL OPR_FILTER_Z(nx,ny,nz, f(3), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
        END DO
     ENDIF
     
  END SELECT
  
  RETURN
END SUBROUTINE OPR_FILTER

