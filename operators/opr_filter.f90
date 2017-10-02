#include "types.h"
#include "dns_const.h"

SUBROUTINE OPR_FILTER(nx,ny,nz, f, u, wrk1d,wrk2d,txc)
        
  USE DNS_TYPES,  ONLY : filter_dt
  USE DNS_GLOBAL, ONLY : isize_txc_field, g
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz
  TYPE(filter_dt), DIMENSION(3),                 INTENT(IN)    :: f
  TREAL,           DIMENSION(nx,ny,nz),          INTENT(INOUT) :: u   ! Inplace operation
  TREAL,           DIMENSION(ny,*),              INTENT(INOUT) :: wrk1d
  TREAL,           DIMENSION(nx,nz,2),           INTENT(INOUT) :: wrk2d
  TREAL,           DIMENSION(isize_txc_field,*), INTENT(INOUT) :: txc ! size 2 if ADM
                                                                      ! size 3 if ALPHA
                                                                      ! size 4 if SPECTRAL
! -------------------------------------------------------------------
  TREAL dummy
  TINTEGER k
  
! ###################################################################
! Global filters
  SELECT CASE( f(1)%type )

  CASE( DNS_FILTER_ALPHA                  )
     DO k = 1,nz
        wrk2d(:,k,1) = u(:,1, k)
        wrk2d(:,k,2) = u(:,ny,k)
     ENDDO
     u = u *f(1)%parameters(2)
     CALL OPR_HELMHOLTZ_FXZ(nx,ny,nz, g, i0, f(1)%parameters(2),&
          u, txc(1,1),txc(1,2), &
          wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,wrk1d(1,5),txc(1,3))
     
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
        CALL OPR_FILTER_X(nx,ny,nz, f(1), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
     ENDIF
     
     IF ( f(2)%type .NE. DNS_FILTER_NONE ) THEN
        CALL OPR_FILTER_Y(nx,ny,nz, f(2), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
     ENDIF
     
     IF ( f(3)%type .NE. DNS_FILTER_NONE ) THEN
        CALL OPR_FILTER_Z(nx,ny,nz, f(3), u, txc(1,2), wrk1d,wrk2d,txc(1,1))
     ENDIF

  END SELECT
  
  RETURN
END SUBROUTINE OPR_FILTER

