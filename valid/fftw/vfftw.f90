#include "types.h"  

PROGRAM VFFTW
  
  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_FFTW
#include "fftw3.f"
#endif

  TREAL,    DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL,    DIMENSION(:,:,:), POINTER :: a, b, c
  TCOMPLEX, DIMENSION(:,:,:), POINTER :: a1, a2, a3
  TREAL,    DIMENSION(:),     POINTER :: wrk1d, wrk2d, wrk3d
  
  TCOMPLEX :: Img
 
  integer*4 i, j, k, ij
  INTEGER(8) fft_plan_fx, fft_plan_fz
  INTEGER(8) fft_plan_bx, fft_plan_bz

!  TREAL fft_data_x, fft_data_z
  TREAL dummy, error

  TREAL, DIMENSION(:,:), POINTER :: dx, dy, dz

! ###################################################################
  CALL DNS_INITIALIZE
  
  CALL DNS_READ_GLOBAL('dns.ini')

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*10))
  ALLOCATE(wrk2d(isize_wrk2d* 5))
  ALLOCATE(wrk3d(imax*jmax*kmax),a(imax,jmax,kmax),b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(a1(imax/2+1,jmax,kmax), a2(kmax,imax/2+1,jmax), a3(kmax,imax/2+1,jmax))

  Img=(-1.0,0.0)
  Img=sqrt(Img)

#include "dns_read_grid.h"

! ###################################################################
!  Define forcing term
! ###################################################################
!  DO k = 1,kmax
!     DO j = 1,jmax
!        DO i = 1,imax
!           a(i,j,k) = sin(C_2_R*C_PI_R/scalex*x(i)*C_2_R) &
!                    * cos(C_2_R*C_PI_R/scalez*z(k)*C_5_R) &
!                    * M_REAL(j-1)/M_REAL(jmax-1)*C_01_R
!!           c(i,j,k) = sin(C_2_R*C_PI_R/scalex*x(i)*C_2_R) &
!!                    * sin(C_2_R*C_PI_R/scalez*z(k)*C_5_R) * (-C_2_R*C_PI_R/scalez*C_5_R)&
!!                    * M_REAL(j-1)/M_REAL(jmax-1)*C_01_R
!           
!        ENDDO
!     ENDDO
!  ENDDO
!  CALL DNS_WRITE_FIELDS('field.inp', imax, jmax, kmax, kmax_total, i1, i1, i1, a, dummy)
  CALL DNS_READ_FIELDS('field.inp', i1, imax,jmax,kmax, i1,i0, i1, a, dummy)

!  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc, dx, a, c, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc, dz, a, c, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL DNS_WRITE_FIELDS('field.ref', i1, imax,jmax,kmax, i1, i1, c, dummy)

! ###################################################################
  CALL dfftw_plan_dft_r2c_1d&
       (fft_plan_fx, imax,       wrk1d, wrk1d, FFTW_ESTIMATE)
  CALL dfftw_plan_dft_1d&
       (fft_plan_fz, kmax_total, wrk1d, wrk1d, FFTW_FORWARD, FFTW_ESTIMATE)

  CALL dfftw_plan_dft_c2r_1d&
       (fft_plan_bx, imax,       wrk1d, wrk1d, FFTW_ESTIMATE)
  CALL dfftw_plan_dft_1d&
       (fft_plan_bz, kmax_total, wrk1d, wrk1d, FFTW_BACKWARD, FFTW_ESTIMATE)

! ###################################################################
! Fourier
! ###################################################################
! Forward Real FFT in x
  DO k = 1,kmax
     DO j = 1,jmax
        CALL dfftw_execute_dft_r2c(fft_plan_fx, a(1,j,k), a1(1,j,k))
     ENDDO
  ENDDO

! make j the last index, k the first
  DO k = 1,kmax_total
     DO ij = 1,(imax/2+1)*jmax
        a2(k,ij,1) = a1(ij,1,k)
     ENDDO
  ENDDO

! Forward complex FFT in z
  DO j = 1,jmax
     DO i = 1,imax/2+1
        CALL dfftw_execute_dft(fft_plan_fz, a2(1,i,j), a3(1,i,j))
     ENDDO
  ENDDO

! ###################################################################
! Manipulation
! ###################################################################
  DO j = 1,jmax
     DO i = 1,imax/2+1
        DO k = 1,kmax_total
!           dummy = C_2_R*C_PI_R*M_REAL(i-1)/M_REAL(imax)
!           dummy = ( C_14_R/C_9_R*sin(dummy) + C_1_R/C_18_R*sin(2*dummy) ) &
!                 / ( C_1_R + C_2_R/C_3_R*cos(dummy) ) 
!           a3(k,i,j) = dummy*Img * a3(k,i,j) * M_REAL(imax)/scalex

           IF ( k .LE. kmax_total/2 ) THEN
              dummy = C_2_R*C_PI_R*M_REAL(k-1)/M_REAL(kmax_total)
           ELSE
              dummy = C_2_R*C_PI_R*M_REAL(k-1-kmax_total)/M_REAL(kmax_total)
           ENDIF
           dummy = ( C_14_R/C_9_R*sin(dummy) + C_1_R/C_18_R*sin(2*dummy) ) &
                 / ( C_1_R + C_2_R/C_3_R*cos(dummy) ) 
           a3(k,i,j) = dummy*Img * a3(k,i,j) * M_REAL(kmax_total)/scalez

        ENDDO
     ENDDO
  ENDDO

! ###################################################################
! Back
! ###################################################################
! backwards complex FFT in z
  DO j = 1,jmax
     DO i = 1,imax/2+1
           CALL dfftw_execute_dft(fft_plan_bz, a3(1,i,j), a2(1,i,j))
     ENDDO
  ENDDO

! make k the last index, kx the first
  DO k = 1,kmax_total
     DO ij = 1,(imax/2+1)*jmax
        a1(ij,1,k)=a2(k,ij,1)/M_REAL(imax*kmax_total) ! normalize
     ENDDO
  ENDDO

! backwards real FFT in z
  DO k = 1,kmax
     DO j = 1,jmax
        CALL dfftw_execute_dft_c2r(fft_plan_bx, a1(1,j,k), b(1,j,k))
     ENDDO
  ENDDO

  CALL DNS_WRITE_FIELDS('field.out', i1, imax,jmax,kmax, i1, i1, b, dummy)

! ###################################################################
! Error
! ###################################################################
  error = C_0_R
  dummy = C_0_R
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
           error = error + (c(i,j,k)-b(i,j,k))*(c(i,j,k)-b(i,j,k))
           dummy = dummy + c(i,j,k)*c(i,j,k)
        ENDDO
     ENDDO
  ENDDO
  WRITE(*,*) 'Relative error ....: ', sqrt(error)/sqrt(dummy)
  
  STOP
END PROGRAM VFFTW
