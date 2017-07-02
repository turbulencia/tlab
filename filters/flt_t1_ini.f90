#include "types.h"

!########################################################################
!# DESCRIPTION
!#
!# Top-hat filter coefficients initialization
!# Trapezoidal rule
!# Free boundary conditions
!#  (Ghost cells w/ linear extrapolation of function from the two nodes next to boundary)
!# 
!########################################################################
SUBROUTINE FLT_T1_INI(scalex, x, f, wrk1d)

  USE DNS_TYPES,  ONLY : filter_dt

  IMPLICIT NONE

  TREAL,                      INTENT(IN)    :: scalex
  TREAL, DIMENSION(*),        INTENT(IN)    :: x
  TYPE(filter_dt),            INTENT(INOUT) :: f
  TREAL, DIMENSION(f%size,3), INTENT(INOUT) :: wrk1d
  
! -----------------------------------------------------------------------
  TREAL dum
  TINTEGER i, ii, im, ic, ip, nx

! #######################################################################
  wrk1d = C_0_R
  nx = f%delta

! #######################################################################
  IF ( f%uniform ) THEN
     IF ( .NOT. f%periodic ) THEN ! I only need info for the two nodes next to the boundary
        DO i = 1,nx/2
           im = nx/2-i+1
           dum = C_0_R
           DO ii=1,im
              dum = dum+M_REAL(ii)
           ENDDO
           ip = (i-1)*2 + 1
           f%coeffs(ip,1) = dum+C_05_R*M_REAL(im+1)
           IF ( nx .EQ. 2 ) THEN
              f%coeffs(ip+1,1) = C_05_R-dum+C_05_R*M_REAL(im)
           ELSE
              f%coeffs(ip+1,1) = C_1_R -dum+C_05_R*M_REAL(im)
           ENDIF
        ENDDO
     ENDIF

! #######################################################################
  ELSE
     ! calculate delta(i)
     DO i = 1,f%size-1
        wrk1d(i,1) = x(i+1)-x(i)
     ENDDO
     IF ( f%periodic ) THEN
        wrk1d(f%size,1) = scalex-(x(f%size)-x(1))
     ELSE
        wrk1d(f%size,1) = wrk1d(f%size-1,1)
     ENDIF
     ! calculate deltasum(i) = delta(i-1)+delta(i)
     DO i = 2,f%size
        wrk1d(i,2) = wrk1d(i-1,1)+wrk1d(i,1)
     ENDDO
     IF ( f%periodic ) THEN
        wrk1d(1,2) = wrk1d(f%size,1)+wrk1d(1,1)
     ELSE
        wrk1d(1,2) = wrk1d(1,1)*C_2_R
     ENDIF
     ! calculate deltaf(i) = delta(i-nx/2) + delta(i-nx/2+1) + ... + delta(i+nx/2-1) 
     IF ( f%periodic ) THEN
        DO i = 1,f%size
           wrk1d(i,3) = C_0_R
           DO ii = i-nx/2,i+nx/2-1
              im = ii+f%size-1
              im = MOD(im,f%size)+1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(im,1)
           ENDDO
        ENDDO
     ELSE
        DO i = 1,nx/2
           wrk1d(i,3) = wrk1d(1,1)*M_REAL(nx/2-i+1)
           DO ii = 1,i+nx/2-1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(ii,1)
           ENDDO
        ENDDO
        DO i = 1+nx/2,f%size-nx/2
           DO ii = i-nx/2,i+nx/2-1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(ii,1)
           ENDDO
        ENDDO
        DO i = f%size-nx/2+1,f%size
           wrk1d(i,3) = wrk1d(f%size,1)*M_REAL(nx/2-(f%size-i))
           DO ii = i-nx/2,f%size-1
              wrk1d(i,3) = wrk1d(i,3)+wrk1d(ii,1)
           ENDDO
        ENDDO
     ENDIF

! -----------------------------------------------------------------------
!    construct the coefficients array as if periodic; corrections for nonperiodic below
     DO i = 1,f%size
        ii = i-nx/2               ! I need to use delta(ii)
        im = ii+f%size-1          ! The index im deals with periodicity
        im = MOD(im,f%size)+1
        ic = ii-i+nx/2+1
        ip = (i-1)*(nx+1) + ic
        f%coeffs(ip,1) = C_05_R*wrk1d(im,1)/wrk1d(i,3)

        DO ii = i-nx/2+1,i+nx/2-1 ! I need to use deltasum(ii)
           im = ii+f%size-1       ! The index im deals with periodicity
           im = MOD(im,f%size)+1
           ic = ii-i+nx/2+1
           ip = (i-1)*(nx+1) + ic
           f%coeffs(ip,1) = C_05_R*wrk1d(im,2)/wrk1d(i,3)
        ENDDO

        ii = i+nx/2               ! I need to use delta(ii-1)
        im = (ii-1)+f%size-1      ! The index im deals with periodicity
        im = MOD(im,f%size)+1
        ic = ii-i+nx/2+1
        ip = (i-1)*(nx+1) + ic
        f%coeffs(ip,1) = C_05_R*wrk1d(im,1)/wrk1d(i,3)

     ENDDO

! -----------------------------------------------------------------------
! modification in case of free boundary
     IF ( .NOT. f%periodic ) THEN
        DO i = 1,nx/2
           im = nx/2-i+1
           dum = C_0_R
           DO ii=1,im
              dum = dum+M_REAL(ii)
           ENDDO
           ic = nx/2-i+2
           ip = (i-1)*(nx+1) + ic
           f%coeffs(ip,1) = f%coeffs(ip,1) + &
                C_05_R*(wrk1d(1,2)*(dum-C_1_R) + &
                wrk1d(1,1)*M_REAL(im+1))/wrk1d(i,3)
           ic = nx/2-i+3
           ip = (i-1)*(nx+1) + ic
           f%coeffs(ip,1) = f%coeffs(ip,1) - &
                C_05_R*(wrk1d(1,2)*(dum-M_REAL(im)) +&
                wrk1d(1,1)*M_REAL(im))/wrk1d(i,3)
           ! pad with ceros
           DO ic=1,nx/2-i+1
              ip = (i-1)*(nx+1) + ic
              f%coeffs(ip,1) = C_0_R
           ENDDO
        ENDDO

        DO i = f%size-nx/2+1,f%size
           im = i-f%size+nx/2
           dum = C_0_R
           DO ii=1,im
              dum = dum+M_REAL(ii)
           ENDDO
           ic = nx+1-im-1
           ip = (i-1)*(nx+1) + ic
           f%coeffs(ip,1) = f%coeffs(ip,1) - &
                C_05_R*(wrk1d(f%size,2)*(dum-M_REAL(im)) +&
                wrk1d(f%size,1)*M_REAL(im))/wrk1d(i,3)
           ic = nx+1-im
           ip = (i-1)*(nx+1) + ic
           f%coeffs(ip,1) = f%coeffs(ip,1) + &
                C_05_R*(wrk1d(f%size,2)*(dum-C_1_R) + &
                wrk1d(f%size,1)*M_REAL(im+1))/wrk1d(i,3)
          ! pad with ceros
           DO ic=nx+1-im+1,nx+1
              ip = (i-1)*(nx+1) + ic
              f%coeffs(ip,1) = C_0_R
           ENDDO
        ENDDO
        
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE FLT_T1_INI
