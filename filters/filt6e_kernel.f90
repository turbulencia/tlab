SUBROUTINE FILT6E_KERNEL(imax, kmax, k1bc, kflt1bc, kfltmxbc, u, uflt)

  IMPLICIT NONE

#include "types.h"

  TINTEGER imax, kmax
  TINTEGER k1bc, kflt1bc, kfltmxbc
  TREAL u(imax,kmax)
  TREAL uflt(imax,kmax)

  TINTEGER i,k
  TINTEGER kstart, kstop
  TREAL b0, b1, b2, b3
  TREAL b_a(7), b_b(7), b_c(7)

  TINTEGER i1, i2, i3, i4, im3, ip3, im2, ip2, im1, ip1

#ifdef SINGLE_PREC
  b0 = 11.0e0/16.0e0
  b1 = 15.0e0/64.0e0
  b2 = -3.0e0/32.0e0
  b3 = 1.0e0/64.0e0

  b_a(1) = 63.0e0/64.0e0
  b_a(2) = 3.0e0/32.0e0
  b_a(3) = -15.0e0/64.0e0
  b_a(4) = 5.0e0/16.0e0
  b_a(5) = -15.0e0/64.0e0
  b_a(6) = 3.0e0/32.0e0
  b_a(7) = -1.0e0/64.0e0

  b_b(1) = 1.0e0/16.0e0
  b_b(2) = 3.0e0/4.0e0
  b_b(3) = 3.0e0/8.0e0
  b_b(4) = -1.0e0/4.0e0
  b_b(5) = 1.0e0/16.0e0
  b_b(6) = 0.0e0
  b_b(7) = 0.0e0

  b_c(1) = -1.0e0/32.0e0
  b_c(2) = 5.0e0/32.0e0
  b_c(3) = 11.0e0/16.0e0
  b_c(4) = 5.0e0/16.0e0
  b_c(5) = -5.0e0/32.0e0
  b_c(6) = 1.0e0/32.0e0
  b_c(7) = 0.0e0

#else
  b0 = 11.0d0/16.0d0
  b1 = 15.0d0/64.0d0
  b2 = -3.0d0/32.0d0
  b3 = 1.0d0/64.0d0

  b_a(1) = 63.0d0/64.0d0
  b_a(2) = 3.0d0/32.0d0
  b_a(3) = -15.0d0/64.0d0
  b_a(4) = 5.0d0/16.0d0
  b_a(5) = -15.0d0/64.0d0
  b_a(6) = 3.0d0/32.0d0
  b_a(7) = -1.0d0/64.0d0

  b_b(1) = 1.0d0/16.0d0
  b_b(2) = 3.0d0/4.0d0
  b_b(3) = 3.0d0/8.0d0
  b_b(4) = -1.0d0/4.0d0
  b_b(5) = 1.0d0/16.0d0
  b_b(6) = 0.0d0
  b_b(7) = 0.0d0

  b_c(1) = -1.0d0/32.0d0
  b_c(2) = 5.0d0/32.0d0
  b_c(3) = 11.0d0/16.0d0
  b_c(4) = 5.0d0/16.0d0
  b_c(5) = -5.0d0/32.0d0
  b_c(6) = 1.0d0/32.0d0
  b_c(7) = 0.0d0
#endif

  kstart = 1
  kstop = kmax

  IF ( k1bc .NE. 0 ) THEN
!    Non periodic case

     k=1
     DO i=1, imax
        uflt(i,k) = u(i,k)
     ENDDO

     IF ( kflt1bc .EQ. 1 ) THEN
!       If the second plane is included do it separately.

        k=2
        DO i=1, imax
           uflt(i,k) = &
                b_b(1) * u(i,k-1)&
                + b_b(2) * u(i,k)&
                + b_b(3) * u(i,k+1) &
                + b_b(4) * u(i,k+2)&
                + b_b(5) * u(i,k+3)&
                + b_b(6) * u(i,k+4)&
                + b_b(7) * u(i,k+5)
        ENDDO

!       If the third plane is included do it separately.

        k=3
        DO i=1, imax
           uflt(i,k) = &
                b_c(1) * u(i,k-2)&
                + b_c(2) * u(i,k-1)&
                + b_c(3) * u(i,k) &
                + b_c(4) * u(i,k+1)&
                + b_c(5) * u(i,k+2)&
                + b_c(6) * u(i,k+3)&
                + b_c(7) * u(i,k+4)
        ENDDO

     ELSE

        DO k=2,3
           DO i=1, imax
              uflt(i,k) = u(i,k)
           ENDDO
        ENDDO

     ENDIF

     kstart = 4

     k=kmax
     DO i=1, imax
        uflt(i,k) = u(i,k)
     ENDDO

     IF ( kfltmxbc .EQ. 1 ) THEN
!       If the second plane is included do it separately.

        k=kmax-1
        DO i=1, imax
           uflt(i,k) = &
                  b_b(1) * u(i,k+1) &
                + b_b(2) * u(i,k)   &
                + b_b(3) * u(i,k-1) &
                + b_b(4) * u(i,k-2) &
                + b_b(5) * u(i,k-3) &
                + b_b(6) * u(i,k-4) &
                + b_b(7) * u(i,k-5)
        ENDDO

!       If the third plane is included do it separately.

        k=kmax-2
        DO i=1, imax
           uflt(i,k) = &
                  b_c(1) * u(i,k+2) &
                + b_c(2) * u(i,k+1) &
                + b_c(3) * u(i,k)   &
                + b_c(4) * u(i,k-1) &
                + b_c(5) * u(i,k-2) &
                + b_c(6) * u(i,k-3) &
                + b_c(7) * u(i,k-4)
        ENDDO
     ELSE

        DO k=kmax-2,kmax-1
           DO i=1, imax
              uflt(i,k) = u(i,k)
           ENDDO
        ENDDO

     ENDIF

     kstop = kmax-3

  ENDIF

  i1 = 1
  i2 = 2 
  i3 = 3
  i4 = 4

  DO k=kstart, kstop
     im3 = MOD(k-i4+kmax,kmax)+1
     im2 = MOD(k-i3+kmax,kmax)+1
     im1 = MOD(k-i2+kmax,kmax)+1
     ip1 = MOD(k+kmax,kmax)+1
     ip2 = MOD(k+i1+kmax,kmax)+1
     ip3 = MOD(k+i2+kmax,kmax)+1
     DO i=1, imax
        uflt(i,k) = b3*(u(i,im3)+u(i,ip3))&
             + b2*(u(i,im2)+u(i,ip2)) &
             + b1*(u(i,im1)+u(i,ip1)) &
             + b0*u(i,k)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FILT6E_KERNEL
