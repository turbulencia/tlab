#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2007/01/01 - H. Foysi
!#              Created
!# 2008/11/25 - J.P. Mellado
!#              Modified
!#
!########################################################################
!# DESCRIPTION
!#
!# Explicit 4th order, nonuniform case.
!#
!# Coincides with filter as described in Stolz's thesis, imported from 
!# Holger Foysi. 
!# Coefficient \alpha_{-2} determined by setting minimum imaginary part,
!# which implies alpha_m2 = alpha_p2
!#
!########################################################################
SUBROUTINE FLT_E4_INI(scalez, z, f)

  USE TLAB_TYPES,  ONLY : filter_dt

  IMPLICIT NONE

  TREAL,               INTENT(IN)    :: scalez
  TREAL, DIMENSION(*), INTENT(IN)    :: z
  TYPE(filter_dt),     INTENT(INOUT) :: f

! -----------------------------------------------------------------------
  TINTEGER k, km2, km1, kp1, kp2, kmin_loc, kmax_loc
  TREAL D0,D1,D2,Dm1, zm2,zm1,zp1,zp2,zp3

#define alpha_m2(k)  f%coeffs(k,1)
#define alpha_m1(k)  f%coeffs(k,2)
#define alpha(k)     f%coeffs(k,3)
#define alpha_p1(k)  f%coeffs(k,4)
#define alpha_p2(k)  f%coeffs(k,5)

  IF ( f%size .LE. 1 ) RETURN

! #######################################################################
! Coefficients near wall
! Vanishing moment of third order
! #######################################################################
  IF ( .NOT. f%periodic ) THEN

! -----------------------------------------------------------------------
! Point 2
! -----------------------------------------------------------------------
     k = 2
     
     zm1 = abs(z(k)  -z(k-1))
     zp1 = abs(z(k+1)-z(k)  )
     zp2 = abs(z(k+2)-z(k)  )
     zp3 = abs(z(k+3)-z(k)  )
     
     D2  = zp2*(-zp1+zp2+zm1)
     D1  =-zp1**2+zm1**2+zp2*zp1+zp2*zm1
     D0  = D2
     Dm1 = D1

     alpha_m2(k) = (zp2**3*zp1*zm1/(C_2_R*D2) - zp1**3*(zm1**2+zp2*zm1)/(C_2_R*D1) + &
          zm1**3*(zp1**2-zp2*zp1)/(C_2_R*Dm1))/                                      &
          ( zp3**3 - zp2**3*(-zp1*zp3-zp1*zm1+zp3**2+zm1*zp3)/D2 + &
          zp1**3*(-zm1**2-zp2*zp3-zp2*zm1+C_2_R*zp3**2)/D1       - &
          zm1**3*(-zp2*zp3-zp1**2+zp3**2+zp2*zp1)/Dm1            )
     
     alpha_m1(k) =-C_05_R*(zp1**2-zp2*zp1+C_2_R*alpha_m2(k)*(-zp2*zp3+zp3**2-zp1**2+zp2*zp1))/Dm1
     alpha(k)    = C_05_R*(-zp2*zp1+zp2**2+zp2*zm1+zp1*zm1-C_2_R*alpha_m2(k)*(zp1*zp3+zp1*zm1-zp3**2-zm1*zp3))/D0
     alpha_p1(k) = C_05_R*(zm1**2+zp2*zm1+C_2_R*alpha_m2(k)*(-zm1**2-zp2*zp3-zp2*zm1+zp3**2))/D1
     alpha_p2(k) =-C_05_R*(C_2_R*alpha_m2(k)*(-zp1*zp3-zp1*zm1+zp3**2+zm1*zp3)+zp1*zm1)/D2
     
! -----------------------------------------------------------------------
! Point N-1
! -----------------------------------------------------------------------
     k = f%size-1
     
     zm1 = abs(z(k+1)-z(k)  )
     zp1 = abs(z(k)  -z(k-1))
     zp2 = abs(z(k)  -z(k-2))
     zp3 = abs(z(k)  -z(k-3))
     
     D2  = zp2*(-zp1+zp2+zm1)
     D1  =-zp1**2+zm1**2+zp2*zp1+zp2*zm1
     D0  = D2
     Dm1 = D1
     
     alpha_p2(k) = (zp2**3*zp1*zm1/(C_2_R*D2) - zp1**3*(zm1**2+zp2*zm1)/(C_2_R*D1) + &
          zm1**3*(zp1**2-zp2*zp1)/(C_2_R*Dm1))/                                      &
          ( zp3**3 - zp2**3*(-zp1*zp3-zp1*zm1+zp3**2+zm1*zp3)/D2 + &
          zp1**3*(-zm1**2-zp2*zp3-zp2*zm1+C_2_R*zp3**2)/D1       - &
          zm1**3*(-zp2*zp3-zp1**2+zp3**2+zp2*zp1)/Dm1            )

     alpha_p1(k) =-C_05_R*(zp1**2-zp2*zp1+C_2_R*alpha_p2(k)*(-zp2*zp3+zp3**2-zp1**2+zp2*zp1))/Dm1
     alpha(k)    = C_05_R*(-zp2*zp1+zp2**2+zp2*zm1+zp1*zm1-C_2_R*alpha_p2(k)*(zp1*zp3+zp1*zm1-zp3**2-zm1*zp3))/D0
     alpha_m1(k) = C_05_R*(zm1**2+zp2*zm1+C_2_R*alpha_p2(k)*(-zm1**2-zp2*zp3-zp2*zm1+zp3**2))/D1
     alpha_m2(k) =-C_05_R*(C_2_R*alpha_p2(k)*(-zp1*zp3-zp1*zm1+zp3**2+zm1*zp3)+zp1*zm1)/D2
     
     kmin_loc = 3
     kmax_loc = f%size-2

  ELSE
     kmin_loc = 1
     kmax_loc = f%size

  ENDIF

! #######################################################################
! inner points
! #######################################################################
  DO k = kmin_loc,kmax_loc
     km2 = k-2; km2 = MOD(km2+f%size-1,f%size) + 1
     km1 = k-1; km1 = MOD(km1+f%size-1,f%size) + 1
     kp1 = k+1; kp1 = MOD(kp1+f%size-1,f%size) + 1
     kp2 = k+2; kp2 = MOD(kp2+f%size-1,f%size) + 1

     zm2 = z(k)  -z(km2); IF ( zm2 .LT. C_0_R ) zm2 = zm2 + scalez
     zm1 = z(k)  -z(km1); IF ( zm1 .LT. C_0_R ) zm1 = zm1 + scalez
     zp1 = z(kp1)-  z(k); IF ( zp1 .LT. C_0_R ) zp1 = zp1 + scalez
     zp2 = z(kp2)-  z(k); IF ( zp2 .LT. C_0_R ) zp2 = zp2 + scalez

     D2  = zp2*(zp1-zp2-zm1)-(zp1*zm2+zm2**2-zm2*zm1)
     D1  = (zp1-zp2-zm1)*(zp1+zm1)
     D0  = zp2*(zp1-zp2-zm1)
     
     alpha_p2(k) = C_05_R*(zp1*zm1)/D2
     alpha_m2(k) = alpha_p2(k)
     alpha_p1(k) =-C_05_R*(C_2_R*alpha_m2(k)*zm2*(zm2+zp2)+zm1*(zp2+zm1))/D1
     alpha_m1(k) = C_05_R*(zp1**2-zp2*zp1+C_2_R*alpha_m2(k)*(zm2**2+zm2*zp2))/D1
     alpha(k)    = C_05_R*(C_1_R-(zp1*zm1+C_2_R*alpha_m2(k)*(zp1*zm2+zm2**2-zm1*zm2))/D0) - alpha_m2(k)
  ENDDO

  RETURN
END SUBROUTINE FLT_E4_INI

