!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2000/07/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Random number generators.
!# From Press, Flannery, Teukolsky & Vettering, "Numerical Recipes"
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"

!########################################################################
! Gaussian pdf according to Box-Muller method
!########################################################################
FUNCTION RANG(mean, sigma, seed)

  IMPLICIT NONE
  
  TREAL mean, sigma
  TINTEGER seed
  TREAL RANG
  
! -----------------------------------------------------------------------
  TREAL v1,v2,r
  TREAL RAN0
      
! #######################################################################
1 v1 = C_2_R*RAN0(seed) - C_1_R
  v2 = C_2_R*RAN0(seed) - C_1_R
  r= v1*v1 + v2*v2
  IF ( r .GE. C_1_R ) GO TO 1
  v2 = v1*SQRT(-C_2_R*LOG(r)/r)
  
  RANG = mean + v2*sigma

  RETURN
END FUNCTION RANG

! #######################################################################
! Uniform pdf
! #######################################################################
FUNCTION RAN0(IDUM)

  IMPLICIT NONE

  TREAL RAN0
  TINTEGER IDUM

  TINTEGER IA, IM, IQ, IR, NTAB, NDIV
  TREAL AM, EPS, RNMX
  PARAMETER (IA=16807, IM=2147483647, AM=C_1_R/IM, &
       IQ=127773, IR=2836)
  PARAMETER (NTAB=32,NDIV=1+(IM-1)/NTAB, EPS=C_12EM7_R, &
       RNMX=C_1_R-EPS)

  TINTEGER j, k, iv(NTAB), iy
  SAVE iv, iy
  DATA iv /NTAB*0/, iy /0/

  IF ( IDUM .LE. 0 .OR. iy .EQ. 0 ) THEN
     IDUM = MAX(-IDUM,1)
     DO j=NTAB+8,1,-1
        k=IDUM/IQ
        IDUM=IA*(IDUM-k*IQ)-IR*k
        IF ( IDUM .LT. 0 ) IDUM = IDUM + IM
        IF ( j .LE. NTAB ) iv(j) = IDUM
     ENDDO
     iy = iv(1)
  ENDIF

  k=IDUM/IQ
  IDUM=IA*(IDUM-k*IQ)-IR*k
  IF ( IDUM .LT. 0 ) IDUM = IDUM + IM
  j = 1 + iy/NDIV
  iy = iv(j)
  iv(j) = IDUM
  RAN0 = MIN(AM*iy,RNMX)

  RETURN
END FUNCTION RAN0

! #######################################################################
! Uniform pdf
! #######################################################################
FUNCTION RAN1(IDUM)

  IMPLICIT NONE

  TREAL RAN1
  TINTEGER IDUM

  TINTEGER M1, M2, M3
  TINTEGER IA1, IA2, IA3
  TINTEGER IC1, IC2, IC3
  TREAL RM1, RM2
  TREAL R(97)

  TINTEGER IX1, IX2, IX3, IFF, J

  PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=C_1_R/M1)
  PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=C_1_R/M2)
  PARAMETER (M3=243000,IA3=4561,IC3=51349)
  DATA IFF /0/

  IF ( IDUM .LT. 0 .OR. IFF .EQ. 0 ) THEN
     IFF=1
     IX1=MOD(IC1-IDUM,M1)
     IX1=MOD(IA1*IX1+IC1,M2)
     IX2=MOD(IX1,M2)
     IX1=MOD(IA1*IX1+IC1,M1)
     IX3=MOD(IX1,M3)
     DO J=1,97
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IA2*IX2+IC2,M2)
        R(J)=(M_REAL(IX1)+M_REAL(IX2)*RM2)*RM1
     ENDDO
     IDUM=1
  ENDIF
  IX1=MOD(IA1*IX1+IC1,M1)
  IX2=MOD(IA2*IX2+IC2,M2)
  IX3=MOD(IA3*IX3+IC2,M3)
  J=1+(97*IX3)/M3
  IF ( J .GT. 97 .OR. J .LT. 1 ) THEN
     WRITE(*,*) 'Error J out of limits in RAN1'
     STOP
  ENDIF
  RAN1=R(J)
  R(J)=(M_REAL(IX1)+M_REAL(IX2)*RM2)*RM1

  RETURN
END FUNCTION RAN1
