#include "types.h"

!########################################################################
!# Tool/Library Math
!#
!########################################################################
!# HISTORY
!#
!# 2002/03/01 - J.P. Mellado
!#              Created
!# 2004/09/01 - J.P. Mellado
!#              Modified
!#
!########################################################################
!# DESCRIPTION
!#
!# (Ln of ) Gamma function of XX>0. 
!# Full accuracy for XX>1; if 0<XX<1 use reflection
!# formula GAMM(1-XX)=PI*XX/(GAMM(1+XX)*SIN(PI*XX))
!# Created from Press, Flannery, Teukolsky, Vetterling, "Numerical Recipes"
!# 
!# Revision  to use Lanczos approximation from
!# Numerical Recipes in C (2nd ed. Cambridge University Press, 1992)
!#
!# Recall that G(z)=(z-1)*G(z-1)
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
FUNCTION GAMMLN(XX)
  
  IMPLICIT NONE

  TREAL XX,GAMMLN

! -----------------------------------------------------------------------
#ifdef SINGLE_PREC
#define C_5P5_L    5.5e+00
#define C_GAMM0_L  0.1000000000190015e+01
#define C_GAMM1_L  0.7618009172947146e+02
#define C_GAMM2_L -0.8650532032941677e+02
#define C_GAMM3_L  0.2401409824083091e+02
#define C_GAMM4_L -0.1231739572450155e+01
#define C_GAMM5_L  0.1208650973866179e-02
#define C_GAMM6_L -0.5395239384953e-05
#define C_GAMM7_L  0.250662827465e+01

#else
#define C_5P5_L    5.5d+00
#define C_GAMM0_L  0.1000000000190015d+01
#define C_GAMM1_L  0.7618009172947146d+02
#define C_GAMM2_L -0.8650532032941677d+02
#define C_GAMM3_L  0.2401409824083091d+02
#define C_GAMM4_L -0.1231739572450155d+01
#define C_GAMM5_L  0.1208650973866179d-02
#define C_GAMM6_L -0.5395239384953d-05
#define C_GAMM7_L  0.250662827465d+01

#endif

  TREAL COF(6),STP,TMP,SER,X
  TINTEGER J, IER
  
  DATA COF,STP / C_GAMM1_L, C_GAMM2_L, C_GAMM3_L, C_GAMM4_L, &
                 C_GAMM5_L, C_GAMM6_L, C_GAMM7_L             /

! #######################################################################
  IER = 0
  IF ( XX .LT. C_0_R ) IER=1
  IF ( IER .NE. C_0_R ) GO TO 999
  
  X = XX-C_1_R
  TMP = X+C_5P5_L
  TMP = (X+C_05_R)*LOG(TMP)-TMP
!  SER = C_1_R
  SER = C_GAMM0_L
  DO J = 1,6
     X = X+C_1_R
     SER = SER+COF(J)/X
  ENDDO
  GAMMLN=TMP+LOG(STP*SER)
  
999 IF ( IER .NE. 0 ) THEN
     WRITE(*,*) 'mathlib : gammln : ier = ',IER
     STOP
  ENDIF
  
  RETURN
END FUNCTION GAMMLN

!########################################################################
!# Tool/Library Math
!#
!########################################################################
!# HISTORY
!#
!# 2002/03/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Incomplete Gamma function, P(A,X). 
!# From Press, Flannery, Teukolsky, Vetterling, "Numerical recipes"
!# 
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################

FUNCTION GAMMP(A,X)

  IMPLICIT NONE

  TREAL GAMMP, A, X
  
! -----------------------------------------------------------------------
  TINTEGER IER
  TREAL GAMMP_SER, GAMMP_CON

! #######################################################################
  IER = 0
  IF ( X .LT. C_0_R ) IER=1
  IF ( A .LT. C_0_R ) IER=1
  IF ( IER .NE. C_0_R ) GO TO 999
  
  IF ( X .LT. A+C_1_R ) THEN
! use series representation
     GAMMP = GAMMP_SER(A,X)
  ELSE
! use continued fraction representation and
! take its complement
     GAMMP = C_1_R - GAMMP_CON(A,X)
  ENDIF
  
999 IF ( IER .NE. 0 ) THEN
     WRITE(*,*) 'mathlib : gammp : ier = ',IER
     STOP
  ENDIF
  
  RETURN
END FUNCTION GAMMP

! #######################################################################
! Incomplete Gamma function, P(A,X), evaluated as a series.
! #######################################################################
FUNCTION GAMMP_SER(A,X)
  
  IMPLICIT NONE
      
  TREAL GAMMP_SER, A, X
  
! -----------------------------------------------------------------------
  TINTEGER ITMAX, IER, N
  TREAL EPS, GLN, AP, SUM, DEL
  TREAL GAMMLN
  PARAMETER(ITMAX=100,EPS=C_1EM8_R)
  
! #######################################################################
  IER=0
  
  IF (X .EQ. C_0_R) THEN
     GAMMP_SER = C_0_R
  ELSE
     GLN=GAMMLN(A)
     AP=A
     SUM=C_1_R/A
     DEL=SUM
     DO N = 1,ITMAX
        AP=AP+C_1_R
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF ( ABS(DEL) .LT. ABS(SUM)*EPS ) EXIT
     ENDDO
     IF ( N .GT. ITMAX ) IER=1
     GAMMP_SER=SUM*EXP(-X+A*LOG(X)-GLN)
  ENDIF
  
999 IF ( IER .NE. 0 ) THEN
     WRITE(*,*) 'mathlib : gammp_ser : ier = ',IER
     STOP
  ENDIF
  
  RETURN
END FUNCTION GAMMP_SER

! #######################################################################
! Incomplete Gamma function, P(A,X), evaluated as a continuous fraction representation. 
! #######################################################################
FUNCTION GAMMP_CON(A,X)
  
  IMPLICIT NONE
  
  TREAL GAMMP_CON, A, X

! -----------------------------------------------------------------------
  TINTEGER ITMAX, N, IER
  TREAL EPS, GLN, G, GOLD, A0, A1, B0, B1, FAC, AN, ANA, ANF
  TREAL GAMMLN
  PARAMETER(ITMAX=100,EPS=C_1EM8_R)
  
! #######################################################################
  IER = 0
  
  GLN=GAMMLN(A)
  GOLD=C_0_R
  A0=C_1_R
  A1=X
  B0=C_0_R
  B1=C_1_R
  FAC=C_1_R
  DO N = 1,ITMAX
     AN=M_REAL(N)
     ANA=AN-A
     A0=(A1+A0*ANA)*FAC
     B0=(B1+B0*ANA)*FAC
     ANF=AN*FAC
     A1=X*A0+ANF*A1
     B1=X*B0+ANF*B1
     IF ( A1 .NE. C_0_R ) THEN
        FAC=C_1_R/A1
        G=B1*FAC
        IF ( ABS((G-GOLD)/G) .LT. EPS ) GO TO 999
        GOLD=G
     ENDIF
  ENDDO
  IF ( N .GT. ITMAX ) IER=1
  
999 IF ( IER .NE. 0 ) THEN
     WRITE(*,*) 'mathlib : gammp_ser : ier = ',IER
     STOP
  ENDIF

  GAMMP_CON=EXP(-X+A*LOG(X)-GLN)*G
  
  RETURN
END FUNCTION GAMMP_CON

!########################################################################
!# Tool/Library Math
!#
!########################################################################
!# HISTORY
!#
!# 2002/03/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Incomplete Beta function, B(A,B,X). 
!# From Press, Flannery, Teukolsky, Vetterling, "Numerical recipes"
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
FUNCTION BETAI(A,B,X)

  IMPLICIT NONE

  TREAL BETAI,A,B,X

! -----------------------------------------------------------------------
  TINTEGER IER
  TREAL BT, GAMMLN, BETAFC, XSYM

! #######################################################################
  IER=0
  IF ( X .LT. C_0_R ) IER=1
  IF ( X .GT. C_1_R ) IER=1
  IF ( IER .NE. 0 ) GO TO 999

  XSYM = C_1_R-X

  IF ( X .EQ. C_0_R .OR. X .EQ. C_1_R ) THEN
     BT = C_0_R
  ELSE
     BT = EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)+ A*LOG(X)+B*LOG(XSYM))
  ENDIF

  IF ( X .LT. (A+C_1_R)/(A+B+C_2_R) ) THEN
     BETAI = BT*BETAFC(A,B,X)/A
  ELSE
     BETAI = C_1_R-BT*BETAFC(B,A,XSYM)/B
  ENDIF

999 IF ( IER .NE. 0 ) THEN
     WRITE(*,*) 'mathlib : betai : ier = ',IER
     STOP
  ENDIF

  RETURN
END FUNCTION BETAI

! #######################################################################
! #######################################################################
FUNCTION BETAFC(A,B,X)

  IMPLICIT NONE

  TREAL BETAFC, A,B,X

  TREAL AM, BM, AZ, BZ, AP, BP, QAB, QAP, QAM
  TREAL EM, TEM, AOLD, D, APP, BPP

  TINTEGER i, n, IER
  TREAL EPS

  n = 100
  EPS = C_1EM8_R
  IER = 0

  AM=C_1_R
  BM=C_1_R
  AZ=C_1_R
  QAB=A+B
  QAP=A+C_1_R
  QAM=A-C_1_R
  BZ=C_1_R-QAB*X/QAP

  DO i = 1,n
     EM=M_REAL(i)
     TEM=EM+EM
     D=EM*(B-EM)*X/((QAM+TEM)*(A+TEM))
     AP=AZ+D*AM
     BP=BZ+D*BM
     D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
     APP=AP+D*AZ
     BPP=BP+D*BZ
     AOLD=AZ
     AM=AP/BPP
     BM=BP/BPP
     AZ=APP/BPP
     BZ=C_1_R
     IF ( ABS(AZ-AOLD) .LT. EPS*ABS(AZ) ) EXIT
  ENDDO

  IF ( i .GT. n ) IER=1
  IF ( IER .NE. 0 ) GO TO 999

  BETAFC=AZ

999 IF ( IER .NE. 0 ) THEN
     WRITE(*,*) 'mathlib : betacf : ier = ',IER
     STOP
  ENDIF

  RETURN
END FUNCTION BETAFC
