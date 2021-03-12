!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2004/11/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Inverse error function.
!# From Abramowitz & Stegun, "Handbook of mathematical functions", page 933
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################

#include "types.h"

FUNCTION ERFINV(X)

  IMPLICIT NONE
  
  TREAL X, ERFINV
  
! -----------------------------------------------------------------------
  TREAL P, T, NUM, DEN
  TINTEGER IER

! #######################################################################
#ifdef SINGLE_PREC
#define C0_LOC 2.515517e0
#define C1_LOC 0.802853e0
#define C2_LOC 0.010328e0

#define D1_LOC 1.432788e0
#define D2_LOC 0.189269e0
#define D3_LOC 0.001308e0

#define FACTOR 0.707106781e0 

#else
#define C0_LOC 2.515517d0
#define C1_LOC 0.802853d0
#define C2_LOC 0.010328d0

#define D1_LOC 1.432788d0
#define D2_LOC 0.189269d0
#define D3_LOC 0.001308d0

#define FACTOR 0.707106781d0 

#endif

  IER = 0
  IF ( ABS(X) .EQ. C_1_R ) IER=1
  IF ( IER .NE. C_0_R ) GO TO 999

  IF ( X .EQ. C_0_R ) THEN
     ERFINV = C_0_R
  ELSE
     P = (C_1_R-ABS(X))*C_05_R
     T = SQRT(-C_2_R*LOG(P))
     NUM = C0_LOC + T*(C1_LOC+ T* C2_LOC             )
     DEN = C_1_R +  T*(D1_LOC+ T*(D2_LOC+ T* D3_LOC ))
     IF ( X .LT. 0 ) THEN; ERFINV =-FACTOR*(T-NUM/DEN)
     ELSE;                 ERFINV = FACTOR*(T-NUM/DEN); ENDIF
  ENDIF

999 IF ( IER .NE. 0 ) THEN
     WRITE(*,*) 'mathlib : erfinv : ier = ',IER
     STOP
  ENDIF

  RETURN
END FUNCTION ERFINV
