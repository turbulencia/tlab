!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/08/09 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Multiplication of pentadiagonal arrays A and B to give the 
!# nonadiadonal matrix C=AB
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE PENTADMT(nmax, a, b, c)

  IMPLICIT NONE

#include "types.h"

  TINTEGER nmax
  TREAL a(nmax,5), b(nmax,5), c(nmax,9)

! -----------------------------------------------------------------------
  TINTEGER n

! #######################################################################
#define A_a(n) a(n,1)
#define A_b(n) a(n,2)
#define A_c(n) a(n,3)
#define A_d(n) a(n,4)
#define A_e(n) a(n,5)

#define B_a(n) b(n,1)
#define B_b(n) b(n,2)
#define B_c(n) b(n,3)
#define B_d(n) b(n,4)
#define B_e(n) b(n,5)

#define C_a(n) c(n,1)
#define C_b(n) c(n,2)
#define C_c(n) c(n,3)
#define C_d(n) c(n,4)
#define C_e(n) c(n,5)
#define C_f(n) c(n,6)
#define C_g(n) c(n,7)
#define C_h(n) c(n,8)
#define C_i(n) c(n,9)

  n=1
  C_a(n) = C_0_R ! padding
  C_b(n) = C_0_R ! padding
  C_c(n) = C_0_R ! padding
  C_d(n) = C_0_R ! padding
  C_e(n) =                                   A_c(n)*B_c(n) +A_d(n)*B_b(n+1) +A_e(n)*B_a(n+2)
  C_f(n) =                                   A_c(n)*B_d(n) +A_d(n)*B_c(n+1) +A_e(n)*B_b(n+2)
  C_g(n) =                                   A_c(n)*B_e(n) +A_d(n)*B_d(n+1) +A_e(n)*B_c(n+2)
  C_h(n) =                                                  A_d(n)*B_e(n+1) +A_e(n)*B_d(n+2)
  C_i(n) =                                                                   A_e(n)*B_e(n+2)

  n=2
  C_a(n) = C_0_R ! padding
  C_b(n) = C_0_R ! padding
  C_c(n) = C_0_R ! padding
  C_d(n) =                  A_b(n)*B_c(n-1) +A_c(n)*B_b(n) +A_d(n)*B_a(n+1) 
  C_e(n) =                  A_b(n)*B_d(n-1) +A_c(n)*B_c(n) +A_d(n)*B_b(n+1) +A_e(n)*B_a(n+2)
  C_f(n) =                  A_b(n)*B_e(n-1) +A_c(n)*B_d(n) +A_d(n)*B_c(n+1) +A_e(n)*B_b(n+2)
  C_g(n) =                                   A_c(n)*B_e(n) +A_d(n)*B_d(n+1) +A_e(n)*B_c(n+2)
  C_h(n) =                                                  A_d(n)*B_e(n+1) +A_e(n)*B_d(n+2)
  C_i(n) =                                                                   A_e(n)*B_e(n+2)

  n=3
  C_a(n) = C_0_R ! padding
  C_b(n) = C_0_R ! padding
  C_c(n) = A_a(n)*B_c(n-2) +A_b(n)*B_b(n-1) +A_c(n)*B_a(n) 
  C_d(n) = A_a(n)*B_d(n-2) +A_b(n)*B_c(n-1) +A_c(n)*B_b(n) +A_d(n)*B_a(n+1) 
  C_e(n) = A_a(n)*B_e(n-2) +A_b(n)*B_d(n-1) +A_c(n)*B_c(n) +A_d(n)*B_b(n+1) +A_e(n)*B_a(n+2)
  C_f(n) =                  A_b(n)*B_e(n-1) +A_c(n)*B_d(n) +A_d(n)*B_c(n+1) +A_e(n)*B_b(n+2)
  C_g(n) =                                   A_c(n)*B_e(n) +A_d(n)*B_d(n+1) +A_e(n)*B_c(n+2)
  C_h(n) =                                                  A_d(n)*B_e(n+1) +A_e(n)*B_d(n+2)
  C_i(n) =                                                                   A_e(n)*B_e(n+2)

  n=4
  C_a(n) = C_0_R ! padding
  C_b(n) = A_a(n)*B_b(n-2) +A_b(n)*B_a(n-1) 
  C_c(n) = A_a(n)*B_c(n-2) +A_b(n)*B_b(n-1) +A_c(n)*B_a(n) 
  C_d(n) = A_a(n)*B_d(n-2) +A_b(n)*B_c(n-1) +A_c(n)*B_b(n) +A_d(n)*B_a(n+1) 
  C_e(n) = A_a(n)*B_e(n-2) +A_b(n)*B_d(n-1) +A_c(n)*B_c(n) +A_d(n)*B_b(n+1) +A_e(n)*B_a(n+2)
  C_f(n) =                  A_b(n)*B_e(n-1) +A_c(n)*B_d(n) +A_d(n)*B_c(n+1) +A_e(n)*B_b(n+2)
  C_g(n) =                                   A_c(n)*B_e(n) +A_d(n)*B_d(n+1) +A_e(n)*B_c(n+2)
  C_h(n) =                                                  A_d(n)*B_e(n+1) +A_e(n)*B_d(n+2)
  C_i(n) =                                                                   A_e(n)*B_e(n+2)

  DO n = 5,nmax-4
     C_a(n) = A_a(n)*B_a(n-2) 
     C_b(n) = A_a(n)*B_b(n-2) +A_b(n)*B_a(n-1) 
     C_c(n) = A_a(n)*B_c(n-2) +A_b(n)*B_b(n-1) +A_c(n)*B_a(n) 
     C_d(n) = A_a(n)*B_d(n-2) +A_b(n)*B_c(n-1) +A_c(n)*B_b(n) +A_d(n)*B_a(n+1) 
     C_e(n) = A_a(n)*B_e(n-2) +A_b(n)*B_d(n-1) +A_c(n)*B_c(n) +A_d(n)*B_b(n+1) +A_e(n)*B_a(n+2)
     C_f(n) =                  A_b(n)*B_e(n-1) +A_c(n)*B_d(n) +A_d(n)*B_c(n+1) +A_e(n)*B_b(n+2)
     C_g(n) =                                   A_c(n)*B_e(n) +A_d(n)*B_d(n+1) +A_e(n)*B_c(n+2)
     C_h(n) =                                                  A_d(n)*B_e(n+1) +A_e(n)*B_d(n+2)
     C_i(n) =                                                                   A_e(n)*B_e(n+2)
  ENDDO

  n=nmax-3
  C_a(n) = A_a(n)*B_a(n-2) 
  C_b(n) = A_a(n)*B_b(n-2) +A_b(n)*B_a(n-1) 
  C_c(n) = A_a(n)*B_c(n-2) +A_b(n)*B_b(n-1) +A_c(n)*B_a(n) 
  C_d(n) = A_a(n)*B_d(n-2) +A_b(n)*B_c(n-1) +A_c(n)*B_b(n) +A_d(n)*B_a(n+1) 
  C_e(n) = A_a(n)*B_e(n-2) +A_b(n)*B_d(n-1) +A_c(n)*B_c(n) +A_d(n)*B_b(n+1) +A_e(n)*B_a(n+2)
  C_f(n) =                  A_b(n)*B_e(n-1) +A_c(n)*B_d(n) +A_d(n)*B_c(n+1) +A_e(n)*B_b(n+2)
  C_g(n) =                                   A_c(n)*B_e(n) +A_d(n)*B_d(n+1) +A_e(n)*B_c(n+2)
  C_h(n) =                                                  A_d(n)*B_e(n+1) +A_e(n)*B_d(n+2)
  C_i(n) = C_0_R ! padding

  n=nmax-2
  C_a(n) = A_a(n)*B_a(n-2) 
  C_b(n) = A_a(n)*B_b(n-2) +A_b(n)*B_a(n-1) 
  C_c(n) = A_a(n)*B_c(n-2) +A_b(n)*B_b(n-1) +A_c(n)*B_a(n) 
  C_d(n) = A_a(n)*B_d(n-2) +A_b(n)*B_c(n-1) +A_c(n)*B_b(n) +A_d(n)*B_a(n+1) 
  C_e(n) = A_a(n)*B_e(n-2) +A_b(n)*B_d(n-1) +A_c(n)*B_c(n) +A_d(n)*B_b(n+1) +A_e(n)*B_a(n+2)
  C_f(n) =                  A_b(n)*B_e(n-1) +A_c(n)*B_d(n) +A_d(n)*B_c(n+1) +A_e(n)*B_b(n+2)
  C_g(n) =                                   A_c(n)*B_e(n) +A_d(n)*B_d(n+1) +A_e(n)*B_c(n+2)
  C_h(n) = C_0_R ! padding
  C_i(n) = C_0_R ! padding

  n=nmax-1
  C_a(n) = A_a(n)*B_a(n-2) 
  C_b(n) = A_a(n)*B_b(n-2) +A_b(n)*B_a(n-1) 
  C_c(n) = A_a(n)*B_c(n-2) +A_b(n)*B_b(n-1) +A_c(n)*B_a(n) 
  C_d(n) = A_a(n)*B_d(n-2) +A_b(n)*B_c(n-1) +A_c(n)*B_b(n) +A_d(n)*B_a(n+1) 
  C_e(n) = A_a(n)*B_e(n-2) +A_b(n)*B_d(n-1) +A_c(n)*B_c(n) +A_d(n)*B_b(n+1) 
  C_f(n) =                  A_b(n)*B_e(n-1) +A_c(n)*B_d(n) +A_d(n)*B_c(n+1) 
  C_g(n) = C_0_R ! padding
  C_h(n) = C_0_R ! padding
  C_i(n) = C_0_R ! padding

  n=nmax
  C_a(n) = A_a(n)*B_a(n-2) 
  C_b(n) = A_a(n)*B_b(n-2) +A_b(n)*B_a(n-1) 
  C_c(n) = A_a(n)*B_c(n-2) +A_b(n)*B_b(n-1) +A_c(n)*B_a(n) 
  C_d(n) = A_a(n)*B_d(n-2) +A_b(n)*B_c(n-1) +A_c(n)*B_b(n)
  C_e(n) = A_a(n)*B_e(n-2) +A_b(n)*B_d(n-1) +A_c(n)*B_c(n)
  C_f(n) = C_0_R ! padding
  C_g(n) = C_0_R ! padding
  C_h(n) = C_0_R ! padding
  C_i(n) = C_0_R ! padding

  RETURN
END SUBROUTINE PENTADMT
