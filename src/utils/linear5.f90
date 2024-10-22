!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/05/21 - J.P. Mellado
!#              Created
!# 2011/11/01 - C. Ansorge
!#              OpenMP Optimization
!#
!########################################################################
!# DESCRIPTION
!#
!# Solve LU of pentadiagonal system
!#
!########################################################################
!# ARGUMENTS 
!#
!# nmax    In     Size of the pentadiagonal system
!# len     In     Number of simultaneous systems to be solved
!# f       In     Array with forcing term
!#         Out    Array with solution
!#
!########################################################################
#include "types.h"

! #######################################################################
! LU factorization stage
! #######################################################################
SUBROUTINE PENTADFS(nmax, a,b,c,d,e)

  IMPLICIT NONE

  TINTEGER nmax
  TREAL, DIMENSION(nmax), INTENT(INOUT) :: a,b,c,d,e

! -----------------------------------------------------------------------
  TINTEGER n

! #######################################################################
!   a(1) = C_1_R ! padding
!   b(1) = C_1_R ! padding

!   a(2) = C_1_R ! padding
  b(2) = b(2)                         /c(1)
  c(2) = c(2)             - b(2)*d(1)
  d(2) = d(2) - b(2)*e(1)

  DO n = 3,nmax-1
     a(n) = a(n)                             /c(n-2)
     b(n) =(b(n)               - a(n)*d(n-2))/c(n-1)
     c(n) = c(n) - b(n)*d(n-1) - a(n)*e(n-2)
     d(n) = d(n) - b(n)*e(n-1)
  ENDDO
!   e(n-1) = C_1_R ! padding
  
  a(nmax) = a(nmax)                                         /c(nmax-2)
  b(nmax) =(b(nmax)                     - a(nmax)*d(nmax-2))/c(nmax-1)
  c(nmax) = c(nmax) - b(nmax)*d(nmax-1) - a(nmax)*e(nmax-2)
!   d(nmax) = C_1_R ! padding
!   e(nmax) = C_1_R ! padding

! Final operations
  a(3:) =-a(3:)
  b(2:) =-b(2:)
  c = C_1_R/c
  d(:nmax-1) =-d(:nmax-1)
  e(:nmax-2) =-e(:nmax-2)

  RETURN
END SUBROUTINE PENTADFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE PENTADSS(nmax,len, a,b,c,d,e, f)

  IMPLICIT NONE

  TINTEGER nmax, len
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: a,b,c,d,e
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: f
! -----------------------------------------------------------------------
  TINTEGER n,l, omp_srt,omp_end,omp_siz

!! !$omp parallel default(none) & 
!! !$omp shared(f,a,b,c,d,e,nmax,len) &
!! !$omp private(l,n,omp_srt,omp_end,omp_siz)   

!  CALL DNS_OMP_PARTITION(len,omp_srt,omp_end,omp_siz)
  omp_srt = 1
  omp_end = len
  omp_siz = len 

! #######################################################################
! -----------------------------------------------------------------------
! Solve Ly=f, forward
! -----------------------------------------------------------------------
  n = 2
  DO l=omp_srt,omp_end
     f(l,n) = f(l,n) + f(l,n-1)*b(2)
  ENDDO
  
  DO n = 3,nmax
     DO l=omp_srt,omp_end
        f(l,n) = f(l,n) + f(l,n-1)*b(n) + f(l,n-2)*a(n)
     ENDDO
  ENDDO

! -----------------------------------------------------------------------
! Solve Ux=y, backward
! -----------------------------------------------------------------------
  n = nmax
  DO l=omp_srt,omp_end
     f(l,n) = f(l,n) *c(n)
  ENDDO

  n = nmax-1
  DO l=omp_srt,omp_end
     f(l,n) =(f(l,n) +f(l,n+1)*d(n))*c(n)
  ENDDO

  DO n = nmax-2,1,-1
     DO l=omp_srt,omp_end
        f(l,n) =(f(l,n) +f(l,n+1)*d(n) +f(l,n+2)*e(n))*c(n)
     ENDDO
  ENDDO
!! !$omp end parallel 
  RETURN
END SUBROUTINE PENTADSS

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/05/21 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Perform LE decomposition of a pentadiagonal system.
!# Reverse ordering.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################

! #######################################################################
! LU factorization stage
! #######################################################################
SUBROUTINE PENTADFS2(nmax, a,b,c,d,e)

  IMPLICIT NONE

  TINTEGER nmax
  TREAL, DIMENSION(nmax), INTENT(INOUT) :: a,b,c,d,e

! -----------------------------------------------------------------------
  TINTEGER n

! #######################################################################
  n = nmax
  e(n) = C_1_R ! padding
  d(n) = C_1_R ! padding
  c(n) = c(n) 
  b(n) = b(n) 

  n = nmax-1
  e(n) = C_1_R ! padding
  d(n) =(d(n)              )/c(n+1)
  c(n) = c(n) - d(n)*b(n+1) 
  b(n) = b(n) - d(n)*a(n+1)

  DO n = nmax-2,3,-1
     e(n) = e(n)                             /c(n+2)
     d(n) =(d(n)               - e(n)*b(n+2))/c(n+1)
     c(n) = c(n) - d(n)*b(n+1) - e(n)*a(n+2)
     b(n) = b(n) - d(n)*a(n+1)
  ENDDO

  n = 2
  e(n) = e(n)                             /c(n+2)
  d(n) =(d(n)               - e(n)*b(n+2))/c(n+1)
  c(n) = c(n) - d(n)*b(n+1) - e(n)*a(n+2)
  b(n) = b(n) - d(n)*a(n+1)
  a(n) = C_1_R ! padding
  
  n = 1
  e(n) = e(n)                             /c(n+2)
  d(n) =(d(n)               - e(n)*b(n+2))/c(n+1)
  c(n) = c(n) - d(n)*b(n+1) - e(n)*a(n+2)
  b(n) = C_1_R ! padding
  a(n) = C_1_R ! padding

  RETURN
END SUBROUTINE PENTADFS2

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE PENTADSS2(nmax,len, a,b,c,d,e, f)

  IMPLICIT NONE

  TINTEGER nmax, len
  TREAL, DIMENSION(nmax),     INTENT(IN)    :: a,b,c,d,e
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: f

! -----------------------------------------------------------------------
  TINTEGER n

! #######################################################################
! -----------------------------------------------------------------------
! Solve Ly=f, forward
! -----------------------------------------------------------------------
  n = nmax-1
  f(:,n) = f(:,n) - f(:,n+1)*d(n)

  DO n = nmax-2,1,-1
     f(:,n) = f(:,n) -f(:,n+1)*d(n) -f(:,n+2)*e(n)
  ENDDO

! -----------------------------------------------------------------------
! Solve Ux=y, backward
! -----------------------------------------------------------------------
  n = 1
  f(:,n) = f(:,n) /c(n)

  n = 2
  f(:,n) =(f(:,n) -f(:,n-1)*b(n))/c(n)

  DO n = 3,nmax
     f(:,n) =(f(:,n) -f(:,n-1)*b(n) -f(:,n-2)*a(n))/c(n)
  ENDDO

  RETURN
END SUBROUTINE PENTADSS2

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/09 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Pentadiagonal Toeplitz solver for a circulant system (periodic bcs),
!# based on the Sherman–Morrison–Woodbury Formula, described in: "A new 
!# algorithm for solving nearly penta-diagonal Toeplitz linear systems"
!# https://doi.org/10.1016/j.camwa.2011.12.044
!# Algorithm 2.2 is implemented here.
!# Ensure that Matrix M is invertible!
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "dns_error.h"
! #######################################################################
! LU factorization stage
! #######################################################################
SUBROUTINE PENTADPFS(nmax, a,b,c,d,e,f,g)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLab_WorkFlow

  IMPLICIT NONE

  TINTEGER,               INTENT(IN   ) :: nmax  
  TREAL, DIMENSION(nmax), INTENT(INOUT) :: a,b,c,d,e ! Diagonals
  TREAL, DIMENSION(nmax), INTENT(INOUT) :: f,g       ! Additional (u1,u2)

! -----------------------------------------------------------------------
  TREAL :: a0, b0, en, dn
  TREAL :: m1, m2, m3, m4

! #######################################################################
! Build regular modified pentadiagonal matrix A1 (eq. 2.8) 

! Save cyclic entries
  a0 = a(1)    ! upper-right corner
  b0 = b(1)  
  en = e(nmax) ! lower-left  corner
  dn = d(nmax) 

! Modified entries of A1
  b(2) = b(2) - d(nmax)
  c(1) = c(1) - e(nmax)
  c(2) = c(2) - e(nmax)
  
  c(nmax-1) = c(nmax-1) - a(1)
  c(nmax  ) = c(nmax  ) - a(1)
  d(nmax-1) = d(nmax-1) - b(1)

! Set off-digonal entries to zero for A1
  a(1) = C_0_R
  a(2) = C_0_R
  b(1) = C_0_R

  d(nmax  ) = C_0_R
  e(nmax  ) = C_0_R 
  e(nmax-1) = C_0_R 

! Regular forward step for A1
  CALL PENTADFS2(nmax,a,b,c,d,e)

! Save cyclic entries again
  a(1)    = a0 ! upper-right corner
  b(1)    = b0  
  e(nmax) = en ! lower-left corner
  d(nmax) = dn  

! Define matrix u [here: u1,u2 stored in additional diagonals f, g] (eq. 2.8)
  f         = C_0_R ! u1
  f(1)      = C_1_R
  f(nmax-1) = C_1_R
  g         = C_0_R ! u2
  g(2)      = C_1_R
  g(nmax)   = C_1_R

! Regular backward step for u1, u2
  CALL PENTADSS2(nmax, 1, a,b,c,d,e, f)
  CALL PENTADSS2(nmax, 1, a,b,c,d,e, g)

! Compute entries of matrix M[2x2] once
  m1 = e(nmax) * f(1) + a(1)    * f(nmax-1) + b(1) * f(nmax) + C_1_R
  m2 = e(nmax) * g(1) + a(1)    * g(nmax-1) + b(1) * g(nmax)
  m3 = d(nmax) * f(1) + e(nmax) * f(2)      + a(1) * f(nmax)
  m4 = d(nmax) * g(1) + e(nmax) * g(2)      + a(1) * g(nmax) + C_1_R
! Check if M is invertible (eq. 2.9)
  IF ( (m1*m4 - m2*m3) < 1e-8 ) THEN
    CALL TLAB_WRITE_ASCII(efile, 'FDM_INITIALIZE. Pendad - matrix M not invertible.')
    CALL TLAB_STOP(DNS_ERROR_PENTADP)
  ENDIF

  RETURN
END SUBROUTINE PENTADPFS

! #######################################################################
! Backward substitution step in the Thomas algorith
! #######################################################################
SUBROUTINE PENTADPSS(nmax,len, a,b,c,d,e,f,g, frc)

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN   ) :: nmax, len
  TREAL, DIMENSION(nmax),     INTENT(IN   ) :: a,b,c,d,e ! Diagonals
  TREAL, DIMENSION(nmax),     INTENT(IN   ) :: f,g       ! Additional (z,w)
  TREAL, DIMENSION(len,nmax), INTENT(INOUT) :: frc       ! Rhs

! -----------------------------------------------------------------------
  TINTEGER :: n, l
  TREAL    :: m1, m2, m3, m4
  TREAL    :: di, d11, d12, d13, d14, d21, d22, d23, d24
  TREAL    :: dummy1, dummy2

! #######################################################################

! Regular backward step for rhs
  CALL PENTADSS2(nmax,len, a,b,c,d,e, frc)

! Compute entries of matrix m[2x2]
  m1 = e(nmax) * f(1) + a(1)    * f(nmax-1) + b(1) * f(nmax) + C_1_R
  m2 = e(nmax) * g(1) + a(1)    * g(nmax-1) + b(1) * g(nmax)
  m3 = d(nmax) * f(1) + e(nmax) * f(2)      + a(1) * f(nmax)
  m4 = d(nmax) * g(1) + e(nmax) * g(2)      + a(1) * g(nmax) + C_1_R

! Compute coefficients
  di  = 1  / (m1*m4      - m2*m3)
  d11 = di * (m4*e(nmax) - m2*d(nmax))
  d12 = di * (m4*b(1)    - m2*a(1))
  d13 = di *  m4*a(1)
  d14 = di *  m2*e(nmax)
  d21 = di * (m1*d(nmax) - m3*e(nmax))
  d22 = di * (m1*a(1)    - m3*b(1))
  d23 = di *  m3*a(1)
  d24 = di *  m1*e(nmax)

! Solve
  DO n = 3,nmax-3,1 ! Main loop
    DO l=1,len,1
      dummy1   = d11*frc(l,1) + d12*frc(l,nmax) + d13*frc(l,nmax-1) - d14*frc(l,2)
      dummy2   = d21*frc(l,1) + d22*frc(l,nmax) - d23*frc(l,nmax-1) + d24*frc(l,2) 
      !
      frc(l,n) = frc(l,n) - dummy1*f(n) - dummy2*g(n)
    ENDDO
  ENDDO
  !
  DO l=1,len,1     ! Boundaries
    dummy1 = d11*frc(l,1) + d12*frc(l,nmax) + d13*frc(l,nmax-1) - d14*frc(l,2)
    dummy2 = d21*frc(l,1) + d22*frc(l,nmax) - d23*frc(l,nmax-1) + d24*frc(l,2) 
    DO n = 1,2,1
      frc(l,n) = frc(l,n) - dummy1*f(n) - dummy2*g(n)
    ENDDO
    DO n = nmax-2,nmax,1
      frc(l,n) = frc(l,n) - dummy1*f(n) - dummy2*g(n)
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE PENTADPSS