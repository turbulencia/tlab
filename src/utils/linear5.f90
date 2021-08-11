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
  a(1) = C_1_R ! padding
  b(1) = C_1_R ! padding

  a(2) = C_1_R ! padding
  b(2) = b(2)                         /c(1)
  c(2) = c(2)             - b(2)*d(1)
  d(2) = d(2) - b(2)*e(1)

  DO n = 3,nmax-1
     a(n) = a(n)                             /c(n-2)
     b(n) =(b(n)               - a(n)*d(n-2))/c(n-1)
     c(n) = c(n) - b(n)*d(n-1) - a(n)*e(n-2)
     d(n) = d(n) - b(n)*e(n-1)
  ENDDO
  e(n-1) = C_1_R ! padding
  
  a(nmax) = a(nmax)                                         /c(nmax-2)
  b(nmax) =(b(nmax)                     - a(nmax)*d(nmax-2))/c(nmax-1)
  c(nmax) = c(nmax) - b(nmax)*d(nmax-1) - a(nmax)*e(nmax-2)
  d(nmax) = C_1_R ! padding
  e(nmax) = C_1_R ! padding

! Final operations
  a =-a
  b =-b
  c = C_1_R/c
  d =-d
  e =-e

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
