      SUBROUTINE GAUSS(n, A, b)

! Adpated from matlab code
! gauss.m
! Solves the system Ax=b for x using Gaussian elimination without
! pivoting.  The matrix A is replaced by the m_ij and U on exit, and
! the vector b is replaced by the solution x of the original system.

      IMPLICIT NONE

#include "types.h"

      TINTEGER n
      TREAL A(n,n), b(n)

      TINTEGER i, j, k
      TREAL tmp

!  -------------- FORWARD SWEEP --------------

      DO j = 1,n-1     ! For each column j<n,
         DO i = j+1,n  ! loop through the elements a_ij below the pivot a_jj.

! Compute m_ij.  Note that we can store m_ij in the location
! (below the diagonal!) that a_ij used to sit without disrupting
! the rest of the algorithm, as a_ij is set to zero by construction
! during this iteration.

            A(i,j)     = - A(i,j) / A(j,j)

! Add m_ij times the upper triangular part of the j'th row of
! the augmented matrix to the i'th row of the augmented matrix.
            
            DO k = j+1,n
               A(i,k) = A(i,k) + A(i,j)*A(j,k)
            ENDDO
            b(i) = b(i) + A(i,j)*b(j)
         ENDDO
      ENDDO

!  ------------ BACK SUBSTITUTION ------------

      b(n) = b(n) / A(n,n)   ! Initialize the backwards march
      DO i = n-1,1,-1

! Note that an inner product is performed at the multiplication
! sign here, accounting for all values of x already determined:

         tmp = C_0_R
         DO j = i+1,n
            tmp = tmp + A(i,j)*b(j)
         ENDDO
         b(i) = (b(i)-tmp) / A(i,i)
!     b(i) = ( b(i) - A(i,i+1:n) * b(i+1:n) ) / A(i,i)
      ENDDO

      RETURN
      END
