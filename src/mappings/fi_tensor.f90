#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate properties of a real symmetric tensor.
!#
!########################################################################

!########################################################################
! eigenvalues: it assumes that no eigenvalue is zero, i.e. no zero determinant
!########################################################################
SUBROUTINE FI_TENSOR_EIGENVALUES(nx,ny,nz, tensor, result)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,6), INTENT(IN)  :: tensor ! a11, a22, a33, a12, a13, a23
  TREAL, DIMENSION(nx*ny*nz,3), INTENT(OUT) :: result ! eigenvalues: lambda1>lambda2>lambda3
  
! -----------------------------------------------------------------------
  TINTEGER ij
  TREAL d, e, lambda1, lambda2, lambda3, r13
  TCOMPLEX root1, root2, root3, unit

! #######################################################################
  r13  = C_1_R/C_3_R
  unit = (-0.5d0,0.866025403784439d0)

! #######################################################################
! Define coefficients of characteristic polynomial
! lambda^3 + a lambda^2 + b lambda + c = 0
! result1 contains a, result2 b and result3 c
! #######################################################################
  DO ij = 1,nx*ny*nz
     result(ij,1) =-(tensor(ij,1)+tensor(ij,2)+tensor(ij,3))
     result(ij,2) = &
            tensor(ij,1)*tensor(ij,2) + tensor(ij,1)*tensor(ij,3) + tensor(ij,2)*tensor(ij,3) &
          - tensor(ij,4)**2 - tensor(ij,5)**2 - tensor(ij,6)**2
     result(ij,3) = &
          - tensor(ij,1)*(tensor(ij,2)*tensor(ij,3)-tensor(ij,6)*tensor(ij,6)) &
          + tensor(ij,4)*(tensor(ij,4)*tensor(ij,3)-tensor(ij,6)*tensor(ij,5)) &
          + tensor(ij,5)*(tensor(ij,4)*tensor(ij,6)-tensor(ij,2)*tensor(ij,5))
  ENDDO

! #######################################################################
! Solve cubic
! #######################################################################
  DO ij = 1,nx*ny*nz
! Depressed eqn (modified) coefficients
     d = result(ij,2) - result(ij,1)*result(ij,1)/C_3_R 
     d = d/C_3_R
     e = result(ij,3) + result(ij,1)*(C_2_R*result(ij,1)*result(ij,1)-C_9_R*result(ij,2))/C_27_R
     e = e/C_2_R

! roots in complex plane
! you could check e**2 + d**3 <= 0, condition for real solns.
     root1 = e**2 + d**3; root1 =-e+sqrt(root1); root1 = root1**r13 ! fundamental root
     root2 = root1*unit
     root3 = root2*unit

     root1 =-d/root1 + root1 - result(ij,1)/C_3_R
     root2 =-d/root2 + root2 - result(ij,1)/C_3_R
     root3 =-d/root3 + root3 - result(ij,1)/C_3_R

! we know they are real
     lambda1 = real(root1)
     lambda2 = real(root2)
     lambda3 = real(root3)

! order from max to min
     d = -result(ij,1) ! trace = lambda1+lambda2+lambda3
     result(ij,1) = max(lambda1,max(lambda2,lambda3))
     result(ij,3) = min(lambda1,min(lambda2,lambda3))
     result(ij,2) = d - result(ij,1) - result(ij,3)

   ENDDO
     
  RETURN
END SUBROUTINE FI_TENSOR_EIGENVALUES

!########################################################################
! Calculate two first eigenvectors; the third can be computed a posteriori with cross product and it is not done here to save space.
! Global orientation is positive Ox for the first one, positive Oy for the second one.
!########################################################################
SUBROUTINE FI_TENSOR_EIGENFRAME(nx,ny,nz, tensor, lambda)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,6), INTENT(INOUT) :: tensor ! in: a11, a22, a33, a12, a13, a23
                                                        ! out: First two normalized eigenvectors (positions 1-3 and 4-6) 
  TREAL, DIMENSION(nx*ny*nz,3), INTENT(IN)    :: lambda ! in: Eigenvalues. The third one is not used.
  
! -----------------------------------------------------------------------
  TINTEGER ij, il
  TREAL e1(2), e2(2), e3(2), d

! #######################################################################
  DO ij = 1,nx*ny*nz
! Loop in number of eigenvalues
     DO il = 1,2

! If two first eqns are l.i.
        d = (tensor(ij,1)-lambda(ij,il))*(tensor(ij,2)-lambda(ij,il))-tensor(ij,4)*tensor(ij,4)
        IF ( d .NE. C_0_R ) THEN
           e1(il) =-( tensor(ij,5)*(tensor(ij,2)-lambda(ij,il)) - tensor(ij,6)*tensor(ij,4) )/d
           e2(il) =-( tensor(ij,6)*(tensor(ij,1)-lambda(ij,il)) - tensor(ij,5)*tensor(ij,4) )/d
           e3(il) = C_1_R
        ELSE
! then two last eqns are l.i.
           d = (tensor(ij,3)-lambda(ij,il))*(tensor(ij,2)-lambda(ij,il))-tensor(ij,6)*tensor(ij,6)
           e1(il) = C_1_R
           e2(il) =-( tensor(ij,4)*(tensor(ij,3)-lambda(ij,il)) - tensor(ij,6)*tensor(ij,5) )/d
           e3(il) =-( tensor(ij,5)*(tensor(ij,2)-lambda(ij,il)) - tensor(ij,6)*tensor(ij,4) )/d
        ENDIF

     ENDDO

! normalize, orient and store in array tensor
     d = sign(C_1_R,e1(1))/SQRT(e1(1)*e1(1)+e2(1)*e2(1)+e3(1)*e3(1)) ! orientation positive Ox
     tensor(ij,3) = e3(1)*d
     tensor(ij,2) = e2(1)*d
     tensor(ij,1) = e1(1)*d
     d = sign(C_1_R,e2(2))/SQRT(e1(2)*e1(2)+e2(2)*e2(2)+e3(2)*e3(2)) ! orientation positive Oy
     tensor(ij,6) = e3(2)*d
     tensor(ij,5) = e2(2)*d
     tensor(ij,4) = e1(2)*d

  ENDDO
     
  RETURN
END SUBROUTINE FI_TENSOR_EIGENFRAME
