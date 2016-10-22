#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2013/01/20 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Compact 2nd-order 6th-order Nonuniform Direct
!#
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme for non-uniform grids
!# according to Shukla and Zhong (2004).
!#
!# Equations (15) for the interior points (beware that there is a typo in paper).
!# Equations (16) for the first/last points. 
!# Table B.2 for the second/second-to-last points.
!#
!# We set up mmax linear system of size nmax of the form
!#                         A up = B u 
!# in this routine. The matrix A is tridiagonal, and B is pentadiagonal,
!# except for the boundary points.
!#
!# We normalize system to get a diagonal 1 in B.
!# 
!########################################################################

!########################################################################
! Create diagonals
!########################################################################
SUBROUTINE FDM_C2N6ND_INITIALIZE(nmax, x, lhs, rhs)

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN) :: nmax
  TREAL,   DIMENSION(nmax),   INTENT(IN) :: x
  TREAL,   DIMENSION(nmax,3), INTENT(OUT):: lhs ! LHS diagonals (#=3)
  TREAL,   DIMENSION(nmax,4), INTENT(OUT):: rhs ! RHS diagonals (#=5-1 because of normalization)

! -------------------------------------------------------------------
  TINTEGER n, j
  TREAL am1, a, ap1                            ! Left-hand side
  TREAL bm2, bm1, b, bp1, bp2, bp3             ! Right-hand side 
  TREAL tmp1,tmp2, tmpa,tmpb,tmpd, tmpp,tmpq   ! Intermediate ops
  TREAL dummy1, dummy2

! #######################################################################
! n = 1
! #######################################################################
  n = 1

! left-hand side
  am1 = C_0_R
  a   = C_1_R
  ap1 = (x(n+1)-x(n))*(x(n)-x(n+2)) + (x(n+1)-x(n))*(x(n)-x(n+3)) - (x(n)-x(n+2))*(x(n)-x(n+3))
  ap1 = ap1 / &
      ( (x(n+1)-x(n))*(x(n+1)-x(n+2)) + (x(n+1)-x(n))*(x(n+1)-x(n+3)) + (x(n+1)-x(n+2))*(x(n+1)-x(n+3)) )

! right-hand side
  bm2 = C_0_R
  bm1 = C_0_R
  b   = (x(n+1)-x(n+2))*(x(n+1)-x(n+3)) + C_2_R*(x(n+1)-x(n))*( (x(n+1)-x(n+2)) + (x(n+1)-x(n+3)) )
  b   = b / &
      ( (x(n+1)-x(n+2))*(x(n+1)-x(n+3)) +       (x(n+1)-x(n))*( (x(n+1)-x(n+2)) + (x(n+1)-x(n+3)) ) )
  b   = b *( (x(n)-x(n+2)) + (x(n)-x(n+3)) ) / &
        (x(n)-x(n+1)) + C_1_R
  b   = b / &
      ( (x(n)-x(n+2))*(x(n)-x(n+3)) ) 
  tmp1= ( (x(n+1)-x(n+2)) + (x(n+1)-x(n+3)) ) / &
      ( (x(n+1)-x(n+2))*(x(n+1)-x(n+3)) +       (x(n+1)-x(n))*( (x(n+1)-x(n+2)) + (x(n+1)-x(n+3)) ) )
  b   = ( b + tmp1 / (x(n+1)-x(n)) ) *C_2_R

  bp1 = (x(n)-x(n+2)) + (x(n)-x(n+3)) + ap1 *( (x(n+1)-x(n)) + (x(n+1)-x(n+2)) + (x(n+1)-x(n+3)) )
  bp1 = bp1 *C_2_R / &
      ( (x(n+1)-x(n))*(x(n+1)-x(n+2))*(x(n+1)-x(n+3)) )
     
  bp2 = (x(n)-x(n+1))*(       (x(n+1)-x(n)) + C_2_R*(x(n+1)-x(n+3)) ) + &
        (x(n)-x(n+3))*( C_2_R*(x(n+1)-x(n)) + C_3_R*(x(n+1)-x(n+3)) )
  bp2 = bp2 / &
      ( (x(n+1)-x(n))*(x(n+1)-x(n+3)) - (x(n+2)-x(n+1))*( (x(n+1)-x(n)) + (x(n+1)-x(n+3)) ) )
  bp2 = bp2 *C_2_R*(x(n+1)-x(n)) / &
      ( (x(n+2)-x(n+1))*(x(n+2)-x(n))*(x(n+2)-x(n+3)) )

  bp3 = (x(n)-x(n+1))*(       (x(n+1)-x(n)) + C_2_R*(x(n+1)-x(n+2)) ) + &
        (x(n)-x(n+2))*( C_2_R*(x(n+1)-x(n)) + C_3_R*(x(n+1)-x(n+2)) )
  bp3 = bp3 / &
      ( (x(n+1)-x(n))*(x(n+1)-x(n+2)) - (x(n+3)-x(n+1))*( (x(n+1)-x(n)) + (x(n+1)-x(n+2)) ) )
  bp3 = bp3 *C_2_R*(x(n+1)-x(n)) / &
      ( (x(n+3)-x(n+1))*(x(n+3)-x(n))*(x(n+3)-x(n+2)) )

! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
  tmp1 = C_1_R/b

  lhs(n,1) = C_0_R
  lhs(n,2) = a   *tmp1 
  lhs(n,3) = ap1 *tmp1 
  
  rhs(n,1) = bp3 *tmp1 ! saving here the last term that goes into the 3. superdiagonal
  rhs(n,2) = C_0_R
  rhs(n,3) = bp1 *tmp1 
  rhs(n,4) = bp2 *tmp1 

! #######################################################################
! n = 2
! #######################################################################
  n = 2

  bm2 = C_0_R; bp2 = C_0_R 

  tmp1 = C_1_R / &
       ( (x(n+1)-x(n))*(x(n)-x(n-1)) + (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) )
  
  bm1 = (x(n+1)-x(n)) / (x(n+1)-x(n-1))
  b   =-C_1_R
  bp1 = (x(n)-x(n-1)) / (x(n+1)-x(n-1))

  am1 = (x(n)-x(n-1))*(x(n)-x(n-1)) + (x(n+1)-x(n))*(x(n)-x(n-1)) - (x(n+1)-x(n))*(x(n+1)-x(n))
  am1 = am1 *bm1 *tmp1
  a   = C_1_R
  ap1 = (x(n+1)-x(n))*(x(n+1)-x(n)) + (x(n+1)-x(n))*(x(n)-x(n-1)) - (x(n)-x(n-1))*(x(n)-x(n-1))
  ap1 = ap1 *bp1 *tmp1

  bm1 = bm1 *C_12_R *tmp1
  b   = b   *C_12_R *tmp1
  bp1 = bp1 *C_12_R *tmp1

! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
  tmp1 = C_1_R/b

  lhs(n,1) = am1 *tmp1
  lhs(n,2) = a   *tmp1 
  lhs(n,3) = ap1 *tmp1 
  
  rhs(n,1) = C_0_R
  rhs(n,2) = bm1 *tmp1
  rhs(n,3) = bp1 *tmp1 
  rhs(n,4) = C_0_R

! #######################################################################
! 2 < n < nmax-1
! #######################################################################
  DO n = 3,nmax-2

  tmpa = (x(n-1)-x(n+1)) * &
       ( C_1_R/(x(n-1)-x(n-2)) + C_1_R/(x(n-1)-x(n+2)) + C_1_R/(x(n-1)-x(n)) )
  tmpb = (x(n+1)-x(n-1)) * &
       ( C_1_R/(x(n+1)-x(n-2)) + C_1_R/(x(n+1)-x(n+2)) + C_1_R/(x(n+1)-x(n)) )
  tmpd = C_1_R / &                                    ! constant 2/D in Table 2 in Shukla&Zhong(2005)
       ( (tmpb+C_2_R) *(tmpa+C_2_R) - C_1_R )
  tmpa = tmpa + C_1_R                                 ! constant a in ...
  tmpb = tmpb + C_1_R                                 ! constant b in ...

  tmpp = (x(n+1)-x(n-1)) *(x(n+1)-x(n-1)) * &         ! constant p
       ( C_1_R/( (x(n-1)-x(n)  ) *(x(n-1)-x(n+2)) ) &
       + C_1_R/( (x(n-1)-x(n-2)) *(x(n-1)-x(n+2)) ) &
       + C_1_R/( (x(n-1)-x(n-2)) *(x(n-1)-x(n)  ) ) )
  tmpq = (x(n+1)-x(n-1)) *(x(n+1)-x(n-1)) * &         ! constant q
       ( C_1_R/( (x(n+1)-x(n)  ) *(x(n+1)-x(n+2)) ) &
       + C_1_R/( (x(n+1)-x(n-2)) *(x(n+1)-x(n+2)) ) &
       + C_1_R/( (x(n+1)-x(n-2)) *(x(n+1)-x(n)  ) ) )

! -------------------------------------------------------------------
  a   = C_1_R

! -------------------------------------------------------------------
  dummy1 = C_1_R + tmpb
  dummy2 =       - tmpb

  tmp1 =   (x(n)-x(n-2)) * (x(n)-x(n+2))   / ( (x(n-1)-x(n-2))*(x(n-1)-x(n+2)) )
  tmp2 = dummy1 *( C_1_R + (x(n)-x(n+1))/(x(n)-x(n-1)) ) &
       + dummy2 *( C_2_R*(x(n)-x(n+1)) + (x(n)-x(n-1)) ) / (x(n+1)-x(n-1))
  am1  = tmp1 *tmp2

  tmp1 = ( (x(n)-x(n-2)) + (x(n)-x(n+2)) ) / ( (x(n-1)-x(n-2))*(x(n-1)-x(n+2)) ) *(x(n)-x(n+1))
  tmp2 = dummy1 &
       + dummy2 *(x(n)-x(n-1))/(x(n+1)-x(n-1)) 
  tmp2 = tmp1 *tmp2

  am1  = ( am1 + tmp2 ) *tmpd

! -------------------------------------------------------------------
  dummy1 = C_4_R + C_2_R*(C_1_R+      tmpb)*(tmpa      -C_2_R) + C_2_R*tmpp*(C_1_R+tmpb) 
  dummy2 = C_1_R +       (C_1_R-C_2_R*tmpa)*(C_2_R*tmpb-C_1_R) - C_2_R*tmpp*tmpb

  tmp1 =   (x(n)-x(n-2)) * (x(n)-x(n+2))   / ( (x(n-1)-x(n-2))*(x(n-1)-x(n+2)) )
  tmp2 = dummy1 *( C_1_R + (x(n)-x(n+1))/(x(n)-x(n-1)) ) &
       + dummy2 *( C_2_R*(x(n)-x(n+1)) + (x(n)-x(n-1)) ) /(x(n+1)-x(n-1)) &
       + C_2_R/tmpd *(x(n+1)-x(n-1)) /(x(n)-x(n-1))
  bm1  = tmp1 *tmp2

  tmp1 = ( (x(n)-x(n-2)) + (x(n)-x(n+2)) ) / ( (x(n-1)-x(n-2))*(x(n-1)-x(n+2)) ) *(x(n)-x(n+1))
  tmp2 = dummy1 &
       + dummy2 *(x(n)-x(n-1))/(x(n+1)-x(n-1)) &
       + C_2_R/tmpd *(x(n+1)-x(n-1)) /(x(n)-x(n-1)) 
  tmp2 = tmp1 *tmp2

  bm1  = ( bm1 + tmp2 ) *tmpd /( (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) )

! -------------------------------------------------------------------
  dummy1 = C_1_R + tmpa
  dummy2 =       - tmpa

  tmp1 =   (x(n)-x(n+2)) * (x(n)-x(n-2))   / ( (x(n+1)-x(n+2))*(x(n+1)-x(n-2)) )
  tmp2 = dummy1 *( C_1_R + (x(n)-x(n-1))/(x(n)-x(n+1)) ) &
       + dummy2 *( C_2_R*(x(n)-x(n-1)) + (x(n)-x(n+1)) ) / (x(n-1)-x(n+1))
  ap1  = tmp1 *tmp2

  tmp1 = ( (x(n)-x(n-2)) + (x(n)-x(n+2)) ) / ( (x(n+1)-x(n-2))*(x(n+1)-x(n+2)) ) *(x(n)-x(n-1))
  tmp2 = dummy1 &
       + dummy2 *(x(n)-x(n+1))/(x(n-1)-x(n+1)) 
  tmp2 = tmp1 *tmp2

  ap1  = ( ap1 + tmp2 ) *tmpd

! -------------------------------------------------------------------
  dummy1 = C_4_R + C_2_R*(C_1_R+      tmpa)*(tmpb      -C_2_R) + C_2_R*tmpq*(C_1_R+tmpa) 
  dummy2 = C_1_R +       (C_1_R-C_2_R*tmpa)*(C_2_R*tmpb-C_1_R) - C_2_R*tmpq*tmpa

  tmp1 =   (x(n)-x(n-2)) * (x(n)-x(n+2))   / ( (x(n+1)-x(n-2))*(x(n+1)-x(n+2)) )
  tmp2 = dummy1 *( C_1_R + (x(n)-x(n-1))/(x(n)-x(n+1)) ) &
       + dummy2 *( C_2_R*(x(n)-x(n-1)) + (x(n)-x(n+1)) ) / (x(n-1)-x(n+1)) &
       + C_2_R/tmpd *(x(n-1)-x(n+1)) /(x(n)-x(n+1))
  bp1  = tmp1 *tmp2

  tmp1 = ( (x(n)-x(n-2)) + (x(n)-x(n+2)) ) / ( (x(n+1)-x(n-2))*(x(n+1)-x(n+2)) ) *(x(n)-x(n-1))
  tmp2 = dummy1 &
       + dummy2 *(x(n)-x(n+1))/(x(n-1)-x(n+1)) &
       + C_2_R/tmpd *(x(n-1)-x(n+1)) /(x(n)-x(n+1)) 
  tmp2 = tmp1 *tmp2

  bp1  = ( bp1 + tmp2 ) *tmpd /( (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) )

! -------------------------------------------------------------------
  j = n

  dummy1 = (x(n+1)-x(j)) / (x(n+1)-x(n-1))
  dummy2 = (x(n-1)-x(j)) / (x(n+1)-x(n-1))

  tmp1 = ( dummy1/dummy2 - dummy2/dummy1 )*(C_1_R-tmpa)*(tmpb-C_1_R)           &
       + (C_1_R-tmpa) *( C_2_R/dummy1 + C_2_R/dummy2 - C_1_R/(dummy1*dummy1) ) &
       - (tmpb-C_1_R) *( C_2_R/dummy1 + C_2_R/dummy2 + C_1_R/(dummy2*dummy2) ) &
       - ( C_1_R/dummy1 + C_1_R/dummy2 )*( C_5_R + C_2_R/(dummy1*dummy2) )
  tmp2 = C_1_R/(dummy1*dummy1) + C_1_R/(dummy2*dummy2) + C_1_R/(dummy1*dummy2) &
       - ( (tmpb-C_1_R)/dummy2 * (C_1_R-tmpa)/dummy1 )                         &
       + ( (tmpb-C_1_R)/dummy2 - (C_1_R-tmpa)/dummy1 ) *( C_1_R/dummy1 + C_1_R/dummy2 )

  b    = (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) / ( (x(n)-x(n-2))*(x(n)-x(n+2)) )     &
       -         (x(n+1)-x(n-1)) *( C_1_R/(x(n)-x(n-2)) + C_1_R/(x(n)-x(n+2)) ) * ( C_1_R/dummy1 + C_1_R/dummy2 ) &
       + C_1_R/(dummy1*dummy2)
  b    = b /tmpd &
       + tmp1 *( (x(n+1)-x(n-1)) *( C_1_R/(x(n)-x(n-2)) + C_1_R/(x(n)-x(n+2)) ) - ( C_1_R/dummy1 + C_1_R/dummy2 ) )
  b    = ( b + tmp2 ) *C_2_R *tmpd /( (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) )

! -------------------------------------------------------------------
  j = n - 2

  dummy1 = (x(n+1)-x(j)) / (x(n+1)-x(n-1))
  dummy2 = (x(n-1)-x(j)) / (x(n+1)-x(n-1))

  tmp1 = ( dummy1/dummy2 - dummy2/dummy1 )*(C_1_R-tmpa)*(tmpb-C_1_R)           &
       + (C_1_R-tmpa) *( C_2_R/dummy1 + C_2_R/dummy2 - C_1_R/(dummy1*dummy1) ) &
       - (tmpb-C_1_R) *( C_2_R/dummy1 + C_2_R/dummy2 + C_1_R/(dummy2*dummy2) ) &
       - ( C_1_R/dummy1 + C_1_R/dummy2 )*( C_5_R + C_2_R/(dummy1*dummy2) )
  tmp2 = C_1_R/(dummy1*dummy1) + C_1_R/(dummy2*dummy2) + C_1_R/(dummy1*dummy2) &
       - ( (tmpb-C_1_R)/dummy2 * (C_1_R-tmpa)/dummy1 )                         &
       + ( (tmpb-C_1_R)/dummy2 - (C_1_R-tmpa)/dummy1 ) *( C_1_R/dummy1 + C_1_R/dummy2 )

  dummy1 = C_1_R/(x(n)-x(n+1)) + C_1_R/(x(n)-x(n-1))
  dummy1 = tmp1       *( C_1_R + (x(n)-x(j))     *dummy1 )                               &
         + tmp2       *( C_2_R + (x(n)-x(j))     *dummy1 ) *(x(n)-x(j)) /(x(n+1)-x(n-1)) &
         + C_1_R/tmpd *(         (x(n+1)-x(n-1)) *dummy1 )
  dummy1 = dummy1 *(x(n+1)-x(n-1)) /(x(n+2)-x(n-2)) *(x(n  )-x(n+2)) /(x(n)-x(n-2)) ! for n-2

  dummy2 = (x(n)-x(j)) /(x(n+1)-x(n-1))
  dummy2 = tmp1 *dummy2         &
         + tmp2 *dummy2 *dummy2 &
         + C_1_R/tmpd 
  dummy2 = dummy2 *(x(n+1)-x(n-1)) /(x(n+2)-x(n-2)) *(x(n+1)-x(n-1)) /(x(n)-x(n-2)) ! for n-2

  bm2 = ( dummy1 + dummy2 ) *tmpd /( (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) ) &
       *C_2_R *(x(n)-x(n+1)) /(x(j)-x(n+1)) *(x(n)-x(n-1)) /(x(j)-x(n-1))

! -------------------------------------------------------------------
  j = n + 2

  dummy1 = (x(n+1)-x(j)) / (x(n+1)-x(n-1))
  dummy2 = (x(n-1)-x(j)) / (x(n+1)-x(n-1))

  tmp1 = ( dummy1/dummy2 - dummy2/dummy1 )*(C_1_R-tmpa)*(tmpb-C_1_R)           &
       + (C_1_R-tmpa) *( C_2_R/dummy1 + C_2_R/dummy2 - C_1_R/(dummy1*dummy1) ) &
       - (tmpb-C_1_R) *( C_2_R/dummy1 + C_2_R/dummy2 + C_1_R/(dummy2*dummy2) ) &
       - ( C_1_R/dummy1 + C_1_R/dummy2 )*( C_5_R + C_2_R/(dummy1*dummy2) )
  tmp2 = C_1_R/(dummy1*dummy1) + C_1_R/(dummy2*dummy2) + C_1_R/(dummy1*dummy2) &
       - ( (tmpb-C_1_R)/dummy2 * (C_1_R-tmpa)/dummy1 )                         &
       + ( (tmpb-C_1_R)/dummy2 - (C_1_R-tmpa)/dummy1 ) *( C_1_R/dummy1 + C_1_R/dummy2 )

  dummy1 = C_1_R/(x(n)-x(n+1)) + C_1_R/(x(n)-x(n-1))
  dummy1 = tmp1       *( C_1_R + (x(n)-x(j))     *dummy1 )                               &
         + tmp2       *( C_2_R + (x(n)-x(j))     *dummy1 ) *(x(n)-x(j)) /(x(n+1)-x(n-1)) &
         + C_1_R/tmpd *(         (x(n+1)-x(n-1)) *dummy1 )
  dummy1 = dummy1 *(x(n+1)-x(n-1)) /(x(n+2)-x(n-2)) *(x(n  )-x(n-2)) /(x(n+2)-x(n)) ! for n+2

  dummy2 = (x(n)-x(j)) /(x(n+1)-x(n-1))
  dummy2 = tmp1 *dummy2         &
         + tmp2 *dummy2 *dummy2 &
         + C_1_R/tmpd 
  dummy2 = dummy2 *(x(n+1)-x(n-1)) /(x(n+2)-x(n-2)) *(x(n+1)-x(n-1)) /(x(n+2)-x(n)) ! for n+2

  bp2 = ( dummy1 + dummy2 ) *tmpd /( (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) ) &
       *C_2_R *(x(n)-x(n+1)) /(x(j)-x(n+1)) *(x(n)-x(n-1)) /(x(j)-x(n-1))

! -------------------------------------------------------------------
! if uniform, we should have ( 2/11 1 2/11 ) and ( 3/44 12/11 -51/22 12/11 3/44 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
  tmp1 = C_1_R /b

  lhs(n,1) = am1 *tmp1
  lhs(n,2) = a   *tmp1 
  lhs(n,3) = ap1 *tmp1 
  
  rhs(n,1) = bm2 *tmp1
  rhs(n,2) = bm1 *tmp1
  rhs(n,3) = bp1 *tmp1 
  rhs(n,4) = bp2 *tmp1

  ENDDO

! #######################################################################
! n = nmax-1; just a copy of case n = 2
! #######################################################################
  n = nmax-1

  bm2 = C_0_R; bp2 = C_0_R 

  tmp1 = C_1_R / &
       ( (x(n+1)-x(n))*(x(n)-x(n-1)) + (x(n+1)-x(n-1))*(x(n+1)-x(n-1)) )
  
  bm1 = (x(n+1)-x(n)) / (x(n+1)-x(n-1))
  b   =-C_1_R
  bp1 = (x(n)-x(n-1)) / (x(n+1)-x(n-1))

  am1 = (x(n)-x(n-1))*(x(n)-x(n-1)) + (x(n+1)-x(n))*(x(n)-x(n-1)) - (x(n+1)-x(n))*(x(n+1)-x(n))
  am1 = am1 *bm1 *tmp1
  a   = C_1_R
  ap1 = (x(n+1)-x(n))*(x(n+1)-x(n)) + (x(n+1)-x(n))*(x(n)-x(n-1)) - (x(n)-x(n-1))*(x(n)-x(n-1))
  ap1 = ap1 *bp1 *tmp1

  bm1 = bm1 *C_12_R *tmp1
  b   = b   *C_12_R *tmp1
  bp1 = bp1 *C_12_R *tmp1

! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
  tmp1 = C_1_R/b

  lhs(n,1) = am1 *tmp1
  lhs(n,2) = a   *tmp1 
  lhs(n,3) = ap1 *tmp1 
  
  rhs(n,1) = C_0_R
  rhs(n,2) = bm1 *tmp1
  rhs(n,3) = bp1 *tmp1 
  rhs(n,4) = C_0_R

! #######################################################################
! n = nmax; same as n = 1, but changing the signs of the increments w.r.t. n
! To understand it, e.g., define a new variable i = -j, where i is the 
! discrete variable moving around n
! #######################################################################
  n = nmax

! left-hand side
  am1 = C_0_R
  a   = C_1_R
  ap1 = (x(n-1)-x(n))*(x(n)-x(n-2)) + (x(n-1)-x(n))*(x(n)-x(n-3)) - (x(n)-x(n-2))*(x(n)-x(n-3))
  ap1 = ap1 / &
      ( (x(n-1)-x(n))*(x(n-1)-x(n-2)) + (x(n-1)-x(n))*(x(n-1)-x(n-3)) + (x(n-1)-x(n-2))*(x(n-1)-x(n-3)) )

! right-hand side
  bm2 = C_0_R
  bm1 = C_0_R
  b   = (x(n-1)-x(n-2))*(x(n-1)-x(n-3)) + C_2_R*(x(n-1)-x(n))*( (x(n-1)-x(n-2)) + (x(n-1)-x(n-3)) )
  b   = b / &
      ( (x(n-1)-x(n-2))*(x(n-1)-x(n-3)) +       (x(n-1)-x(n))*( (x(n-1)-x(n-2)) + (x(n-1)-x(n-3)) ) )
  b   = b *( (x(n)-x(n-2)) + (x(n)-x(n-3)) ) / &
        (x(n)-x(n-1)) + C_1_R
  b   = b / &
      ( (x(n)-x(n-2))*(x(n)-x(n-3)) ) 
  tmp1= ( (x(n-1)-x(n-2)) + (x(n-1)-x(n-3)) ) / &
      ( (x(n-1)-x(n-2))*(x(n-1)-x(n-3)) +       (x(n-1)-x(n))*( (x(n-1)-x(n-2)) + (x(n-1)-x(n-3)) ) )
  b   = ( b + tmp1 / (x(n-1)-x(n)) ) *C_2_R

  bp1 = (x(n)-x(n-2)) + (x(n)-x(n-3)) + ap1 *( (x(n-1)-x(n)) + (x(n-1)-x(n-2)) + (x(n-1)-x(n-3)) )
  bp1 = bp1 *C_2_R / &
      ( (x(n-1)-x(n))*(x(n-1)-x(n-2))*(x(n-1)-x(n-3)) )
     
  bp2 = (x(n)-x(n-1))*(       (x(n-1)-x(n)) + C_2_R*(x(n-1)-x(n-3)) ) + &
        (x(n)-x(n-3))*( C_2_R*(x(n-1)-x(n)) + C_3_R*(x(n-1)-x(n-3)) )
  bp2 = bp2 / &
      ( (x(n-1)-x(n))*(x(n-1)-x(n-3)) - (x(n-2)-x(n-1))*( (x(n-1)-x(n)) + (x(n-1)-x(n-3)) ) )
  bp2 = bp2 *C_2_R*(x(n-1)-x(n)) / &
      ( (x(n-2)-x(n-1))*(x(n-2)-x(n))*(x(n-2)-x(n-3)) )

  bp3 = (x(n)-x(n-1))*(       (x(n-1)-x(n)) + C_2_R*(x(n-1)-x(n-2)) ) + &
        (x(n)-x(n-2))*( C_2_R*(x(n-1)-x(n)) + C_3_R*(x(n-1)-x(n-2)) )
  bp3 = bp3 / &
      ( (x(n-1)-x(n))*(x(n-1)-x(n-2)) - (x(n-3)-x(n-1))*( (x(n-1)-x(n)) + (x(n-1)-x(n-2)) ) )
  bp3 = bp3 *C_2_R*(x(n-1)-x(n)) / &
      ( (x(n-3)-x(n-1))*(x(n-3)-x(n))*(x(n-3)-x(n-2)) )

! if uniform, we should have ( 11 1 0 ) and ( -1 15 -27 13 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
  tmp1 = C_1_R/b

  lhs(n,3) = C_0_R
  lhs(n,2) = a   *tmp1 
  lhs(n,1) = ap1 *tmp1 
  
  rhs(n,4) = bp3 *tmp1 ! saving here the last term that goes into the 3. subdiagonal
  rhs(n,3) = C_0_R
  rhs(n,2) = bp1 *tmp1 
  rhs(n,1) = bp2 *tmp1 

  RETURN
END SUBROUTINE FDM_C2N6ND_INITIALIZE

! #######################################################################
! Constructing forcing term
! #######################################################################
SUBROUTINE FDM_C2N6ND_RHS(nmax,mmax, rhs, u, d)

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: nmax,mmax ! m linear systems or size n
  TREAL,   DIMENSION(nmax,4),    INTENT(IN) :: rhs       ! RHS diagonals (#=5-1 because of normalization)
  TREAL,   DIMENSION(mmax,nmax), INTENT(IN) :: u         ! function
  TREAL,   DIMENSION(mmax,nmax), INTENT(OUT):: d         ! RHS

! -------------------------------------------------------------------
  TINTEGER n

! #######################################################################
! Boundary
  n = 1 ! rhs(1,1) contains 3. superdiagonal
     d(:,n) =                                         &
            + u(:,n)                                  &
            + u(:,n+1) *rhs(n,3) + u(:,n+2) *rhs(n,4) + u(:,n+3) *rhs(n,1)

  n = 2
     d(:,n) =                      u(:,n-1) *rhs(n,2) &
            + u(:,n)                                  &
            + u(:,n+1) *rhs(n,3)                     

! Interior points
  DO n = 3,nmax-2
     d(:,n) = u(:,n-2) *rhs(n,1) + u(:,n-1) *rhs(n,2) &
            + u(:,n)                                  &
            + u(:,n+1) *rhs(n,3) + u(:,n+2) *rhs(n,4)
  ENDDO
  
! Boundary
  n = nmax-1
     d(:,n) =                      u(:,n-1) *rhs(n,2) &
            + u(:,n)                                  &
            + u(:,n+1) *rhs(n,3)                     

  n = nmax ! rhs(1,4) contains 3. subdiagonal
     d(:,n) = u(:,n-3) *rhs(n,4) + u(:,n-2) *rhs(n,1) + u(:,n-1) *rhs(n,2) &
            + u(:,n)

  RETURN
END SUBROUTINE FDM_C2N6ND_RHS
