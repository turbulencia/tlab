#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2010/10/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Solving (u')' = f  (General case with \lambda=0)
!#
!########################################################################
!# ARGUMENTS 
!#
!# bcs     In    BCs data
!# u       Out   Solution
!# f       In    Forcing
!# tmp1    Out   First derivative
!#
!########################################################################

!########################################################################
!Dirichlet/Neumann boundary conditions at imin/imax
!########################################################################
SUBROUTINE FDE_BVP_SINGULAR_DN(imode_fdm, imax,jkmax, dx, u,f,bcs, tmp1, wrk1d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, imax, jkmax
  TREAL, DIMENSION(imax)           :: dx
  TREAL, DIMENSION(imax,7), TARGET :: wrk1d
  TREAL, DIMENSION(jkmax,imax)     :: u, f, tmp1
  TREAL, DIMENSION(jkmax,2)        :: bcs

! -----------------------------------------------------------------------
  TINTEGER i
  TREAL dummy
  TREAL, DIMENSION(:), POINTER :: a,b,c,d,e

! #######################################################################
  a => wrk1d(:,1)
  b => wrk1d(:,2)
  c => wrk1d(:,3)
  d => wrk1d(:,4)
  e => wrk1d(:,5)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_imax given
! -----------------------------------------------------------------------
  f(:,1) = C_0_R
  IF ( imode_fdm .EQ. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT .OR. imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6_LHS(imax,       i2,     a,b,c,d,e)
     CALL INT_C1N6_RHS(imax,jkmax, i2, dx, f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(1,6) = dx(1); wrk1d(2,6) = a(1)*dx(1) ! for v^1
  ENDIF
  CALL PENTADFS(imax-1,       a,b,c,d,e)
  
! obtain v^0, array tmp1
  CALL PENTADSS(imax-1,jkmax, a,b,c,d,e, tmp1)
  tmp1(:,imax) = C_0_R
  DO i = 1,imax
     tmp1(:,i) = tmp1(:,i) + bcs(:,2) ! add v_N to free array bcs(:,2)
  ENDDO

! obtain v^1, array wrk1d(:,6)
  CALL PENTADSS(imax-1,i1,    a,b,c,d,e, wrk1d(1,6))
  wrk1d(imax,6) = C_0_R

! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT .OR. imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6_LHS(imax,       i1,     a,b,c,d,e)
     CALL INT_C1N6_RHS(imax,jkmax, i1, dx, tmp1,u)
  ENDIF
  CALL PENTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2))
     
!obtain u^0
  CALL PENTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2), u(1,2))
  bcs(:,2) = u(:,1); u(:,1) = C_0_R
  bcs(:,2) =(bcs(:,2)+ c(1)*u(:,1)+ d(1)*u(:,2)+ e(1)*u(:,3))/dx(1) !u^(0)'_1
  
!obtain u^1, array wrk1d(:,7)
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT .OR. imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6_RHS(imax,i1, i1, dx, wrk1d(1,6),wrk1d(1,7))
  ENDIF
  CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), wrk1d(2,7))
  dummy = wrk1d(1,7); wrk1d(1,7) = C_0_R
  dummy =(dummy+ c(1)*wrk1d(1,7)+ d(1)*wrk1d(2,7)+ e(1)*wrk1d(3,7))/dx(1) ! u^(1)'_1

! Constraint
  dummy = C_1_R/(dummy-wrk1d(1,6))
  bcs(:,2) = (tmp1(:,1) - bcs(:,2))*dummy

! Result
  DO i = 1,imax
     u(:,i)    = u(:,i)    + bcs(:,2)*wrk1d(i,7) + bcs(:,1)
     tmp1(:,i) = tmp1(:,i) + bcs(:,2)*wrk1d(i,6)
  ENDDO

  RETURN
END SUBROUTINE FDE_BVP_SINGULAR_DN

!########################################################################
!Neumann/Dirichlet boundary conditions at imin/imax
!########################################################################
SUBROUTINE FDE_BVP_SINGULAR_ND(imode_fdm, imax,jkmax, dx, u,f,bcs, tmp1, wrk1d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, imax, jkmax
  TREAL, DIMENSION(imax)           :: dx
  TREAL, DIMENSION(imax,7), TARGET :: wrk1d
  TREAL, DIMENSION(jkmax,imax)     :: u, f, tmp1
  TREAL, DIMENSION(jkmax,2)        :: bcs

! -----------------------------------------------------------------------
  TINTEGER i
  TREAL dummy
  TREAL, DIMENSION(:), POINTER :: a,b,c,d,e

! #######################################################################
  a => wrk1d(:,1)
  b => wrk1d(:,2)
  c => wrk1d(:,3)
  d => wrk1d(:,4)
  e => wrk1d(:,5)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
  f(:,imax) = C_0_R
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT .OR. imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6_LHS(imax,       i1,     a,b,c,d,e)
     CALL INT_C1N6_RHS(imax,jkmax, i1, dx, f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(imax,6) = dx(imax); wrk1d(imax-1,6) = e(imax)*dx(imax) ! for v^1
  ENDIF
  CALL PENTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2))
  
! obtain v^0, array tmp1
  CALL PENTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2), tmp1(1,2))
  tmp1(:,1) = C_0_R
  DO i = 1,imax
     tmp1(:,i) = tmp1(:,i) + bcs(:,1) ! add v_1 to free array bcs(:,1)
  ENDDO

! obtain v^1, array wrk1d(:,6)
  CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), wrk1d(2,6))
  wrk1d(1,6) = C_0_R

! -----------------------------------------------------------------------
! solve for u in u' = v, u_N given
! -----------------------------------------------------------------------
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT .OR. imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6_LHS(imax,       i2,     a,b,c,d,e)
     CALL INT_C1N6_RHS(imax,jkmax, i2, dx, tmp1,u)
  ENDIF
  CALL PENTADFS(imax-1,       a,b,c,d,e)
     
!obtain u^0
  CALL PENTADSS(imax-1,jkmax, a,b,c,d,e, u)
  bcs(:,1) = u(:,imax); u(:,imax) = C_0_R
  bcs(:,1) =(bcs(:,1)+ a(imax)*u(:,imax-2)+ b(imax)*u(:,imax-1)+ c(imax)*u(:,imax))/dx(imax) !u^(0)'_imax
  
!obtain u^1, array wrk1d(:,7)
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT .OR. imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6_RHS(imax,i1, i2, dx, wrk1d(1,6),wrk1d(1,7))
  ENDIF
  CALL PENTADSS(imax-1,i1,    a,b,c,d,e, wrk1d(1,7))
  dummy = wrk1d(imax,7); wrk1d(imax,7) = C_0_R
  dummy =(dummy+ a(imax)*wrk1d(imax-2,7)+ b(imax)*wrk1d(imax-1,7)+ c(imax)*wrk1d(imax,7))/dx(imax) ! u^(1)'_imax

! Constraint
  dummy = C_1_R/(dummy-wrk1d(imax,6))
  bcs(:,1) = (tmp1(:,imax) - bcs(:,1))*dummy

! Result
  DO i = 1,imax
     u(:,i)    = u(:,i)    + bcs(:,1)*wrk1d(i,7) + bcs(:,2)
     tmp1(:,i) = tmp1(:,i) + bcs(:,1)*wrk1d(i,6)
  ENDDO

  RETURN
END SUBROUTINE FDE_BVP_SINGULAR_ND

!########################################################################
!Dirichlet/Dirichlet boundary conditions at imin/imax
!########################################################################
SUBROUTINE FDE_BVP_SINGULAR_DD(imode_fdm, imax,jkmax, x,dx, u,f,bcs, tmp1, wrk1d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, imax, jkmax
  TREAL, DIMENSION(imax)           :: dx, x
  TREAL, DIMENSION(imax,9), TARGET :: wrk1d
  TREAL, DIMENSION(jkmax,imax)     :: u, f, tmp1
  TREAL, DIMENSION(jkmax,3)        :: bcs

! -----------------------------------------------------------------------
  TINTEGER i
  TREAL dummy
  TREAL, DIMENSION(:), POINTER :: a,b,c,d,e,g,h

! #######################################################################
  a => wrk1d(:,1)
  b => wrk1d(:,2)
  c => wrk1d(:,3)
  d => wrk1d(:,4)
  e => wrk1d(:,5)
  ! additional diagonals
  g => wrk1d(:,8)
  h => wrk1d(:,9)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v = u' in (u')' = f , u'_imax given
! -----------------------------------------------------------------------
  f(:,1) = C_0_R
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_LHS(imax,       i2,     a,b,c,d,e)
     CALL INT_C1N6_RHS(imax,jkmax, i2, dx, f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(1,6) = dx(1); wrk1d(2,6) = a(1)*dx(1) ! for v^1 
 
     CALL PENTADFS(imax-1,       a,b,c,d,e)
  
! obtain v^0, array tmp1
     CALL PENTADSS(imax-1,jkmax, a,b,c,d,e, tmp1)
     tmp1(:,imax) = C_0_R

! obtain v^1, array wrk1d(:,6)
     CALL PENTADSS(imax-1,i1,    a,b,c,d,e, wrk1d(1,6))
     wrk1d(imax,6) = C_0_R

  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_LHS(imax,       i2,     a,b,c,d,e,g,h)
     CALL INT_C1N6M_RHS(imax,jkmax, i2, dx, f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(1,6) = dx(1); wrk1d(2,6) = b(1)*dx(1) ! for v^1 

     CALL HEPTADFS(imax-1,       a,b,c,d,e,g,h)
 
! obtain v^0, array tmp1
     CALL HEPTADSS(imax-1,jkmax, a,b,c,d,e,g,h, tmp1)
     tmp1(:,imax) = C_0_R

! obtain v^1, array wrk1d(:,6)
     CALL HEPTADSS(imax-1,i1,    a,b,c,d,e,g,h, wrk1d(1,6))
     wrk1d(imax,6) = C_0_R

  ENDIF
! -----------------------------------------------------------------------
! solve for u in u' v f, u_1 given
! -----------------------------------------------------------------------
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_LHS(imax,       i1,     a,b,c,d,e)
     CALL INT_C1N6_RHS(imax,jkmax, i1, dx, tmp1,u)
  
     CALL PENTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2))
     
!obtain u^0
     CALL PENTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2), u(1,2))
     bcs(:,3) = u(:,1); u(:,1) = C_0_R
     bcs(:,3) =(bcs(:,3)+ c(1)*u(:,1)+ d(1)*u(:,2)+ e(1)*u(:,3))/dx(1) !u^(0)'_1
  
  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_LHS(imax,       i1,     a,b,c,d,e,g,h)
     CALL INT_C1N6M_RHS(imax,jkmax, i1, dx, tmp1,u)
  
     CALL HEPTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2),g(2),h(2))
    
!obtain u^0
     CALL HEPTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2),g(2),h(2), u(1,2))
     bcs(:,3) = u(:,1); u(:,1) = C_0_R
     bcs(:,3) =(bcs(:,3)+ d(1)*u(:,1)+ e(1)*u(:,2)+ g(1)*u(:,3))/dx(1) !u^(0)'_1
  ENDIF
  
!obtain u^1, array wrk1d(:,7)
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_RHS(imax,i1, i1, dx, wrk1d(1,6),wrk1d(1,7))
  
     CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), wrk1d(2,7))
     dummy = wrk1d(1,7); wrk1d(1,7) = C_0_R
     dummy =(dummy+ c(1)*wrk1d(1,7)+ d(1)*wrk1d(2,7)+ e(1)*wrk1d(3,7))/dx(1) ! u^(1)'_1

  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_RHS(imax,i1, i1, dx, wrk1d(1,6),wrk1d(1,7))
   
     CALL HEPTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2),g(2),h(2), wrk1d(2,7))
     dummy = wrk1d(1,7); wrk1d(1,7) = C_0_R
     dummy =(dummy+ d(1)*wrk1d(1,7)+ e(1)*wrk1d(2,7)+ g(1)*wrk1d(3,7))/dx(1) ! u^(1)'_1

  ENDIF
! Constraint
  dummy = C_1_R/(dummy-wrk1d(1,6))
  bcs(:,3) = (tmp1(:,1) - bcs(:,3))*dummy

! BCs
  dummy = C_1_R/(x(imax)-x(1))
  bcs(:,2) =(bcs(:,2) - u(:,imax) - bcs(:,3)*wrk1d(imax,7) - bcs(:,1))*dummy

! Result
  DO i = 1,imax
     u(:,i)    = u(:,i)    + bcs(:,3)*wrk1d(i,7) + bcs(:,2)*(x(i)-x(1)) + bcs(:,1)
     tmp1(:,i) = tmp1(:,i) + bcs(:,3)*wrk1d(i,6) + bcs(:,2)
  ENDDO

  RETURN
END SUBROUTINE FDE_BVP_SINGULAR_DD

!########################################################################
!Neumann/Neumann boundary conditions at imin/imax; must be compatible!
!########################################################################
SUBROUTINE FDE_BVP_SINGULAR_NN(imode_fdm, imax,jkmax, dx, u,f,bcs, tmp1, wrk1d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, imax, jkmax
  TREAL, DIMENSION(imax)           :: dx
  TREAL, DIMENSION(imax,7), TARGET :: wrk1d
  TREAL, DIMENSION(jkmax,imax)     :: u, f, tmp1
  TREAL, DIMENSION(jkmax,2)        :: bcs

! -----------------------------------------------------------------------
  TINTEGER i
  TREAL, DIMENSION(:), POINTER :: a,b,c,d,e,g,h

! #######################################################################
  a => wrk1d(:,1)
  b => wrk1d(:,2)
  c => wrk1d(:,3)
  d => wrk1d(:,4)
  e => wrk1d(:,5)
  ! additional diagonals
  g => wrk1d(:,6)
  h => wrk1d(:,7)

! #######################################################################
! -----------------------------------------------------------------------
! solve for v in v' = f , v_1 given
! -----------------------------------------------------------------------
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_LHS(imax,       i1,     a,b,c,d,e)
     CALL INT_C1N6_RHS(imax,jkmax, i1, dx, f,tmp1)
     !
     CALL PENTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2))
     CALL PENTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2), tmp1(1,2))
  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_LHS(imax,       i1,     a,b,c,d,e,g,h)
     CALL INT_C1N6M_RHS(imax,jkmax, i1, dx, f,tmp1)
     !
     CALL HEPTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2),g(2),h(2))
     CALL HEPTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2),g(2),h(2), tmp1(1,2))
  ENDIF

  tmp1(:,1) = C_0_R
  DO i = 1,imax
     tmp1(:,i) = tmp1(:,i) + bcs(:,1) ! this step assumes compatible problem
  ENDDO

! -----------------------------------------------------------------------
! solve for u in u' = v, u_1 given
! -----------------------------------------------------------------------
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
!    same l.h.s. as before
     CALL INT_C1N6_RHS(imax,jkmax, i1, dx, tmp1,u)
! same l.h.s. as before
     CALL PENTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2), u(1,2))
  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_RHS(imax,jkmax, i1, dx, tmp1,u)
     CALL HEPTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2),g(2),h(2), u(1,2))
  ENDIF

  u(:,1) = C_0_R ! this integration constant is free and set to zero

  RETURN
END SUBROUTINE FDE_BVP_SINGULAR_NN

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2010/10/11 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Solving (u')' - \lambda^2 u = f
!#
!########################################################################
!# ARGUMENTS 
!#
!# bcs     In    BCs data
!# u       Out   Solution
!# f       In    Forcing
!# tmp1    Out   First derivative
!#
!########################################################################

!########################################################################
!Neumann/Neumann boundary conditions at imin/imax
!########################################################################
SUBROUTINE FDE_BVP_REGULAR_NN(imode_fdm, imax,jkmax, cst, dx, u,f,bcs, tmp1, wrk1d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, imax, jkmax
  TREAL cst
  TREAL, DIMENSION(imax)            :: dx
  TREAL, DIMENSION(imax,12), TARGET :: wrk1d
  TREAL, DIMENSION(jkmax,imax)      :: u, f, tmp1
  TREAL, DIMENSION(jkmax,3)         :: bcs

! -----------------------------------------------------------------------
  TINTEGER i
  TREAL lambda, dummy, g_1, g_2, a_1, b_1, deti
  TREAL, DIMENSION(:), POINTER :: a,b,c,d,e,g,h, ep,em

! #######################################################################
  a => wrk1d(:,1)
  b => wrk1d(:,2)
  c => wrk1d(:,3)
  d => wrk1d(:,4)
  e => wrk1d(:,5)
  ! additional diagonals 
  g => wrk1d(:,11)
  h => wrk1d(:,12)

  ep=> wrk1d(:,9)
  em=> wrk1d(:,10)

  lambda = SQRT(cst)

! #######################################################################
! -----------------------------------------------------------------------
! 1st step; solve for v^(0) and v^(1)
! -----------------------------------------------------------------------
  dummy =-lambda
  f(:,1) = C_0_R
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_LHS_E(imax,       i2, dx, dummy, a,b,c,d,e, ep)
     CALL INT_C1N6_RHS  (imax,jkmax, i2, dx,        f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(1,6) = dx(1); wrk1d(2,6) = a(1)*dx(1) ! for v^1

     CALL PENTADFS(imax-1,       a,b,c,d,e)
     
! obtain e^(+), array ep
     CALL PENTADSS(imax-1,i1,    a,b,c,d,e, ep)
   
! obtain v^(0), array tmp1
     CALL PENTADSS(imax-1,jkmax, a,b,c,d,e, tmp1)
     tmp1(:,imax) = C_0_R
   
! obtain v^(1), array wrk1d(:,6)
     CALL PENTADSS(imax-1,i1,    a,b,c,d,e, wrk1d(1,6))
     wrk1d(imax,6) = C_0_R

  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_LHS_E(imax,       i2, dx, dummy, a,b,c,d,e,g,h, ep)
     CALL INT_C1N6M_RHS  (imax,jkmax, i2, dx,        f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(1,6) = dx(1); wrk1d(2,6) = b(1)*dx(1) ! for v^1

     CALL HEPTADFS(imax-1,       a,b,c,d,e,g,h)
    
! obtain e^(+), array ep
     CALL HEPTADSS(imax-1,i1,    a,b,c,d,e,g,h, ep)
  
! obtain v^(0), array tmp1
     CALL HEPTADSS(imax-1,jkmax, a,b,c,d,e,g,h, tmp1)
     tmp1(:,imax) = C_0_R
  
! obtain v^(1), array wrk1d(:,6)
     CALL HEPTADSS(imax-1,i1,    a,b,c,d,e,g,h, wrk1d(1,6))
     wrk1d(imax,6) = C_0_R
  
  ENDIF

! -----------------------------------------------------------------------
! 2nd step; solve for u^(0) and u^(1) and u^(2)
! -----------------------------------------------------------------------
  dummy = lambda
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_LHS_E(imax,       i1, dx, dummy, a,b,c,d,e, em)
     CALL INT_C1N6_RHS  (imax,jkmax, i1, dx,        tmp1,u)
  
     CALL PENTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2))
        
! obtain e^(m), array em
     CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), em(2))
     g_1 =(c(1)*em(1)+ d(1)*em(2)+ e(1)*em(3))/dx(1)/lambda + C_1_R ! e^(-)'_1/\lambda + 1
   
! obtain u^(2), array wrk1d(:,8)
     CALL INT_C1N6_RHS(imax,i1, i1, dx, ep,wrk1d(1,8))
     CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), wrk1d(2,8))
     g_2 = wrk1d(1,8); wrk1d(1,8) = C_0_R
     g_2 =(g_2+ c(1)*wrk1d(1,8)+ d(1)*wrk1d(2,8)+ e(1)*wrk1d(3,8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1 
   
! obtain u^(0), array u
     CALL PENTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2), u(1,2))
   
! BCs; intermediate step to save memory space
     dummy = C_1_R - wrk1d(imax,8)*lambda
     bcs(:,3) = tmp1(:,1)        - bcs(:,1)
     bcs(:,2) = u(:,imax)*lambda + bcs(:,2)
     bcs(:,1) = bcs(:,2)       + bcs(:,3)*em(imax)! a_0 *det
     bcs(:,2) = bcs(:,2)*ep(1) + bcs(:,3)*dummy   ! b_0 *det
   
     bcs(:,3) = u(:,1); u(:,1) = C_0_R
     bcs(:,3) =(bcs(:,3)+ c(1)*u(:,1)+ d(1)*u(:,2)+ e(1)*u(:,3))/dx(1) !u^(0)'_1
     
!obtain u^(1), array wrk1d(:,7)
     CALL INT_C1N6_RHS(imax,i1, i1, dx, wrk1d(1,6),wrk1d(1,7))
     CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), wrk1d(2,7))
     dummy = wrk1d(1,7); wrk1d(1,7) = C_0_R
     dummy =(dummy+ c(1)*wrk1d(1,7)+ d(1)*wrk1d(2,7)+ e(1)*wrk1d(3,7))/dx(1) ! u^(1)'_1

  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_LHS_E(imax,       i1, dx, dummy, a,b,c,d,e,g,h, em)
     CALL INT_C1N6M_RHS  (imax,jkmax, i1, dx,        tmp1,u)
  
     CALL HEPTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2),g(2),h(2))
        
! obtain e^(m), array em
     CALL HEPTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2),g(2),h(2), em(2))
     g_1 =(d(1)*em(1)+ e(1)*em(2)+ g(1)*em(3))/dx(1)/lambda + C_1_R ! e^(-)'_1/\lambda + 1
  
! obtain u^(2), array wrk1d(:,8)
     CALL INT_C1N6M_RHS(imax,i1, i1, dx, ep,wrk1d(1,8))
     CALL HEPTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2),g(2),h(2), wrk1d(2,8))
     g_2 = wrk1d(1,8); wrk1d(1,8) = C_0_R
     g_2 =(g_2+ c(1)*wrk1d(1,8)+ e(1)*wrk1d(2,8)+ g(1)*wrk1d(3,8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1 
  
! obtain u^(0), array u
     CALL HEPTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2),g(2),h(2), u(1,2))
   
! BCs; intermediate step to save memory space
     dummy = C_1_R - wrk1d(imax,8)*lambda
     bcs(:,3) = tmp1(:,1)        - bcs(:,1)
     bcs(:,2) = u(:,imax)*lambda + bcs(:,2)
     bcs(:,1) = bcs(:,2)       + bcs(:,3)*em(imax)! a_0 *det
     bcs(:,2) = bcs(:,2)*ep(1) + bcs(:,3)*dummy   ! b_0 *det
   
     bcs(:,3) = u(:,1); u(:,1) = C_0_R
     bcs(:,3) =(bcs(:,3)+ d(1)*u(:,1)+ e(1)*u(:,2)+ g(1)*u(:,3))/dx(1) !u^(0)'_1
    
! obtain u^(1), array wrk1d(:,7)
     CALL INT_C1N6M_RHS(imax,i1, i1, dx, wrk1d(1,6),wrk1d(1,7))
     CALL HEPTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2),g(2),h(2), wrk1d(2,7))
     dummy = wrk1d(1,7); wrk1d(1,7) = C_0_R
     dummy =(dummy+ d(1)*wrk1d(1,7)+ e(1)*wrk1d(2,7)+ g(1)*wrk1d(3,7))/dx(1) ! u^(1)'_1

  ENDIF

! BCs; final step
  deti = C_1_R/(C_1_R-ep(1)*em(imax)-wrk1d(imax,8)*lambda) ! inverse of determinant
  bcs(:,1) = bcs(:,1) *deti ! a_0
  bcs(:,2) = bcs(:,2) *deti ! b_0
  a_1 =(lambda*wrk1d(imax,7)       + wrk1d(1,6)*em(imax)                    )*deti
  b_1 =(lambda*wrk1d(imax,7)*ep(1) + wrk1d(1,6)*(C_1_R-wrk1d(imax,8)*lambda))*deti

! Constraint
  dummy = C_1_R/( dummy + b_1*g_1 - wrk1d(1,6) + a_1*g_2 )
  bcs(:,3) = (tmp1(:,1) - bcs(:,3) - bcs(:,2)*g_1 - bcs(:,1)*g_2)*dummy

  dummy = C_1_R/lambda
  bcs(:,1) = bcs(:,1) + bcs(:,3)*a_1
  bcs(:,2) =(bcs(:,2) + bcs(:,3)*b_1)*dummy

! Result
  DO i = 1,imax
     u(:,i)    = u(:,i)    + bcs(:,3)*wrk1d(i,7) + bcs(:,1)*wrk1d(i,8) + bcs(:,2)*em(i)
     tmp1(:,i) = tmp1(:,i) + bcs(:,3)*wrk1d(i,6) + bcs(:,1)*ep(i) - lambda*u(:,i)
  ENDDO

  RETURN
END SUBROUTINE FDE_BVP_REGULAR_NN

!########################################################################
!Dirichlet/Dirichlet boundary conditions at imin/imax
!########################################################################
SUBROUTINE FDE_BVP_REGULAR_DD(imode_fdm, imax,jkmax, cst, dx, u,f,bcs, tmp1, wrk1d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, imax, jkmax
  TREAL cst
  TREAL, DIMENSION(imax)            :: dx
  TREAL, DIMENSION(imax,12), TARGET :: wrk1d
  TREAL, DIMENSION(jkmax,imax)      :: u, f, tmp1
  TREAL, DIMENSION(jkmax,3)         :: bcs

! -----------------------------------------------------------------------
  TINTEGER i
  TREAL lambda, dummy, g_1, g_2, a_1, deti
  TREAL, DIMENSION(:), POINTER :: a,b,c,d,e,g,h, ep,em

! #######################################################################
  a => wrk1d(:,1)
  b => wrk1d(:,2)
  c => wrk1d(:,3)
  d => wrk1d(:,4)
  e => wrk1d(:,5)
  ! additional diagonals 
  g => wrk1d(:,11)
  h => wrk1d(:,12)

  ep=> wrk1d(:,9)
  em=> wrk1d(:,10)

  lambda = SQRT(cst)

! #######################################################################
! -----------------------------------------------------------------------
! 1st step; solve for v^(0) and v^(1)
! -----------------------------------------------------------------------
  dummy =-lambda
  f(:,1) = C_0_R
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_LHS_E(imax,       i2, dx, dummy, a,b,c,d,e, ep)
     CALL INT_C1N6_RHS  (imax,jkmax, i2, dx,        f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(1,6) = dx(1); wrk1d(2,6) = a(1)*dx(1) ! for v^1
  
     CALL PENTADFS(imax-1,       a,b,c,d,e)
     
!    obtain e^(+), array ep
     CALL PENTADSS(imax-1,i1,    a,b,c,d,e, ep)
   
!    obtain v^(0), array tmp1
     CALL PENTADSS(imax-1,jkmax, a,b,c,d,e, tmp1)
     tmp1(:,imax) = C_0_R
   
!    obtain v^(1), array wrk1d(:,6)
     CALL PENTADSS(imax-1,i1,    a,b,c,d,e, wrk1d(1,6))
     wrk1d(imax,6) = C_0_R

  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN

     CALL INT_C1N6M_LHS_E(imax,       i2, dx, dummy, a,b,c,d,e,g,h, ep)
     CALL INT_C1N6M_RHS  (imax,jkmax, i2, dx,        f,tmp1)
     wrk1d(:,6) = C_0_R; wrk1d(1,6) = dx(1); wrk1d(2,6) = b(1)*dx(1) ! for v^1
  
     CALL HEPTADFS(imax-1,       a,b,c,d,e,g,h)
     
!     obtain e^(+), array ep
     CALL HEPTADSS(imax-1,i1,    a,b,c,d,e,g,h, ep)
   
!     obtain v^(0), array tmp1
     CALL HEPTADSS(imax-1,jkmax, a,b,c,d,e,g,h, tmp1)
     tmp1(:,imax) = C_0_R
   
!     obtain v^(1), array wrk1d(:,6)
     CALL HEPTADSS(imax-1,i1,    a,b,c,d,e,g,h, wrk1d(1,6))
     wrk1d(imax,6) = C_0_R

  ENDIF

! -----------------------------------------------------------------------
! 2nd step; solve for u^(0) and u^(1) and u^(2)
! -----------------------------------------------------------------------
  dummy = lambda
  IF ( imode_fdm .eq. FDM_COM6_JACOBIAN .OR. imode_fdm .EQ. FDM_COM6_DIRECT ) THEN
     CALL INT_C1N6_LHS_E(imax,       i1, dx, dummy, a,b,c,d,e, em)
     CALL INT_C1N6_RHS  (imax,jkmax, i1, dx,        tmp1,u)
     
     CALL PENTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2))
        
!    obtain e^(m), array em
     CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), em(2))
     g_1 =(c(1)*em(1)+ d(1)*em(2)+ e(1)*em(3))/dx(1)/lambda + C_1_R ! e^(-)'_1/\lambda + 1
   
!    obtain u^(2), array wrk1d(:,8)
     CALL INT_C1N6_RHS(imax,i1, i1, dx, ep,wrk1d(1,8))
     CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), wrk1d(2,8))
     g_2 = wrk1d(1,8); wrk1d(1,8) = C_0_R
     g_2 =(g_2+ c(1)*wrk1d(1,8)+ d(1)*wrk1d(2,8)+ e(1)*wrk1d(3,8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1 
   
!    obtain u^(0), array u
     CALL PENTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2), u(1,2))
     bcs(:,3) = u(:,1); u(:,1) = C_0_R
     bcs(:,3) =(bcs(:,3)+ c(1)*u(:,1)+ d(1)*u(:,2)+ e(1)*u(:,3))/dx(1) !u^(0)'_1
   
!    obtain u^(1), array wrk1d(:,7)
     CALL INT_C1N6_RHS(imax,i1, i1, dx, wrk1d(1,6),wrk1d(1,7))
     CALL PENTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2), wrk1d(2,7))
     dummy = wrk1d(1,7); wrk1d(1,7) = C_0_R
     dummy =(dummy+ c(1)*wrk1d(1,7)+ d(1)*wrk1d(2,7)+ e(1)*wrk1d(3,7))/dx(1) ! u^(1)'_1

  ELSEIF ( imode_fdm .EQ. FDM_COM6_JACPENTA ) THEN
     CALL INT_C1N6M_LHS_E(imax,       i1, dx, dummy, a,b,c,d,e,g,h, em)
     CALL INT_C1N6M_RHS  (imax,jkmax, i1, dx,        tmp1,u)
     
     CALL HEPTADFS(imax-1,       a(2),b(2),c(2),d(2),e(2),g(2),h(2))
        
!     obtain e^(m), array em
     CALL HEPTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2),g(2),h(2), em(2))
     g_1 =(d(1)*em(1)+ e(1)*em(2)+ g(1)*em(3))/dx(1)/lambda + C_1_R ! e^(-)'_1/\lambda + 1
   
!     obtain u^(2), array wrk1d(:,8)
     CALL INT_C1N6M_RHS(imax,i1, i1, dx, ep,wrk1d(1,8))
     CALL HEPTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2),g(2),h(2), wrk1d(2,8))
     g_2 = wrk1d(1,8); wrk1d(1,8) = C_0_R
     g_2 =(g_2+ d(1)*wrk1d(1,8)+ e(1)*wrk1d(2,8)+ g(1)*wrk1d(3,8))/dx(1) - ep(1)! u^(2)'_1 - e^(+)|_1 
   
!     obtain u^(0), array u
     CALL HEPTADSS(imax-1,jkmax, a(2),b(2),c(2),d(2),e(2),g(2),h(2), u(1,2))
     bcs(:,3) = u(:,1); u(:,1) = C_0_R
     bcs(:,3) =(bcs(:,3)+ d(1)*u(:,1)+ e(1)*u(:,2)+ g(1)*u(:,3))/dx(1) !u^(0)'_1
   
!     obtain u^(1), array wrk1d(:,7)
     CALL INT_C1N6M_RHS(imax,i1, i1, dx, wrk1d(1,6),wrk1d(1,7))
     CALL HEPTADSS(imax-1,i1,    a(2),b(2),c(2),d(2),e(2),g(2),h(2), wrk1d(2,7))
     dummy = wrk1d(1,7); wrk1d(1,7) = C_0_R
     dummy =(dummy+ d(1)*wrk1d(1,7)+ e(1)*wrk1d(2,7)+ g(1)*wrk1d(3,7))/dx(1) ! u^(1)'_1

  ENDIF


! BCs; final step
  deti = C_1_R/wrk1d(imax,8) ! inverse of determinant
  bcs(:,2) =(bcs(:,2) - u(:,imax) - bcs(:,1)*em(imax))*deti ! a_0
! bcs(:,1) = bcs(:,1)                                       ! b_0/lambda
  a_1 =-wrk1d(imax,7)*deti
! b_1 = C_0_R

! Constraint
  g_1 = g_1*lambda
  dummy = C_1_R/( dummy - wrk1d(1,6) + a_1*g_2 )
  bcs(:,3) = (tmp1(:,1) - bcs(:,3) - bcs(:,1)*g_1 - bcs(:,2)*g_2)*dummy

  bcs(:,2) = bcs(:,2) + bcs(:,3)*a_1 ! a = v_N
! bcs(:,1) = bcs(:,1)                ! b = u_1

! Result
  DO i = 1,imax
     u(:,i)    = u(:,i)    + bcs(:,3)*wrk1d(i,7) + bcs(:,2)*wrk1d(i,8) + bcs(:,1)*em(i)
     tmp1(:,i) = tmp1(:,i) + bcs(:,3)*wrk1d(i,6) + bcs(:,2)*ep(i) - lambda*u(:,i)
  ENDDO

  RETURN
END SUBROUTINE FDE_BVP_REGULAR_DD