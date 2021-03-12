PROGRAM VNONAD

  IMPLICIT NONE
  
#include "types.h"
#include "integers.h"

  TINTEGER,               parameter :: nmax = 32, len = 5
  TREAL, dimension(nmax,5)          :: A, B 
  TREAL, dimension(nmax,9)          :: C
  TREAL, dimension(len,nmax)        :: x, frc9, frc5

  TINTEGER n, ij, seed, nmax_loc
  TREAL RAN0, error, sol, diff

! ###################################################################
#define A_a(n) A(n,1)
#define A_b(n) A(n,2)
#define A_c(n) A(n,3)
#define A_d(n) A(n,4)
#define A_e(n) A(n,5)

#define B_a(n) B(n,1)
#define B_b(n) B(n,2)
#define B_c(n) B(n,3)
#define B_d(n) B(n,4)
#define B_e(n) B(n,5)

#define C_a(n) C(n,1)
#define C_b(n) C(n,2)
#define C_c(n) C(n,3)
#define C_d(n) C(n,4)
#define C_e(n) C(n,5)
#define C_f(n) C(n,6)
#define C_g(n) C(n,7)
#define C_h(n) C(n,8)
#define C_i(n) C(n,9)

  seed = 256

! create diagonals
  DO n = 1,nmax
     A_a(n) = RAN0(seed)
     A_b(n) = RAN0(seed)
     A_c(n) = RAN0(seed)
     A_d(n) = RAN0(seed)
     A_e(n) = RAN0(seed)
  ENDDO
  DO n = 1,nmax
     B_a(n) = RAN0(seed)
     B_b(n) = RAN0(seed)
     B_c(n) = RAN0(seed)
     B_d(n) = RAN0(seed)
     B_e(n) = RAN0(seed)
  ENDDO

! padding
  A_a(1) = C_0_R; A_b(1) = C_0_R
                  A_a(2) = C_0_R
  A_e(nmax-1) = C_0_R
  A_d(nmax  ) = C_0_R; A_e(nmax) = C_0_R

  B_a(1) = C_0_R; B_b(1) = C_0_R
                  B_a(2) = C_0_R
  B_e(nmax-1) = C_0_R
  B_d(nmax  ) = C_0_R; B_e(nmax) = C_0_R

! check [(N-2)xN] X [Nx(N-2)]
  A_c(1) = C_0_R;    A_d(1) = C_0_R;    A_e(1) = C_0_R
  A_a(nmax) = C_0_R; A_b(nmax) = C_0_R; A_c(nmax) = C_0_R

  B_c(1) = C_0_R; B_e(nmax-2) = C_0_R
  B_b(2) = C_0_R; B_d(nmax-1) = C_0_R
  B_a(3) = C_0_R; B_c(nmax  ) = C_0_R
!

  CALL PENTADMT(nmax, A, B, C)

! create solution
  DO n = 1,nmax
     DO ij = 1,len
        x(ij,n) = RAN0(seed)
     ENDDO
  ENDDO
! check [(N-2)xN] X [Nx(N-2)]
  n = 1
  DO ij = 1,len
     x(ij,n) = C_0_R
  ENDDO
  n = nmax
  DO ij = 1,len
     x(ij,n) = C_0_R
  ENDDO
!

! compute forcing term
  DO n = 1,nmax
     DO ij = 1,len
        frc9(ij,n) = x(ij,n-4)*C_a(n) + x(ij,n-3)*C_b(n) + x(ij,n-2)*C_c(n) + x(ij,n-1)*C_d(n) + x(ij,n)*C_e(n) &
                   + x(ij,n+1)*C_f(n) + x(ij,n+2)*C_g(n) + x(ij,n+3)*C_h(n) + x(ij,n+4)*C_i(n)
        frc5(ij,n) = frc9(ij,n)
     ENDDO
     write(*,'(9f10.3)') (C(n,ij),ij=1,9)
  ENDDO

! ###################################################################
! solve nonadiagonal system
  nmax_loc = nmax-2
  CALL NONADFS(nmax_loc, i1,  C_a(2),C_b(2),C_c(2),C_d(2),C_e(2),C_f(2),C_g(2),C_h(2),C_i(2))
  CALL NONADSS(nmax_loc, len, C_a(2),C_b(2),C_c(2),C_d(2),C_e(2),C_f(2),C_g(2),C_h(2),C_i(2), frc9(1,2))

! solve two pentadiagonal systems
!  CALL PENTADFS(nmax,      A_a(1), A_b(1), A_c(1), A_d(1), A_e(1))
!  CALL PENTADSS(nmax, len, A_a(1), A_b(1), A_c(1), A_d(1), A_e(1), frc5)
!  CALL PENTADFS(nmax,      B_a(1), B_b(1), B_c(1), B_d(1), B_e(1))
!  CALL PENTADSS(nmax, len, B_a(1), B_b(1), B_c(1), B_d(1), B_e(1), frc5)

! error
  diff  = C_0_R
  error = C_0_R
  sol   = C_0_R
  DO n = 1,nmax
     DO ij = 1,len
        diff  = diff  + (frc9(ij,n)-frc5(ij,n))*(frc9(ij,n)-frc5(ij,n))
        error = error + (frc9(ij,n)-x(ij,n))*(frc9(ij,n)-x(ij,n))
        sol   = sol   + x(ij,n)*x(ij,n)
     ENDDO
  ENDDO
  WRITE(*,*) 'Solution L2-norm ..:', sqrt(sol)
  WRITE(*,*) 'Relative error ....:', sqrt(error)/sqrt(sol)
  WRITE(*,*) 'Diff 9 and 5 ......:', sqrt(diff)/sqrt(sol)

  STOP
END PROGRAM VNONAD
