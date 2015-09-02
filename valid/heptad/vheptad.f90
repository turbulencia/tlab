PROGRAM VHEPTAD

  IMPLICIT NONE
  
#include "types.h"
#include "integers.h"

  TINTEGER,               parameter :: nmax = 32, len = 5
  TREAL, dimension(nmax)            :: a, b, c, d, e, f, g, h, i
  TREAL, dimension(nmax)            :: a5, b5, c5, d5, e5
  TREAL, dimension(len,nmax)        :: x, frc9, frc7, frc5

  TINTEGER n, ij, seed
  TREAL RAN0, error, sol, diff

! ###################################################################
  seed = 256

! create diagonals
  DO n = 1,nmax
!     a(n) = RAN0(seed)
!     b(n) = RAN0(seed)
     a(n) = C_0_R
     b(n) = C_0_R
     c(n) = RAN0(seed); a5(n) = c(n)
     d(n) = RAN0(seed); b5(n) = d(n)
     e(n) = RAN0(seed); c5(n) = e(n)
     f(n) = RAN0(seed); d5(n) = f(n)
     g(n) = RAN0(seed); e5(n) = g(n)
!     h(n) = RAN0(seed)
!     i(n) = RAN0(seed)
     h(n) = C_0_R
     i(n) = C_0_R
  ENDDO

! padding
  a(1) = C_0_R; a(2) = C_0_R; a(3) = C_0_R; a(4) = C_0_R
  b(1) = C_0_R; b(2) = C_0_R; b(3) = C_0_R
  c(1) = C_0_R; c(2) = C_0_R
  d(1) = C_0_R
  i(nmax) = C_0_R; i(nmax-1) = C_0_R; i(nmax-2) = C_0_R; i(nmax-3) = C_0_R
  h(nmax) = C_0_R; h(nmax-1) = C_0_R; h(nmax-2) = C_0_R
  g(nmax) = C_0_R; g(nmax-1) = C_0_R 
  f(nmax) = C_0_R

  a5(1) = C_0_R; a5(2) = C_0_R
  b5(1) = C_0_R
  e5(nmax) = C_0_R; e5(nmax-1) = C_0_R
  d5(nmax) = C_0_R

! create solution
  DO n = 1,nmax
     DO ij = 1,len
        x(ij,n) = RAN0(seed)
     ENDDO
  ENDDO

! compute forcing term
  DO n = 1,nmax
     DO ij = 1,len
        frc9(ij,n) = x(ij,n-4)*a(n) + x(ij,n-3)*b(n) + x(ij,n-2)*c(n) + x(ij,n-1)*d(n) + x(ij,n)*e(n)&
                   + x(ij,n+1)*f(n) + x(ij,n+2)*g(n) + x(ij,n+3)*h(n) + x(ij,n+4)*i(n)
!        frc7(ij,n) = x(ij,n-3)*a(n) + x(ij,n-2)*b(n) + x(ij,n-1)*c(n) + x(ij,n)*d(n)&
!                   + x(ij,n+1)*e(n) + x(ij,n+2)*f(n) + x(ij,n+3)*g(n)
        frc5(ij,n) = frc9(ij,n)
     ENDDO
  ENDDO

! ###################################################################
! solve system
  CALL NONADFS(nmax, i1,  a, b, c, d, e, f, g, h, i)
  CALL NONADSS(nmax, len, a, b, c, d, e, f, g, h, i, frc9)

! solve system
!  CALL HEPTADFS(nmax, a, b, c, d, e, f, g)
!  CALL HEPTADSS(nmax, len, a, b, c, d, e, f, g, frc7)

! solve system pentadiagonal
  CALL PENTADFS(nmax,      a5, b5, c5, d5, e5)
  CALL PENTADSS(nmax, len, a5, b5, c5, d5, e5, frc5)

! error
  diff  = C_0_R
  error = C_0_R
  sol   = C_0_R
  DO n = 1,nmax
     DO ij = 1,len
!        diff  = diff  + (frc7(ij,n)-frc5(ij,n))*(frc7(ij,n)-frc5(ij,n))
!        error = error + (frc7(ij,n)-x(ij,n))*(frc7(ij,n)-x(ij,n))
        diff  = diff  + (frc9(ij,n)-frc5(ij,n))*(frc9(ij,n)-frc5(ij,n))
        error = error + (frc9(ij,n)-x(ij,n))*(frc9(ij,n)-x(ij,n))
        sol   = sol   + x(ij,n)*x(ij,n)
     ENDDO
  ENDDO
  WRITE(*,*) 'Solution L2-norm ..:', sqrt(sol)
  WRITE(*,*) 'Relative error ....:', sqrt(error)/sqrt(sol)
  WRITE(*,*) 'Diff ? and 5 ......:', sqrt(diff)/sqrt(sol)

  STOP
END PROGRAM VHEPTAD
