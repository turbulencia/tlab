program VMATRIX

    implicit none

#include "types.h"

    TINTEGER, parameter :: nmax = 5, len = 1
    TREAL, dimension(nmax, 5) :: a
    TREAL, dimension(nmax, nmax) :: c, d, wrk2d
!  TREAL, dimension(len,nmax)        :: x, f3, f5

    TINTEGER n, ij, seed
    TREAL RAN0 !, error, sol

! ###################################################################
#define a_a(n) a(n,1)
#define a_b(n) a(n,2)
#define a_c(n) a(n,3)
#define a_d(n) a(n,4)
#define a_e(n) a(n,5)

#define b_a(n) b(n,1)
#define b_b(n) b(n,2)
#define b_c(n) b(n,3)
#define b_d(n) b(n,4)
#define b_e(n) b(n,5)

    seed = 256

! create diagonals
    do n = 1, nmax
        a_a(n) = RAN0(seed)
        a_b(n) = RAN0(seed)
        a_c(n) = RAN0(seed)
        a_d(n) = RAN0(seed)
        a_e(n) = RAN0(seed)
    end do

! padding
    a_a(1) = C_0_R; a_a(2) = C_0_R; 
    a_b(1) = C_0_R
    a_e(nmax - 1) = C_0_R; a_e(nmax) = C_0_R; 
    a_d(nmax) = C_0_R

! ###################################################################
    call TRIDFS(nmax, a_b(1), a_c(1), a_d(1))
    do n = 1, nmax
        write (*, *) a_b(n), a_c(n), a_d(n)
    end do

    print *, 'Inverse of L'
    call TRIDINV(nmax, a_b(2), c)
    do n = 1, nmax
        write (*, '(6F10.5)') (c(n, ij), ij=1, nmax)
    end do

    print *, 'Inverse of U'
    call TRIDINV(nmax, a_d(1), wrk2d)
    call TLab_Transpose(wrk2d, nmax, nmax, nmax, d, nmax)
    do n = 1, nmax
        write (*, '(6F10.5)') (d(n, ij), ij=1, nmax)
    end do

! ###################################################################
!! create solution
!  DO n = 1,nmax
!     DO ij = 1,len
!        x(ij,n) = RAN0(seed)
!     ENDDO
!  ENDDO
!
!! compute forcing term
!  DO n = 1,nmax
!     DO ij = 1,len
!        f5(ij,n) = x(ij,n-2)*a_a(n) + x(ij,n-1)*a_b(n) + x(ij,n)*a_c(n) + x(ij,n+1)*a_d(n) + x(ij,n+2)*a_e(n)
!     ENDDO
!  ENDDO
!
!!solve system
!  CALL PENTADFS(nmax,      a_a(1), a_b(1), a_c(1), a_d(1), a_e(1))
!  CALL PENTADSS(nmax, len, a_a(1), a_b(1), a_c(1), a_d(1), a_e(1), f5)
!
!  error = C_0_R
!  sol   = C_0_R
!  DO n = 1,nmax
!     DO ij = 1,len
!        error = error + (f5(ij,n)-x(ij,n))*(f5(ij,n)-x(ij,n))
!        sol   = sol   + x(ij,n)*x(ij,n)
!     ENDDO
!  ENDDO
!  WRITE(*,*) 'Solution L2-norm ..:', sqrt(sol)
!  WRITE(*,*) 'Relative error ....:', sqrt(error)/sqrt(sol)

    stop
end program VMATRIX
