#include "types.h"

PROGRAM VHELMHOLTZ_FDE

  USE DNS_GLOBAL, ONLY : jmax_total

  IMPLICIT NONE
  
#include "integers.h"
  
  TINTEGER imode_fdm, imax, jmax, kmax
  TINTEGER len, i, l, wk, i1bc, iunif, inb_grid
  PARAMETER(imax=2048,len=1,inb_grid=2+4*3+4*3)
  TREAL scalex
  TREAL x(imax), dx(imax,inb_grid)
  TREAL, DIMENSION(len,imax) :: u, du1_n, du2_n, u_n, a, tmp1
  TREAL wrk1d(imax,5+2+2+1), wrk2d(imax*len*5), wrk3d(len,imax*5)
  TREAL sol, error, w2, bcs(len,2+1)

! ###################################################################
  scalex = C_2_R*C_PI_R
  jmax  = 1
  kmax  = 1
  jmax_total = imax
  imode_fdm = 6
  iunif = 0
  i1bc  = 1

  WRITE(*,*) 'Wavenumber ?'
  READ(*,*) wk

  WRITE(*,*) 'Eigenvalue ?'
  READ(*,*) w2

  DO i = 1,imax
     x(i) = M_REAL(i-1)/M_REAL(imax-1)*scalex
  ENDDO
!   OPEN(21,file='grid.dat')
!   DO i = 1,imax
!      READ(21,*) x(i)
!   ENDDO
!   scalex=x(imax)-x(1)
!   CLOSE(21)

  CALL DNS_INITIALIZE
 ! TO BE REVIEWED 
!  CALL FDM_INITIALIZE(iunif, imode_fdm, imax, i1bc, scalex, x, dx, wrk1d)

! ###################################################################
! Define the function
  DO i = 1,imax
     DO l = 1,len
! single-mode
      u(l,i) = SIN(C_2_R*C_PI_R/scalex*M_REAL(wk)*x(i)+C_PI_R/C_4_R)
!        u(l,i) =-COS(C_2_R*C_PI_R*x(i))*C_PI_R*C_PI_R
!        u(l,i) = C_025_R*(COS(C_2_R*C_PI_R*x(i))-C_1_R)
! Gaussian
!      u(l,i) = EXP(-(x(i)-C_05_R*scalex)**2/(C_2_R*(scalex/M_REAL(wk/l))**2))
! Exponential
!      u(l,i) = EXP((x(i)-x(imax))*M_REAL(wk/l)/scalex)
! linear
!      u(l,i) = -x(i)
! 2 delta wave
!      u(l,i) = MOD(i,2)
! zero
!        u(l,i) = C_0_R
     ENDDO
  ENDDO

!   OPEN(22,file='field.dat')
!   l = 1
!   DO i = 1,imax
!      READ(22,*) du1_n(l,i)
!   ENDDO
!   CLOSE(22)

! ###################################################################
  CALL PARTIAL_YY(i1,iunif,imode_fdm, len,imax,i1, i1bc, dx, u,     du2_n, &
       i0,i0, i0,i0, du1_n, wrk1d,wrk2d,wrk3d)
!  CALL PARTIAL_Y(imode_fdm, len,imax,i1, i1bc, dx, u,     du1_n, i0,i0, wrk1d,wrk2d,wrk3d)
!  CALL PARTIAL_Y(imode_fdm, len,imax,i1, i1bc, dx, du1_n, du2_n, i0,i0, wrk1d,wrk2d,wrk3d)
!  a = du1_n
  a = du2_n + w2*u 

  bcs(:,1) = u(:,1); bcs(:,2) = u(:,imax)
!  bcs(:,1) = C_0_R; bcs(:,2) = u(:,imax)

!  CALL FDE_BVP_SINGULAR_DN(imode_fdm, imax,len,       dx, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_SINGULAR_ND(imode_fdm, imax,len,       dx, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_SINGULAR_NN(imode_fdm, imax,len,       dx, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_SINGULAR_DD(imode_fdm, imax,len,     x,dx, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_REGULAR_NN(imode_fdm, imax,len, w2,    dx, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_REGULAR_DD(imode_fdm, imax,len, w2,    dx, u_n,a,bcs, tmp1, wrk1d)

!  CALL PARTIAL_Y(imode_fdm, len,imax,i1, i1bc, dx, u_n,   du1_n, i0,i0, wrk1d,wrk2d,wrk3d)
!  CALL PARTIAL_Y(imode_fdm, len,imax,i1, i1bc, dx, du1_n, du2_n, i0,i0, wrk1d,wrk2d,wrk3d)
!  u_n = du1_n !- w2*u_n

  u_n(:,1)=bcs(:,1); u_n(:,imax) = bcs(:,2)
  ! CALL LIN_C2N6(i0, imax, len, w2, dx, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),&
  !      a, u_n)
  CALL PENTADFS(imax,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
  CALL PENTADSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_n)

! ###################################################################
  OPEN(20,file='helmholtz.dat')
  error = C_0_R; sol = C_0_R
  DO i = 1,imax!-1
     DO l = 1,len
        WRITE(20,1000) x(i), u(l,i), u_n(l,i), u(l,i)-u_n(l,i), tmp1(l,i)-du1_n(l,i)
        error = error + (u(l,i)-u_n(l,i))*(u(l,i)-u_n(l,i))
        sol   = sol   +  u(l,i)          * u(l,i)
     ENDDO
  ENDDO
  CLOSE(20)
  IF ( sol .GT. C_0_R ) WRITE(*,*) 'Relative error solution ........:', sqrt(error)/sqrt(sol)
                        WRITE(*,*) 'L2-norm ........................:', sqrt(sol)

  error = C_0_R; sol = C_0_R
  DO i = 1,imax!-1
     DO l = 1,len
        error = error + (tmp1(l,i)-du1_n(l,i))*(tmp1(l,i)-du1_n(l,i))
        sol   = sol   +  du1_n(l,i)           * du1_n(l,i)
     ENDDO
  ENDDO

  IF ( sol .GT. C_0_R ) WRITE(*,*) 'Relative error 1st-helmholtz....:', sqrt(error)/sqrt(sol)
                        WRITE(*,*) 'L2-norm ........................:', sqrt(sol)

  STOP

1000 FORMAT(5(1x,e18.10e3)) 

END PROGRAM VHELMHOLTZ_FDE
