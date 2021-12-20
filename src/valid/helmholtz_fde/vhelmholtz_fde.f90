#include "types.h"
#include "dns_const.h"


PROGRAM VHELMHOLTZ_FDE

  USE TLAB_TYPES, ONLY : grid_dt

  IMPLICIT NONE
  
#include "integers.h"
  
  TYPE(grid_dt) :: g
  TINTEGER imode_fdm, imax, jmax, kmax
  TINTEGER len, i, l, wk, inb_grid
  PARAMETER(imax=2048,len=1,inb_grid=44)!3+4*3+4*3)
  TREAL, DIMENSION(imax,inb_grid) :: x
  TREAL, DIMENSION(len,imax) :: u, du1_n, du2_n, u_n, a, tmp1
  TREAL wrk1d(imax,5+2+2+1), wrk2d(imax*len*5), wrk3d(len,imax*5)
  TREAL sol, error, w2, bcs(len,2+1)

! ###################################################################
  bcs = 0
  
  g%size     = imax
  g%scale    = C_2_R*C_PI_R
  g%mode_fdm = FDM_COM6_JACOBIAN
  g%uniform  = .TRUE.
  g%periodic = .FALSE.

  jmax = 1
  kmax = 1

  WRITE(*,*) 'Wavenumber ?'
  READ(*,*) wk

  WRITE(*,*) 'Eigenvalue ?'
  READ(*,*) w2

  DO i = 1,imax
     x(i,1) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
  ENDDO
!   OPEN(21,file='grid.dat')
!   DO i = 1,imax
!      READ(21,*) x(i)
!   ENDDO
!   g%scale=x(imax)-x(1)
!   CLOSE(21)

  CALL FDM_INITIALIZE(x, g, wrk1d)

! ###################################################################
! Define the function
  DO i = 1,imax
     DO l = 1,len
! single-mode
      u(l,i) = SIN(C_2_R*C_PI_R/g%scale*M_REAL(wk)*g%nodes(i)+C_PI_R/C_4_R)
!        u(l,i) =-COS(C_2_R*C_PI_R*g%nodes(i))*C_PI_R*C_PI_R
!        u(l,i) = C_025_R*(COS(C_2_R*C_PI_R*g%nodes(i))-C_1_R)
! Gaussian
!      u(l,i) = EXP(-(g%nodes(i)-C_05_R*g%scale)**2/(C_2_R*(g%scale/M_REAL(wk/l))**2))
! Exponential
!      u(l,i) = EXP((g%nodes(i)-g%nodes(imax))*M_REAL(wk/l)/g%scale)
! linear
!      u(l,i) = -g%nodes(i)
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
  CALL OPR_PARTIAL_Y(OPR_P2_P1, len,imax,i1, bcs, g, u, du2_n, du1_n, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_Y(OPR_P1, len,imax,i1, bcs, g, u,     du1_n, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_Y(OPR_P1, len,imax,i1, bcs, g, du1_n, du2_n, wrk3d, wrk2d,wrk3d)
!  a = du1_n
  a = du2_n + w2*u 

  bcs(:,1) = u(:,1); bcs(:,2) = u(:,imax)
!  bcs(:,1) = C_0_R; bcs(:,2) = u(:,imax)

!  CALL FDE_BVP_SINGULAR_DN(OPR_P1, imax,len,       g%jac, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_SINGULAR_ND(OPR_P1, imax,len,       g%jac, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_SINGULAR_NN(OPR_P1, imax,len,       g%jac, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_SINGULAR_DD(OPR_P1, imax,len,     g%nodes,g%jac, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_REGULAR_NN(OPR_P1, imax,len, w2,    g%jac, u_n,a,bcs, tmp1, wrk1d)
!  CALL FDE_BVP_REGULAR_DD(OPR_P1, imax,len, w2,    g%jac, u_n,a,bcs, tmp1, wrk1d)

!  CALL OPR_PARTIAL_Y(OPR_P1, len,imax,i1, bcs, g, u_n,   du1_n, wrk3d, wrk2d,wrk3d)
!  CALL OPR_PARTIAL_Y(OPR_P1, len,imax,i1, bcs, g, du1_n, du2_n, wrk3d, wrk2d,wrk3d)
!  u_n = du1_n !- w2*u_n

  u_n(:,1)=bcs(:,1); u_n(:,imax) = bcs(:,2)
  ! CALL LIN_C2N6(i0, imax, len, w2, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5),&
  !      a, u_n)
  CALL PENTADFS(imax,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
  CALL PENTADSS(imax,len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_n)

! ###################################################################
  OPEN(20,file='helmholtz.dat')
  error = C_0_R; sol = C_0_R
  DO i = 1,imax!-1
     DO l = 1,len
        WRITE(20,1000) g%nodes(i), u(l,i), u_n(l,i), u(l,i)-u_n(l,i), tmp1(l,i)-du1_n(l,i)
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
