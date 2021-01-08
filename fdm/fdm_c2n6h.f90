#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2020/08/22 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme.
!# Interior points according to JCP, Lamballais et al. 2011, JCP 230:3270-3275
!# Eqs. 1,3 with kc = pi**2.
!# It adds one term in the RHS to Lele's Eq. 2.2.7 scheme to better match the
!# exact transfer function (hyperdiffusive instead of hypodiffusive).
!# Rest according to JCP Lele 1992, nonperiodic.
!# The first point at the boundary from Eq. 4.3.1, third-order
!# The second point from Eq. 2.2.6 forth-order (b=0).
!# The linear system is divided by the f_i coefficient in the right
!# hand side to reduce number of operations.
!#
!########################################################################
#define C_22D51_L .4313725490196078d+0
#define C_04D51_L .7843137254901960d-1
#define C_24D51_L .4705882352941176d+0
#define C_3D102_L .2941176470588235d-1

#define C_242D51_L .4745098039215686d+1

#define C_01D24_L .4166666666666666d-1
#define C_10D24_L .4166666666666666d+0
#define C_05_L    .5000000000000000d+0

#define C_01D30_L .3333333333333333d-1
#define C_01D03_L .3333333333333333d+0
#define C_07D60_L .1166666666666666d+0

!########################################################################
!Left-hand side; tridiagonal matrix of the linear system
!########################################################################
#define C_LHS0_L .522941253359481d+0
#define C_LHS1_L .190603035365364d+0

SUBROUTINE FDM_C2N6H_LHS(imax, imin_set_zero,imax_set_zero, dx, a,b,c)

  IMPLICIT NONE

  TINTEGER,                  INTENT(IN) :: imax
  TINTEGER,                  INTENT(IN) :: imin_set_zero, imax_set_zero ! Set endpoint derivative to zero
  TREAL,   DIMENSION(imax,2),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),  INTENT(OUT):: a,b,c

  ! -------------------------------------------------------------------
  TINTEGER i

  ! #######################################################################
  a(1) = C_0_R
  b(1) =         dx(1,1) *dx(1,1)
  c(1) = C_11_R *dx(2,1) *dx(2,1)

  a(2) = C_01D24_L *dx(1,1) *dx(1,1)  ! Compact centered 4th order
  b(2) = C_10D24_L *dx(2,1) *dx(2,1)
  c(2) = C_01D24_L *dx(3,1) *dx(3,1)
  ! a(2) = C_01D30_L *dx(1,1) *dx(1,1)  ! Compact biased 5th order
  ! b(2) = C_01D03_L *dx(2,1) *dx(2,1)
  ! c(2) =-C_07D60_L *dx(3,1) *dx(3,1)

  a(3) = C_04D51_L *dx(2,1)*dx(2,1)
  b(3) = C_22D51_L *dx(3,1)*dx(3,1)
  c(3) = C_04D51_L *dx(4,1)*dx(4,1)

  DO i = 4,imax-3
    a(i) = C_LHS1_L *dx(i-1,1)*dx(i-1,1)
    b(i) = C_LHS0_L *dx(i,  1)*dx(i,  1)
    c(i) = C_LHS1_L *dx(i+1,1)*dx(i+1,1)
  ENDDO

  a(imax-2) = C_04D51_L *dx(imax-3,1) *dx(imax-3,1)
  b(imax-2) = C_22D51_L *dx(imax-2,1) *dx(imax-2,1)
  c(imax-2) = C_04D51_L *dx(imax-1,1) *dx(imax-1,1)

  a(imax-1) = C_01D24_L *dx(imax-2,1) *dx(imax-2,1)  ! Compact centered 4th order
  b(imax-1) = C_10D24_L *dx(imax-1,1) *dx(imax-1,1)
  c(imax-1) = C_01D24_L *dx(imax,  1) *dx(imax,  1)
  ! a(imax-1) =-C_07D60_L *dx(imax-2,1) *dx(imax-2,1)  ! Compact biased 5th order
  ! b(imax-1) = C_01D03_L *dx(imax-1,1) *dx(imax-1,1)
  ! c(imax-1) = C_01D30_L *dx(imax,  1) *dx(imax,  1)

  a(imax  ) = C_11_R *dx(imax-1,1) *dx(imax-1,1)
  b(imax  ) =         dx(imax,  1) *dx(imax,  1)
  c(imax  ) = C_0_R

  ! Set endpoint derivative to zero
  IF ( imin_set_zero .EQ. 1 ) c(1)    = C_0_R
  IF ( imax_set_zero .EQ. 1 ) a(imax) = C_0_R

  RETURN
END SUBROUTINE FDM_C2N6H_LHS

! #######################################################################
! Right-hand side; forcing term, uniform case
! #######################################################################
#define C_RHS1_L .355555050467152d+0
#define C_RHS2_L .150282454434515d+0
#define C_RHS3_L .583750490166691d-2

#define C_33D80_L .4125d+0
#define C_62D80_L .7750d+0
#define C_01D05_L .2000d+0
#define C_01D80_L .1250d-1

SUBROUTINE FDM_C2N6H_RHS(imax,jkmax, imin_set_zero,imax_set_zero, u,d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TINTEGER,                      INTENT(IN) :: imin_set_zero, imax_set_zero ! Set endpoint derivative to zero
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u                            ! function to be diferentiated
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d                            ! right-hand side vector of the linear system

  ! -------------------------------------------------------------------
  TINTEGER i, jk, srt,end,siz
#ifdef USE_OPENMP
  TINTEGER thread
#endif

  ! #######################################################################
  !$omp parallel default ( none ) &
  !$omp private( jk,i, srt,end,siz, thread ) &
  !$omp shared(jkmax,d,u,imax,imin_set_zero,imax_set_zero)

  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz)

  IF ( imin_set_zero .EQ. 1 ) THEN
    DO jk = srt,end
      d(jk,1)     = C_0_R
      d(jk,2)     = C_05_L   *( u(jk,1)+u(jk,3) )                                 -u(jk,2) ! Compact centered 4th order
      ! d(jk,2)     = C_33D80_L*u(jk,1) +C_62D80_L*u(jk,3) -C_01D05_L*u(jk,4) +C_01D80_L*u(jk,5) -u(jk,2)  ! Compact biased 5th order
      d(jk,3)     = C_24D51_L*( u(jk,4)+u(jk,2) ) + C_3D102_L*( u(jk,5)+u(jk,1) ) -u(jk,3)
    ENDDO
  ELSE
    DO jk = srt,end
      d(jk,1)     = C_13_R *u(jk,1)    -C_27_R *u(jk,2)      +C_15_R *u(jk,3)      -u(jk,4)
      d(jk,2)     = C_05_L   *( u(jk,1)+u(jk,3) )                                 -u(jk,2) ! Compact centered 4th order
      ! d(jk,2)     = C_33D80_L*u(jk,1) +C_62D80_L*u(jk,3) -C_01D05_L*u(jk,4) +C_01D80_L*u(jk,5) -u(jk,2)  ! Compact biased 5th order
      d(jk,3)     = C_24D51_L*( u(jk,4)+u(jk,2) ) + C_3D102_L*( u(jk,5)+u(jk,1) )              -u(jk,3)
    ENDDO
  ENDIF

  IF ( imax_set_zero .EQ. 1 ) THEN
    DO jk = srt,end
      d(jk,imax-2)= C_24D51_L*( u(jk,imax-1)+u(jk,imax-3) ) +C_3D102_L *( u(jk,imax)+u(jk,imax-4) ) -u(jk,imax-2)
      d(jk,imax-1)= C_05_L   *( u(jk,imax  )+u(jk,imax-2) )                                         -u(jk,imax-1) ! Compact centered 4th order
      ! d(jk,imax-1)= C_33D80_L*u(jk,imax  ) +C_62D80_L*u(jk,imax-2) -C_01D05_L*u(jk,imax-3) +C_01D80_L*u(jk,imax-4) -u(jk,imax-1) ! Compact centered 4th order
      d(jk,imax  )= C_0_R
    ENDDO
  ELSE
    DO jk = srt,end
      d(jk,imax-2)= C_24D51_L*( u(jk,imax-1)+u(jk,imax-3) ) +C_3D102_L*( u(jk,imax)+u(jk,imax-4) ) -u(jk,imax-2)
      d(jk,imax-1)= C_05_L   *( u(jk,imax  )+u(jk,imax-2) )                                        -u(jk,imax-1) ! Compact centered 4th order
      ! d(jk,imax-1)= C_33D80_L*u(jk,imax  ) +C_62D80_L*u(jk,imax-2) -C_01D05_L*u(jk,imax-3) +C_01D80_L*u(jk,imax-4) -u(jk,imax-1) ! Compact centered 4th order
      d(jk,imax  )= C_13_R *u(jk,imax) -C_27_R *u(jk,imax-1) +C_15_R *u(jk,imax-2) -u(jk,imax-3)
    ENDDO
  ENDIF

  DO i = 4,imax-3 ! Interior points
    DO jk = srt,end
      d(jk,i) = C_RHS1_L*( u(jk,i+1)+u(jk,i-1) ) +C_RHS2_L*( u(jk,i+2)+u(jk,i-2) ) -C_RHS3_L*( u(jk,i+3)+u(jk,i-3) ) - u(jk,i)
    ENDDO
  ENDDO

  !$omp end parallel

  RETURN
END SUBROUTINE FDM_C2N6H_RHS

! #######################################################################
! Right-hand side; forcing term, non-uniform case using Jacobian
! #######################################################################
SUBROUTINE FDM_C2N6HNJ_RHS(imax,jkmax, imin_set_zero,imax_set_zero, dx, u,up,d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TINTEGER,                      INTENT(IN) :: imin_set_zero, imax_set_zero ! Set endpoint derivative to zero
  TREAL,   DIMENSION(imax,2),    INTENT(IN) :: dx
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u                            ! function to be diferentiated
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: up                           ! first derivative, needed in non-uniform case
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d                            ! right-hand side vector of the linear system

  ! -------------------------------------------------------------------
  TINTEGER i, jk, srt,end,siz
#ifdef USE_OPENMP
  TINTEGER thread
#endif

  ! #######################################################################
  !$omp parallel default ( none ) &
  !$omp private( jk,i, srt,end,siz, thread ) &
  !$omp shared(jkmax,d,u,imax,imin_set_zero,imax_set_zero,dx,up)

  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz)

  IF ( imin_set_zero .EQ. 1 ) THEN
    DO jk = srt,end
      d(jk,1)     = C_0_R
      d(jk,2)     = C_05_L   *( u(jk,1)+u(jk,3) )                                -u(jk,2)&
      -C_01D24_L*( dx(1,2)*up(jk,1) +dx(3,2)*up(jk,3) ) -C_10D24_L *dx(2,2)*up(jk,2)
      d(jk,3)     = C_24D51_L*( u(jk,4)+u(jk,2) ) +C_3D102_L*( u(jk,5)+u(jk,1) ) -u(jk,3)&
      -C_04D51_L*( dx(2,2)*up(jk,2) +dx(4,2)*up(jk,4) ) -C_22D51_L *dx(3,2)*up(jk,3)
    ENDDO
  ELSE
    DO jk = srt,end
      d(jk,1)     = C_13_R *u(jk,1)    -C_27_R *u(jk,2)      +C_15_R *u(jk,3)      -u(jk,4)      &
      - ( dx(1,2)*up(jk,1) +C_11_R*dx(2,2)*up(jk,2) )
      d(jk,2)     = C_05_L   *( u(jk,1)+u(jk,3) )                                -u(jk,2)&
      -C_01D24_L*( dx(1,2)*up(jk,1) +dx(3,2)*up(jk,3) ) -C_10D24_L *dx(2,2)*up(jk,2)
      d(jk,3)     = C_24D51_L*( u(jk,4)+u(jk,2) ) +C_3D102_L*( u(jk,5)+u(jk,1) ) -u(jk,3)&
      -C_04D51_L*( dx(2,2)*up(jk,2) +dx(4,2)*up(jk,4) ) -C_22D51_L *dx(3,2)*up(jk,3)
    ENDDO
  ENDIF

  IF ( imax_set_zero .EQ. 1 ) THEN
    DO jk = srt,end
      d(jk,imax-2)= C_24D51_L*( u(jk,imax-1)+u(jk,imax-3) ) + C_3D102_L*( u(jk,imax)+u(jk,imax-4) ) -u(jk,imax-2)&
      -C_04D51_L*( dx(imax-3,2)*up(jk,imax-3) +dx(imax-1,2)*up(jk,imax-1) ) -C_22D51_L *dx(imax-2,2)*up(jk,imax-2)
      d(jk,imax-1)= C_05_L   *( u(jk,imax  )+u(jk,imax-2) )                                         -u(jk,imax-1)&
      -C_01D24_L*( dx(imax-2,2)*up(jk,imax-2) +dx(imax,  2)*up(jk,imax  ) ) -C_10D24_L *dx(imax-1,2)*up(jk,imax-1)
      d(jk,imax)  = C_0_R
    ENDDO
  ELSE
    DO jk = srt,end
      d(jk,imax-2)= C_24D51_L*( u(jk,imax-1)+u(jk,imax-3) ) + C_3D102_L*( u(jk,imax)+u(jk,imax-4) ) -u(jk,imax-2)&
      -C_04D51_L*( dx(imax-3,2)*up(jk,imax-3) +dx(imax-1,2)*up(jk,imax-1) ) -C_22D51_L *dx(imax-2,2)*up(jk,imax-2)
      d(jk,imax-1)= C_05_L   *( u(jk,imax  )+u(jk,imax-2) )                                         -u(jk,imax-1)&
      -C_01D24_L*( dx(imax-2,2)*up(jk,imax-2) +dx(imax,  2)*up(jk,imax  ) ) -C_10D24_L *dx(imax-1,2)*up(jk,imax-1)
      d(jk,imax  )= C_13_R *u(jk,imax) -C_27_R *u(jk,imax-1) +C_15_R *u(jk,imax-2) -u(jk,imax-3) &
      - ( C_11_R*dx(imax-1,2)*up(jk,imax-1) +dx(imax,  2)*up(jk,imax  ) )
    ENDDO
  ENDIF

  DO i = 4,imax-3 ! Interior points
    DO jk = srt,end
      d(jk,i) = C_RHS1_L*( u(jk,i+1)+u(jk,i-1) ) +C_RHS2_L*( u(jk,i+2)+u(jk,i-2) ) -C_RHS3_L*( u(jk,i+3)+u(jk,i-3) ) - u(jk,i) &
      -C_LHS1_L*( dx(i-1,2)*up(jk,i-1) +dx(i+1,2)*up(jk,i+1) ) -C_LHS0_L*dx(i,  2)*up(jk,i)
    ENDDO
  ENDDO

  !$omp end parallel

  RETURN
END SUBROUTINE FDM_C2N6HNJ_RHS
