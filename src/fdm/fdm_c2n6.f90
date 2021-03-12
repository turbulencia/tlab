#include "types.h"

!########################################################################
!# Tool/Library PADE
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/09/05 - J.P. Mellado
!#              Optimized
!# 2007/12/02 - J.P. Mellado
!#              Nonuniform case completed
!# 2009/12/09 - J.P. Mellado
!#              Splitting into two routines
!#
!########################################################################
!# DESCRIPTION
!#
!# Implementation of the second derivative finite difference with
!# 6th order tridiagonal compact scheme by JCP Lele 1992, nonperiodic.
!# Interior points according to Eq. 2.2.7
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

#define C_2P2D51_L .4313725490196078d-1
#define C_242D51_L .4745098039215686d+1
#define C_44D85_L  .5176470588235294d+0

!########################################################################
!Left-hand side; tridiagonal matrix of the linear system
!########################################################################
SUBROUTINE FDM_C2N6_LHS(imax, imin_set_zero,imax_set_zero, dx, a,b,c)

  IMPLICIT NONE

  TINTEGER,                  INTENT(IN) :: imax
  TINTEGER,                  INTENT(IN) :: imin_set_zero, imax_set_zero ! Set endpoint derivative to zero
  TREAL,   DIMENSION(imax,2),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),  INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i

! #######################################################################
  a(1) = C_0_R
  b(1) = C_22D51_L  *dx(1,1) *dx(1,1)
  c(1) = C_242D51_L *dx(2,1) *dx(2,1)

  a(2) = C_2P2D51_L *dx(1,1) *dx(1,1)
  b(2) = C_22D51_L  *dx(2,1) *dx(2,1)
  c(2) = C_2P2D51_L *dx(3,1) *dx(3,1)

  DO i = 3,imax-2
     a(i) = C_04D51_L *dx(i-1,1)*dx(i-1,1)
     b(i) = C_22D51_L *dx(i,  1)*dx(i,  1)
     c(i) = C_04D51_L *dx(i+1,1)*dx(i+1,1)
  ENDDO

  a(imax-1) = C_2P2D51_L *dx(imax-2,1) *dx(imax-2,1)
  b(imax-1) = C_22D51_L  *dx(imax-1,1) *dx(imax-1,1)
  c(imax-1) = C_2P2D51_L *dx(imax,  1) *dx(imax,  1)
  
  a(imax  ) = C_242D51_L *dx(imax-1,1) *dx(imax-1,1)
  b(imax  ) = C_22D51_L  *dx(imax,  1) *dx(imax,  1)
  c(imax  ) = C_0_R

! Set endpoint derivative to zero
  IF ( imin_set_zero .EQ. 1 ) c(1)    = C_0_R
  IF ( imax_set_zero .EQ. 1 ) a(imax) = C_0_R

  RETURN
END SUBROUTINE FDM_C2N6_LHS

! #######################################################################
! Right-hand side; forcing term, uniform case
! #######################################################################
SUBROUTINE FDM_C2N6_RHS(imax,jkmax, imin_set_zero,imax_set_zero, u,d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TINTEGER,                      INTENT(IN) :: imin_set_zero, imax_set_zero ! Set endpoint derivative to zero
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u ! function to be diferentiated
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d ! right-hand side vector of the linear system

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
        d(jk,2)     = C_44D85_L*(         u(jk,1)    -C_2_R  *u(jk,2)      +        u(jk,3)                    )
     ENDDO
  ELSE
     DO jk = srt,end
        d(jk,1)     = C_22D51_L*( C_13_R *u(jk,1)    -C_27_R *u(jk,2)      +C_15_R *u(jk,3)      -u(jk,4)      )
        d(jk,2)     = C_44D85_L*(         u(jk,1)    -C_2_R  *u(jk,2)      +        u(jk,3)                    )
     ENDDO
  ENDIF

  IF ( imax_set_zero .EQ. 1 ) THEN
     DO jk = srt,end
        d(jk,imax-1)= C_44D85_L*(         u(jk,imax) -C_2_R  *u(jk,imax-1) +        u(jk,imax-2)               )
        d(jk,imax  )= C_0_R
     ENDDO
  ELSE
     DO jk = srt,end
        d(jk,imax-1)= C_44D85_L*(         u(jk,imax) -C_2_R  *u(jk,imax-1) +        u(jk,imax-2)               )
        d(jk,imax  )= C_22D51_L*( C_13_R *u(jk,imax) -C_27_R *u(jk,imax-1) +C_15_R *u(jk,imax-2) -u(jk,imax-3) )
     ENDDO
  ENDIF

  DO i = 3,imax-2 ! Interior points
     DO jk = srt,end
        d(jk,i) = C_24D51_L*(u(jk,i+1)+u(jk,i-1)) + C_3D102_L*(u(jk,i+2)+u(jk,i-2)) - u(jk,i)
     ENDDO
  ENDDO

!$omp end parallel

  RETURN
END SUBROUTINE FDM_C2N6_RHS

! #######################################################################
! Right-hand side; forcing term, non-uniform case using Jacobian
! #######################################################################
SUBROUTINE FDM_C2N6NJ_RHS(imax,jkmax, imin_set_zero,imax_set_zero, dx, u,up,d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: imax, jkmax
  TINTEGER,                      INTENT(IN) :: imin_set_zero, imax_set_zero ! Set endpoint derivative to zero
  TREAL,   DIMENSION(imax,2),    INTENT(IN) :: dx
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u  ! function to be diferentiated
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: up ! first derivative, needed in non-uniform case
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d  ! right-hand side vector of the linear system

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
        d(jk,2)     = C_44D85_L*(         u(jk,1)    -C_2_R  *u(jk,2)      +        u(jk,3)                    )&
                    - ( C_2P2D51_L*dx(1,2)*up(jk,1) +C_22D51_L *dx(2,2)*up(jk,2) +C_2P2D51_L*dx(3,2)*up(jk,3) )
     ENDDO
  ELSE
     DO jk = srt,end
        d(jk,1)     = C_22D51_L*( C_13_R *u(jk,1)    -C_27_R *u(jk,2)      +C_15_R *u(jk,3)      -u(jk,4)      )&
                    - ( C_22D51_L *dx(1,2)*up(jk,1) +C_242D51_L*dx(2,2)*up(jk,2) )
        d(jk,2)     = C_44D85_L*(         u(jk,1)    -C_2_R  *u(jk,2)      +        u(jk,3)                    )&
                    - ( C_2P2D51_L*dx(1,2)*up(jk,1) +C_22D51_L *dx(2,2)*up(jk,2) +C_2P2D51_L*dx(3,2)*up(jk,3) )
     ENDDO
  ENDIF

  IF ( imax_set_zero .EQ. 1 ) THEN
     DO jk = srt,end
        d(jk,imax-1)= C_44D85_L*(         u(jk,imax) -C_2_R  *u(jk,imax-1) +        u(jk,imax-2)               )&
                    - ( C_2P2D51_L*dx(imax-2,2)*up(jk,imax-2) +C_22D51_L *dx(imax-1,2)*up(jk,imax-1) +C_2P2D51_L*dx(imax,  2)*up(jk,imax  ) )           
        d(jk,imax)  = C_0_R
     ENDDO
  ELSE
     DO jk = srt,end
        d(jk,imax-1)= C_44D85_L*(         u(jk,imax) -C_2_R  *u(jk,imax-1) +        u(jk,imax-2)               )&
                    - ( C_2P2D51_L*dx(imax-2,2)*up(jk,imax-2) +C_22D51_L *dx(imax-1,2)*up(jk,imax-1) +C_2P2D51_L*dx(imax,  2)*up(jk,imax  ) )
        d(jk,imax  )= C_22D51_L*( C_13_R *u(jk,imax) -C_27_R *u(jk,imax-1) +C_15_R *u(jk,imax-2) -u(jk,imax-3) )&
                    - ( C_242D51_L*dx(imax-1,2)*up(jk,imax-1) +C_22D51_L *dx(imax,  2)*up(jk,imax  ) )
     ENDDO
  ENDIF

  DO i = 3,imax-2 ! Interior points
     DO jk = srt,end
        d(jk,i) = C_24D51_L*(u(jk,i+1)+u(jk,i-1)) + C_3D102_L*(u(jk,i+2)+u(jk,i-2)) - u(jk,i)&
                - ( C_04D51_L*dx(i-1,2)*up(jk,i-1) + C_22D51_L*dx(i,  2)*up(jk,i  ) + C_04D51_L*dx(i+1,2)*up(jk,i+1) )
     ENDDO
  ENDDO

!$omp end parallel

  RETURN
END SUBROUTINE FDM_C2N6NJ_RHS
