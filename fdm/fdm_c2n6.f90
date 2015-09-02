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
!# ARGUMENTS 
!#
!# u    In    function to be diferentiated
!# up   In    du/dx, needed in nonuniform case
!# d    Out   right-hand side vector of the linear system
!#
!########################################################################
#include "types.h"

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
SUBROUTINE FDM_C2N6_LHS(imax, imin_zet_zero,imax_set_zero, dx, a,b,c)

  IMPLICIT NONE

  TINTEGER,                  INTENT(IN) :: imax, imin_zet_zero, imax_set_zero
  TREAL,   DIMENSION(imax,2),INTENT(IN) :: dx
  TREAL,   DIMENSION(imax),  INTENT(OUT):: a,b,c

! -------------------------------------------------------------------
  TINTEGER i
  TREAL dx1, dx2, dxm, dxn
  TREAL vmult_imin, vmult_imax

! #######################################################################
! Set Endpoint Derivative Zeroing Multiplier
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_zet_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! #######################################################################
  dx1 = dx(1,1)*dx(1,1)
  dx2 = dx(2,1)*dx(2,1)
  dxm = dx(imax-1,1)*dx(imax-1,1)
  dxn = dx(imax,  1)*dx(imax,  1)

  a(1) = C_0_R
  b(1) = C_22D51_L  *dx1
  c(1) = C_242D51_L *dx2 *vmult_imin
  a(2) = C_2P2D51_L *dx1
  b(2) = C_22D51_L  *dx2
  c(2) = C_2P2D51_L *dx(3,1)*dx(3,1)

  DO i = 3,imax-2
     a(i) = C_04D51_L *dx(i-1,1)*dx(i-1,1)
     b(i) = C_22D51_L *dx(i,  1)*dx(i,  1)
     c(i) = C_04D51_L *dx(i+1,1)*dx(i+1,1)
  ENDDO

  a(imax-1) = C_2P2D51_L *dx(imax-2,1)*dx(imax-2,1)
  b(imax-1) = C_22D51_L  *dxm
  c(imax-1) = C_2P2D51_L *dxn
  a(imax)   = C_242D51_L *dxm *vmult_imax
  b(imax)   = C_22D51_L  *dxn
  c(imax)   = C_0_R

  RETURN
END SUBROUTINE FDM_C2N6_LHS

! #######################################################################
! Right-hand side; forcing term
! #######################################################################
SUBROUTINE FDM_C2N6_RHS(iunif, imax,jkmax, imin_zet_zero,imax_set_zero, dx, u,up,d)
#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  IMPLICIT NONE

  TINTEGER,                      INTENT(IN) :: iunif, imax, jkmax, imin_zet_zero, imax_set_zero
  TREAL,   DIMENSION(imax,2),    INTENT(IN) :: dx
  TREAL,   DIMENSION(jkmax,imax),INTENT(IN) :: u, up
  TREAL,   DIMENSION(jkmax,imax),INTENT(OUT):: d

! -------------------------------------------------------------------
  TINTEGER i, jk, srt,end,siz
  TREAL vmult_imin, vmult_imax
  TREAL cl01, cl13, cl27, cl15, cr01, cr13, cr27, cr15
#ifdef USE_OPENMP
  TINTEGER thread
#endif

! #######################################################################
! Set Endpoint Derivative Zeroing Multiplier
  vmult_imin = C_1_R
  vmult_imax = C_1_R
  IF (imin_zet_zero .EQ. 1) vmult_imin = C_0_R
  IF (imax_set_zero .EQ. 1) vmult_imax = C_0_R

! #######################################################################
!$omp parallel default ( none ) &
!$omp private( jk,i, cl01,cl13,cl27,cl15, cr01,cr13,cr27,cr15,srt,end,siz, thread ) &
!$omp shared(jkmax,iunif,d,u,imax,vmult_imin,vmult_imax,dx,up)

  CALL DNS_OMP_PARTITION(jkmax,srt,end,siz)

! Left (i=1) boundary, third-order
  cl01 = C_22D51_L   *vmult_imin
  cl13 = C_13_R*cl01 
  cl27 = C_27_R*cl01 
  cl15 = C_15_R*cl01 
! Right (i=imax) boundary, third-order
  cr01 = C_22D51_L   *vmult_imax 
  cr13 = C_13_R*cr01
  cr27 = C_27_R*cr01
  cr15 = C_15_R*cr01

! -------------------------------------------------------------------
! Uniform case
! -------------------------------------------------------------------
  IF ( iunif .EQ. 0 ) THEN

     DO jk = srt,end
        d(jk,1)     = cl13*u(jk,1) - cl27*u(jk,2) + cl15*u(jk,3) - cl01*u(jk,4)
        d(jk,2)     = C_44D85_L*( u(jk,1)    - C_2_R*u(jk,2)      + u(jk,3)     )
        d(jk,imax-1)= C_44D85_L*( u(jk,imax) - C_2_R*u(jk,imax-1) + u(jk,imax-2))
        d(jk,imax)  = cr13*u(jk,imax)   - cr27*u(jk,imax-1) &
             + cr15*u(jk,imax-2) - cr01*u(jk,imax-3)
     ENDDO

! Interior points
     DO i = 3,imax-2
        DO jk = srt,end
           d(jk,i) = C_24D51_L*(u(jk,i+1)+u(jk,i-1)) + C_3D102_L*(u(jk,i+2)+u(jk,i-2)) - u(jk,i)
        ENDDO
     ENDDO

! -------------------------------------------------------------------
! Nonuniform case
! -------------------------------------------------------------------
  ELSE
     DO jk = srt,end
        d(jk,1)     = cl13*u(jk,1) - cl27*u(jk,2) &
             + cl15*u(jk,3) - cl01*u(jk,4)&
             - ( C_22D51_L* dx(1,2)*up(jk,1) + C_242D51_L*dx(2,2)*up(jk,2) )* vmult_imin

        d(jk,2)     = C_44D85_L*( u(jk,1)-C_2_R*u(jk,2)+u(jk,3) )&
             - ( C_2P2D51_L*dx(1,2)*up(jk,1) + C_22D51_L *dx(2,2)*up(jk,2)&
             + C_2P2D51_L*dx(3,2)*up(jk,3) )

        d(jk,imax-1)= C_44D85_L*( u(jk,imax) - C_2_R*u(jk,imax-1) + u(jk,imax-2) )&
             - ( C_2P2D51_L*dx(imax-2,2)*up(jk,imax-2) + C_22D51_L *dx(imax-1,2)*up(jk,imax-1)&
             + C_2P2D51_L*dx(imax,  2)*up(jk,imax  ) )

        d(jk,imax)  = cr13*u(jk,imax  ) - cr27*u(jk,imax-1) &
             + cr15*u(jk,imax-2) - cr01*u(jk,imax-3)&
             - ( C_242D51_L*dx(imax-1,2)*up(jk,imax-1) &
             + C_22D51_L *dx(imax,  2)*up(jk,imax  ) )* vmult_imax 
     ENDDO

! Interior points
     DO i = 3,imax-2
        DO jk = srt,end
           d(jk,i) = C_24D51_L*(u(jk,i+1)+u(jk,i-1)) + C_3D102_L*(u(jk,i+2)+u(jk,i-2)) - u(jk,i)&
                - ( C_04D51_L*dx(i-1,2)*up(jk,i-1) + &
                C_22D51_L*dx(i,  2)*up(jk,i  ) + C_04D51_L*dx(i+1,2)*up(jk,i+1) )
        ENDDO
     ENDDO

  ENDIF

!$omp end parallel

  RETURN
END SUBROUTINE FDM_C2N6_RHS
