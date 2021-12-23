!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/23 - J. Kostelecky
!#              Created
!#              
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate compact interpolation schemes.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

PROGRAM INTERPOL

  USE TLAB_TYPES, ONLY : grid_dt
  
  IMPLICIT NONE
 
#include "integers.h"
  
  TYPE(grid_dt)                      :: g
   
  TINTEGER                           :: jmax,kmax, i, l, test_type
  
  TREAL                              :: lambda, error, sol

  TINTEGER, PARAMETER                :: imax=100, len=1, inb_grid=44
  
  TREAL,    DIMENSION(imax,inb_grid) :: x
  TREAL,    DIMENSION(imax)          :: x_int, x_aux
  TREAL,    DIMENSION(len,imax)      :: u, u_int, u_aux, u_a, u_b, u_c
  TREAL,    DIMENSION(len,imax)      :: dudx, dudx_int, dudx_aux
  TREAL,    DIMENSION(imax,5)        :: wrk1d
  TREAL,    DIMENSION(len)           :: wrk2d

! ###################################################################
! Initialize
  g%size     = imax 
  g%scale    = C_1_R
  g%uniform  = .TRUE.
  jmax       = 1
  kmax       = 1
  g%mode_fdm = FDM_COM6_JACOBIAN 
   
! Valid stettings
  g%periodic = .TRUE. ! non-periodic not implmented yet !
  lambda     = 1               
  test_type  = 4
 
! ###################################################################
! Initialize grid    
  IF ( g%periodic ) THEN
    DO i = 1,imax
      x(i,1) = M_REAL(i-1)/M_REAL(imax)*g%scale
    ENDDO
  ELSE
    DO i = 1,imax
      x(i,1) = M_REAL(i-1)/M_REAL(imax-1)*g%scale
    ENDDO
  ENDIF

  CALL FDM_INITIALIZE(x, g, wrk1d)

! Interpolation grid -- midpoints
  IF ( g%periodic ) THEN
    DO i = 1,imax
      x_int(i) = g%nodes(i) + 0.5 * g%jac(i,1)
    ENDDO
  ELSE
    DO i = 1,imax-1
      x_int(i) = g%nodes(i) + 0.5 * g%jac(i,1)
    ENDDO
    x_int(imax) = C_0_R
  ENDIF
  x_aux(:) = x_int(:)

! Define the function on both grids
  DO i = 1,imax
    DO l = 1,len
      u(l,i)        =                                              &
                       SIN(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))
      u_int(l,i)    =                                              &
                       SIN(C_2_R*C_PI_R/g%scale*lambda*x_int(i)  )
      dudx(l,i)     = (C_2_R*C_PI_R/g%scale*lambda)                &
                      *COS(C_2_R*C_PI_R/g%scale*lambda*g%nodes(i))
      dudx_int(l,i) = (C_2_R*C_PI_R/g%scale*lambda)                &
                      *COS(C_2_R*C_PI_R/g%scale*lambda*x_int(i)  )
      u_aux(l,i)    = u_int(l,i)
      dudx_aux(l,i) = dudx_int(l,i)
    ENDDO
  ENDDO

! Switch grid and u (for test_type==2)
  IF ( test_type .EQ. 2 .OR. test_type .EQ. 4) THEN 
    DO i = 1,imax
      x_int(i) = x(i,1) 
      x(i,1)   = x_aux(i) 
      DO l = 1,len
        u_int(l,i)    = u(l,i) 
        u(l,i)        = u_aux(l,i)
        ! 
        dudx_int(l,i) = dudx(l,i) 
        dudx(l,i)     = dudx_aux(l,i) 
      ENDDO
    ENDDO
  ENDIF
! ###################################################################
! Testing interpolation and interpolatory 1st derivative
  IF (     test_type .EQ. 1 .AND. g%periodic .EQV. .TRUE. ) THEN 
    WRITE(*,*) '.......... Interpolation to the right .......... '
    CALL FDM_CINT6P_LHS( imax,  wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
    CALL FDM_CINT6PR_RHS(imax,len, u, u_a)
    CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
    CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
  ELSEIF ( test_type .EQ. 2 .AND. g%periodic .EQV. .TRUE. ) THEN 
    WRITE(*,*) '.......... Interpolation to the left ........... '
    CALL FDM_CINT6P_LHS( imax,  wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
    CALL FDM_CINT6PL_RHS(imax,len, u, u_a)
    CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
    CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
  ELSEIF ( test_type .EQ. 3 .AND. g%periodic .EQV. .TRUE. ) THEN 
    WRITE(*,*) '.... 1st Interpolatory deriv to the right ..... '
    CALL FDM_C1INT6P_LHS( imax, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
    CALL FDM_C1INT6PR_RHS(imax,len, u, u_a)
    CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
    CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
    DO i = 1,imax
      DO l = 1,len
        u_int(l,i) = dudx_int(l,i) 
      ENDDO
    ENDDO
  ELSEIF ( test_type .EQ. 4 .AND. g%periodic .EQV. .TRUE. ) THEN 
    WRITE(*,*) '.... 1st Interpolatory deriv to the left ..... '
    CALL FDM_C1INT6P_LHS( imax, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
    CALL FDM_C1INT6PL_RHS(imax,len, u, u_a)
    CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
    CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
    DO i = 1,imax
      DO l = 1,len
        u_int(l,i) = dudx_int(l,i) 
      ENDDO
    ENDDO
  ENDIF
! ###################################################################
! IO - Error and function values
  OPEN(20,file='interpol.dat')
  error = C_0_R
  sol   = C_0_R
  DO i = 1,imax
    DO l = 1,len
      WRITE(20,1000) x(i,1), x_int(i), u(l,i), u_int(l,i), u_a(l,i), u_a(l,i) - u_int(l,i)
      u_c(l,i) = ABS(u_a(l,i) - u_int(l,i))
      error = error + u_c(l,i)   * u_c(l,i)
      sol   = sol   + u_int(l,i) * u_int(l,i)
    ENDDO
  ENDDO
  CLOSE(20)

  WRITE(*,2000) 'Solution L2-norm ...........:', SQRT(g%jac(1,1)*sol) / M_REAL(len)
  IF ( sol .EQ. C_0_R ) STOP
  WRITE(*,2000) 'Error L2-norm ..............:', SQRT(g%jac(1,1)*error) 
  WRITE(*,2000) 'Error Linf-norm ............:', MAXVAL(u_c(1,1:imax))
  WRITE(*,2000) 'Relative error .............:', sqrt(error)/sqrt(sol)

  STOP

1000 FORMAT(6(1x,e16.10))
2000 FORMAT(1(1x,A,3x,e16.10))

END PROGRAM INTERPOL