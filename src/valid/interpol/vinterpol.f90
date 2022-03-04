!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2021/12/23 - J. Kostelecky
!#              Created           
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
  USE TLAB_PROCS
  
  IMPLICIT NONE
 
#include "integers.h"
  
  TYPE(grid_dt)                       :: g, g_pre
  TINTEGER                            :: jmax,kmax, i, l, test_type, periodic
  TREAL                               :: lambda, error, sol
 
  TINTEGER, PARAMETER                 :: imax=32, len=10, inb_grid=57
  TINTEGER, PARAMETER                 :: imaxp=imax-1 
 
  TREAL,    DIMENSION(imax, inb_grid) :: x
  TREAL,    DIMENSION(imaxp,inb_grid) :: x_pre ! pressure grid (for non-periodic case)
 
  TREAL,    DIMENSION(imax)           :: x_int, x_aux
  TREAL,    DIMENSION(len,imax)       :: u, u_int, u_aux, u_a, u_b, u_c
  TREAL,    DIMENSION(len,imax)       :: dudx, dudx_int, dudx_aux
  TREAL,    DIMENSION(imax,5)         :: wrk1d
  TREAL,    DIMENSION(len)            :: wrk2d

! ###################################################################
! Initialize
  g%size     = imax 
  g%scale    = C_1_R
  g%uniform  = .TRUE.
  jmax       = 1
  kmax       = 1
  g%mode_fdm = FDM_COM6_JACOBIAN 
  lambda     = 1               
   
! Input
  WRITE(*,*) '........... Test type (1/2/3/4)? ...............'
  WRITE(*,*) 'Interpolation from vel. to pre. grid        = 1'
  WRITE(*,*) 'Interpolation from pre. to vel. grid        = 2'
  WRITE(*,*) '1st Interpol. deriv. from vel. to pre. grid = 3'
  WRITE(*,*) '1st Interpol. deriv. from pre. to vel. grid = 4'
  WRITE(*,*) '................................................'
  READ(*,*) test_type
  IF ( test_type .LT. 1 .OR. test_type .GT. 4 ) THEN 
    WRITE(*,*) 'ERROR, STOP. Chose value between 1 and 4.'
    CALL TLAB_STOP(i0)
  ENDIF
  WRITE(*,*) 'Testing     periodic schemes = 0'
  WRITE(*,*) 'Testing non-periodic schemes = 1'
  WRITE(*,*) '................................................'
  READ(*,*) periodic
  IF ( periodic .NE. 0 .AND. periodic .NE. 1) THEN 
    WRITE(*,*) 'ERROR, STOP. Chose value between 0 and 1.'
    CALL TLAB_STOP(i0)
  ENDIF
  IF (     periodic .EQ. 0 ) THEN 
    g%periodic = .TRUE.
    WRITE(*,*) '.......... Testing periodic schemes ............'
  ELSEIF ( periodic .EQ. 1 ) THEN 
    g%periodic = .FALSE.
    WRITE(*,*) '........ Testing non-periodic schemes ..........'
  ENDIF
  WRITE(*,*) '................................................'
! ###################################################################
! Initialize grid 
  IF ( g%periodic ) THEN
    DO i = 1,imax
      x(i,1)   = M_REAL(i-1)/M_REAL(imax)*g%scale
    ENDDO
  ELSE
    DO i = 1,imax
      x(i,1)   = M_REAL(i-1)/M_REAL(imax-1)*g%scale
    ENDDO
  ENDIF

! Velocity grid
  CALL FDM_INITIALIZE(x, g, wrk1d) 

! Initialize grids (interpolation grid on midpoints)  
  IF ( g%periodic ) THEN
    DO i = 1,imax
      x_int(i) = g%nodes(i) + 0.5 * g%jac(i,1)
    ENDDO
  ELSE
    DO i = 1,imaxp
      x_int(i) = g%nodes(i) + 0.5 * g%jac(i,1)
    ENDDO
    x_int(imax) = C_0_R
  ENDIF
  x_aux(:) = x_int(:)

! Initialize pressure grid (only needed for non-periodic case)
! (here: periodic case implies the usage of uniform grids!)
  DO i = 1,imaxp; x_pre(i,1) = x_int(i); ENDDO
  g_pre%size=imaxp; g_pre%scale=x_int(imaxp); g_pre%uniform=g%uniform
  g_pre%mode_fdm=g%mode_fdm; g_pre%periodic=.FALSE. 
  CALL FDM_INITIALIZE(x_pre, g_pre, wrk1d) 

! Define the function + deriv. on both grids
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

! Switch grids and functions (according to interpolation direction [vp <--> pv])
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
! Testing interpolation
! -------------------------------------------------------------------
  IF (       test_type .EQ. 1 ) THEN 
    WRITE(*,*) '..... Interpolation from vel. to pre. grid .....'
    IF (     g%periodic .EQV. .TRUE. ) THEN
      CALL FDM_C0INT6P_LHS(  imax,  wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C0INTVP6P_RHS(imax,len, u, u_a)
      CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
      CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
    ElSEIF ( g%periodic .NEQV. .TRUE.) THEN
      CALL FDM_C0INTVP6_LHS( imaxp, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C0INTVP6_RHS( imax, imaxp, len, u, u_a)
      CALL TRIDFS(imaxp,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL TRIDSS(imaxp, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), u_a(1,1))
    ENDIF
! -------------------------------------------------------------------
  ELSEIF (   test_type .EQ. 2 ) THEN 
    WRITE(*,*) '..... Interpolation from pre. to vel. grid .....'
    IF (     g%periodic .EQV. .TRUE. ) THEN
      CALL FDM_C0INT6P_LHS(  imax,  wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C0INTPV6P_RHS(imax,len, u, u_a)
      CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
      CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
    ElSEIF ( g%periodic .NEQV. .TRUE.) THEN
      CALL FDM_C0INTPV6_LHS( imax, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C0INTPV6_RHS( imax, imaxp, len, u, u_a)
      CALL TRIDFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL TRIDSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), u_a(1,1))
    ENDIF
! -------------------------------------------------------------------
! Testing interpolatory 1st derivative
! -------------------------------------------------------------------
  ELSEIF (   test_type .EQ. 3 ) THEN 
    WRITE(*,*) '1st Interpol. deriv. from vel. to pre. grid ....'
    IF (     g%periodic .EQV. .TRUE. ) THEN
      CALL FDM_C1INT6P_LHS(  imax, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C1INTVP6P_RHS(imax,len, u, u_a)
      CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
      CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
      DO i = 1,imax
        DO l = 1,len
          u_int(l,i) = dudx_int(l,i) 
        ENDDO
      ENDDO
    ElSEIF ( g%periodic .NEQV. .TRUE. ) THEN
      CALL FDM_C1INTVP6_LHS( imaxp, g_pre%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C1INTVP6_RHS( imax, imaxp, len, u, u_a)
      CALL TRIDFS(imaxp,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL TRIDSS(imaxp, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), u_a(1,1))
      DO i = 1,imax
        DO l = 1,len
          u_int(l,i) = dudx_int(l,i) 
        ENDDO
      ENDDO
    ENDIF
! -------------------------------------------------------------------
  ELSEIF ( test_type .EQ. 4 ) THEN 
    WRITE(*,*) '1st Interpol. deriv. from pre. to vel. grid ....'
    IF (     g%periodic .EQV. .TRUE. ) THEN
      CALL FDM_C1INT6P_LHS(  imax, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C1INTPV6P_RHS(imax,len, u, u_a)
      CALL TRIDPFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
      CALL TRIDPSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), u_a(1,1),wrk2d(1))
      DO i = 1,imax
        DO l = 1,len
          u_int(l,i) = dudx_int(l,i) 
        ENDDO
      ENDDO
    ElSEIF ( g%periodic .NEQV. .TRUE. ) THEN
      CALL FDM_C1INTPV6_LHS( imax, g%jac, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL FDM_C1INTPV6_RHS( imax, imaxp, len, u, u_a)
      CALL TRIDFS(imax,      wrk1d(1,1),wrk1d(1,2),wrk1d(1,3))
      CALL TRIDSS(imax, len, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3), u_a(1,1))
      DO i = 1,imax
        DO l = 1,len
          u_int(l,i) = dudx_int(l,i) 
        ENDDO
      ENDDO
    ENDIF    
  ENDIF
! ###################################################################
! IO - Error and function values
  IF ( g%periodic .NEQV. .TRUE. ) THEN
    IF ( test_type .EQ. 1 .OR. test_type .EQ. 3 ) THEN
      u_a(:,imax)   = C_0_R
      u_int(:,imax) = C_0_R
      u(:,imax)     = C_0_R
    ENDIF
  ENDIF

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

  WRITE(*,2000) 'Solution L2-norm ...........:', SQRT(g%jac(1,1)*sol   / M_REAL(len))
  IF ( sol .EQ. C_0_R ) STOP
  WRITE(*,2000) 'Error L2-norm ..............:', SQRT(g%jac(1,1)*error / M_REAL(len)) 
  WRITE(*,2000) 'Error Linf-norm ............:', MAXVAL(u_c(1,1:imax))
  WRITE(*,2000) 'Relative error .............:', sqrt(error)/sqrt(sol)

  STOP

1000 FORMAT(6(1x,e16.10))
2000 FORMAT(1(1x,A,3x,e16.10))

END PROGRAM INTERPOL