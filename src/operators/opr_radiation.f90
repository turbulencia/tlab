#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2011/08/29 - A de Lozar
!#              Created
!# 2013/12/20 - J.P. Mellado
!#              Cleaned
!# 2016/01/05 - J.P. Mellado
!#              Further cleaned and vectorized
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the radiation term assuming the 1D approximation and
!# using compact schemes to calculate the integral term.
!# 
!########################################################################
SUBROUTINE OPR_RADIATION(radiation, nx,ny,nz, g, s, r, wrk1d,wrk3d)

  USE TLAB_TYPES, ONLY : term_dt, grid_dt

  IMPLICIT NONE

#include "integers.h"

  TYPE(term_dt),       INTENT(IN)           :: radiation
  TINTEGER,            INTENT(IN)           :: nx,ny,nz
  TYPE(grid_dt),       INTENT(IN)           :: g
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s        ! Radiatively active scalar
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: r        ! Radiative heating rate
  TREAL, DIMENSION(ny,*),     INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: wrk3d

  TARGET s, r, wrk3d
  
! -----------------------------------------------------------------------
  TINTEGER j, ip,ip2, nxy,nxz, ibc
  TREAL delta_inv, f0, f1
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst
  
! #######################################################################
  nxy = nx*ny ! For transposition to make y direction the last one
  ibc = 2     ! Boundary condition at the top for integral calulation 

! Prepare the pentadiagonal system
  CALL INT_C1N6_LHS(ny,    ibc,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
  CALL PENTADFS(ny-1,               wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))

  delta_inv = C_1_R /radiation%parameters(2)

! ###################################################################
  IF      ( radiation%type .EQ. EQNS_RAD_BULK1D_GLOBAL ) THEN
     nxz = 1

     IF ( nz .EQ. 1 ) THEN
        p_org => wrk3d
        p_dst => r
     ELSE
        p_org => r
        p_dst => wrk3d
     ENDIF
     
     CALL AVG1V2D_V(nx,ny,nz, i1, s, p_org, p_dst) ! Calculate averaged scalar into p_org; p_dst is auxiliar
     
! ###################################################################
  ELSE IF ( radiation%type .EQ. EQNS_RAD_BULK1D_LOCAL  ) THEN
     nxz = nx*nz

     IF ( nz .EQ. 1 ) THEN
        p_org => s
        p_dst => r
     ELSE 
        p_org => r
        p_dst => wrk3d

#ifdef USE_ESSL
        CALL DGETMO(s, nxy, nxy, nz, p_org, nz)
#else
        CALL DNS_TRANSPOSE(s, nxy, nz, nxy, p_org, nz)
#endif
     ENDIF

  ENDIF

! ###################################################################
! Calculate (negative) integral path.
  CALL INT_C1N6_RHS(ny,nxz, ibc, g%jac, p_org, p_dst)
  CALL PENTADSS(ny-1,nxz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), p_dst)
  ip = nxz *(ny-1) +1; p_dst(ip:ip+nxz-1) = C_0_R ! boundary condition

! Calculate radiative heating rate
  f0 = radiation%parameters(1)
  IF ( ABS(radiation%parameters(3)) .GT. C_0_R ) THEN
     f1 = radiation%parameters(3)
     DO j = ny,1,-1
        ip = nx*nz *(j-1) +1
        ip2= ip + nx*nz -1
        p_dst(ip:ip2) = p_org(ip:ip2) *( &
               EXP( p_dst(ip:ip2)                *delta_inv ) *f0 &
             + EXP((p_dst(1:nx*nz)-p_dst(ip:ip2))*delta_inv ) *f1 )
     ENDDO
  ELSE
     DO j = 1,ny*nxz
        p_dst(j) = p_org(j) *EXP( p_dst(j) *delta_inv ) *f0
     ENDDO
!  p_dst(1:ny*nxz) = radiation%parameters(1) *p_org(1:ny*nxz) *DEXP( p_dst(1:ny*nxz) *delta_inv ) seg-fault; need ulimit -u unlimited
  ENDIF
  
! ###################################################################
  IF      ( radiation%type .EQ. EQNS_RAD_BULK1D_GLOBAL ) THEN
     DO j = ny,1,-1
        ip = nx*nz *(j-1) +1; p_dst(ip:ip+nx*nz-1) = p_dst(j)
     ENDDO
  ENDIF
  
  IF ( nz .GT. 1 ) THEN ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
     CALL DGETMO(p_dst, nz, nz, nxy, r, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst, nz, nxy, nz, r, nxy)
#endif
  ENDIF
  
  NULLIFY(p_org,p_dst)

  RETURN
END SUBROUTINE OPR_RADIATION

!########################################################################
!########################################################################
SUBROUTINE OPR_RADIATION_FLUX(radiation, nx,ny,nz, g, s, r, wrk1d,wrk3d)

  USE TLAB_TYPES, ONLY : term_dt, grid_dt

  IMPLICIT NONE

#include "integers.h"

  TYPE(term_dt),              INTENT(IN)    :: radiation
  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TYPE(grid_dt),              INTENT(IN)    :: g
  TREAL, DIMENSION(nx*ny*nz), INTENT(IN)    :: s        ! Radiatively active scalar
  TREAL, DIMENSION(nx*ny*nz), INTENT(OUT)   :: r        ! Radiative flux
  TREAL, DIMENSION(ny,*),     INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(nx*ny*nz), INTENT(INOUT) :: wrk3d

  TARGET s, r, wrk3d
  
! -----------------------------------------------------------------------
  TINTEGER j, ip,ip2, nxy,nxz, ibc
  TREAL delta_inv, f0,f1
  TREAL, DIMENSION(:), POINTER :: p_org, p_dst
  
! #######################################################################
  nxy = nx*ny ! For transposition to make y direction the last one
  ibc = 2     ! Boundary condition at the top for integral calulation 

! Prepare the pentadiagonal system
  CALL INT_C1N6_LHS(ny,    ibc,     wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
  CALL PENTADFS(ny-1,               wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))

  delta_inv = C_1_R /radiation%parameters(2)

! ###################################################################
  IF      ( radiation%type .EQ. EQNS_RAD_BULK1D_GLOBAL ) THEN
     nxz = 1

     IF ( nz .EQ. 1 ) THEN
        p_org => wrk3d
        p_dst => r
     ELSE
        p_org => r
        p_dst => wrk3d
     ENDIF
     
     CALL AVG1V2D_V(nx,ny,nz, i1, s, p_org, p_dst) ! Calculate averaged scalar into p_org; p_dst is auxiliar
     
! ###################################################################
  ELSE IF ( radiation%type .EQ. EQNS_RAD_BULK1D_LOCAL  ) THEN
     nxz = nx*nz

     IF ( nz .EQ. 1 ) THEN
        p_org => s
        p_dst => r
     ELSE 
        p_org => r
        p_dst => wrk3d

#ifdef USE_ESSL
        CALL DGETMO(s, nxy, nxy, nz, p_org, nz)
#else
        CALL DNS_TRANSPOSE(s, nxy, nz, nxy, p_org, nz)
#endif
     ENDIF

  ENDIF

! ###################################################################
! Calculate (negative) integral path.
  CALL INT_C1N6_RHS(ny,nxz, ibc, g%jac, p_org, p_dst)
  CALL PENTADSS(ny-1,nxz, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5), p_dst)
  ip = nxz *(ny-1) +1; p_dst(ip:ip+nxz-1) = C_0_R ! boundary condition

! Calculate radiative flux
  f0 = -radiation%parameters(1) *radiation%parameters(2)
  IF ( ABS(radiation%parameters(3)) .GT. C_0_R ) THEN
     f1 = -radiation%parameters(3) *radiation%parameters(2)
     DO j = ny,1,-1
        ip = nx*nz *(j-1) +1
        ip2= ip + nx*nz -1
        p_dst(ip:ip2) = EXP( p_dst(ip:ip2)                *delta_inv ) *f0 &
                      - EXP((p_dst(1:nx*nz)-p_dst(ip:ip2))*delta_inv ) *f1
     ENDDO
  ELSE
     DO j = 1,ny*nxz
        p_dst(j) =           EXP( p_dst(j) *delta_inv ) *f0
     ENDDO
!  p_dst(1:ny*nxz) = radiation%parameters(1) *p_org(1:ny*nxz) *DEXP( p_dst(1:ny*nxz) *delta_inv ) seg-fault; need ulimit -u unlimited
  ENDIF
  
! ###################################################################
  IF      ( radiation%type .EQ. EQNS_RAD_BULK1D_GLOBAL ) THEN
     DO j = ny,1,-1
        ip = nx*nz *(j-1) +1; p_dst(ip:ip+nx*nz-1) = p_dst(j)
     ENDDO
  ENDIF
  
  IF ( nz .GT. 1 ) THEN ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
     CALL DGETMO(p_dst, nz, nz, nxy, r, nxy)
#else
     CALL DNS_TRANSPOSE(p_dst, nz, nxy, nz, r, nxy)
#endif
  ENDIF
  
  NULLIFY(p_org,p_dst)

  RETURN
END SUBROUTINE OPR_RADIATION_FLUX
