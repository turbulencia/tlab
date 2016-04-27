#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/06/18 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Determine the buoyancy term (in particular, it gives the density difference \rho-\rho_0)
!# to be added the the RHS of momentum equation when it is a function of a scalar Z.
!# A reference profile is passed through array wrk1d
!#
!########################################################################
!# ARGUMENTS 
!#
!# In the piecewise-linear case (buoyancy reversal)
!# param1    In    b2-b1
!# param2    In    bm-b1
!# param3    In    Sm
!# param4    In    delta
!#
!########################################################################
SUBROUTINE FI_BUOYANCY(ibodyforce, nx,ny,nz, param, s, b, wrk1d)

  USE DNS_GLOBAL, ONLY : inb_scal_array

  IMPLICIT NONE
  
  TINTEGER,                     INTENT(IN) :: ibodyforce, nx, ny, nz
  TREAL, DIMENSION(*),          INTENT(IN) :: param
  TREAL, DIMENSION(nx,ny,nz,*), INTENT(IN) :: s
  TREAL, DIMENSION(nx,ny,nz),   INTENT(OUT):: b
  TREAL, DIMENSION(ny),         INTENT(IN) :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER i,j,k, is
  TREAL c0_loc,c1_loc,c2_loc,c3_loc, dummy
  TREAL delta, delta_inv

! #######################################################################
! Constant homogeneous value
! #######################################################################
  IF      ( ibodyforce .EQ. EQNS_BOD_HOMOGENEOUS  ) THEN
     b = param(1)

! #######################################################################
! Linear relation
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_LINEAR       ) THEN
     c1_loc = param(1); c2_loc = param(2); c3_loc = param(3) ! proportionality factors
     c0_loc = param(inb_scal_array+1)                        ! independent term

     IF      ( inb_scal_array .EQ. 1 .OR. &
             ( inb_scal_array .EQ. 2 .AND. ABS(param(2)) .LT. C_SMALL_R ) ) THEN ! avoid mem call for one common case is=2
        DO k = 1,nz; DO j = 1,ny
           dummy = wrk1d(j) - c0_loc
           b(1:nx,j,k) = c1_loc *s(1:nx,j,k,1) - dummy
        ENDDO; ENDDO
        
     ELSE IF ( inb_scal_array .EQ. 2 ) THEN
        DO k = 1,nz; DO j = 1,ny
           dummy = wrk1d(j) - c0_loc
           b(1:nx,j,k) = c1_loc *s(1:nx,j,k,1) + c2_loc *s(1:nx,j,k,2) - dummy
        ENDDO; ENDDO
        
     ELSE IF ( inb_scal_array .EQ. 3 ) THEN
        DO k = 1,nz; DO j = 1,ny
           dummy = wrk1d(j) - c0_loc
           b(1:nx,j,k) = c1_loc *s(1:nx,j,k,1) + c2_loc *s(1:nx,j,k,2) + c3_loc *s(1:nx,j,k,3) - dummy
        ENDDO; ENDDO

     ELSE ! general
        DO k = 1,nz; DO j = 1,ny
           dummy = wrk1d(j)
           b(1:nx,j,k) = c0_loc -dummy
        ENDDO; ENDDO

        DO is = 1,inb_scal_array
           IF ( ABS(param(is)) .GT. C_SMALL_R ) b(:,:,:)= b(:,:,:) + param(is) *s(:,:,:,is)         
        ENDDO

     ENDIF

! #######################################################################
! Bilinear
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_BILINEAR     ) THEN
     c0_loc = param(1); c1_loc = param(2); c2_loc = param(3)

     DO k = 1,nz; DO j = 1,ny
        dummy = wrk1d(j)
        b(1:nx,j,k) = c0_loc*s(1:nx,j,k,1) + c1_loc*s(1:nx,j,k,2) + c2_loc*s(1:nx,j,k,1)*s(1:nx,j,k,2) - dummy
     ENDDO; ENDDO

! #######################################################################
! Quadratic relation
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_QUADRATIC    ) THEN
     c0_loc =-param(1)/(param(2)/C_2_R)**2 
     c1_loc = param(2)

     DO k = 1,nz; DO j = 1,ny
        dummy = wrk1d(j)
        b(1:nx,j,k) = c0_loc*s(1:nx,j,k,1)*(s(1:nx,j,k,1)-c1_loc) - dummy
     ENDDO; ENDDO

  ENDIF

  RETURN
END SUBROUTINE FI_BUOYANCY

!########################################################################
!########################################################################
SUBROUTINE FI_BUOYANCY_SOURCE(ibodyforce, isize_field, param, s, gradient, b_source)

  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: ibodyforce, isize_field
  TREAL, DIMENSION(*),             INTENT(IN)    :: param
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: s
  TREAL, DIMENSION(isize_field,*), INTENT(INOUT) :: gradient
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: b_source

! -----------------------------------------------------------------------
  TINTEGER ij
  TREAL c0_loc, c1_loc, c2_loc, c3_loc, delta, delta_inv

! #######################################################################
! Constant homogeneous value
! #######################################################################
  IF      ( ibodyforce .EQ. EQNS_BOD_HOMOGENEOUS  ) THEN
     b_source = C_0_R

! #######################################################################
! Linear relation
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_LINEAR       ) THEN
     b_source = C_0_R

! #######################################################################
! Bilinear (2 scalars)
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_BILINEAR     ) THEN
     b_source = C_0_R

! #######################################################################
! Quadratic relation
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_QUADRATIC    ) THEN
     c0_loc =-param(1)/(param(2)/C_2_R)**2; c0_loc = c0_loc*C_2_R
     
     b_source = c0_loc* gradient(:,1)

  ENDIF

  RETURN
END SUBROUTINE FI_BUOYANCY_SOURCE
