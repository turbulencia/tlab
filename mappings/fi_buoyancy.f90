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
SUBROUTINE FI_BUOYANCY(buoyancy, nx,ny,nz, s, b, wrk1d)

  USE DNS_TYPES,  ONLY : term_dt
  USE DNS_GLOBAL, ONLY : inb_scal_array

  IMPLICIT NONE
  
  TYPE(term_dt),         INTENT(IN) :: buoyancy
  TINTEGER,                     INTENT(IN) :: nx,ny,nz
  TREAL, DIMENSION(nx,ny,nz,*), INTENT(IN) :: s
  TREAL, DIMENSION(nx,ny,nz),   INTENT(OUT):: b
  TREAL, DIMENSION(ny),         INTENT(IN) :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER j,k, is
  TREAL c0_loc,c1_loc,c2_loc,c3_loc, dummy

! #######################################################################
  SELECT CASE( buoyancy%type )

! -----------------------------------------------------------------------
  CASE( EQNS_BOD_HOMOGENEOUS )
     b = buoyancy%parameters(1)

! -----------------------------------------------------------------------
  CASE( EQNS_BOD_LINEAR )
     c1_loc = buoyancy%parameters(1); c2_loc = buoyancy%parameters(2); c3_loc = buoyancy%parameters(3) ! proportionality factors
     c0_loc = buoyancy%parameters(inb_scal_array+1)                                                    ! independent term

     IF      ( inb_scal_array .EQ. 1 .OR. &
             ( inb_scal_array .EQ. 2 .AND. ABS(buoyancy%parameters(2)) .LT. C_SMALL_R ) ) THEN ! avoid mem call for one common case is=2
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
           IF ( ABS(buoyancy%parameters(is)) .GT. C_SMALL_R ) b(:,:,:)= b(:,:,:) + buoyancy%parameters(is) *s(:,:,:,is)         
        ENDDO

     ENDIF

! -----------------------------------------------------------------------
  CASE ( EQNS_BOD_BILINEAR )
     c0_loc = buoyancy%parameters(1); c1_loc = buoyancy%parameters(2); c2_loc = buoyancy%parameters(3)

     DO k = 1,nz; DO j = 1,ny
        dummy = wrk1d(j)
        b(1:nx,j,k) = c0_loc*s(1:nx,j,k,1) + c1_loc*s(1:nx,j,k,2) + c2_loc*s(1:nx,j,k,1)*s(1:nx,j,k,2) - dummy
     ENDDO; ENDDO

! -----------------------------------------------------------------------
  CASE( EQNS_BOD_QUADRATIC )
     c0_loc =-buoyancy%parameters(1)/(buoyancy%parameters(2)/C_2_R)**2 
     c1_loc = buoyancy%parameters(2)

     DO k = 1,nz; DO j = 1,ny
        dummy = wrk1d(j)
        b(1:nx,j,k) = c0_loc*s(1:nx,j,k,1)*(s(1:nx,j,k,1)-c1_loc) - dummy
     ENDDO; ENDDO

  END SELECT

  RETURN
END SUBROUTINE FI_BUOYANCY

!########################################################################
!########################################################################
SUBROUTINE FI_BUOYANCY_SOURCE(buoyancy, isize_field, s, gradient, b_source)

  USE DNS_TYPES, ONLY : term_dt

  IMPLICIT NONE

  TYPE(term_dt),                   INTENT(IN)    :: buoyancy
  TINTEGER,                        INTENT(IN)    :: isize_field
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: s
  TREAL, DIMENSION(isize_field,*), INTENT(INOUT) :: gradient
  TREAL, DIMENSION(isize_field),   INTENT(OUT)   :: b_source

! -----------------------------------------------------------------------
  TREAL c0_loc

! #######################################################################
  SELECT CASE( buoyancy%type )

! -----------------------------------------------------------------------
  CASE( EQNS_BOD_HOMOGENEOUS )
     b_source = C_0_R

! -----------------------------------------------------------------------
  CASE( EQNS_BOD_LINEAR )
     b_source = C_0_R

! -----------------------------------------------------------------------
  CASE( EQNS_BOD_BILINEAR )
     b_source = C_0_R

! -----------------------------------------------------------------------
  CASE( EQNS_BOD_QUADRATIC )
     c0_loc =-buoyancy%parameters(1)/(buoyancy%parameters(2)/C_2_R)**2; c0_loc = c0_loc*C_2_R
     
     b_source = c0_loc* gradient(:,1)

  END SELECT

  RETURN
END SUBROUTINE FI_BUOYANCY_SOURCE
