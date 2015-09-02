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
!# 2013/06/12 - A. de Lozar
!#              Adding FI_LIQUIDWATER
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

  IMPLICIT NONE
  
  TINTEGER,                     INTENT(IN) :: ibodyforce, nx, ny, nz
  TREAL, DIMENSION(*),          INTENT(IN) :: param
  TREAL, DIMENSION(nx,ny,nz,*), INTENT(IN) :: s
  TREAL, DIMENSION(nx,ny,nz),   INTENT(OUT):: b
  TREAL, DIMENSION(ny),         INTENT(IN) :: wrk1d

! -----------------------------------------------------------------------
  TINTEGER i, j, k
  TREAL c0_loc, c1_loc, c2_loc, c3_loc, delta, delta_inv, dummy

! #######################################################################
! Constant homogeneous value
! #######################################################################
  IF      ( ibodyforce .EQ. EQNS_BOD_HOMOGENEOUS  ) THEN
     b = param(1)

! #######################################################################
! Linear relation
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_LINEAR       ) THEN
     c0_loc = param(1)

     DO k = 1,nz; DO j = 1,ny
        dummy = wrk1d(j)
        b(1:nx,j,k) = c0_loc*s(1:nx,j,k,1) - dummy
     ENDDO; ENDDO

! #######################################################################
! Bilinear (2 scalars)
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

! #######################################################################
! Piecewise linear (buoyancy reversal)
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_PIECEWISE_LINEAR     ) THEN
     c0_loc = param(2)/param(3)
     c2_loc = (param(1)-param(2))/(C_1_R-param(3))
     delta  = param(4)

     IF ( delta .EQ. C_0_R ) THEN
        c1_loc = param(1) - c2_loc
        IF ( c0_loc .GT. 0 ) THEN
           DO k = 1,nz; DO j = 1,ny
              DO i = 1,nx
                 b(i,j,k) = MIN(c0_loc*s(i,j,k,1),c1_loc+s(i,j,k,1)*c2_loc) - wrk1d(j)
              ENDDO
           ENDDO; ENDDO
        ELSE
           DO k = 1,nz; DO j = 1,ny
              DO i = 1,nx
                 b(i,j,k) = MAX(c0_loc*s(i,j,k,1),c1_loc+s(i,j,k,1)*c2_loc) - wrk1d(j)
              ENDDO
           ENDDO; ENDDO
        ENDIF

     ELSE
        delta_inv = C_1_R/delta
        c2_loc    = delta*(c2_loc-c0_loc)
        c1_loc    = param(3)
        DO k = 1,nz; DO j = 1,ny
           dummy = wrk1d(j)
           b(1:nx,j,k) = c0_loc*s(1:nx,j,k,1) + &
                c2_loc*LOG(C_1_R+EXP((s(1:nx,j,k,1)-c1_loc)*delta_inv)) - dummy
        ENDDO; ENDDO

     ENDIF

! #######################################################################
! Piecewise bilinear (buoyancy reversal)
! In this case of the mixture the first scale stores the mixed fraction , the 
! second scalar the deviations from the mixed fraction and the third scalar stores 
! the liquid water (calculated in FI_LIQUIDWATER, possibly smoothed)
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_PIECEWISE_BILINEAR  ) THEN
     c0_loc = (param(1)-param(2))         /( C_1_R-param(3) )
     c2_loc = (param(1)*param(3)-param(2))/( C_1_R-param(3) )
     c3_loc =  param(5)

     DO k = 1,nz; DO j = 1,ny
        DO i = 1,nx
           b(i,j,k) = c0_loc*s(i,j,k,1) + c3_loc*s(i,j,k,2) + c2_loc*( s(i,j,k,3) -C_1_R )  - wrk1d(j)
        ENDDO
     ENDDO; ENDDO

! #######################################################################
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_PIECEWISE_TRILINEAR ) THEN
     c0_loc = (param(1)-param(2))         /( C_1_R-param(3) )
     c2_loc = (param(1)*param(3)-param(2))/( C_1_R-param(3) )
     c3_loc =  param(5)

     DO k = 1,nz; DO j = 1,ny
        DO i = 1,nx
           b(i,j,k) = c0_loc*s(i,j,k,1) + c3_loc*s(i,j,k,2) + c2_loc*( s(i,j,k,4)  -C_1_R ) + s(i,j,k,3)  - wrk1d(j)
        ENDDO
     ENDDO; ENDDO

  ENDIF

  RETURN
END SUBROUTINE FI_BUOYANCY

!########################################################################
!########################################################################
SUBROUTINE FI_LIQUIDWATER(ibodyforce, nx,ny,nz, param, s,sout)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: ibodyforce, nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)  :: param
  TREAL, DIMENSION(nx*ny*nz,2), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: sout
!--------------------------------------------------------------------------
  TINTEGER ij
  TREAL c1_loc, c2_loc, c4_loc, delta, delta_inv, dummy

! #######################################################################
! Piecewise linear (buoyancy reversal)
! #######################################################################
  IF       ( ibodyforce .EQ. EQNS_BOD_PIECEWISE_LINEAR   ) THEN
     c1_loc = param(3)
     c2_loc = param(4)/param(3)
     delta  = param(4)
     
     IF ( delta .GT. C_SMALL_R ) THEN
        delta_inv = C_1_R/delta
        DO ij = 1,nx*ny*nz
           sout(ij) = c2_loc*LOG( C_1_R + EXP( ( c1_loc - s(ij,1) )*delta_inv) )
        ENDDO
     ELSE
        DO ij = 1,nx*ny*nz
           sout(ij) = MAX(c1_loc - s(ij,1),C_0_R)
        ENDDO
     ENDIF

! #######################################################################
! Piecewise bilinear (buoyancy reversal)
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_PIECEWISE_BILINEAR .OR. ibodyforce .EQ. EQNS_BOD_PIECEWISE_TRILINEAR ) THEN
     c1_loc = param(3)
     c2_loc = param(4)/param(3)
     delta  = param(4)

     dummy = ( param(1)*param(3)-param(2) ) /( C_1_R-param(3) )/param(3)
     IF ( dummy .GT. C_SMALL_R ) THEN; c4_loc = param(5)*param(6)/dummy
     ELSE;                             c4_loc = C_0_R; ENDIF ! Radiation only
     
     IF ( delta .GT. C_SMALL_R ) THEN
        delta_inv = C_1_R/delta
        DO ij = 1,nx*ny*nz
           sout(ij) = c2_loc*LOG( C_1_R + EXP( ( c1_loc - s(ij,1) - c4_loc*s(ij,2) )*delta_inv) )
        ENDDO
     ELSE
        DO ij = 1,nx*ny*nz
           sout(ij) = MAX(c1_loc - s(ij,1) - c4_loc*s(ij,2),C_0_R) /c1_loc
        ENDDO
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE FI_LIQUIDWATER

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

! #######################################################################
! Piecewise linear (buoyancy reversal)
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_PIECEWISE_LINEAR ) THEN
     c1_loc = param(3)
     c2_loc = (param(1)*param(3)-param(2))/( C_1_R-param(3) )/param(3)
     delta  = param(4)
     
     IF ( delta .EQ. C_0_R ) THEN
        b_source = C_BIG_R
        
     ELSE
        delta_inv = C_05_R/delta
        c2_loc    = c2_loc*C_05_R*delta_inv

        DO ij = 1,isize_field
           b_source(ij) =-gradient(ij,1)*c2_loc / (COSH((s(ij,1)-c1_loc)*delta_inv)**2)
        ENDDO

     ENDIF

! #######################################################################
! Piecewise bilinear (buoyancy reversal)
! #######################################################################
  ELSE IF ( ibodyforce .EQ. EQNS_BOD_PIECEWISE_BILINEAR .OR. ibodyforce .EQ. EQNS_BOD_PIECEWISE_TRILINEAR ) THEN
     c1_loc = param(3)
     c2_loc = (param(1)*param(3)-param(2))/( C_1_R-param(3) )/param(3)
     c3_loc = param(5)*param(6)
     delta  = param(4)
     
     IF ( delta .EQ. C_0_R ) THEN
        b_source = C_BIG_R
        
     ELSE
        delta_inv = C_05_R/delta
        c2_loc    = c2_loc*C_05_R*delta_inv

        DO ij = 1,isize_field
           b_source(ij) =-gradient(ij,1)*c2_loc / (COSH(s(ij,1)*delta_inv)**2) ! s is scalar xi
        ENDDO 

        IF ( c3_loc .GT. C_0_R ) THEN ! Couplig field of radiation and evaporative cooling active
           c3_loc    = -c3_loc
           delta_inv = -C_1_R/delta
           
           DO ij = 1,isize_field
              gradient(ij,1) = c3_loc/(C_1_R + EXP(s(ij,1)*delta_inv))
           ENDDO
        ELSE
           gradient(:,1) = C_0_R
        ENDIF

     ENDIF

  ENDIF

  RETURN
END SUBROUTINE FI_BUOYANCY_SOURCE
