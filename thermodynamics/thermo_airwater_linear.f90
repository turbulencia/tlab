#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2016/02/16 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculating the equilibrium liquid from the mixture fraction and the
!# enthalpy deviation according to the linearized equilibrium thermodynamics
!#
!########################################################################
SUBROUTINE THERMO_AIRWATER_LINEAR(nx,ny,nz, s, l, wrk3d)
  
  USE DNS_GLOBAL,    ONLY : inb_scal
  USE THERMO_GLOBAL, ONLY : thermo_param
  
  IMPLICIT NONE
  
#include "integers.h"
  
  TINTEGER,                            INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,inb_scal), INTENT(IN)    :: s     ! chi, psi
  TREAL, DIMENSION(nx*ny*nz),          INTENT(OUT)   :: l     ! normalized liquid
  TREAL, DIMENSION(nx*ny*nz),          INTENT(INOUT) :: wrk3d ! xi
  
! -------------------------------------------------------------------
  TINTEGER ij
  TREAL dummy, dummy2

! ###################################################################
  IF ( inb_scal .EQ. 1 ) THEN
     wrk3d = C_1_R + thermo_param(1)*s(:,1)
  ELSE
     wrk3d = C_1_R + thermo_param(1)*s(:,1) + thermo_param(2)*s(:,2)
  ENDIF

  IF ( ABS(thermo_param(inb_scal+1)) .LT. C_SMALL_R ) THEN
     DO ij = 1,nx*ny*nz
        l(ij) = MAX( wrk3d(ij), C_0_R )
     ENDDO

  ELSE
     dummy  = thermo_param(inb_scal+1)
     dummy2 = C_1_R /dummy
!     l = dummy *LOG( EXP(dummy2 *wrk3d) +C_1_R )
     DO ij = 1,nx*ny*nz
        l(ij) = dummy *LOG( EXP(dummy2 *wrk3d(ij)) +C_1_R )
     ENDDO

  ENDIF
  
  RETURN
END SUBROUTINE THERMO_AIRWATER_LINEAR

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_LINEAR_SOURCE(nx,ny,nz, s, der1,der2, wrk3d)
  
  USE DNS_GLOBAL, ONLY : inb_scal
  USE THERMO_GLOBAL, ONLY : thermo_param
  
  IMPLICIT NONE
  
#include "integers.h"
  
  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)    :: s          ! chi, psi
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT)   :: der1,der2
  TREAL, DIMENSION(nx*ny*nz),   INTENT(INOUT) :: wrk3d      ! xi
  
! -------------------------------------------------------------------
  TINTEGER ij
  TREAL dummy

! ###################################################################
   IF ( inb_scal .EQ. 1 ) THEN
     wrk3d = C_1_R + thermo_param(1)*s(:,1)
  ELSE
     wrk3d = C_1_R + thermo_param(1)*s(:,1) + thermo_param(2)*s(:,2)
  ENDIF
  
  IF ( ABS(thermo_param(inb_scal+1)) .LT. C_SMALL_R ) THEN
     der2 = C_BIG_R

     DO ij = 1,nx*ny*nz
        IF ( wrk3d(ij) .LE. C_0_R ) THEN; der1(ij) = C_0_R;
        ELSE;                             der1(ij) = C_1_R; ENDIF
     ENDDO

  ELSE
     dummy  =-C_1_R/thermo_param(inb_scal+1)

     der1 = C_1_R / ( C_1_R + EXP(dummy *wrk3d) )

     der2 = (der1-C_1_R) *der1 *dummy

  ENDIF
  
  RETURN
END SUBROUTINE THERMO_AIRWATER_LINEAR_SOURCE
