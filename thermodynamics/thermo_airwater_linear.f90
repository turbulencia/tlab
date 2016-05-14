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
SUBROUTINE THERMO_AIRWATER_LINEAR(nx,ny,nz, s, l)
  
  USE DNS_GLOBAL,    ONLY : inb_scal
  USE THERMO_GLOBAL, ONLY : thermo_param
  
  IMPLICIT NONE
  
#include "integers.h"
  
  TINTEGER,                            INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,inb_scal), INTENT(IN)    :: s     ! chi, psi
  TREAL, DIMENSION(nx*ny*nz),          INTENT(OUT)   :: l     ! normalized liquid
  
! -------------------------------------------------------------------
  TINTEGER ij
  TREAL dummy, dummy2

! ###################################################################
! Calculating \xi
  IF ( inb_scal .EQ. 1 ) THEN
     l = C_1_R + thermo_param(1)*s(:,1)
  ELSE
     l = C_1_R + thermo_param(1)*s(:,1) + thermo_param(2)*s(:,2)
  ENDIF

! Calculating liquid; \xi is overwritten in this routine
  IF ( ABS(thermo_param(inb_scal+1)) .LT. C_SMALL_R ) THEN
     DO ij = 1,nx*ny*nz
        l(ij) = MAX( l(ij), C_0_R )
     ENDDO

  ELSE
     dummy  = thermo_param(inb_scal+1)
     dummy2 = C_1_R /dummy
!     l = dummy *LOG( EXP(dummy2 *l) +C_1_R )
     DO ij = 1,nx*ny*nz
        l(ij) = dummy *LOG( EXP(dummy2 *l(ij)) +C_1_R )
     ENDDO

  ENDIF
  
  RETURN
END SUBROUTINE THERMO_AIRWATER_LINEAR

!########################################################################
!########################################################################
SUBROUTINE THERMO_AIRWATER_LINEAR_SOURCE(nx,ny,nz, s, xi,der1,der2)
  
  USE DNS_GLOBAL, ONLY : inb_scal
  USE THERMO_GLOBAL, ONLY : thermo_param
  
  IMPLICIT NONE
  
#include "integers.h"
  
  TINTEGER,                            INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,inb_scal), INTENT(IN)  :: s          ! chi, psi
  TREAL, DIMENSION(nx*ny*nz),          INTENT(OUT) :: xi,der1,der2
  
! -------------------------------------------------------------------
  TINTEGER ij
  TREAL dummy

! ###################################################################
   IF ( inb_scal .EQ. 1 ) THEN
     xi(:) = C_1_R + thermo_param(1)*s(:,1)
  ELSE
     xi(:) = C_1_R + thermo_param(1)*s(:,1) + thermo_param(2)*s(:,2)
  ENDIF
  
  IF ( ABS(thermo_param(inb_scal+1)) .LT. C_SMALL_R ) THEN
     der2 = C_BIG_R

     DO ij = 1,nx*ny*nz
        IF ( xi(ij) .LE. C_0_R ) THEN; der1(ij) = C_0_R;
        ELSE;                          der1(ij) = C_1_R; ENDIF
     ENDDO

  ELSE
! Formulation in terms of the exponential might lead to NAN because of too large numbers!
!     dummy  =-C_1_R/thermo_param(inb_scal+1)
!     der1 = C_1_R / ( C_1_R + EXP(dummy *xi) ) 
     
     dummy= C_05_R /thermo_param(inb_scal+1)
     der1 = C_05_R *( C_1_R + TANH( dummy *xi ) )
     
     der2 = (der1-C_1_R) *der1 *dummy

  ENDIF
  
  RETURN
END SUBROUTINE THERMO_AIRWATER_LINEAR_SOURCE
