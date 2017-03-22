#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2017/03/22 - J.P. Mellado
!#              Extracted from airwater_pool
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate thermodynamic properties from h, p and composition in the
!# incompressible formulation, when the thermodynamic pressure is
!# a given profile
!#
!########################################################################

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_TEMPERATURE(nx,ny,nz, s, e, T)

  USE THERMO_GLOBAL, ONLY : imixture, THERMO_AI

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: T

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL E_LOC
  
! ###################################################################
  IF      ( imixture .EQ. 0 ) THEN
! s1 is specific static energy
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T = s(ij,1) - E_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
! s1 is specific static energy
! s2 is total water specific humidity = water vapor specific humidity
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T = (s(ij,1) - E_LOC ) / &
                ( (C_1_R-s(ij,2))*THERMO_AI(1,1,2) + s(ij,2)*THERMO_AI(1,1,1) )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! s1 is specific static energy
! s2 is total water specific humidity
! s3 is liquid water specific humidity
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T = (s(ij,1) - E_LOC - s(ij,3)*THERMO_AI(6,1,3) ) / &
                ( (C_1_R-s(ij,2))*THERMO_AI(1,1,2) + (s(ij,2)-s(ij,3))*THERMO_AI(1,1,1) &
                + s(ij,3)* THERMO_AI(1,1,3) )
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_TEMPERATURE

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_DENSITY(nx,ny,nz, s, e,p, rho)

  USE THERMO_GLOBAL, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL WMEAN_INV, P_LOC, E_LOC, T, dummy
  
! ###################################################################
  IF      ( imixture .EQ. 0 ) THEN
! s1 is specific static energy
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T = s(ij,1) - E_LOC
           rho(ij) = P_LOC /T
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
! s1 is specific static energy
! s2 is total water specific humidity = water vapor specific humidity
     dummy = WGHT_INV(1)-WGHT_INV(2)
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T = (s(ij,1) - E_LOC ) / &
                ( (C_1_R-s(ij,2))*THERMO_AI(1,1,2) + s(ij,2)*THERMO_AI(1,1,1) )
           WMEAN_INV = WGHT_INV(2) + s(ij,2)*dummy
           rho(ij) = P_LOC /(WMEAN_INV*T)
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! s1 is specific static energy
! s2 is total water specific humidity
! s3 is liquid water specific humidity
     dummy = WGHT_INV(1)-WGHT_INV(2)
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T = (s(ij,1) - E_LOC - s(ij,3)*THERMO_AI(6,1,3) ) / &
                ( (C_1_R-s(ij,2))*THERMO_AI(1,1,2) + (s(ij,2)-s(ij,3))*THERMO_AI(1,1,1) &
                + s(ij,3)* THERMO_AI(1,1,3) )
           WMEAN_INV = WGHT_INV(2) + s(ij,2)*dummy - s(ij,3)*WGHT_INV(1)
           rho(ij) = P_LOC /(WMEAN_INV*T)
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_DENSITY

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_BUOYANCY(nx,ny,nz, s, e,p,r, b)

  USE THERMO_GLOBAL, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p,r
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: b

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL WMEAN_INV, P_LOC, E_LOC, R_LOC, R_LOC_INV, T, dummy
  
! ###################################################################
  IF      ( imixture .EQ. 0 ) THEN
! s1 is specific static energy
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T = s(ij,1) - E_LOC
           b(ij) = R_LOC_INV *( P_LOC /T -R_LOC )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
! s1 is specific static energy
! s2 is total water specific humidity = water vapor specific humidity
     dummy = WGHT_INV(1)-WGHT_INV(2)
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T = (s(ij,1) - E_LOC ) / &
                ( (C_1_R-s(ij,2))*THERMO_AI(1,1,2) + s(ij,2)*THERMO_AI(1,1,1) )
           WMEAN_INV = WGHT_INV(2) + s(ij,2)*dummy
           b(ij) = R_LOC_INV *( P_LOC /(WMEAN_INV*T) -R_LOC )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
! s1 is specific static energy
! s2 is total water specific humidity
! s3 is liquid water specific humidity
     dummy = WGHT_INV(1)-WGHT_INV(2)
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T = (s(ij,1) - E_LOC - s(ij,3)*THERMO_AI(6,1,3) ) / &
                ( (C_1_R-s(ij,2))*THERMO_AI(1,1,2) + (s(ij,2)-s(ij,3))*THERMO_AI(1,1,1) &
                + s(ij,3)* THERMO_AI(1,1,3) )
           WMEAN_INV = WGHT_INV(2) + s(ij,2)*dummy - s(ij,3)*WGHT_INV(1)
           b(ij) = R_LOC_INV *( P_LOC /(WMEAN_INV*T) -R_LOC )
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_BUOYANCY
