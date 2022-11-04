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
!# s1 is specific static energy
!# s2 is total water specific humidity
!# s3 is liquid water specific humidity, if any
!#
!########################################################################

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_TEMPERATURE(nx,ny,nz, s, e, T)

  USE THERMO_VARS, ONLY : imixture, THERMO_AI

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: T

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL E_LOC

  TREAL Cd, Cdv, Lv0, Cvl

! ###################################################################
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = s(ij,1) - E_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )

        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_TEMPERATURE

!########################################################################
! Calculating h_l - h; very similar to the temperature routine
!########################################################################
SUBROUTINE THERMO_ANELASTIC_STATIC_L(nx,ny,nz, s, e, result)

  USE THERMO_VARS, ONLY : imixture, THERMO_AI

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: result

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL E_LOC

  TREAL Cd, Cdv, Lv0, Cvl, Cl

! ###################################################################
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)
  Cl = THERMO_AI(1,1,3)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           result(ij) = s(ij,1) - E_LOC
           
           result(ij) = Cl *result(ij) +E_LOC -Lv0 -s(ij,1)

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           result(ij) = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           
           result(ij) = Cl *result(ij) +E_LOC -Lv0 -s(ij,1)

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           result(ij) = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )

           result(ij) = Cl *result(ij) +E_LOC -Lv0 -s(ij,1)
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_STATIC_L

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_DENSITY(nx,ny,nz, s, e,p, rho)

  USE THERMO_VARS, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rho

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL WMEAN_INV, P_LOC, E_LOC, T_LOC
  
  TREAL Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl

! ###################################################################
  Rv = WGHT_INV(1)
  Rd = WGHT_INV(2)
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           rho(ij) = P_LOC /T_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           WMEAN_INV = Rd + s(ij,2) *Rdv
           rho(ij) = P_LOC /( WMEAN_INV *T_LOC )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           WMEAN_INV = Rd + s(ij,2) *Rdv - s(ij,3) *Rv
           rho(ij) = P_LOC /(WMEAN_INV*T_LOC)
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_DENSITY

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_BUOYANCY(nx,ny,nz, s, e,p,r, b)

  USE THERMO_VARS, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p,r
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: b

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL WMEAN_INV, P_LOC, E_LOC, R_LOC, R_LOC_INV, T_LOC
  
  TREAL Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl

! ###################################################################
  Rv = WGHT_INV(1)
  Rd = WGHT_INV(2)
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           b(ij) = R_LOC_INV *(R_LOC -P_LOC /T_LOC )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           WMEAN_INV = Rd + s(ij,2) *Rdv
           b(ij) = R_LOC_INV *( R_LOC -P_LOC /(WMEAN_INV*T_LOC) )
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        R_LOC_INV = C_1_R /R_LOC
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3) *Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           WMEAN_INV = Rd + s(ij,2) *Rdv - s(ij,3) *Rv
           b(ij) = R_LOC_INV *( R_LOC -P_LOC /(WMEAN_INV*T_LOC) )
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_BUOYANCY

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_QVEQU(nx,ny,nz, s, e,p, T,qvequ)

  USE THERMO_VARS, ONLY : imixture, THERMO_AI, WGHT_INV, MRATIO, THERMO_PSAT, NPSAT

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p 
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: qvequ, T
  
! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is, ipsat
  TREAL psat, E_LOC, P_LOC
  
  TREAL Cd, Cdv, Lv0, Cvl, rd_ov_rv
  
! ###################################################################
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  rd_ov_rv = WGHT_INV(2) /WGHT_INV(1)
  
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO*p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           
           psat = C_0_R
           DO ipsat = NPSAT,1,-1
              psat = psat*T(ij) + THERMO_PSAT(ipsat)
           ENDDO
           qvequ(ij) = rd_ov_rv *( C_1_R -s(ij,2) ) /( P_LOC/psat -C_1_R )
           
        ENDDO

     ENDDO
  ENDIF
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_QVEQU

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_RELATIVEHUMIDITY(nx,ny,nz, s, e,p, T,rh)

  USE THERMO_VARS, ONLY : imixture, THERMO_AI, THERMO_PSAT, NPSAT, WGHT_INV, MRATIO

  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p 
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: rh, T
  
! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is, ipsat
  TREAL psat, E_LOC, P_LOC
  
  TREAL Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl
  
! ###################################################################
  Rv = WGHT_INV(1)
  Rd = WGHT_INV(2)
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO*p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T(ij) = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           
           psat = C_0_R
           DO ipsat = NPSAT,1,-1
              psat = psat*T(ij) + THERMO_PSAT(ipsat)
           ENDDO
           rh(ij) = (s(ij,2) -s(ij,3)) *P_LOC/psat *Rv /( Rd + s(ij,2) *Rdv - s(ij,3) *Rv )
           rh(ij) = rh(ij) *C_100_R
           
        ENDDO

     ENDDO
  ENDIF
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_RELATIVEHUMIDITY

!########################################################################
! Potential temperature
!########################################################################
SUBROUTINE THERMO_ANELASTIC_THETA(nx,ny,nz, s, e,p, theta)

  USE THERMO_VARS, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO, GRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: theta

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL P_LOC, E_LOC, T_LOC

  TREAL Rd, Cd, Cdv, Lv0, Cvl, kappa
  
! ###################################################################
  Rd = WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  kappa = Rd *GRATIO /Cd
          
  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           
           theta(ij) = T_LOC / P_LOC**kappa

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC )  / ( Cd + s(ij,2) *Cdv )

           theta(ij) = T_LOC / P_LOC**kappa

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )

           theta(ij) = T_LOC / P_LOC**kappa

        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_THETA

!########################################################################
! Virtual Potential temperature
!########################################################################
SUBROUTINE THERMO_ANELASTIC_THETA_V(nx,ny,nz, s, e,p, theta)

  USE THERMO_VARS, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO, GRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: theta

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL P_LOC, E_LOC, T_LOC

  TREAL Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl, kappa
  
! ###################################################################
  Rv = WGHT_INV(1)
  Rd = WGHT_INV(2)
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  kappa = Rd *GRATIO /Cd
          
  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           
           theta(ij) = T_LOC / P_LOC**kappa

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC )  / ( Cd + s(ij,2) *Cdv )

           theta(ij) = T_LOC *( Rd + s(ij,2) *Rdv )/ P_LOC**kappa

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )

           theta(ij) = T_LOC *( Rd + s(ij,2) *Rdv - s(ij,3) *Rv )/ P_LOC**kappa

        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_THETA_V

!########################################################################
! Liquid water potential temperature
! Gas constants are multiplied by GRATIO because they always enter as ratios wrt Cps
!########################################################################
SUBROUTINE THERMO_ANELASTIC_THETA_L(nx,ny,nz, s, e,p, theta)

  USE THERMO_VARS, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO, GRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: theta

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL P_LOC, E_LOC, T_LOC

  TREAL Rd, Rdv, Rl, Cd, Cdv, Cl, Lv0, Cvl, Lv
  
! ###################################################################
  Rd = WGHT_INV(2)               *GRATIO
  Rdv=(WGHT_INV(1) - WGHT_INV(2))*GRATIO
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           
           Rl = Rd /Cd
          
           theta(ij) = T_LOC / P_LOC**Rl

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC )  / ( Cd + s(ij,2) *Cdv )
           Cl = C_1_R/( Cd  + s(ij,2) *Cdv )
           Rl =(Rd  + s(ij,2) *Rdv) *Cl
           theta(ij) = T_LOC / P_LOC**Rl

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           Cl = C_1_R/( Cd  + s(ij,2) *Cdv )
           Rl =(Rd  + s(ij,2) *Rdv) *Cl
           Lv = Lv0 - T_LOC *Cvl
           theta(ij) = T_LOC / P_LOC**Rl *EXP(-Lv *s(ij,3) *Cl / T_LOC )

        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_THETA_L

!########################################################################
! Equivalent potential temperature
! Gas constants are multiplied by GRATIO because they always enter as ratios wrt Cps
!########################################################################
SUBROUTINE THERMO_ANELASTIC_THETA_E(nx,ny,nz, s, e,p, theta)

  USE THERMO_VARS, ONLY : imixture, WGHT_INV, THERMO_AI, MRATIO, GRATIO

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: theta

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL P_LOC, E_LOC, T_LOC

  TREAL Rd, Rdv, Cd, Cdv, Cdl, Lv0, Cvl, Lv
  TREAL Re, Ce
  
! ###################################################################
  Rd = WGHT_INV(2)               *GRATIO
  Rdv=(WGHT_INV(1) - WGHT_INV(2))*GRATIO
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Cdl= THERMO_AI(1,1,3) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           
           Re = Rd /Cd
          
           theta(ij) = T_LOC / P_LOC**Re

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv )
           Ce = C_1_R/( Cd  + s(ij,2) *Cdl )
           Re = Rd  *( C_1_R -s(ij,2) )*Ce
           Lv = Lv0 - T_LOC *Cvl
           theta(ij) = T_LOC / P_LOC**Re *EXP( Lv *s(ij,2) *Ce /T_LOC )

        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           Ce = C_1_R/( Cd  + s(ij,2) *Cdl )
           Re = Rd  *( C_1_R -s(ij,2) )*Ce
           Lv = Lv0 - T_LOC *Cvl
!           Lv = Lv0
           theta(ij) = T_LOC / P_LOC**Re *EXP( Lv *(s(ij,2)-s(ij,3)) *Ce /T_LOC )

        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_THETA_E

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_LAPSE_FR(nx,ny,nz, s, dTdy, e, lapse, frequency)

  USE THERMO_VARS, ONLY : imixture, THERMO_AI, GRATIO, scaleheight
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: dTdy
  TREAL, DIMENSION(*),          INTENT(IN)  :: e
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: lapse, frequency

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL E_LOC, T_LOC

  TREAL Cd, Cdv, Cvl, Lv0
  
! ###################################################################
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)
  Lv0=-THERMO_AI(6,1,3)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     lapse = GRATIO /scaleheight

     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           frequency(ij) = ( lapse(ij) + dTdy(ij) )/ T_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN 
     lapse = GRATIO /scaleheight /( Cd  + s(:,2) *Cdv )

     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           frequency(ij) = ( lapse(ij) + dTdy(ij) )/ T_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     lapse = GRATIO /scaleheight /( Cd + s(:,2) *Cdv + s(:,3) *Cvl )
     
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )
           frequency(ij) = ( lapse(ij) + dTdy(ij) )/ T_LOC
           
        ENDDO
     
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE THERMO_ANELASTIC_LAPSE_FR

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_LAPSE_EQU(nx,ny,nz, s, dTdy,dqldy, e,p,r, lapse, frequency)

  USE THERMO_VARS
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: dTdy,dqldy 
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p,r
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: lapse, frequency

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is, ipsat
  TREAL P_LOC, E_LOC, T_LOC, R_LOC, RT_INV, psat, dpsat

  TREAL Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl, Lv, C

  TREAL scaleheightinv, one_p_eps, rd_ov_rv, qvequ, qsat, dummy
  
! ###################################################################
  Rv = WGHT_INV(1)
  Rd = WGHT_INV(2)              
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)
  Lv0=-THERMO_AI(6,1,3)

  scaleheightinv = GRATIO /scaleheight
  rd_ov_rv = WGHT_INV(2) /WGHT_INV(1)
    
  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     lapse     = scaleheightinv

     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = s(ij,1) - E_LOC
           frequency(ij) = ( lapse(ij) + dTdy(ij) )/ T_LOC
           
        ENDDO
     
     ENDDO
     
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     lapse = scaleheightinv /( Cd  + s(:,2) *Cdv )

     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           frequency(ij) = ( lapse(ij) + dTdy(ij) )/ T_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        
        RT_INV = R_LOC /P_LOC /scaleheight
        
        DO i = 1,nx
           ij = ij +1

           C = Cd + s(ij,2) *Cdv + s(ij,3) *Cvl ! I need it below
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / C
           
           psat = THERMO_PSAT(NPSAT); dpsat = C_0_R
           DO ipsat = NPSAT-1,1,-1
              psat  = psat *T_LOC + THERMO_PSAT(ipsat)
              dpsat = dpsat*T_LOC + THERMO_PSAT(ipsat+1) *M_REAL(ipsat)
           ENDDO
           dummy = rd_ov_rv /( P_LOC/psat -C_1_R )
           qsat = dummy /( C_1_R +dummy )

! We cannot use ql directly (s(:,3)) because the smoothing function imposes
! an exponentially small value, but nonzero           
           IF ( qsat .GE. s(ij,2) ) THEN
              lapse(ij)     = scaleheightinv /( Cd  + s(ij,2) *Cdv )
              frequency(ij) = ( lapse(ij) + dTdy(ij) )/ T_LOC
              
           ELSE
              qvequ = dummy *( C_1_R -s(ij,2) )
              
              one_p_eps = C_1_R /( C_1_R - psat /P_LOC )
              Lv = Lv0 - T_LOC *Cvl
              
              lapse(ij) = scaleheightinv + qvequ *one_p_eps *Lv *RT_INV 
              lapse(ij) = lapse(ij) /( C + qvequ *one_p_eps *Lv *dpsat /psat )

              frequency(ij) = qvequ *( lapse(ij) *dpsat /psat -RT_INV )* one_p_eps - dqldy(ij)
              frequency(ij) = ( lapse(ij) + dTdy(ij) )/ T_LOC &
                            - frequency(ij) *Rv /( Rd + s(ij,2) *Rdv - s(ij,3) *Rv )

           ENDIF
           
        ENDDO
        
     ENDDO
     
  ENDIF
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_LAPSE_EQU

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_DEWPOINT(nx,ny,nz, s, e,p,r, Td,Lapse)

  USE THERMO_VARS
  
  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e,p,r
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: Td,Lapse

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is, ipsat
  TINTEGER inr, nrmax
  TREAL P_LOC, E_LOC, T_LOC, R_LOC, psat, dpsat, dummy, qsat

  TREAL Rv, Rd, Rdv, Cd, Cdv, Lv0, Cvl, rd_ov_rv 
!  TREAL NEWTONRAPHSON_ERROR
  TREAL scaleheightinv
  
! ###################################################################
! maximum number of iterations in Newton-Raphson
  nrmax = 5

! initialize homogeneous data
  Rv = WGHT_INV(1)
  Rd = WGHT_INV(2)              
  Rdv= WGHT_INV(1) - WGHT_INV(2)
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)
  Lv0=-THERMO_AI(6,1,3)

  rd_ov_rv = WGHT_INV(2) /WGHT_INV(1)
  scaleheightinv = C_1_R /scaleheight

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        
        DO i = 1,nx
           ij = ij +1
           
! Using actual temperature as initial condition
           T_LOC = ( s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           
! executing Newton-Raphson 
           DO inr = 1,nrmax
              psat = THERMO_PSAT(NPSAT); dpsat = C_0_R
              DO ipsat = NPSAT-1,1,-1
                 psat  = psat *T_LOC + THERMO_PSAT(ipsat)
                 dpsat = dpsat*T_LOC + THERMO_PSAT(ipsat+1) *M_REAL(ipsat)
              ENDDO
! we seek root of function:    psat -P_LOC *s(ij,2) *Rv /( Rd +s(ij,2) *Rdv )                               
              T_LOC = T_LOC -( psat -P_LOC *s(ij,2) *Rv /( Rd +s(ij,2) *Rdv ) )/dpsat
           ENDDO
!           NEWTONRAPHSON_ERROR = MAX(NEWTONRAPHSON_ERROR,ABS(psat/dpsat)/T_LOC)
           Td(ij) = T_LOC
           Lapse(ij) = scaleheightinv *R_LOC/ P_LOC *psat /dpsat
           
        ENDDO
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        P_LOC = MRATIO *p(is)
        E_LOC = e(is)
        R_LOC = r(is)
        
        DO i = 1,nx
           ij = ij +1

! Using actual temperature as initial condition
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 ) / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )

! We cannot use ql directly (s(:,3)) because the smoothing function imposes
! and exponentially small value, but nonzero           
           psat = THERMO_PSAT(NPSAT)
           DO ipsat = NPSAT-1,1,-1
              psat  = psat *T_LOC + THERMO_PSAT(ipsat)
           ENDDO
           dummy = rd_ov_rv /( P_LOC/psat -C_1_R )
           qsat = dummy /( C_1_R +dummy )

           IF ( qsat .LE. s(ij,2) ) THEN
              Td(ij) = T_LOC

           ELSE
! executing Newton-Raphson 
              DO inr = 1,nrmax
                 psat = THERMO_PSAT(NPSAT); dpsat = C_0_R
                 DO ipsat = NPSAT-1,1,-1
                    psat  = psat *T_LOC + THERMO_PSAT(ipsat)
                    dpsat = dpsat*T_LOC + THERMO_PSAT(ipsat+1) *M_REAL(ipsat)
                 ENDDO
!                                 psat -P_LOC *(s(ij,2)-s(ij,3)) *Rv /( Rd +s(ij,2) *Rdv -s(ij,3) *Rv )
! we seek root of function:       psat -P_LOC *s(ij,2) *Rv /( Rd +s(ij,2) *Rdv )                 
                 T_LOC = T_LOC -( psat -P_LOC *s(ij,2) *Rv /( Rd +s(ij,2) *Rdv ) )/dpsat
              ENDDO
              ! NEWTONRAPHSON_ERROR = MAX(NEWTONRAPHSON_ERROR,ABS(psat/dpsat)/T_LOC)
              ! print*,NEWTONRAPHSON_ERROR
              Td(ij) = T_LOC
              Lapse(ij) = scaleheightinv *R_LOC/ P_LOC *psat /dpsat
              
           ENDIF
           
        ENDDO
     ENDDO
     
  ENDIF
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_DEWPOINT

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_WEIGHT_INPLACE(nx,ny,nz, weight, a)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)    :: weight
  TREAL, DIMENSION(nx,ny*nz),   INTENT(INOUT) :: a

! -------------------------------------------------------------------
  TINTEGER jk, j

! ###################################################################
  DO jk = 1,ny*nz
     j = MOD(jk-1,ny) +1
     
     a(1:nx,jk) = a(1:nx,jk) *weight(j)
     
  ENDDO
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_WEIGHT_INPLACE

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx,ny,nz, weight, a, b)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)  :: weight
  TREAL, DIMENSION(nx,ny*nz),   INTENT(IN)  :: a
  TREAL, DIMENSION(nx,ny*nz),   INTENT(OUT) :: b

! -------------------------------------------------------------------
  TINTEGER jk, j

! ###################################################################
  DO jk = 1,ny*nz
     j = MOD(jk-1,ny) +1
     
     b(1:nx,jk) = a(1:nx,jk) *weight(j)
     
  ENDDO

  RETURN
END SUBROUTINE THERMO_ANELASTIC_WEIGHT_OUTPLACE

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_WEIGHT_ADD(nx,ny,nz, weight, a, b)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)    :: weight
  TREAL, DIMENSION(nx,ny*nz),   INTENT(IN)    :: a
  TREAL, DIMENSION(nx,ny*nz),   INTENT(INOUT) :: b

! -------------------------------------------------------------------
  TINTEGER jk, j

! ###################################################################
  DO jk = 1,ny*nz
     j = MOD(jk-1,ny) +1
     
     b(1:nx,jk) = b(1:nx,jk) +a(1:nx,jk) *weight(j)
     
  ENDDO

  RETURN
END SUBROUTINE THERMO_ANELASTIC_WEIGHT_ADD

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_WEIGHT_SUBSTRACT(nx,ny,nz, weight, a, b)

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(*),          INTENT(IN)    :: weight
  TREAL, DIMENSION(nx,ny*nz),   INTENT(IN)    :: a
  TREAL, DIMENSION(nx,ny*nz),   INTENT(INOUT) :: b

! -------------------------------------------------------------------
  TINTEGER jk, j

! ###################################################################
  DO jk = 1,ny*nz
     j = MOD(jk-1,ny) +1
     
     b(1:nx,jk) = b(1:nx,jk) -a(1:nx,jk) *weight(j)
     
  ENDDO

  RETURN
END SUBROUTINE THERMO_ANELASTIC_WEIGHT_SUBSTRACT

!########################################################################
!########################################################################
SUBROUTINE THERMO_ANELASTIC_LWP(nx,ny,nz, g, r, ql, lwp, wrk1d,wrk3d)

  USE TLAB_TYPES, ONLY : grid_dt

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TYPE(grid_dt),              INTENT(IN)    :: g
  TREAL, DIMENSION(*),        INTENT(IN)    :: r
  TREAL, DIMENSION(nx*nz,ny), INTENT(IN)    :: ql
  TREAL, DIMENSION(nx,nz),    INTENT(OUT)   :: lwp
  TREAL, DIMENSION(ny),       INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(nx,ny,nz), INTENT(INOUT) :: wrk3d

! -------------------------------------------------------------------
  TINTEGER i,j,k
  TREAL SIMPSON_NU

! ###################################################################
  CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx,ny,nz, r, ql, wrk3d)

  DO k = 1,nz
     DO i = 1,nx
        DO j = 1,ny
           wrk1d(j) = wrk3d(i,j,k)
        ENDDO
        lwp(i,k) = SIMPSON_NU(ny, wrk1d, g%nodes)
     ENDDO
  ENDDO
     
  RETURN
END SUBROUTINE THERMO_ANELASTIC_LWP

!########################################################################
!########################################################################
! Just to check what the effect of using a wrong cp would be
SUBROUTINE THERMO_ANELASTIC_STATIC_CONSTANTCP(nx,ny,nz, s, e, result)

  USE THERMO_VARS, ONLY : imixture, THERMO_AI

  IMPLICIT NONE

  TINTEGER,                     INTENT(IN)  :: nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: s
  TREAL, DIMENSION(*),          INTENT(IN)  :: e
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: result

! -------------------------------------------------------------------
  TINTEGER ij, i, jk, is
  TREAL E_LOC, T_LOC

  TREAL Cd, Cdv, Lv0, Cvl, Lv

! ###################################################################
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Lv0=-THERMO_AI(6,1,3)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  IF      ( imixture .EQ. MIXT_TYPE_AIR      ) THEN
     result = s(:,1)

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRVAPOR ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
        
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC ) / ( Cd + s(ij,2) *Cdv )
           result(ij) = T_LOC *Cd +E_LOC
           
        ENDDO
     
     ENDDO

  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     ij = 0
     DO jk = 0,ny*nz-1
        is = MOD(jk,ny) +1
        E_LOC = e(is)
       
        DO i = 1,nx
           ij = ij +1
           
           T_LOC = (s(ij,1) - E_LOC + s(ij,3)*Lv0 )  / ( Cd + s(ij,2) *Cdv + s(ij,3) *Cvl )

           Lv = Lv0 - T_LOC *Cvl
           result(ij) = T_LOC *Cd +E_LOC - s(ij,3) *Lv
           
        ENDDO
     
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE THERMO_ANELASTIC_STATIC_CONSTANTCP

