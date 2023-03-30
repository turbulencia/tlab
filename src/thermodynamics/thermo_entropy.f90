#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool DNS
!#
!########################################################################
!# HISTORY
!#
!# 2007/05/22 - J.P. Mellado
!#              Created
!# 2007/07/04 - J.P. Mellado
!#              Including case AIRWATER.
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate entropy from T, p and composition. The reference state 
!# is T_0, which is set to 298 K in the case of multispecies, and 
!# pbg%mean. Nondimensional with C_{p,0}.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_ENTROPY(nx,ny,nz, z1,T,p, s)

  USE TLAB_VARS, ONLY : pbg

  USE THERMO_VARS, ONLY : imixture, GRATIO
  USE THERMO_VARS, ONLY : NSP, NCP, WGHT_INV, THERMO_AI, THERMO_TLIM

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz),   INTENT(IN)  :: T,p
  TREAL, DIMENSION(nx*ny*nz,*), INTENT(IN)  :: z1
  TREAL, DIMENSION(nx*ny*nz),   INTENT(OUT) :: s

! -------------------------------------------------------------------
  TINTEGER ij, is, im, icp
  TREAL WMEAN_INV, ENTROPY_I, XMOL_I
  TREAL, dimension(NSP) :: YMASS

! ###################################################################
! Single species
! ###################################################################
  IF ( imixture .EQ. 0 ) THEN
     DO ij = 1,nx*ny*nz
        s(ij) = log(T(ij)/(p(ij)/pbg%mean)**GRATIO)
     ENDDO

! ###################################################################
! mixture MIXT_TYPE_AIRWATER
! ###################################################################
  ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     DO ij = 1,nx*ny*nz
        YMASS(1) = z1(ij,1)-z1(ij,2) !q_v=q_t-q_l
        YMASS(2) = C_1_R-z1(ij,1)    !q_d=1-q_t
! calculate temperature part
        WMEAN_INV = C_0_R
        DO is = 1,2
           IF ( T(ij) .LT. THERMO_TLIM(3,is) ) THEN; im = 2
           ELSE;                                     im = 1; ENDIF
           ENTROPY_I = C_0_R
           DO icp = NCP,2,-1
              ENTROPY_I = ENTROPY_I*T(ij) + THERMO_AI(icp,im,is)/M_REAL(icp-1)
           ENDDO
           ENTROPY_I = ENTROPY_I*T(ij) + THERMO_AI(7,im,is) + &
                THERMO_AI(1,im,is)*log(T(ij))
           s(ij)      = s(ij)     + YMASS(is)*ENTROPY_I
           WMEAN_INV  = WMEAN_INV + YMASS(is)*WGHT_INV(is)
        ENDDO
! calculate pressure part
        DO is = 1,2
           XMOL_I = YMASS(is)*WGHT_INV(is)/WMEAN_INV
           IF ( XMOL_I .GT. C_0_R ) THEN
              s(ij) = s(ij) - GRATIO*YMASS(is)*WGHT_INV(is)*log(XMOL_I)
           ENDIF
        ENDDO
        s(ij) = s(ij) - GRATIO*WMEAN_INV*log(p(ij)/pbg%mean)

        s(ij) = s(ij) + z1(ij,2)*(THERMO_AI(7,im,3)+THERMO_AI(1,1,3)*log(T(ij)))

     ENDDO

! ###################################################################
! General mixture of NSP species. z1 contains NSP-1 mass fractions
! ###################################################################
  ELSE 
     s(1:nx*ny*nz) = C_0_R
     DO ij = 1,nx*ny*nz
! pass species to YMASS vector
        YMASS(NSP) = C_1_R
        DO is = 1, NSP-1
           YMASS(is) = z1(ij,is)
           YMASS(NSP)= YMASS(NSP)-YMASS(is)
        ENDDO
! calculate temperature part
        WMEAN_INV = C_0_R
        DO is = 1,NSP
           IF ( T(ij) .LT. THERMO_TLIM(3,is) ) THEN; im = 2
           ELSE;                                     im = 1; ENDIF
           ENTROPY_I = C_0_R
           DO icp = NCP,2,-1
              ENTROPY_I = ENTROPY_I*T(ij) + THERMO_AI(icp,im,is)/M_REAL(icp-1)
           ENDDO
           ENTROPY_I = ENTROPY_I*T(ij) + THERMO_AI(7,im,is) + &
                THERMO_AI(1,im,is)*log(T(ij))
           s(ij)     = s(ij)     + YMASS(is)*ENTROPY_I
           WMEAN_INV = WMEAN_INV + YMASS(is)*WGHT_INV(is)
        ENDDO
! calculate pressure part
        DO is = 1,NSP
           XMOL_I = YMASS(is)*WGHT_INV(is)/WMEAN_INV
           IF ( XMOL_I .GT. C_0_R ) THEN
              s(ij) = s(ij) - GRATIO*YMASS(is)*WGHT_INV(is)*log(XMOL_I)
           ENDIF
        ENDDO
        s(ij) = s(ij) - GRATIO*WMEAN_INV*log(p(ij)/pbg%mean)
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE THERMO_ENTROPY
