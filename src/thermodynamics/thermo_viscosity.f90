#include "types.h"
#include "dns_const.h"
  
!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the nondimensional dynamic viscosity \mu=\mu(T) depending
!# on the given flag itransport:
!#
!# 0 None
!# 1 Powerlaw
!# 2 Sutherland
!# 
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE THERMO_VISCOSITY(nx,ny,nz, T, mu)
  
  USE TLAB_VARS, ONLY : itransport

  IMPLICIT NONE
  
  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz) :: T, mu

! -----------------------------------------------------------------------
  TREAL dummy

! #######################################################################
  IF      ( itransport .EQ. EQNS_TRANS_SUTHERLAND ) THEN ! not yet implemented
     mu = C_1_R  

  ELSE IF ( itransport .EQ. EQNS_TRANS_POWERLAW   ) THEN
     dummy = C_7_R/C_10_R
     mu = T**dummy

  ELSE
     mu = C_1_R 
     
  ENDIF
          
  RETURN
END SUBROUTINE THERMO_VISCOSITY
