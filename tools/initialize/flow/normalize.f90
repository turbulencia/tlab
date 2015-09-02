!########################################################################
!# Tool/Library INIT/FLOW
!#
!########################################################################
!# HISTORY
!#
!# 2006/10/07 - J.P. Mellado
!#              Created
!# 2007/03/10 - J.P. Mellado
!#              Uniform treatment of different flows based on K
!#
!########################################################################
!# DESCRIPTION
!#
!# Originally extracted from brdb.f to use also in discr.f
!# It scales the velocity field to get a given value of turbulent kinetic
!# energy, K. It assumes fields with zero mean.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

SUBROUTINE NORMALIZE(nx,ny,nz, u,v,w, tke1)

  USE DNS_GLOBAL, ONLY : imode_flow

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(nx,ny,nz) :: u, v, w
  TREAL tke1

! -------------------------------------------------------------------
  TINTEGER j
  TREAL q, tke0, factor, AVG1V2D, AVG1V3D 

! ###################################################################
! Get actual reference value of K
! ###################################################################
  IF      ( imode_flow .EQ. DNS_FLOW_SHEAR .OR. imode_flow .EQ. DNS_FLOW_JET ) THEN
! maximum across the layer
     tke0= C_0_R
     DO j = 1,ny
        q = AVG1V2D(nx,ny,nz, j, i2, u) +&
            AVG1V2D(nx,ny,nz, j, i2, v) +&
            AVG1V2D(nx,ny,nz, j, i2, w)
        tke0 = MAX(tke0,q)
     ENDDO
     tke0 = C_05_R*tke0

  ELSE IF ( imode_flow .EQ. DNS_FLOW_ISOTROPIC ) THEN
     tke0 = C_05_R*( AVG1V3D(nx,ny,nz, i2, u) +&
                     AVG1V3D(nx,ny,nz, i2, v) +&
                     AVG1V3D(nx,ny,nz, i2, w) )
  ENDIF

  factor = SQRT(tke1/tke0)

! ###################################################################
! Scale flow field
! ###################################################################
  u = u*factor
  v = v*factor
  w = w*factor

  RETURN
END SUBROUTINE NORMALIZE
