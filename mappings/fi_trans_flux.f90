#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2015/06/25 - A de Lozar
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the transport terms due to settling
!# There are two options for the settling model.The first corresponds to a 
!# constant settling velocity while in the second the settling velocity 
!# assumes that the droplet spectrum is represented by a log-normal distribution.
!# For each case we consider a simplified calculation and a more exact one (differences of 10^-4)
!#
!########################################################################
SUBROUTINE FI_TRANS_FLUX(itransport, nx,ny,nz, is, is_trans, param, settling, &
     dy, s,trans, tmp, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imode_fdm, j1bc
  
  IMPLICIT NONE

#include "integers.h"

  TINTEGER,                          INTENT(IN)    :: itransport, nx,ny,nz
  TINTEGER,                          INTENT(IN)    :: is         ! Scalar for which the transport term is calculated
  TINTEGER,                          INTENT(IN)    :: is_trans   ! Scalar used to calculate the transport
  TREAL, DIMENSION(*),               INTENT(IN)    :: param      ! Transport Parameters
  TREAL,                             INTENT(IN)    :: settling   ! Settling Parameter
  TREAL, DIMENSION(nx*ny*nz,*),      INTENT(IN)    :: s          ! Array with all  scalars
  TREAL, DIMENSION(nx*ny*nz,*),      INTENT(INOUT) :: tmp        ! It saves some calculations to include some calculations outside the function in some cases
  TREAL, DIMENSION(nx*ny*nz,*),      INTENT(OUT)   :: trans      ! Transport component. It can have three directions

  TREAL, DIMENSION(*)         :: dy
  TREAL, DIMENSION(nx*ny*nz)  :: wrk3d
  TREAL, DIMENSION(*)         :: wrk1d 
  TREAL, DIMENSION(nx*nz)     :: wrk2d

! -----------------------------------------------------------------------
  TREAL dummy
  TINTEGER ij

!########################################################################
  IF     (itransport .EQ. EQNS_TRANS_C_SET ) THEN ! Constant Settling formulation in simplified way
     IF ( is .EQ. 1 ) THEN
        CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, &
             dy, s(:,is_trans), tmp, i0,i0, wrk1d,wrk2d,wrk3d) !Derivative of the liquid field on tmp
     ENDIF
     
     dummy = settling*param(is)
     DO ij = 1,nx*ny*nz    
        trans(ij,1) =  dummy*tmp(ij,1) ! Linear formulation (all sizes with the same velocity). Assume ds/dy in tmp
     ENDDO

  ELSEIF (itransport .EQ. EQNS_TRANS_LN_SET ) THEN
     IF ( is .EQ. 1 ) THEN
        CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, &
             dy, s(:,is_trans), tmp, i0,i0, wrk1d,wrk2d,wrk3d) !Derivative of the liquid field on tmp
        dummy = C_2_R/C_3_R
        DO ij=1,nx*ny*nz
           tmp(ij,1) = tmp(ij,1)*(s(ij,is_trans)**dummy)
        ENDDO        
     ENDIF

     dummy = settling*param(is)*C_5_R/C_3_R
     DO ij = 1,nx*ny*nz
        trans(ij,1) = dummy*tmp(ij,1) ! Linear formulation (all sizes with the same velocity). Assume ds/dy*s^(2/3) in tmp 
     ENDDO

  ELSEIF (itransport .EQ. EQNS_TRANS_C_SET_FULL ) THEN
     DO ij = 1,nx*ny*nz ! Loop to calculate the vector flux
        tmp(ij,1) = settling*(param(is) - param(3)*s(ij,is))*s(ij,is_trans) ! constant settling velocity
     ENDDO
     CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, &
          dy, tmp(1,1), trans(1,1), i0,i0, wrk1d,wrk2d,wrk3d) !Gradient of the vector flux

  ELSEIF (itransport .EQ. EQNS_TRANS_LN_SET_FULL ) THEN
     dummy = C_5_R/C_3_R
     DO ij = 1,nx*ny*nz
        trans(ij,1) = settling*(param(is) - param(3)*s(ij,is))*(s(ij,is_trans)**dummy) ! constant settling velocity
     ENDDO
     CALL PARTIAL_Y(imode_fdm, nx,ny,nz, j1bc, &
          dy, tmp(1,1), trans(1,1), i0,i0, wrk1d,wrk2d,wrk3d) !Gradient of the vector flux

  ENDIF

  RETURN
END SUBROUTINE FI_TRANS_FLUX
