#include "types.h"
#include "dns_const.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE BOUNDARY_SURFACE_J
! Calculates and updates interactive surface boundary condition
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_SURFACE_J(is,bcs,s,hs,tmp1,tmp2,aux,wrk1d,wrk2d,wrk3d)
#ifdef TRACE_ON
  USE DNS_CONSTANTS,ONLY : tfile
#endif
  USE DNS_CONSTANTS,ONLY : lfile
  USE DNS_GLOBAL,   ONLY : imax,jmax,kmax,inb_flow,inb_vars,inb_scal,g
  USE DNS_GLOBAL,   ONLY : isize_field,isize_wrk1d 
  USE DNS_GLOBAL,   ONLY : visc,itime,schmidt
  USE DNS_GLOBAL,   ONLY : imode_fdm 
  USE DNS_LOCAL,    ONLY : dtime, rkm_substep, rkm_endstep
  USE BOUNDARY_BCS, ONLY : BcsScalJmin, BcsScalJmax 

  IMPLICIT NONE  

#include "integers.h"

  TINTEGER is
  TINTEGER, DIMENSION(2,2), INTENT(IN) :: bcs          ! Boundary conditions from derivative operator
  TREAL, DIMENSION(isize_field,*)      :: s,hs
  TREAL, DIMENSION(isize_field)        :: tmp1,tmp2
  TREAL, DIMENSION(imax,kmax,6),TARGET :: aux 
  TREAL, DIMENSION(isize_wrk1d,*)      :: wrk1d
  TREAL, DIMENSION(*)                  :: wrk2d,wrk3d

  TINTEGER nxy,ip,k
  TREAL, DIMENSION(:,:), POINTER       :: hfx,hfx_anom
  TREAL :: diff,hfx_avg
  TREAL AVG1V2D

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile,'ENTERING SUBROUTINE BOUNDARY_SURFACE_J')
#endif
  diff = visc/schmidt(is)
  nxy = imax*jmax

  ! vertical derivative of scalar for flux at the boundaries
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(:,is), tmp1,wrk3d,wrk2d,wrk3d)

  ! ------------------------------------------------------------
  ! Bottom Boundary
  ! ------------------------------------------------------------
  IF ( BcsScalJmin%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
     hfx =>      aux(:,:,1)
     hfx_anom => aux(:,:,2)
     ip=1
     DO k=1,kmax    ! Calculate the surface flux
        hfx(:,k) = diff*tmp1(ip:ip+imax-1); ip=ip+nxy
     ENDDO
     hfx_avg = diff*AVG1V2D(imax,jmax,kmax,1,1,tmp1)
     hfx_anom = hfx - hfx_avg
     BcsScalJmin%ref(:,:,is) = BcsScalJmin%ref(:,:,is) + BcsScalJmin%cpl(is)*hfx_anom
  ENDIF


  ! ------------------------------------------------------------
  ! Top Boundary
  ! ------------------------------------------------------------
  IF ( BcsScalJmax%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
     hfx =>      aux(:,:,3)
     hfx_anom => aux(:,:,4)
     ip = imax*(jmax-1) + 1
     DO k=1,kmax;     ! Calculate the surface flux
        hfx(:,k) = -diff*tmp1(ip:ip+imax-1); ip=ip+nxy;
     ENDDO
     hfx_avg = diff*AVG1V2D(imax,jmax,kmax,1,1,tmp1)
     hfx_anom = hfx - hfx_avg
     BcsScalJmax%ref(:,:,is) = BcsScalJmax%ref(:,:,is) + BcsScalJmax%cpl(is)*hfx_anom
  ENDIF


#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(TFILE,'LEAVING SUBROUTINE BOUNDAR_SURFACE_J')
#endif

  RETURN

END SUBROUTINE BOUNDARY_SURFACE_J
