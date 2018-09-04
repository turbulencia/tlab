#include "types.h"
#include "dns_const.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE BOUNDARY_SURFACE_J
! Calculates and updates interactive surface boundary condition
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUNDARY_SURFACE_J(is,bcs,q,hq,s,hs,tmp1,tmp2,aux,wrk1d,wrk2d,wrk3d)
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
  TREAL, DIMENSION(isize_field,*)      :: q,hq,s,hs 
  TREAL, DIMENSION(isize_field)        :: tmp1,tmp2
  TREAL, DIMENSION(imax,kmax,6),TARGET :: aux 
  TREAL, DIMENSION(isize_wrk1d,*)      :: wrk1d
  TREAL, DIMENSION(*)                  :: wrk2d,wrk3d

  TINTEGER nxy,ip,j,k
  TREAL, DIMENSION(:,:), POINTER       :: hfx,hfx_anom
  TREAL :: diff,hfx_avg
  TREAL AVG1V2D 
  TREAL var,var2,avg,avg2,avg_anom
  CHARACTER lstr*256

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile,'ENTERING SUBROUTINE BOUNDARY_SURFACE_J')
#endif
  diff = visc/schmidt(is)
  nxy = imax*jmax

  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(:,is), tmp1,wrk3d,wrk2d,wrk3d)
  ! ------------------------------------------------------------
  ! Bottom Boundary
  ! ------------------------------------------------------------
  IF ( BcsScalJmin%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
     hfx =>      aux(:,:,1)
     hfx_anom => aux(:,:,2)
     ip=1
     DO k=1,kmax
        hfx(:,k) = diff*tmp1(ip:ip+imax-1); ip=ip+nxy
     ENDDO
     hfx_avg = diff*AVG1V2D(imax,jmax,kmax,1,1,tmp1)
     hfx_anom = hfx - hfx_avg
     BcsScalJmin%ref(:,:,is) = BcsScalJmin%ref(:,:,is) + BcsScalJmin%cpl(is)*hfx_anom

     avg_anom=AVG1V2D(imax,1,kmax,1,1,hfx_anom)
     avg=AVG1V2D(imax,jmax,kmax,1,1,s)
     var=AVG1V2D(imax,jmax,kmax,1,2,s)
     avg2=AVG1V2D(imax,jmax,kmax,2,1,s)
     var2=AVG1V2D(imax,jmax,kmax,2,2,s) 
     WRITE(lstr,*) itime,'Jmin: Heat flux:', hfx_avg, avg_anom,MINVAL(hfx_anom),MAXVAL(hfx_anom),avg,var-avg**2,var2-avg2**2
     CALL IO_WRITE_ASCII(lfile,lstr)

     ! TESTING: solve ds/dt=s => s(t) = exp(t) on the the surface boundary
     ! WRITE(*,*) 'Jmin-sfc BEF', itime, MINVAL(BcsScalJmin%ref(:,:,is)), MAXVAL(BcsScalJmin%ref(:,:,is)), s(1,is), BcsScalJmin%cpl(is)
     ! ip=1
     ! DO k=1,kmax
     !    BcsScalJmin%ref(:,k,is) = BcsScalJmin%ref(:,k,is) + s(ip:ip+imax-1,is)*BcsScalJmin%cpl(is)
     !    ip = ip+nxy
     ! ENDDO
  ENDIF


  ! ------------------------------------------------------------
  ! Top Boundary
  ! ------------------------------------------------------------
  IF ( BcsScalJmax%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
     hfx =>      aux(:,:,3)
     hfx_anom => aux(:,:,4)
     ip = imax*(jmax-1) + 1
     DO k=1,kmax;
        hfx(:,k) = -diff*tmp1(ip:ip+imax-1); ip=ip+nxy;
     ENDDO

     hfx_avg = diff*AVG1V2D(imax,jmax,kmax,1,1,tmp1)
     hfx_anom = hfx - hfx_avg
     BcsScalJmax%ref(:,:,is) = BcsScalJmax%ref(:,:,is) + BcsScalJmax%cpl(is)*hfx_anom
     IF ( rkm_substep .EQ. rkm_endstep ) THEN
        avg=AVG1V2D(imax,jmax,kmax,g(2)%size,1,s)
        var=AVG1V2D(imax,jmax,kmax,g(2)%size,2,s) 
        avg2=AVG1V2D(imax,jmax,kmax,g(2)%size-1,1,s) 
        var2=AVG1V2D(imax,jmax,kmax,g(2)%size-1,2,s)
        WRITE(lstr,*) itime,'Jmax: Heat flux:', hfx_avg, MINVAL(hfx_anom),MAXVAL(hfx_anom),avg,var-avg**2,var2-avg2**2

        CALL IO_WRITE_ASCII(lfile,lstr)
     ENDIF

     ! TESTING: solve ds/dt=s => s(t) = exp(t) on the the surface boundary
     ! WRITE(*,*) 'Jmax-sfc BEF', itime, MINVAL(BcsScalJmax%ref(:,:,is)), MAXVAL(BcsScalJmax%ref(:,:,is)), s(imax*(jmax-1)+1,is)
     ! ip = imax*(jmax-1) + 1
     ! DO k=1,kmax
     !    BcsScalJmax%ref(:,k,is) = BcsScalJmax%ref(:,k,is) + s(ip:ip+imax-1,is)*BcsScalJmax%cpl(is)
     !    ip = ip+nxy
     ! ENDDO
  ENDIF


#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(TFILE,'LEAVING SUBROUTINE BOUNDAR_SURFACE_J')
#endif

  RETURN

END SUBROUTINE BOUNDARY_SURFACE_J
