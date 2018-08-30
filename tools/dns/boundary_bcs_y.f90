#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2003/06/11 - J.P. Mellado
!#              Created
!# 2007/06/21 - J.P. Mellado
!#              AIRWATER case included.
!#
!########################################################################
!# DESCRIPTION
!# 
!# Non-periodic characteristic BCs at ymin and ymax.
!# The flunctuating inflow forcing has not yet been implemented like
!# in BOUNDARY_BCS_X.
!#
!########################################################################
SUBROUTINE BOUNDARY_BCS_Y(iaux, M2_max, rho,u,v,w,p,gama,z1, &
     h0,h1,h2,h3,h4,zh1, tmp1,tmp2,tmp3,tmp4,tmp5, aux2d, wrk2d,wrk3d)

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture, gama0, THERMO_AI
  USE DNS_LOCAL
  USE BOUNDARY_BCS
  
  IMPLICIT NONE

#include "integers.h"

  TINTEGER iaux
  TREAL M2_max

  TREAL, DIMENSION(imax,jmax,kmax)   :: rho, u, v, w, p, gama, h0, h1, h2, h3, h4
  TREAL, DIMENSION(imax,jmax,kmax)   :: tmp1, tmp2, tmp3, tmp4, tmp5
  TREAL, DIMENSION(imax,jmax,kmax,*) :: z1, zh1
  TREAL, DIMENSION(imax,kmax,*)      :: aux2d
  TREAL, DIMENSION(*)                :: wrk2d, wrk3d

  TARGET aux2d

! -------------------------------------------------------------------
  TINTEGER i, k, is, nt, inb_scal_loc, iflag_min, iflag_max, idir, ip0, bcs(2,1)
  TINTEGER imin_loc, imax_loc
  TREAL prefactor, pl_out, pl_inf

  TREAL, DIMENSION(:,:,:), POINTER :: tmin, lmin, tmax, lmax, inf_rhs

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_Y' )
#endif

#define hr_loc(i,k)  aux2d(i,k,1)
#define hu_loc(i,k)  aux2d(i,k,2)
#define hv_loc(i,k)  aux2d(i,k,3)
#define hw_loc(i,k)  aux2d(i,k,4)
#define he_loc(i,k)  aux2d(i,k,5)
#define hz1_loc(i,k) aux2d(i,k,6)

#define r_loc(i,k)   aux2d(i,k,7)
#define u_loc(i,k)   aux2d(i,k,8)
#define v_loc(i,k)   aux2d(i,k,9)
#define w_loc(i,k)   aux2d(i,k,10)
#define p_loc(i,k)   aux2d(i,k,11)
#define g_loc(i,k)   aux2d(i,k,12)
#define z1_loc(i,k)  aux2d(i,k,13)

#define drdn_loc(i,k)  aux2d(i,k,14)
#define dudn_loc(i,k)  aux2d(i,k,15)
#define dvdn_loc(i,k)  aux2d(i,k,16)
#define dwdn_loc(i,k)  aux2d(i,k,17)
#define dpdn_loc(i,k)  aux2d(i,k,18)
#define dz1dn_loc(i,k) aux2d(i,k,19)

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero
    
  ip0 = 19

  nt = imax*kmax
  prefactor = (gama0-C_1_R)*mach*mach

  IF ( iaux .LT. nt*(19+5*(inb_flow+inb_scal_array)) ) THEN
     CALL IO_WRITE_ASCII(efile, 'RHS_BCS_Y. Not enough space.')
     CALL DNS_STOP(DNS_ERROR_JBC)
  ENDIF

! Define pointers
  inf_rhs => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  tmin    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  lmin    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  tmax    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  lmax    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )

! -------------------------------------------------------------------
! Type of characteristic BCs
! 1. only nonreflective
! 2. add fluctuation
! 3. add mean
! 4. add fluctuation+mean
! 
! Relaxation towards a mean profile (Poinsot & Lele term) 
! The local value of c is added later at the boundary
! Note that pl_??? has dimensions of 1/length
! -------------------------------------------------------------------
  idir = 2
  IF ( bcs_euler_drift .EQ. 1 ) THEN
     pl_out = bcs_sigma_out*(C_1_R-M2_max)/g(2)%scale
     pl_inf = bcs_sigma_inf_j/g(2)%scale
     iflag_min =-3
     iflag_max = 3
  ELSE
     pl_out = C_0_R
     pl_inf = C_0_R
     iflag_min =-1
     iflag_max = 1
  ENDIF

  ! pl_out_min = C_0_R ! default is only nonreflective
  ! pl_inf_min = C_0_R
  ! iflag_min =-1      
  ! IF ( BcsFlowJmin%cinf .GT. 0 ) THEN
  !    pl_inf_min = BcsFlowJmin%cinf /g(2)%scale
  !    iflag_min =-3
  ! ENDIF
  ! IF ( BcsFlowJmin%cout .GT. 0 ) THEN
  !    pl_out_min = BcsFlowJmin%cout *(C_1_R-M2_max) /g(2)%scale
  !    iflag_min =-3
  ! ENDIF
  
  IF      ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN; imin_loc  = 1; imax_loc = imax
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL  ) THEN; imin_loc  = 2; imax_loc = imax-1; ENDIF

! ###################################################################
! Transverse terms
! ###################################################################
  CALL BOUNDARY_BCS_TRANSVERSE_Y(u,v,w,p,rho,gama, z1, &
       tmin,lmin,tmax,lmax, tmp1,tmp2,tmp3, wrk2d,wrk3d)

! ###################################################################
! Flow 
! ###################################################################
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), rho, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), u,   tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v,   tmp3, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), w,   tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), p,   tmp5, wrk3d, wrk2d,wrk3d)

! -------------------------------------------------------------------
! BCs at ymin
! -------------------------------------------------------------------
  DO k = 1,kmax; DO i = 1,imax
     r_loc(i,k) =  rho(i,1,k)
     u_loc(i,k) =    v(i,1,k)
     v_loc(i,k) =    u(i,1,k)
     w_loc(i,k) =    w(i,1,k)
     p_loc(i,k) =    p(i,1,k)
     g_loc(i,k) = gama(i,1,k)
     drdn_loc(i,k) = tmp1(i,1,k)
     dudn_loc(i,k) = tmp3(i,1,k)
     dvdn_loc(i,k) = tmp2(i,1,k)
     dwdn_loc(i,k) = tmp4(i,1,k)
     dpdn_loc(i,k) = tmp5(i,1,k)
  ENDDO; ENDDO
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_2(i0, nt, pl_out, BcsFlowJmin%ref(1,1,5), &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1),&
          buoyancy%vector(2),hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_3(iflag_min, idir, nt, pl_out, pl_inf, inf_rhs, BcsFlowJmin%ref, &
          BcsFlowJmin%ref(1,1,inb_flow+1), &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1),&
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1), &
          buoyancy%vector(2), hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
! add transverse terms
     CALL BOUNDARY_BCS_FLOW_NR_4(iflag_min, idir, nt, BcsFlowJmin%ctan, &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          tmin(1,1,1), tmin(1,1,3), tmin(1,1,2), tmin(1,1,4), tmin(1,1,5), &
          lmin(1,1,1), lmin(1,1,5), &
          hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ENDIF
  DO k = 1,kmax; DO i = imin_loc,imax_loc
     h0(i,1,k) = h0(i,1,k) + hr_loc(i,k)
     h1(i,1,k) = h1(i,1,k) + hv_loc(i,k)
     h2(i,1,k) = h2(i,1,k) + hu_loc(i,k)
     h3(i,1,k) = h3(i,1,k) + hw_loc(i,k)
     h4(i,1,k) = h4(i,1,k) + he_loc(i,k)*prefactor
  ENDDO; ENDDO
  IF ( imixture .GT. 0 ) THEN
     DO k = 1,kmax; DO i = imin_loc,imax_loc
!        h4(i,1,k) = h4(i,1,k) + hr_loc(i,k)*THERMO_AI(6,1,NSP)
        h4(i,1,k) = h4(i,1,k) + hr_loc(i,k)*THERMO_AI(6,1,inb_scal+1)
     ENDDO; ENDDO
  ENDIF

! -------------------------------------------------------------------
! BCs at ymax
! -------------------------------------------------------------------
  DO k = 1,kmax; DO i = 1,imax
     r_loc(i,k) =  rho(i,jmax,k)
     u_loc(i,k) =    v(i,jmax,k)
     v_loc(i,k) =    u(i,jmax,k)
     w_loc(i,k) =    w(i,jmax,k)
     p_loc(i,k) =    p(i,jmax,k)
     g_loc(i,k) = gama(i,jmax,k)
     drdn_loc(i,k) = tmp1(i,jmax,k)
     dudn_loc(i,k) = tmp3(i,jmax,k)
     dvdn_loc(i,k) = tmp2(i,jmax,k)
     dwdn_loc(i,k) = tmp4(i,jmax,k)
     dpdn_loc(i,k) = tmp5(i,jmax,k)
  ENDDO; ENDDO
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_2(i1, nt, pl_out, BcsFlowJmax%ref(1,1,5), &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1),&
          buoyancy%vector(2),hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_3(iflag_max, idir, nt, pl_out, pl_inf, inf_rhs, BcsFlowJmax%ref, & 
          BcsFlowJmax%ref(1,1,inb_flow+1), &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1),&
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1), &
          buoyancy%vector(2),hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
! add transverse terms
     CALL BOUNDARY_BCS_FLOW_NR_4(iflag_max, idir, nt, BcsFlowJmax%ctan, &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          tmax(1,1,1), tmax(1,1,3), tmax(1,1,2), tmax(1,1,4), tmax(1,1,5), &
          lmax(1,1,1), lmax(1,1,5), &
          hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ENDIF
  DO k = 1,kmax; DO i = imin_loc,imax_loc
     h0(i,jmax,k) = h0(i,jmax,k) + hr_loc(i,k)
     h1(i,jmax,k) = h1(i,jmax,k) + hv_loc(i,k)
     h2(i,jmax,k) = h2(i,jmax,k) + hu_loc(i,k)
     h3(i,jmax,k) = h3(i,jmax,k) + hw_loc(i,k)
     h4(i,jmax,k) = h4(i,jmax,k) + he_loc(i,k)*prefactor
  ENDDO; ENDDO
  IF ( imixture .GT. 0 ) THEN
     DO k = 1,kmax; DO i = imin_loc,imax_loc
!        h4(i,jmax,k) = h4(i,jmax,k) + hr_loc(i,k)*THERMO_AI(6,1,NSP)
        h4(i,jmax,k) = h4(i,jmax,k) + hr_loc(i,k)*THERMO_AI(6,1,inb_scal+1)
     ENDDO; ENDDO
  ENDIF

! ###################################################################
! Scalar 
! ###################################################################
  IF ( icalc_scal .EQ. 1 ) THEN
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN; inb_scal_loc = inb_scal + 1
     ELSE;                                         inb_scal_loc = inb_scal     ; ENDIF

     DO is = 1,inb_scal_loc
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), z1(1,1,1,is), tmp2, wrk3d, wrk2d,wrk3d)

! -------------------------------------------------------------------
! BCs at ymin
! -------------------------------------------------------------------
        DO k = 1,kmax; DO i = 1,imax
           r_loc(i,k) =  rho(i,1,k)
           u_loc(i,k) =    v(i,1,k)
           z1_loc(i,k)=   z1(i,1,k,is)
           p_loc(i,k) =    p(i,1,k)
           g_loc(i,k) = gama(i,1,k)
           drdn_loc(i,k) = tmp1(i,1,k)
           dudn_loc(i,k) = tmp3(i,1,k)
           dz1dn_loc(i,k)= tmp2(i,1,k)
           dpdn_loc(i,k) = tmp5(i,1,k)
        ENDDO; ENDDO
        CALL BOUNDARY_BCS_SCAL_NR_3(iflag_min, idir, nt, pl_out, pl_inf, &
             inf_rhs, inf_rhs(1,1,5+is), BcsFlowJmin%ref, BcsScalJmin%ref, BcsScalJmin%ref(1,1,inb_scal+1), &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1),&
             drdn_loc(1,1), dudn_loc(1,1), dz1dn_loc(1,1), dpdn_loc(1,1),&
             buoyancy%vector(2), hz1_loc(1,1))
! add transverse terms
        CALL BOUNDARY_BCS_SCAL_NR_4(iflag_min, nt, BcsScalJmin%ctan, &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1), &
             tmin(1,1,1), tmin(1,1,3), tmin(1,1,5), tmin(1,1,5+is), &
             hz1_loc(1,1))
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
        ELSE
           DO k = 1,kmax; DO i = imin_loc,imax_loc
              zh1(i,1,k,is) = zh1(i,1,k,is) + hz1_loc(i,k)
           ENDDO; ENDDO
        ENDIF
        IF ( imixture .GT. 0 ) THEN
! special case
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
              DO k = 1,kmax; DO i = imin_loc,imax_loc                    
                 h4(i,1,k) = h4(i,1,k) + hz1_loc(i,k)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,1))
              ENDDO; ENDDO
! general case
           ELSE
              DO k = 1,kmax; DO i = imin_loc,imax_loc
!                 h4(i,1,k) = h4(i,1,k) + hz1_loc(i,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                 h4(i,1,k) = h4(i,1,k) + hz1_loc(i,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,inb_scal+1))
              ENDDO; ENDDO
           ENDIF
        ENDIF

! -------------------------------------------------------------------
! BCs at ymax
! -------------------------------------------------------------------
        DO k = 1,kmax; DO i = 1,imax
           r_loc(i,k) =  rho(i,jmax,k)
           u_loc(i,k) =    v(i,jmax,k)
           z1_loc(i,k)=   z1(i,jmax,k,is)
           p_loc(i,k) =    p(i,jmax,k)
           g_loc(i,k) = gama(i,jmax,k)
           drdn_loc(i,k) = tmp1(i,jmax,k)
           dudn_loc(i,k) = tmp3(i,jmax,k)
           dz1dn_loc(i,k)= tmp2(i,jmax,k)
           dpdn_loc(i,k) = tmp5(i,jmax,k)
        ENDDO; ENDDO
        CALL BOUNDARY_BCS_SCAL_NR_3(iflag_max, idir, nt, pl_out, pl_inf, &
             inf_rhs, inf_rhs(1,1,5+is), BcsFlowJmax%ref, BcsScalJmax%ref, BcsScalJmax%ref(1,1,inb_scal+1), &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1),&
             drdn_loc(1,1), dudn_loc(1,1), dz1dn_loc(1,1), dpdn_loc(1,1),&
             buoyancy%vector(2), hz1_loc(1,1))
! add transverse terms
        CALL BOUNDARY_BCS_SCAL_NR_4(iflag_max, nt, BcsScalJmax%ctan, &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1), &
             tmax(1,1,1), tmax(1,1,3), tmax(1,1,5), tmax(1,1,5+is), &
             hz1_loc(1,1))
! special case affects only energy equation
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
        ELSE
           DO k = 1,kmax; DO i = imin_loc,imax_loc
              zh1(i,jmax,k,is) = zh1(i,jmax,k,is) + hz1_loc(i,k)
           ENDDO; ENDDO
        ENDIF
        IF ( imixture .GT. 0 ) THEN
! special case
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
              DO k = 1,kmax; DO i = imin_loc,imax_loc
                 h4(i,jmax,k) = h4(i,jmax,k) + hz1_loc(i,k)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,1))
              ENDDO; ENDDO
! general case
           ELSE
              DO k = 1,kmax; DO i = imin_loc,imax_loc
!                 h4(i,jmax,k) = h4(i,jmax,k) + hz1_loc(i,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                 h4(i,jmax,k) = h4(i,jmax,k) + hz1_loc(i,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,inb_scal+1))
              ENDDO; ENDDO
           ENDIF
        ENDIF
     
     ENDDO

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_Y' )
#endif

  RETURN
END SUBROUTINE BOUNDARY_BCS_Y
