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
!#
!########################################################################
!# DESCRIPTION
!# 
!# Nonperiodic characteristic BCs at xmin and xmax
!#
!########################################################################
SUBROUTINE BOUNDARY_BCS_X(itxc, M2_max, etime, rho,u,v,w,p,gama,z1, &
     x_inf,y_inf,z_inf, q_inf,z1_inf, bcs_vi,bcs_vo, h0,h1,h2,h3,h4,zh1,&
     txc, aux2d, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture, gama0, THERMO_AI
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER itxc

  TREAL M2_max, etime

  TREAL, DIMENSION(imax,jmax,kmax)   :: rho, u, v, w, p, gama, h0, h1, h2, h3, h4
  TREAL, DIMENSION(imax,jmax,kmax,*) :: z1, zh1, txc
  TREAL, DIMENSION(*)                :: x_inf, y_inf, z_inf, q_inf, z1_inf
  TREAL, DIMENSION(jmax,kmax,*)      :: aux2d, bcs_vi, bcs_vo
  TREAL, DIMENSION(*)                :: wrk1d, wrk2d, wrk3d

  TARGET aux2d

! -------------------------------------------------------------------
  TINTEGER j, k, is, nt, inb_scal_loc, isize, iflag_min, iflag_max, idir, ip0
  TREAL prefactor, pl_out, pl_inf1, pl_inf2, pl_inf3, dummy

  TREAL, DIMENSION(:,:,:), POINTER :: tmin, mmin, tmax, mmax, inf_rhs

  TREAL dx(1) ! To use old wrappers to calculate derivatives

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_BCS_X' )
#endif

#define hr_loc(j,k)  aux2d(j,k,1)
#define hu_loc(j,k)  aux2d(j,k,2)
#define hv_loc(j,k)  aux2d(j,k,3)
#define hw_loc(j,k)  aux2d(j,k,4)
#define he_loc(j,k)  aux2d(j,k,5)
#define hz1_loc(j,k) aux2d(j,k,6)

#define r_loc(j,k)   aux2d(j,k,7)
#define u_loc(j,k)   aux2d(j,k,8)
#define v_loc(j,k)   aux2d(j,k,9)
#define w_loc(j,k)   aux2d(j,k,10)
#define p_loc(j,k)   aux2d(j,k,11)
#define g_loc(j,k)   aux2d(j,k,12)
#define z1_loc(j,k)  aux2d(j,k,13)

#define drdn_loc(j,k)  aux2d(j,k,14)
#define dudn_loc(j,k)  aux2d(j,k,15)
#define dvdn_loc(j,k)  aux2d(j,k,16)
#define dwdn_loc(j,k)  aux2d(j,k,17)
#define dpdn_loc(j,k)  aux2d(j,k,18)
#define dz1dn_loc(j,k) aux2d(j,k,19)

  ip0 = 19

  nt = jmax*kmax
  prefactor = (gama0-C_1_R)*mach*mach

  IF ( itxc .LT. nt*(19+5*(inb_flow+inb_scal_array)) .OR. &
       itxc .LT. imax_inf*jmax_inf*kmax_inf ) THEN
     CALL IO_WRITE_ASCII(efile, 'BOUNDARY_BCS_X. Not enough space in txc.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF

! Define pointers
  inf_rhs => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  tmin    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  mmin    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  tmax    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )
  ip0 = ip0 + inb_flow + inb_scal_array
  mmax    => aux2d(:,:,ip0+1:ip0+ inb_flow + inb_scal_array )

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
  idir = 1
  IF      ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN ! not used
  ELSE IF ( imode_sim .EQ. DNS_MODE_SPATIAL  ) THEN; iflag_min =-4; iflag_max = 3; ENDIF

  IF ( bcs_euler_drift .EQ. 1) THEN
     pl_out = bcs_sigma_out*(C_1_R-M2_max)/g(1)%scale
     pl_inf1 = bcs_sigma_inf_imin *qbg(1)%mean/qbg(1)%diam ! jet inflow region (dimensions 1/time)
     pl_inf2 = bcs_sigma_inf_imax/g(1)%scale        ! Lx plane
     pl_inf3 = bcs_sigma_inf_j  /g(2)%scale        ! x0 plane far from jet inflow region
  ELSE
     pl_out  = C_0_R
     pl_inf1 = C_0_R
     pl_inf2 = C_0_R
     pl_inf3 = C_0_R
  ENDIF

! ###################################################################
! forcing terms in array inf_rhs
  IF ( bcs_euler_imin .EQ. DNS_BCS_INFLOW ) THEN
     isize = inb_flow + inb_scal_array
     inf_rhs(:,:,isize) = C_0_R

     IF     ( ifrc_mode .EQ. 1 ) THEN
        CALL BOUNDARY_INFLOW_DISCRETE(etime, inf_rhs)
     ELSEIF ( ifrc_mode .EQ. 2 .OR. ifrc_mode .EQ. 3 ) THEN
        CALL BOUNDARY_INFLOW_BROADBAND&
             (etime, inf_rhs, x_inf,y_inf,z_inf, q_inf,z1_inf, txc, wrk1d,wrk2d,wrk3d)
! not yet developed
!     ELSEIF ( ifrc_mode .EQ. 4 ) THEN
!        CALL BOUNDARY_INFLOW_DISCRETE(etime, inf_rhs)
!        CALL BOUNDARY_INFLOW_BROADBAND&
!             (etime, inf_rhs, x_inf, q_inf, z1_inf, tmp1, wrk1d, wrk2d, wrk3d)
     ENDIF
  ENDIF

! ###################################################################
! Transverse terms
! ###################################################################
  CALL BOUNDARY_BCS_TRANSVERSE_X(u,v,w,p,rho,gama, z1, &
       tmin, mmin, tmax, mmax, txc(1,1,1,1), txc(1,1,1,2), txc(1,1,1,3), wrk1d, wrk2d, wrk3d)

! ###################################################################
! Flow 
! ###################################################################
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u,   txc(1,1,1,2), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, p,   txc(1,1,1,5), i0,i0, wrk1d,wrk2d,wrk3d)         
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, rho, txc(1,1,1,1), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v,   txc(1,1,1,3), i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, w,   txc(1,1,1,4), i0,i0, wrk1d,wrk2d,wrk3d)

! -------------------------------------------------------------------
! Nonreflective BCs at xmin
! -------------------------------------------------------------------
  DO k = 1,kmax
     DO j = 1,jmax
        r_loc(j,k) =  rho(1,j,k)
        u_loc(j,k) =    u(1,j,k)
        v_loc(j,k) =    v(1,j,k)
        w_loc(j,k) =    w(1,j,k)
        p_loc(j,k) =    p(1,j,k)
        g_loc(j,k) = gama(1,j,k)
        drdn_loc(j,k) = txc(1,j,k,1)
        dudn_loc(j,k) = txc(1,j,k,2)
        dvdn_loc(j,k) = txc(1,j,k,3)
        dwdn_loc(j,k) = txc(1,j,k,4)
        dpdn_loc(j,k) = txc(1,j,k,5)
     ENDDO
  ENDDO
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_2(i0, nt, pl_out, bcs_p_imin,&
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1),&
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1), &
          buoyancy%vector(1),hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_3(iflag_min, idir, nt, pl_inf3,pl_inf1, inf_rhs, bcs_vi,&
          bcs_vi(1,1,inb_vars+1), &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1),&
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1), &
          buoyancy%vector(1),hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
! add transverse terms
     CALL BOUNDARY_BCS_FLOW_NR_4(iflag_min, idir, nt, bcs_sigma_trans, &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          tmin(1,1,1), tmin(1,1,2), tmin(1,1,3), tmin(1,1,4), tmin(1,1,5), &
          mmin(1,1,1), mmin(1,1,5), &
          hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
! edge corrections
     CALL BOUNDARY_BCS_FLOW_NR_EDGE(iflag_min, jmax, kmax, bcs_sigma_trans, &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          mmin(1,1,1), mmin(1,1,2), mmin(1,1,3), mmin(1,1,4), mmin(1,1,5), &
          hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ENDIF
  DO k = 1,kmax
     DO j = 1,jmax
        h0(1,j,k) = h0(1,j,k) + hr_loc(j,k)
        h1(1,j,k) = h1(1,j,k) + hu_loc(j,k)
        h2(1,j,k) = h2(1,j,k) + hv_loc(j,k)
        h3(1,j,k) = h3(1,j,k) + hw_loc(j,k)
        h4(1,j,k) = h4(1,j,k) + he_loc(j,k)*prefactor
     ENDDO
  ENDDO
  IF ( imixture .GT. 0 ) THEN
     DO k = 1,kmax
        DO j = 1,jmax
!           h4(1,j,k) = h4(1,j,k) + hr_loc(j,k)*THERMO_AI(6,1,NSP)
           h4(1,j,k) = h4(1,j,k) + hr_loc(j,k)*THERMO_AI(6,1,inb_scal+1)
        ENDDO
     ENDDO
  ENDIF

! -------------------------------------------------------------------
! Nonreflective BCs at xmax
! -------------------------------------------------------------------
  DO k = 1,kmax
     DO j = 1,jmax
        r_loc(j,k) =  rho(imax,j,k)
        u_loc(j,k) =    u(imax,j,k)
        v_loc(j,k) =    v(imax,j,k)
        w_loc(j,k) =    w(imax,j,k)
        p_loc(j,k) =    p(imax,j,k)
        g_loc(j,k) = gama(imax,j,k)
        drdn_loc(j,k) = txc(imax,j,k,1)
        dudn_loc(j,k) = txc(imax,j,k,2)
        dvdn_loc(j,k) = txc(imax,j,k,3)
        dwdn_loc(j,k) = txc(imax,j,k,4)
        dpdn_loc(j,k) = txc(imax,j,k,5)
     ENDDO
  ENDDO
  IF      ( imode_eqns .EQ. DNS_EQNS_TOTAL    ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_2(i1, nt, pl_out, bcs_p_imax,&
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1),&
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1), &
          buoyancy%vector(1), hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     CALL BOUNDARY_BCS_FLOW_NR_3(iflag_max, idir, nt, pl_out, pl_inf2, inf_rhs, bcs_vo,&
          dummy, &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1),&
          drdn_loc(1,1), dudn_loc(1,1), dvdn_loc(1,1), dwdn_loc(1,1), dpdn_loc(1,1), &
          buoyancy%vector(1),hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
! add transverse terms
     CALL BOUNDARY_BCS_FLOW_NR_4(iflag_max, idir, nt, bcs_sigma_trans, &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          tmax(1,1,1), tmax(1,1,2), tmax(1,1,3), tmax(1,1,4), tmax(1,1,5), &
          mmax(1,1,1), mmax(1,1,5), &
          hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
! edge corrections
     CALL BOUNDARY_BCS_FLOW_NR_EDGE(iflag_max, jmax, kmax, bcs_sigma_trans, &
          r_loc(1,1), u_loc(1,1), v_loc(1,1), w_loc(1,1), p_loc(1,1), g_loc(1,1), &
          mmax(1,1,1), mmax(1,1,2), mmax(1,1,3), mmax(1,1,4), mmax(1,1,5), &
          hr_loc(1,1), hu_loc(1,1), hv_loc(1,1), hw_loc(1,1), he_loc(1,1))
  ENDIF
  DO k = 1,kmax
     DO j = 1,jmax
        h0(imax,j,k) = h0(imax,j,k) + hr_loc(j,k)
        h1(imax,j,k) = h1(imax,j,k) + hu_loc(j,k)
        h2(imax,j,k) = h2(imax,j,k) + hv_loc(j,k)
        h3(imax,j,k) = h3(imax,j,k) + hw_loc(j,k)
        h4(imax,j,k) = h4(imax,j,k) + he_loc(j,k)*prefactor
     ENDDO
  ENDDO
  IF ( imixture .GT. 0 ) THEN
     DO k = 1,kmax
        DO j = 1,jmax
!           h4(imax,j,k) = h4(imax,j,k) + hr_loc(j,k)*THERMO_AI(6,1,NSP)
           h4(imax,j,k) = h4(imax,j,k) + hr_loc(j,k)*THERMO_AI(6,1,inb_scal+1)
        ENDDO
     ENDDO
  ENDIF

! ###################################################################
! Scalar 
! ###################################################################
  IF ( icalc_scal .EQ. 1 ) THEN
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN; inb_scal_loc = inb_scal + 1
     ELSE;                                         inb_scal_loc = inb_scal;    ENDIF
     DO is = 1,inb_scal_loc
        CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, z1(1,1,1,is), txc(1,1,1,3), i0,i0, wrk1d,wrk2d,wrk3d)

! -------------------------------------------------------------------
! Nonreflective BCs at xmin
! -------------------------------------------------------------------
        DO k = 1,kmax
           DO j = 1,jmax
              r_loc(j,k) =  rho(1,j,k)
              u_loc(j,k) =    u(1,j,k)
              z1_loc(j,k)=   z1(1,j,k,is)
              p_loc(j,k) =    p(1,j,k)
              g_loc(j,k) = gama(1,j,k)
              drdn_loc(j,k) = txc(1,j,k,1)
              dudn_loc(j,k) = txc(1,j,k,2)
              dz1dn_loc(j,k)= txc(1,j,k,3)
              dpdn_loc(j,k) = txc(1,j,k,5)
           ENDDO
        ENDDO
        CALL BOUNDARY_BCS_SCAL_NR_3(iflag_min, idir, nt, pl_inf3, pl_inf1, &
             inf_rhs, inf_rhs(1,1,5+is), bcs_vi, bcs_vi(1,1,5+is), bcs_vi(1,1,inb_vars+1), &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1),&
             drdn_loc(1,1), dudn_loc(1,1), dz1dn_loc(1,1), dpdn_loc(1,1),&
             buoyancy%vector(1), hz1_loc(1,1))
! add transverse terms
        CALL BOUNDARY_BCS_SCAL_NR_4(iflag_min, nt, bcs_sigma_trans, &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1), &
             tmin(1,1,1), tmin(1,1,2), tmin(1,1,5), tmin(1,1,5+is), &
             hz1_loc(1,1))
! edge corrections
        CALL BOUNDARY_BCS_SCAL_NR_EDGE(iflag_min, jmax, kmax, bcs_sigma_trans, &
             r_loc(1,1), u_loc(1,1), v_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1), &
             mmin(1,1,1), mmin(1,1,2), mmin(1,1,3), mmin(1,1,5), mmin(1,1,5+is), hz1_loc(1,1))
! special case affects only energy equation
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
        ELSE
           DO k = 1,kmax
              DO j = 1,jmax
                 zh1(1,j,k,is) = zh1(1,j,k,is) + hz1_loc(j,k)
              ENDDO
           ENDDO
        ENDIF
        IF ( imixture .GT. 0 ) THEN
! special case
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
              DO k = 1,kmax
                 DO j = 1,jmax
                    h4(1,j,k) = h4(1,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,1))
                 ENDDO
              ENDDO
! general case
           ELSE
              DO k = 1,kmax
                 DO j = 1,jmax
!                    h4(1,j,k) = h4(1,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                    h4(1,j,k) = h4(1,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,inb_scal+1))
                 ENDDO
              ENDDO
           ENDIF
        ENDIF

! -------------------------------------------------------------------
! Nonreflective BCs at xmax
! -------------------------------------------------------------------
        DO k = 1,kmax
           DO j = 1,jmax
              r_loc(j,k) =  rho(imax,j,k)
              u_loc(j,k) =    u(imax,j,k)
              z1_loc(j,k) =  z1(imax,j,k,1)
              p_loc(j,k) =    p(imax,j,k)
              g_loc(j,k) = gama(imax,j,k)
              drdn_loc(j,k) = txc(imax,j,k,1)
              dudn_loc(j,k) = txc(imax,j,k,2)
              dz1dn_loc(j,k)= txc(imax,j,k,3)
              dpdn_loc(j,k) = txc(imax,j,k,5)
           ENDDO
        ENDDO
        CALL BOUNDARY_BCS_SCAL_NR_3(iflag_max, idir, nt, pl_out, pl_inf2, &
             inf_rhs, inf_rhs(1,1,5+is), bcs_vo, bcs_vo(1,1,5+is), dummy, &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1),&
             drdn_loc(1,1), dudn_loc(1,1), dz1dn_loc(1,1), dpdn_loc(1,1),&
             buoyancy%vector(1), hz1_loc(1,1))
! add transverse terms
        CALL BOUNDARY_BCS_SCAL_NR_4(iflag_max, nt, bcs_sigma_trans, &
             r_loc(1,1), u_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1), &
             tmax(1,1,1), tmax(1,1,2), tmax(1,1,5), tmax(1,1,5+is), &
             hz1_loc(1,1))
! edge corrections
        CALL BOUNDARY_BCS_SCAL_NR_EDGE(iflag_max, jmax, kmax, bcs_sigma_trans, &
             r_loc(1,1), u_loc(1,1), v_loc(1,1), z1_loc(1,1), p_loc(1,1), g_loc(1,1), &
             mmax(1,1,1), mmax(1,1,2), mmax(1,1,3), mmax(1,1,5), mmax(1,1,5+is), hz1_loc(1,1))
! special case affects only energy equation
        IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
        ELSE
           DO k = 1,kmax
              DO j = 1,jmax
                 zh1(imax,j,k,is) = zh1(imax,j,k,is) + hz1_loc(j,k)
              ENDDO
           ENDDO
        ENDIF
        IF ( imixture .GT. 0 ) THEN
! special case
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. is .EQ. 2 ) THEN
              DO k = 1,kmax
                 DO j = 1,jmax
                    h4(imax,j,k) = h4(imax,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,3)-THERMO_AI(6,1,1))
                 ENDDO
              ENDDO
! general case
           ELSE
              DO k = 1,kmax
                 DO j = 1,jmax
!                    h4(imax,j,k) = h4(imax,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,NSP))
                    h4(imax,j,k) = h4(imax,j,k) + hz1_loc(j,k)*(THERMO_AI(6,1,is)-THERMO_AI(6,1,inb_scal+1))
                 ENDDO
              ENDDO
           ENDIF
        ENDIF

     ENDDO

  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_BCS_X' )
#endif

  RETURN
END SUBROUTINE BOUNDARY_BCS_X
