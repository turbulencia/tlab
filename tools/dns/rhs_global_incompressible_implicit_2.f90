#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2012/07/10 - C. Ansorge
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Derived from RHS_GLOBAL_INCOMPRESSIBLE_1
!#
!# Momentum equations, nonlinear term in convective form and the 
!# viscous term semi-implicit: 9 1st order derivatives.
!# The explicit laplacian is eliminated by augmenting the array on which 
!# the Helmholtz equation is solved. 
!# 
!# Pressure term requires 3 1st order derivatives
!# BCs need 2 1st order derivatives in Oy
!# Scalar needed for the buoyancy term
!#
!########################################################################
SUBROUTINE  RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2&
     (dte, kex,kim,kco, q,hq, u,v,w,h1,h2,h3, s,hs,&
     tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7, &
     bcs_hb,bcs_ht, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : isize_field, isize_txc_field, isize_wrk1d,  inb_flow,inb_scal
  USE DNS_GLOBAL, ONLY : icalc_scal
  USE DNS_GLOBAL, ONLY : visc, schmidt
  USE DNS_LOCAL,  ONLY : bcs_flow_jmin, bcs_flow_jmax
  USE BOUNDARY_BUFFER

  IMPLICIT NONE

#include "integers.h"

  TREAL dte, kex, kim, kco 
  TREAL, DIMENSION(isize_field,*)                     :: q,hq
  TREAL, DIMENSION(isize_field),         INTENT(INOUT):: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,inb_scal),INTENT(INOUT):: s,hs 
  TREAL, DIMENSION(isize_txc_field),     INTENT(OUT)  :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7
  TREAL, DIMENSION(isize_wrk1d,*),       INTENT(OUT)  :: wrk1d
  TREAL, DIMENSION(*),                   INTENT(OUT)  :: wrk2d,wrk3d
  TREAL, DIMENSION(imax,kmax,*),         INTENT(OUT)  :: bcs_hb, bcs_ht

  TARGET v

! -----------------------------------------------------------------------
  TINTEGER is,ij, k, nxy, ip_b, ip_t 
  TINTEGER ibc, bcs(2,2)
  TREAL dummy, visc_exp, visc_imp, visc_tot, alpha, beta, kef, aug, vdummy!, u_geo, w_geo

  TREAL, DIMENSION(:), POINTER :: p_bcs
  TREAL, DIMENSION(imax,kmax,4):: bcs_locb, bcs_loct    

! #######################################################################
  nxy   = imax*jmax
  
  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero
  
  visc_exp = kex*visc
  visc_imp = kim*visc 
  visc_tot = visc 

  kef= kex/kim
  aug = C_1_R + kef 
  alpha= dte*visc_imp

!  vdummy= -alpha*visc_tot*dte
!  vdummy= -alpha*visc_exp*dte
  vdummy = C_0_R 

! ################################################################################
! SAVE VALUES AT BOUNDARIES. Recovered at the end of the routine
! ################################################################################
  ip_b =                 1 ! bottom
  DO k = 1,kmax
     bcs_hb(1:imax,k,1) = u(ip_b:ip_b+imax-1)
     bcs_hb(1:imax,k,2) = w(ip_b:ip_b+imax-1)
     IF ( icalc_scal .NE. 0 ) THEN 
        DO is=1,inb_scal 
           bcs_hb(1:imax,k,inb_flow+is) = s(ip_b:ip_b+imax-1,is)
        ENDDO
     ENDIF
     ip_b = ip_b + nxy
  ENDDO

  ip_t = imax*(jmax-1) + 1 ! top
  DO k = 1,kmax
     bcs_ht(1:imax,k,1) = u(ip_t:ip_t+imax-1)
     bcs_ht(1:imax,k,2) = w(ip_t:ip_t+imax-1)
     IF ( icalc_scal .NE. 0 ) THEN 
        DO is=1,inb_scal 
           bcs_ht(1:imax,k,inb_flow+is) = s(ip_t:ip_t+imax-1,is)
        ENDDO
     ENDIF
     ip_t = ip_t + nxy
  ENDDO

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  tmp6(1:isize_field) = h3(1:isize_field) ! SAVE old tendencies until end of implicit substep

! h3 contains explicit nonlinear w-tendency (nonlinear operator N(u_n))
  IF ( g(3)%size .GT. 1 ) THEN
     CALL OPR_PARTIAL_X(OPR_P1,    imax,jmax,kmax, bcs, g(1), w, h3,   wrk3d, wrk2d,wrk3d) 
     CALL OPR_PARTIAL_Z(OPR_P1,    imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d) 
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), w, tmp1, tmp2,  wrk2d,wrk3d) ! tmp1 used for BCs below

     DO ij=1,isize_field  
        h3(ij) = -u(ij)*h3(ij) -v(ij)*tmp2(ij) - w(ij)*tmp3(ij)
     ENDDO

! -----------------------------------------------------------------------
! Diffusion part of the BCs for intermediate velocity
! -----------------------------------------------------------------------
     ip_b =                 1
     ip_t = imax*(jmax-1) + 1
     DO k = 1,kmax 
        bcs_locb(:,k,3) = vdummy*tmp1(ip_b:); ip_b = ip_b + nxy ! bottom
        bcs_loct(:,k,3) = vdummy*tmp1(ip_t:); ip_t = ip_t + nxy ! top
     ENDDO

  ENDIF

! #######################################################################

! #######################################################################
  tmp4(1:isize_field) = h1(1:isize_field) ! SAVE old tendencies until end of implicit substep

! h1 contains explicit nonlinear u-tendency  
  CALL OPR_PARTIAL_X(OPR_P1,    imax,jmax,kmax, bcs, g(1), u, h1,   wrk3d, wrk2d,wrk3d) 
  CALL OPR_PARTIAL_Z(OPR_P1,    imax,jmax,kmax, bcs, g(3), u, tmp3, wrk3d, wrk2d,wrk3d) 
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), u, tmp1, tmp2,  wrk2d,wrk3d) ! tmp1 used for BCs below

  DO ij=1,isize_field 
     h1(ij) = -u(ij)*h1(ij) -v(ij)*tmp2(ij) -w(ij)*tmp3(ij) 
  ENDDO

! -----------------------------------------------------------------------
! Diffusion part of the BCs for intermediate velocity
! -----------------------------------------------------------------------
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax 
     bcs_locb(:,k,1) = vdummy*tmp1(ip_b:); ip_b = ip_b + nxy ! bottom
     bcs_loct(:,k,1) = vdummy*tmp1(ip_t:); ip_t = ip_t + nxy ! top
  ENDDO
  
! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
  tmp5(1:isize_field) = h2(1:isize_field) ! SAVE old tendencies until end of implicit substep

! h2 contains explicit nonlinear v-tendency 
  CALL OPR_PARTIAL_X(OPR_P1,    imax,jmax,kmax, bcs, g(1), v, h2,   wrk3d, wrk2d,wrk3d) 
  CALL OPR_PARTIAL_Z(OPR_P1,    imax,jmax,kmax, bcs, g(3), v, tmp3, wrk3d, wrk2d,wrk3d) 
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), v, tmp1, tmp2,  wrk2d,wrk3d)   

  DO ij = 1,isize_field
     h2(ij) = -u(ij)*h2(ij) -v(ij)*tmp2(ij) -w(ij)*tmp3(ij)
  ENDDO

! -----------------------------------------------------------------------
! Diffusion part of the BCs for intermediate velocity
! -----------------------------------------------------------------------
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax 
     bcs_locb(:,k,2) = vdummy*tmp1(ip_b:); ip_b = ip_b + nxy ! bottom
     bcs_loct(:,k,2) = vdummy*tmp1(ip_t:); ip_t = ip_t + nxy ! top
  ENDDO

! #######################################################################
  CALL FI_SOURCES_FLOW(q,s, hq, tmp1, wrk1d,wrk2d,wrk3d)
  
! #######################################################################
! Impose buffer zone as relaxation terms (Flow)
! #######################################################################
  IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(q,hq) 
  ENDIF

! #######################################################################
! Explicit part of time integration for velocities
! #######################################################################
! tmp4-6 contain tendencies, h1-3 contain old tendencies
  DO ij=1,isize_field 
     tmp1(ij) = u(ij)*aug + dte*( h1(ij) + kco*tmp4(ij) );  
     tmp2(ij) = v(ij)*aug + dte*( h2(ij) + kco*tmp5(ij) ); 
     tmp3(ij) = w(ij)*aug + dte*( h3(ij) + kco*tmp6(ij) );
  ENDDO
! tmp1-tmp3:  u,v,w updated with explicit scheme 
! h1-h3:      explicit tendencies from this substep 
! u-w:        old velocities 

! -----------------------------------------------------------------------
! Remaining part of the BCs for intermediate velocity
! -----------------------------------------------------------------------
  ! This block implements the BCs retaining a dt correction according to the
  ! notes but it leads to first-order convergence rate instead of second-order 
  ! and even instability if the diffusion terms are retained.
  !
  ! ip_b =                 1
  ! ip_t = imax*(jmax-1) + 1
  ! DO k = 1,kmax 
  !    bcs_locb(:,k,1) = bcs_locb(:,k,1) - alpha*tmp1(ip_b:) ! bottom
  !    bcs_loct(:,k,1) = bcs_loct(:,k,1) - alpha*tmp1(ip_t:) ! top

  !    bcs_locb(:,k,2) = bcs_locb(:,k,2) - alpha*tmp2(ip_b:) ! bottom
  !    bcs_loct(:,k,2) = bcs_loct(:,k,2) - alpha*tmp2(ip_t:) ! top

  !    bcs_locb(:,k,3) = bcs_locb(:,k,3) - alpha*tmp3(ip_b:) ! bottom
  !    bcs_loct(:,k,3) = bcs_loct(:,k,3) - alpha*tmp3(ip_t:) ! top

  !    ip_b = ip_b + nxy
  !    ip_t = ip_t + nxy
  ! ENDDO

  ! This block retains only the leading order terms and provides second-order
  ! convergence rate.
  bcs_locb(:,:,1) = -alpha*aug*bcs_hb(:,:,1) ! bottom
  bcs_loct(:,:,1) = -alpha*aug*bcs_ht(:,:,1) ! top
  bcs_locb(:,:,2) = C_0_R                    ! bottom
  bcs_loct(:,:,2) = C_0_R                    ! top
  bcs_locb(:,:,3) = -alpha*aug*bcs_hb(:,:,2) ! bottom
  bcs_loct(:,:,3) = -alpha*aug*bcs_ht(:,:,2) ! top

! ################################################################################
! ADVECTION-DIFFUSION FOR SCALAR
! We need this code segment here because we need the old velocities  
! ################################################################################
! hs  -> explicit tendency without diffusion 
  IF ( icalc_scal .NE. 0 ) THEN 

     DO is = 1,inb_scal 

        tmp7(1:isize_field) = hs(1:isize_field,is) ! save old values

        CALL OPR_PARTIAL_X(OPR_P1,    imax,jmax,kmax, bcs, g(1), s(1,is), hs(1,is), wrk3d, wrk2d,wrk3d) 
        CALL OPR_PARTIAL_Z(OPR_P1,    imax,jmax,kmax, bcs, g(3), s(1,is), tmp6,     wrk3d, wrk2d,wrk3d) 
        CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,is), tmp4,     tmp5,  wrk2d,wrk3d)

        DO ij = 1,isize_field; 
          hs(ij,is) = -u(ij)*hs(ij,is) -v(ij)*tmp5(ij) -w(ij)*tmp6(ij) 
        ENDDO

! -----------------------------------------------------------------------
! Diffusion part of the BCs for intermediate scalar
! -----------------------------------------------------------------------
        ip_b =                 1
        ip_t = imax*(jmax-1) + 1
        DO k = 1,kmax 
           bcs_locb(:,k,4) = (vdummy/schmidt(is)) *tmp4(ip_b:); ip_b = ip_b + nxy ! bottom
           bcs_loct(:,k,4) = (vdummy/schmidt(is)) *tmp4(ip_t:); ip_t = ip_t + nxy ! top
        ENDDO
        
! #######################################################################
! Impose buffer zone as relaxation terms (Scalar #is) 
! #######################################################################
        IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
           CALL BOUNDARY_BUFFER_RELAXATION_SCAL(is, q,s(1,is), hs(1,is) ) 
        ENDIF
        
! #######################################################################
! Explicit Time Stepping for scalar #is        
! #######################################################################
        DO ij = 1,isize_field
           tmp4(ij) = s(ij,is)*aug + dte*( hs(ij,is) + kco*tmp7(ij) );  
        ENDDO

! -----------------------------------------------------------------------
! Remaining part of the BCs for intermediate scalar
! -----------------------------------------------------------------------
        ! See velocity part
        ! ip_b =                 1
        ! ip_t = imax*(jmax-1) + 1
        ! DO k = 1,kmax 
        !    bcs_locb(:,k,4) = bcs_locb(:,k,4) - (alpha/schmidt(is))*tmp4(ip_b:); ip_b = ip_b + nxy ! bottom
        !    bcs_loct(:,k,4) = bcs_loct(:,k,4) - (alpha/schmidt(is))*tmp4(ip_t:); ip_t = ip_t + nxy ! top
        ! ENDDO
        
        bcs_locb(:,:,4) = -(alpha/schmidt(is))*aug*bcs_hb(:,:,inb_flow+is)
        bcs_loct(:,:,4) = -(alpha/schmidt(is))*aug*bcs_ht(:,:,inb_flow+is)

! #######################################################################
! semi-Implicit Diffusion for Scalar #is 
! #######################################################################
        beta =-C_1_R/(alpha/schmidt(is))
 
        CALL OPR_HELMHOLTZ_FXZ_2(imax,jmax,kmax, g, i0, beta, &
             tmp4, tmp6,tmp7, bcs_locb(1,1,4), bcs_loct(1,1,4), wrk1d, wrk1d(1,5),wrk3d )

        DO ij = 1,isize_field  
           s(ij,is) = beta*tmp4(ij) - kef*s(ij,is)
        ENDDO

     ENDDO
  ENDIF

! ################################################################################
! Implicit part of time integration for velocities
! ################################################################################ 
  beta =-C_1_R/alpha 

  CALL OPR_HELMHOLTZ_FXZ_2(imax,jmax,kmax, g, i0, beta, &
       tmp1, tmp5,tmp6, bcs_locb(1,1,1), bcs_loct(1,1,1), wrk1d, wrk1d(1,5),wrk3d)
  CALL OPR_HELMHOLTZ_FXZ_2(imax,jmax,kmax, g, i0, beta, &
       tmp2, tmp5,tmp6, bcs_locb(1,1,2), bcs_loct(1,1,2), wrk1d, wrk1d(1,5),wrk3d)
  CALL OPR_HELMHOLTZ_FXZ_2(imax,jmax,kmax, g, i0, beta, &
       tmp3, tmp5,tmp6, bcs_locb(1,1,3), bcs_loct(1,1,3), wrk1d, wrk1d(1,5),wrk3d) 

  DO ij=1,isize_field  
     u(ij) = beta*tmp1(ij) -kef*u(ij) 
     v(ij) = beta*tmp2(ij) -kef*v(ij) 
     w(ij) = beta*tmp3(ij) -kef*w(ij)
  ENDDO

! #######################################################################
! Pressure term (actually solving for dte*p) 
! #######################################################################
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
     
! -----------------------------------------------------------------------
  DO ij = 1,isize_field
     tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
  ENDDO

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     p_bcs => v(ip_b:); bcs_hb(1:imax,k,3) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => v(ip_t:); bcs_ht(1:imax,k,3) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! pressure in tmp1, Oy derivative in tmp3
  ibc = 3
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, ibc, &
       tmp1,tmp3, tmp2,tmp4, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp1, tmp4, wrk3d, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient 
! -----------------------------------------------------------------------
  dummy = C_1_R/dte
  DO ij = 1,isize_field
     u(ij) = u(ij) - tmp2(ij);      h1(ij) = h1(ij) - dummy*tmp2(ij) 
     v(ij) = v(ij) - tmp3(ij);      h2(ij) = h2(ij) - dummy*tmp3(ij)
     w(ij) = w(ij) - tmp4(ij);      h3(ij) = h3(ij) - dummy*tmp4(ij) 
  ENDDO

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
  ibc = 0
! Default is Dirichlet -> values at boundary are kept const;
! bcs_hb and bcs_ht initialized to old BC values above 
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     ! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux) 
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), u, &
          bcs_hb(1,1,1),bcs_ht(1,1,1), wrk1d,tmp1,wrk3d)
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), w, &
          bcs_hb(1,1,2),bcs_ht(1,1,2), wrk1d,tmp1,wrk3d)
  ENDIF

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin 
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     u(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
     v(ip_b:ip_b+imax-1) = C_0_R               ! no penetration
     w(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2); 
     IF ( icalc_scal .NE. 0 ) THEN 
        DO is=1,inb_scal 
           s(ip_b:ip_b+imax-1,is) = bcs_hb(1:imax,k,inb_flow+is) 
        ENDDO
     ENDIF
     ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     u(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
     v(ip_t:ip_t+imax-1) = C_0_R               ! no penetration
     w(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2); 
     IF ( icalc_scal .NE. 0 ) THEN 
        DO is=1,inb_scal 
           s(ip_t:ip_t+imax-1,is) = bcs_ht(1:imax,k,inb_flow+is) 
        ENDDO
     ENDIF
     ip_t = ip_t + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_2
