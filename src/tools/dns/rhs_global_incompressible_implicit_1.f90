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
!# Derived from RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1
!#
!# Semi-implicit time stepping (implicit for viscous term)
!#
!#
!########################################################################
SUBROUTINE  RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_1&
     (dte, kex,kim,kco, &
     q,hq,u,v,w,h1,h2,h3,s,hs,&
     tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9, &
     wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_CONSTANTS,ONLY: efile
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : imax,jmax,kmax
  USE TLAB_VARS, ONLY : isize_field, isize_txc_field, isize_wrk1d, inb_scal, inb_flow
  USE TLAB_VARS, ONLY : icalc_scal
  USE TLAB_VARS, ONLY : visc, schmidt, rossby
  USE TLAB_VARS, ONLY : buoyancy, coriolis
  USE TLAB_VARS, ONLY : bbackground
  USE TLAB_PROCS
  USE DNS_LOCAL,  ONLY : idivergence
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS

  IMPLICIT NONE

#include "integers.h"

  TREAL dte, kex, kim, kco
  TREAL, DIMENSION(isize_field,*)                     :: q,hq
  TREAL, DIMENSION(isize_field),         INTENT(INOUT):: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,inb_scal),INTENT(INOUT):: s,hs
  TREAL, DIMENSION(isize_txc_field),     INTENT(OUT)  :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  TREAL, DIMENSION(isize_wrk1d,*),       INTENT(OUT)  :: wrk1d
  TREAL, DIMENSION(*),                   INTENT(OUT)  :: wrk2d,wrk3d

  TARGET tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,h2,wrk2d,u,v,w

! -----------------------------------------------------------------------
  TINTEGER iq,is,ij, k, nxy, ip, ip_b, ip_t
  TINTEGER ibc, bcs(2,2)
  TREAL dummy, visc_exp, visc_imp, visc_tot, diff, alpha, beta

  TREAL, DIMENSION(:),        POINTER :: p_bcs

#ifdef USE_BLAS
  INTEGER ilen
#endif

! #######################################################################
  nxy    = imax*jmax

#ifdef USE_BLAS
  ilen = isize_field
#endif

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

  visc_exp = kex*visc
  visc_imp = kim*visc
  visc_tot = visc

! SAVE old tendencies until end of implicit substep of integration

  tmp4(1:isize_field) = hq(:,1)
  tmp5(1:isize_field) = hq(:,2)
  tmp6(1:isize_field) = hq(:,3)

! ################################################################################
! SAVE VALUES AT BOUNDARIES
! ################################################################################

  ip = 1
  DO k=1,kmax
     ip_b = (k-1)*imax*jmax +1;
     ip_t = ip_b + (jmax-1)*imax
     IF ( BcsFlowJmin%type(1) .EQ. DNS_BCS_DIRICHLET ) BcsFlowJmin%ref(1:imax,k,1) = u(ip_b:ip_b+imax-1)
     IF ( BcsFlowJmin%type(3) .EQ. DNS_BCS_DIRICHLET ) BcsFlowJmin%ref(1:imax,k,3) = w(ip_b:ip_b+imax-1)
     IF ( BcsFlowJmax%type(1) .EQ. DNS_BCS_DIRICHLET ) BcsFlowJmax%ref(1:imax,k,1) = u(ip_t:ip_t+imax-1)
     IF ( BcsFlowJmax%type(3) .EQ. DNS_BCS_DIRICHLET ) BcsFlowJmax%ref(1:imax,k,3) = w(ip_t:ip_t+imax-1)
     DO is=1,inb_scal
        IF ( BcsScalJmin%type(is) .EQ. DNS_BCS_DIRICHLET .AND. &
             BcsScalJmax%type(is) .EQ. DNS_BCS_DIRICHLET ) THEN
           BcsScalJmin%ref(1:imax,k,is) = s(ip_b:ip_b+imax-1,is)
           BcsScalJmax%ref(1:imax,k,is) = s(ip_t:ip_t+imax-1,is)
        ELSE  ! Only Dirichlet BCs implemented for scalar
           CALL TLAB_WRITE_ASCII(efile,'Only Dirichlet BCs implemented for scalar in implicit mode')
           CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
        ENDIF

        IF ( BcsScalJmin%SfcType(is) .EQ. DNS_SFC_STATIC .AND. &
             BcsScalJmax%SfcType(is) .EQ. DNS_SFC_STATIC ) THEN
           ! Nothing to do
        ELSE
           CALL TLAB_WRITE_ASCII(efile,'Only static surface implemented in implicit mode')
        ENDIF
     ENDDO
  ENDDO

! ################################################################################
! EXPLICIT PART OF SUBSTEP - change viscosity accordingly
! ################################################################################
  visc = visc_exp

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  ! h3   contains explicit nonlinear w-tendency (nonlinear operator N(u_n))
  ! tmp9 contains explicit diffusive tendency for w
  IF ( g(3)%size .GT. 1 ) THEN
     !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),w,    tmp2,wrk3d, wrk2d,wrk3d)
     !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),tmp2, tmp1,wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), w, tmp1, tmp2, wrk2d,wrk3d)   ! tmp2 used for BCs below
     DO ij = 1,isize_field;
          h3(ij) =           -v(ij)*tmp2(ij)
        tmp9(ij) =             visc*tmp1(ij)
     ENDDO

     !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),w,    tmp3,wrk3d, wrk2d,wrk3d)
     !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),tmp3, tmp1,wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), w, tmp1, tmp3, wrk2d,wrk3d)
     DO ij = 1,isize_field;
          h3(ij) =   h3(ij) -u(ij)*tmp3(ij)
        tmp9(ij) = tmp9(ij) + visc*tmp1(ij)
     ENDDO

     !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),w,    tmp3,wrk3d, wrk2d,wrk3d)
     !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),tmp3, tmp1,wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), w, tmp1, tmp3, wrk2d,wrk3d)
! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
     IF ( buoyancy%active(3) ) THEN
        wrk1d(:,1) = C_0_R
        CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
        dummy = buoyancy%vector(3)
        DO ij = 1,isize_field
             h3(ij) =   h3(ij) + dummy*wrk3d(ij)
        ENDDO
     ENDIF

! -----------------------------------------------------------------------
! Coriolis (so far, rotation only in the Oy direction)
! -----------------------------------------------------------------------
     IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
        dummy = C_1_R/rossby
        DO ij = 1,isize_field
             h3(ij) =   h3(ij) - w(ij)*tmp3(ij) + dummy*(u(ij)-C_1_R)
           tmp9(ij) = tmp9(ij) + visc*tmp1(ij)
        ENDDO
! -----------------------------------------------------------------------
     ELSE
        DO ij = 1,isize_field
             h3(ij) =   h3(ij) - w(ij)*tmp3(ij)
           tmp9(ij) = tmp9(ij) + visc*tmp1(ij)
        ENDDO
     ENDIF
  ENDIF

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! #######################################################################
!   h1 contains explicit nonlinear u-tendency
! tmp7 contains explicit diffusive u-tendency
  !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),u,    tmp2,wrk3d, wrk2d,wrk3d)
  !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),tmp2, tmp1,wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), u, tmp1, tmp2, wrk2d,wrk3d)   ! tmp2 used for BCs below
  DO ij=1,isize_field;
       h1(ij)  =           - v(ij)*tmp2(ij)
     tmp7(ij)  =              visc*tmp1(ij)
  ENDDO

  !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),u,    tmp3,wrk3d, wrk2d,wrk3d)
  !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),tmp3, tmp1,wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, tmp3, wrk2d,wrk3d)
  DO ij=1,isize_field;
       h1(ij) =   h1(ij) - u(ij)*tmp3(ij);
     tmp7(ij) = tmp7(ij) + visc*tmp1(ij)
  ENDDO

  !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),u,    tmp3,wrk3d, wrk2d,wrk3d)
  !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),tmp3, tmp1,wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), u, tmp1, tmp3, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( buoyancy%active(1) ) THEN
     wrk1d(:,1) = C_0_R
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
     dummy = buoyancy%vector(1)
     DO ij = 1,isize_field
          h1(ij) =   h1(ij) + dummy*wrk3d(ij)
     ENDDO
  ENDIF

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction.
! -----------------------------------------------------------------------

  IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
     dummy = C_1_R/rossby
     DO ij = 1, isize_field
          h1(ij) =   h1(ij) - w(ij)*tmp3(ij) - dummy*w(ij)
        tmp7(ij) = tmp7(ij) + visc*tmp1(ij)
     ENDDO
! -----------------------------------------------------------------------
  ELSE
     DO ij = 1, isize_field
          h1(ij) =   h1(ij) - w(ij)*tmp3(ij)
        tmp7(ij) = tmp7(ij) + visc*tmp1(ij)
     ENDDO
  ENDIF

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
!   h2 contains explicit nonlinear v-tendency
! tmp8 contains explicit diffusive v-tendency
  !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),v,    tmp2,wrk3d, wrk2d,wrk3d)
  !CALL OPR_PARTIAL_Y(OPR_P1,imax,jmax,kmax,bcs,g(2),tmp2, tmp1,wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), v, tmp1, tmp2, wrk2d,wrk3d)
  DO ij=1,isize_field
       h2(ij) =          - v(ij)*tmp2(ij)
     tmp8(ij) =             visc*tmp1(ij)
  ENDDO


  !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),v,    tmp3,wrk3d, wrk2d,wrk3d)
  !CALL OPR_PARTIAL_X(OPR_P1,imax,jmax,kmax,bcs,g(1),tmp3, tmp1,wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), v, tmp1, tmp3, wrk2d,wrk3d)
  DO ij=1,isize_field
       h2(ij) =   h2(ij) - u(ij)*tmp3(ij)
     tmp8(ij) = tmp8(ij) + visc*tmp1(ij)
  ENDDO


  !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),v,    tmp3,wrk3d, wrk2d,wrk3d)
  !CALL OPR_PARTIAL_Z(OPR_P1,imax,jmax,kmax,bcs,g(3),tmp3, tmp1,wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), v, tmp1, tmp3, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( buoyancy%active(2) ) THEN
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, bbackground)
     dummy = buoyancy%vector(2)
     DO ij = 1,isize_field
          h2(ij) =   h2(ij) - w(ij)*tmp3(ij) + dummy*wrk3d(ij)
        tmp8(ij) = tmp8(ij) + visc*tmp1(ij)
     ENDDO
  ELSE
     DO ij = 1,isize_field
          h2(ij) =   h2(ij) - w(ij)*tmp3(ij)
        tmp8(ij) = tmp8(ij) +  visc*tmp1(ij)
     ENDDO
  ENDIF


! #######################################################################
! Impose buffer zone as relaxation terms (Flow)
! #######################################################################
  IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
     CALL BOUNDARY_BUFFER_RELAX_FLOW(q,hq)
  ENDIF

! Explicit part of time integration
! tmp4-6 contain tendencies, h1-3 contain old tendencies
  DO ij=1,isize_field
     tmp1(ij) = u(ij) + dte*( h1(ij) + tmp7(ij) + kco*tmp4(ij) );
     tmp2(ij) = v(ij) + dte*( h2(ij) + tmp8(ij) + kco*tmp5(ij) );
     tmp3(ij) = w(ij) + dte*( h3(ij) + tmp9(ij) + kco*tmp6(ij) );
  ENDDO
! tmp1-tmp3:  u,v,w updated with explicit scheme
! h1-h3:      explicit tendencies from this substep
! u-w:        old velocities
! tmp7-tmp9:  vacant


! ################################################################################
! ADVECTION-DIFFUSION FOR SCALAR
! ################################################################################
! hs  -> explicit tendency without diffusion
! tmp7-> diffusion contribution to explicit tendency
  IF ( icalc_scal .NE. i0 ) THEN
     DO is=1,inb_scal
        diff = visc_exp / schmidt(is)
        tmp6 = hs(:,is)
        !CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s(1,is), tmp4, tmp5, wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), s(1,is), tmp5, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp5,    tmp4, wrk3d, wrk2d,wrk3d)
        DO ij=1,isize_field;
          hs(ij,is) =           -u(ij)*tmp5(ij)
           tmp7(ij) =             diff*tmp4(ij)
        ENDDO

        !CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s(1,is), tmp4, tmp5, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), s(1,is), tmp5, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp5,    tmp4, wrk3d, wrk2d,wrk3d)
        DO ij=1,isize_field;
          hs(ij,is) =hs(ij,is) -v(ij)*tmp5(ij)
           tmp7(ij) = tmp7(ij)  +diff*tmp4(ij)
        ENDDO

        !CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s(1,is), tmp4, tmp5, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), s(1,is), tmp5, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp5,    tmp4, wrk3d, wrk2d,wrk3d)

        DO ij=1,isize_field;
          hs(ij,is) =hs(ij,is) -w(ij)*tmp5(ij)
           tmp7(ij) = tmp7(ij)  +diff*tmp4(ij)
        ENDDO

! #######################################################################
! Impose buffer zone as relaxation terms (Scalar #is)
! #######################################################################
        IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
           CALL BOUNDARY_BUFFER_RELAX_SCAL_I(is, s(1,is), hs(1,is) )
        ENDIF


        DO ij=1,isize_field
           s(ij,is) = s(ij,is) + dte*( tmp7(ij) + hs(ij,is) + kco*tmp6(ij) );
        ENDDO
        diff= visc_imp/schmidt(is)
        alpha= dte*diff
        beta =-C_1_R/alpha

        ip_b = 1     ; ip_t = 1+ imax*kmax
        DO k=1,kmax
           ip=(k-1)*imax*jmax+1; wrk2d(ip_b:ip_b+imax-1)=-alpha*BcsScalJmin%ref(1:imax,k,is); ip_b=ip_b+imax
           ip=ip+(jmax-1)*imax;  wrk2d(ip_t:ip_t+imax-1)=-alpha*BcsScalJmax%ref(1:imax,k,is); ip_t=ip_t+imax
        ENDDO
        ip_b = 1     ; ip_t = 1+ imax*kmax
        CALL OPR_HELMHOLTZ_FXZ(imax,jmax,kmax, g, i0, beta, &
             s(1,is), tmp5, tmp6, wrk2d(ip_b), wrk2d(ip_t), wrk1d, wrk1d(1,5),wrk3d )

        DO ij=1,isize_field
           s(ij,is) = s(ij,is) * beta
        ENDDO

     ENDDO

  ENDIF

! ################################################################################
! IMPLICIT PART OF SUBSTEP - change viscosity accordingly
! ################################################################################
! TO BE CLEANED: in this section wrk2d is used for the BCS of the implicit solver

  visc = visc_imp
  alpha= dte*visc_imp
  beta =-C_1_R/alpha

  ip_b = 1; ip_t = 1+ imax*kmax
  DO k=1,kmax
     ip=(k-1)*imax*jmax+1; wrk2d(ip_b:ip_b+imax-1)=-alpha*u(ip:ip+imax-1); ip_b=ip_b+imax
     ip=ip+(jmax-1)*imax;  wrk2d(ip_t:ip_t+imax-1)=-alpha*u(ip:ip+imax-1); ip_t=ip_t+imax
  ENDDO

  ip_b = 1; ip_t = 1+imax*kmax
  CALL OPR_HELMHOLTZ_FXZ(imax,jmax,kmax, g, i0, beta, &
       tmp1, tmp5, tmp6, wrk2d(ip_b), wrk2d(ip_t), wrk1d, wrk1d(:,5),wrk3d )

  DO k=1,kmax
     ip=(k-1)*imax*jmax+1; wrk2d(ip_b:ip_b+imax-1)=-alpha*v(ip:ip+imax-1); ip_b=ip_b+imax
     ip=ip+(jmax-1)*imax;  wrk2d(ip_t:ip_t+imax-1)=-alpha*v(ip:ip+imax-1); ip_t=ip_t+imax
  ENDDO
  ip_b = 1     ; ip_t = 1+ imax*kmax
  CALL OPR_HELMHOLTZ_FXZ(imax,jmax,kmax, g, i0, beta, &
       tmp2, tmp5, tmp6, wrk2d(ip_b), wrk2d(ip_t), wrk1d, wrk1d(:,5),wrk3d )

  DO k=1,kmax
     ip=(k-1)*imax*jmax+1; wrk2d(ip_b:ip_b+imax-1)=-alpha*w(ip:ip+imax-1); ip_b=ip_b+imax
     ip=ip+(jmax-1)*imax;  wrk2d(ip_t:ip_t+imax-1)=-alpha*w(ip:ip+imax-1); ip_t=ip_t+imax
  ENDDO
  ip_b = 1     ; ip_t = 1+ imax*kmax
  CALL OPR_HELMHOLTZ_FXZ(imax,jmax,kmax, g, i0, beta, &
       tmp3, tmp5, tmp6, wrk2d(ip_b), wrk2d(ip_t), wrk1d, wrk1d(:,5),wrk3d )
  DO ij=1,isize_field
     u(ij) = tmp1(ij) * beta
     v(ij) = tmp2(ij) * beta
     w(ij) = tmp3(ij) * beta
  ENDDO

! ################################################################################
! END OF IMPLICIT PART OF SUBSTEP - set viscosity back to actual value
! ################################################################################
  visc = visc_tot


! #######################################################################
! Pressure term
! #######################################################################
  IF ( idivergence .EQ. EQNS_DIVERGENCE ) THEN ! remove residual divergence
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, tmp1, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), v, tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), w, tmp3, wrk3d, wrk2d,wrk3d)
  ENDIF

! -----------------------------------------------------------------------
  DO ij = 1,isize_field
     tmp1(ij) = (tmp1(ij) + tmp2(ij) + tmp3(ij))/dte ! forcing term in tmp1
  ENDDO

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     p_bcs => v(ip_b:); BcsFlowJmin%ref(1:imax,k,2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => v(ip_t:); BcsFlowJmax%ref(1:imax,k,2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! pressure in tmp1, Oy derivative in tmp3
  ibc = 3
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, ibc, &
       tmp1,tmp3, tmp2,tmp4, BcsFlowJmin%ref(1,1,2),BcsFlowJmax%ref(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp1, tmp4, wrk3d, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
  DO ij = 1,isize_field
     u(ij) = u(ij) - dte*tmp2(ij);    h1(ij) = h1(ij)-tmp2(ij)
     v(ij) = v(ij) - dte*tmp3(ij);    h2(ij) = h2(ij)-tmp3(ij)
     w(ij) = w(ij) - dte*tmp4(ij);    h3(ij) = h3(ij)-tmp4(ij)
  ENDDO

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
! Default is Dirichlet -> values at boundary are kept const;
! BcsFlowJmin/Jmax for u and w initialized to old BC values above
  BcsFlowJmin%ref(:,:,2) = C_0_R
  BcsFlowJmax%ref(:,:,2) = C_0_R
  DO iq = 1,inb_flow
     ibc = 0
     IF ( BcsFlowJmin%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( BcsFlowJmax%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( ibc .GT. 0 ) THEN
! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), q(1,iq), &
             BcsFlowJmin%ref(1,1,iq),BcsFlowJmax%ref(1,1,iq), wrk1d,tmp1,wrk3d)
     ENDIF
  ENDDO

  DO is = 1,inb_scal
     ibc = 0
     IF ( BcsScalJmin%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( BcsScalJmax%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( ibc .GT. 0 ) THEN
! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), s(1,is), &
             BcsScalJmin%ref(1,1,is),BcsScalJmax%ref(1,1,is), wrk1d,tmp1,wrk3d)
     ENDIF
  ENDDO

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     u(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,1)
     v(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,2)
     w(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,3)
     DO is = 1,inb_scal
        s(ip_b:ip_b+imax-1,is) = BcsScalJmin%ref(1:imax,k,is)
     ENDDO
     ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     u(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,1)
     v(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,2)
     w(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,3)
     DO is = 1,inb_scal
        s(ip_t:ip_t+imax-1,is) = BcsScalJmax%ref(1:imax,k,is)
     ENDDO
     ip_t = ip_t + nxy
  ENDDO

! ! -----------------------------------------------------------------------
! ! Preliminaries
! ! -----------------------------------------------------------------------

!   ibc = 0
! ! Default is Dirichlet -> values at boundary are kept const;
! ! bcs_hb and bcs_ht initialized to old BC values above
!   IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
!   IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
!   IF ( ibc .GT. 0 ) THEN
!      ! WATCH OUT - THIS IS REALLY GOING TO GIVE NO-FLUX (instead of constant flux)
!      CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), u, &
!           bcs_hb(1,1,1),bcs_ht(1,1,1), wrk1d,tmp1,wrk3d)
!      CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), w, &
!           bcs_hb(1,1,2),bcs_ht(1,1,2), wrk1d,tmp1,wrk3d)
!   ENDIF



! ! -----------------------------------------------------------------------
! ! Impose bottom BCs at Jmin
! ! -----------------------------------------------------------------------
!   ip_b =                 1
!   DO k = 1,kmax
!      u(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
!      v(ip_b:ip_b+imax-1) = C_0_R               ! no penetration
!      w(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2);
!      IF ( icalc_scal .NE. 0 ) THEN
!         DO is=1,inb_scal
!            s(ip_b:ip_b+imax-1,is) = BcsScalJmin%ref(1:imax,k,is)
!         ENDDO
!      ENDIF
!      ip_b = ip_b + nxy
!   ENDDO

! ! -----------------------------------------------------------------------
! ! Impose top BCs at Jmax
! ! -----------------------------------------------------------------------
!   ip_t = imax*(jmax-1) + 1
!   DO k = 1,kmax
!      u(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
!      v(ip_t:ip_t+imax-1) = C_0_R               ! no penetration
!      w(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2);
!      IF ( icalc_scal .NE. 0 ) THEN
!         DO is=1,inb_scal
!            s(ip_t:ip_t+imax-1,is) = BcsScalJmax%ref(1:imax,k,is)
!         ENDDO
!      ENDIF
!      ip_t = ip_t + nxy
!   ENDDO


  RETURN
END SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_1
