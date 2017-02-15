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
     bcs_hb,bcs_ht, vaux, &
     wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL
  USE DNS_LOCAL,  ONLY : bcs_flow_jmin, bcs_flow_jmax, bcs_scal_jmin, bcs_scal_jmax
  USE DNS_LOCAL,  ONLY : idivergence
  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : buff_type

  IMPLICIT NONE

#include "integers.h"

  TREAL dte, kex, kim, kco 
  TREAL, DIMENSION(isize_field,*)                     :: q,hq
  TREAL, DIMENSION(isize_field),         INTENT(INOUT):: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,inb_scal),INTENT(INOUT):: s,hs 
  TREAL, DIMENSION(isize_txc_field),     INTENT(OUT)  :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
  TREAL, DIMENSION(isize_wrk1d,*),       INTENT(OUT)  :: wrk1d
  TREAL, DIMENSION(*),                   INTENT(OUT)  :: wrk2d,wrk3d, vaux
  TREAL, DIMENSION(imax,kmax,*),         INTENT(OUT)  :: bcs_hb, bcs_ht

  TARGET tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,h2,wrk2d,u,v,w

! -----------------------------------------------------------------------
  TINTEGER is,ij, k, nxy, ip, ip_b, ip_t 
  TINTEGER ibc
  TREAL dummy, visc_exp, visc_imp, visc_tot, diff, alpha, beta

  TREAL, DIMENSION(:),        POINTER :: p_bcs
  TREAL, DIMENSION(imax,kmax,inb_scal):: bcs_sb, bcs_st 

#ifdef USE_BLAS
  INTEGER ilen
#endif

  TREAL dx(1), dy(1), dz(1) ! To use old wrappers to calculate derivatives

! #######################################################################
  nxy    = imax*jmax

#ifdef USE_BLAS
  ilen = isize_field
#endif

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
     IF ( bcs_flow_jmin .EQ. DNS_BCS_DIRICHLET ) THEN 
        bcs_hb(:,k,1) = u(ip_b:ip_b+imax-1)  
        bcs_hb(:,k,2) = w(ip_b:ip_b+imax-1) 
     ENDIF
     IF ( bcs_flow_jmax .EQ. DNS_BCS_DIRICHLET ) THEN
        bcs_ht(:,k,1) = u(ip_t:ip_t+imax-1) 
        bcs_ht(:,k,2) = w(ip_t:ip_t+imax-1)  
     ENDIF
     DO is=1,inb_scal
        IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_DIRICHLET .AND. &
             bcs_scal_jmax(is) .EQ. DNS_BCS_DIRICHLET ) THEN 
           bcs_sb(1:imax,k,is) = s(ip_b:ip_b+imax-1,is) 
           bcs_st(1:imax,k,is) = s(ip_t:ip_t+imax-1,is) 
        ELSE  ! Only Dirichlet BCs implemented for scalar
           CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
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
  IF ( kmax_total .GT. 1 ) THEN
     !CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,w,    tmp2,i0,i0,wrk1d,wrk2d,wrk3d) 
     !CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,tmp2, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
     CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, w, tmp1, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)   ! tmp2 used for BCs below
     DO ij = 1,isize_field; 
          h3(ij) =           -v(ij)*tmp2(ij)     
        tmp9(ij) =             visc*tmp1(ij)
     ENDDO

     !CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,w,    tmp3,i0,i0,wrk1d,wrk2d,wrk3d) 
     !CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,tmp3, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
     CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, w, tmp1, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
     DO ij = 1,isize_field; 
          h3(ij) =   h3(ij) -u(ij)*tmp3(ij)
        tmp9(ij) = tmp9(ij) + visc*tmp1(ij)
     ENDDO

     !CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,w,    tmp3,i0,i0,wrk1d,wrk2d,wrk3d) 
     !CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,tmp3, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
     CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc, &
          dz, w, tmp1, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)     
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
  !CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,u,    tmp2,i0,i0,wrk1d,wrk2d,wrk3d) 
  !CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,tmp2, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, u, tmp1, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)   ! tmp2 used for BCs below
  DO ij=1,isize_field; 
       h1(ij)  =           - v(ij)*tmp2(ij)  
     tmp7(ij)  =              visc*tmp1(ij)
  ENDDO

  !CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,u,    tmp3,i0,i0,wrk1d,wrk2d,wrk3d) 
  !CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,tmp3, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, u, tmp1, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  DO ij=1,isize_field; 
       h1(ij) =   h1(ij) - u(ij)*tmp3(ij); 
     tmp7(ij) = tmp7(ij) + visc*tmp1(ij)
  ENDDO

  !CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,u,    tmp3,i0,i0,wrk1d,wrk2d,wrk3d) 
  !CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,tmp3, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, u, tmp1, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d) 

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
  !CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,v,    tmp2,i0,i0,wrk1d,wrk2d,wrk3d) 
  !CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,tmp2, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, v, tmp1, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  DO ij=1,isize_field
       h2(ij) =          - v(ij)*tmp2(ij)
     tmp8(ij) =             visc*tmp1(ij) 
  ENDDO


  !CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,v,    tmp3,i0,i0,wrk1d,wrk2d,wrk3d) 
  !CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,tmp3, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, v, tmp1, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  DO ij=1,isize_field
       h2(ij) =   h2(ij) - u(ij)*tmp3(ij)
     tmp8(ij) = tmp8(ij) + visc*tmp1(ij) 
  ENDDO


  !CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,v,    tmp3,i0,i0,wrk1d,wrk2d,wrk3d) 
  !CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,tmp3, tmp1,i0,i0,wrk1d,wrk2d,wrk3d) 
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, v, tmp1, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d) 

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
  IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(&
          vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
          vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), q,hq ) 
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
        !CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
        !     dx, s(1,is), tmp4, i0,i0,i0,i0, tmp5,wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,s(1,is),tmp5,i0,i0,wrk1d,wrk2d,wrk3d) 
        CALL PARTIAL_X(imode_fdm,imax,jmax,kmax,i1bc,dx,tmp5,   tmp4,i0,i0,wrk1d,wrk2d,wrk3d) 
        DO ij=1,isize_field; 
          hs(ij,is) =           -u(ij)*tmp5(ij) 
           tmp7(ij) =             diff*tmp4(ij)  
        ENDDO
        
        !CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
        !     dy, s(1,is), tmp4, i0,i0,i0,i0, tmp5,wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,s(1,is),tmp5,i0,i0,wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Y(imode_fdm,imax,jmax,kmax,j1bc,dy,tmp5,   tmp4,i0,i0,wrk1d,wrk2d,wrk3d) 
        DO ij=1,isize_field; 
          hs(ij,is) =hs(ij,is) -v(ij)*tmp5(ij)
           tmp7(ij) = tmp7(ij)  +diff*tmp4(ij)
        ENDDO
        
        !CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
        !     dz, s(1,is), tmp4, i0,i0,i0,i0, tmp5,wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,s(1,is),tmp5,i0,i0,wrk1d,wrk2d,wrk3d) 
        CALL PARTIAL_Z(imode_fdm,imax,jmax,kmax,k1bc,dz,tmp5,   tmp4,i0,i0,wrk1d,wrk2d,wrk3d) 

        DO ij=1,isize_field; 
          hs(ij,is) =hs(ij,is) -w(ij)*tmp5(ij) 
           tmp7(ij) = tmp7(ij)  +diff*tmp4(ij)
        ENDDO

! #######################################################################
! Impose buffer zone as relaxation terms (Scalar #is) 
! #######################################################################
        IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
           CALL BOUNDARY_BUFFER_RELAXATION_SCAL(is,&
                vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
                vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), q,s(1,is), hs(1,is) ) 
        ENDIF
        

        DO ij=1,isize_field
           s(ij,is) = s(ij,is) + dte*( tmp7(ij) + hs(ij,is) + kco*tmp6(ij) );  
        ENDDO
        diff= visc_imp/schmidt(is) 
        alpha= dte*diff
        beta =-C_1_R/alpha 

        ip_b = 1     ; ip_t = 1+ imax*kmax
        DO k=1,kmax 
           ip=(k-1)*imax*jmax+1; wrk2d(ip_b:ip_b+imax-1)=-alpha*bcs_sb(1:imax,k,is); ip_b=ip_b+imax
           ip=ip+(jmax-1)*imax;  wrk2d(ip_t:ip_t+imax-1)=-alpha*bcs_st(1:imax,k,is); ip_t=ip_t+imax
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
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, w, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
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
     p_bcs => v(ip_b:); bcs_hb(1:imax,k,3) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => v(ip_t:); bcs_ht(1:imax,k,3) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! pressure in tmp1, Oy derivative in tmp3
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, i3, &
       tmp1,tmp3, tmp2,tmp4, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i0, dx, tmp1, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, i0, dz, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)

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
           s(ip_b:ip_b+imax-1,is) = bcs_sb(1:imax,k,is) 
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
           s(ip_t:ip_t+imax-1,is) = bcs_st(1:imax,k,is) 
        ENDDO
     ENDIF
     ip_t = ip_t + nxy
  ENDDO


  RETURN
END SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_IMPLICIT_1
