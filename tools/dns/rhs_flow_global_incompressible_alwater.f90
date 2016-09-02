!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!# 2011/05/01    A. de Lozar
!#               Modified for the cloud mixture
!########################################################################
!# DESCRIPTION
!#
!# Momentum equations, nonlinear term in convective form and the 
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# BCs need 2 1st order derivatives in Oy
!# Scalar needed for the buoyancy term
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

SUBROUTINE  RHS_FLOW_GLOBAL_INCOMPRESSIBLE_ALWATER&
     (dte,etime, x,y,z,dx,dy,dz, u,v,w,h1,h2,h3, s, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     bcs_hb,bcs_ht,b_ref, wrk1d,wrk2d,wrk3d)
  
  USE OMP_LIB
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL, ONLY : bcs_flow_jmin, bcs_flow_jmax

  IMPLICIT NONE

#include "integers.h"

  TREAL dte,etime
  TREAL, DIMENSION(*)             :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field)   :: u,v,w, h1,h2,h3 !, s Alberto, now s it is a multidimensional array
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(*)             :: wrk2d, wrk3d
  TREAL, DIMENSION(jmax)          :: b_ref
  TREAL, DIMENSION(imax,kmax,*)   :: bcs_hb, bcs_ht

  TREAL, DIMENSION(isize_field,*), TARGET   :: s !Alberto. Now mulidimensinal array (including h,qt,ql)

  TARGET tmp2, h2

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip_b, ip_t
  TINTEGER ibc
  TREAL dummy, u_geo, w_geo

  TREAL, DIMENSION(:), POINTER :: p_bcs

!Alberto.Ponter to the qt,ql (2dim array) and enthalpy
  TREAL, DIMENSION(:,:),  POINTER :: al_q
  TREAL, DIMENSION(:),  POINTER ::al_h


#ifdef USE_BLAS
  INTEGER ilen
#endif

! #######################################################################
  nxy   = imax*jmax
  u_geo = COS(coriolis%parameters(1))
  w_geo =-SIN(coriolis%parameters(1))

#ifdef USE_BLAS
  ilen = isize_field
#endif

! #######################################################################
! Validation
! #######################################################################
!  CALL FI_FORCING_1(iunifx,iunify, imode_fdm, imax,jmax,kmax, i1bc,j1bc, &
!       etime,visc, x,y,dx,dy, h1,h2, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)
!  CALL FI_FORCING_0(imax,jmax,kmax, etime,visc, x,y, u,v, h1,h2)

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! #######################################################################
! CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_MPI&
!    (visc, dx,dy,dz, u,v,w,u,h1, tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, u, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, u, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, u, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction. 
! -----------------------------------------------------------------------
  IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
!$omp parallel default( shared ) private( ij, dummy )
     dummy = coriolis%vector(2)
!$omp do
     DO ij = 1,isize_field
        h1(ij) = h1(ij) + dummy*( w_geo-w(ij) ) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
             - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
     ENDDO
!$omp end do
!$omp end parallel

! -----------------------------------------------------------------------
  ELSE
!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        h1(ij) = h1(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
             - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
     ENDDO
!$omp end do
!$omp end parallel

  ENDIF

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( buoyancy%active(1) ) THEN
     wrk1d(:,1) = C_0_R
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)

#ifdef USE_BLAS
     dummy = buoyancy%vector(1)
     CALL DAXPY(ilen, dummy, wrk3d, 1, h1, 1)

#else
!$omp parallel default( shared ) private( ij, dummy )
     dummy = buoyancy%vector(1)
!$omp do
     DO ij = 1,isize_field
        h1(ij) = h1(ij) + dummy*wrk3d(ij)
     ENDDO
!$omp end do
!$omp end parallel
#endif

  ENDIF

! -----------------------------------------------------------------------
! BCs s.t. derivative d/dy(u) is zero
! -----------------------------------------------------------------------
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) THEN
     ip_b =                 1
     DO k = 1,kmax
        p_bcs => tmp2(ip_b:); bcs_hb(1:imax,k,1) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     ENDDO
  ENDIF
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) THEN
     ip_t = imax*(jmax-1) + 1
     DO k = 1,kmax
        p_bcs => tmp2(ip_t:); bcs_ht(1:imax,k,1) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
     ENDDO
  ENDIF

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  IF ( kmax_total .GT. 1 ) THEN
! CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_MPI&
!    (visc, dx,dy,dz, u,v,w,w,h3, tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

     CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, w, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, w, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, w, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction. 
! -----------------------------------------------------------------------
     IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
!$omp parallel default( shared ) private( ij, dummy )
        dummy = coriolis%vector(2)
!$omp do
        DO ij = 1,isize_field
           h3(ij) = h3(ij) + dummy*( u(ij)-u_geo ) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
                - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
        ENDDO
!$omp end do
!$omp end parallel

! -----------------------------------------------------------------------
     ELSE
!$omp parallel default( shared ) private( ij )
!$omp do
        DO ij = 1,isize_field
           h3(ij) = h3(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
                - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
        ENDDO
!$omp end do
!$omp end parallel
     ENDIF

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
     IF ( buoyancy%active(3) ) THEN
        wrk1d(:,1) = C_0_R
        CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)

#ifdef USE_BLAS
        dummy = buoyancy%vector(3)
        CALL DAXPY(ilen, dummy, wrk3d, 1, h3, 1)

#else
!$omp parallel default( shared ) private( ij, dummy )
        dummy = buoyancy%vector(3)
!$omp do
        DO ij = 1,isize_field
           h3(ij) = h3(ij) + dummy*wrk3d(ij)
        ENDDO
!$omp end do
!$omp end parallel
#endif

     ENDIF

! -----------------------------------------------------------------------
! BCs s.t. derivative d/dy(w) is zero
! -----------------------------------------------------------------------
     IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) THEN
        ip_b =                 1
        DO k = 1,kmax
           p_bcs => tmp2(ip_b:); bcs_hb(1:imax,k,2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
        ENDDO
     ENDIF
     IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) THEN
        ip_t = imax*(jmax-1) + 1
        DO k = 1,kmax
           p_bcs => tmp2(ip_t:); bcs_ht(1:imax,k,2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
        ENDDO
     ENDIF

  ENDIF

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
! CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_MPI&
!    (visc, dx,dy,dz, u,v,w,v,h2, tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, v, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, v, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, v, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( buoyancy%active(2) ) THEN
!ALberto: Changed to calculate the density before performing any operation
     IF ( (imixture .EQ. MIXT_TYPE_AIRWATER) .OR. (imixture .EQ. MIXT_TYPE_SUPSAT) ) THEN !Alberto: everything equally valid for Super Saturation
        al_h=>s(:,1)  !We point to the scalar array which contains the enthalpy
        al_q=> s(:,2:3) ! We point to the scalar array which contains qt,ql
!Calculate the density and store it in wrk3d
        CALL THERMO_THERMAL_DENSITY_HP_ALWATER(imax,jmax,kmax, al_q,al_h,p_init,wrk3d)   
!Alberto: Now we take care of the bouyancy term	
!$omp parallel default( shared ) private( ij, dummy)
        dummy =  buoyancy%vector(2)/mean_rho
!$omp do
        DO ij = 1,isize_field
           h2(ij) = h2(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
                - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) ) &
                + dummy*(wrk3d(ij)-mean_rho)
        ENDDO
!$omp end do
!$omp end parallel

!************* END changes Alberto     
     ELSE
        CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, b_ref) 
!$omp parallel default( shared ) private( ij, dummy)
        dummy = buoyancy%vector(2)
!$omp do
        DO ij = 1,isize_field
           h2(ij) = h2(ij) + dummy*wrk3d(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
                - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
        ENDDO
!$omp end do
!$omp end parallel
     ENDIF

  ELSE
!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        h2(ij) = h2(ij) + visc*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
             - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
     ENDDO
!$omp end do
!$omp end parallel

  ENDIF

! #######################################################################
! Pressure term
! #######################################################################
#ifdef USE_BLAS
  dummy=C_1_R/dte
  CALL DZAXPY(ilen, dummy, v, 1, h2, 1, tmp2, 1)
  CALL DZAXPY(ilen, dummy, u, 1, h1, 1, tmp3, 1)
  CALL DZAXPY(ilen, dummy, w, 1, h3, 1, tmp4, 1)

#else
!$omp parallel default( shared ) private( ij, dummy )
  dummy=C_1_R/dte
!$omp do
  DO ij = 1,isize_field
     tmp2(ij) = h2(ij) + v(ij)*dummy
     tmp3(ij) = h1(ij) + u(ij)*dummy
     tmp4(ij) = h3(ij) + w(ij)*dummy
  ENDDO
!$omp end do
!$omp end parallel

#endif

  CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp2, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp3, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp4, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,isize_field
     tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
  ENDDO
!$omp end do
!$omp end parallel

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     p_bcs => h2(ip_b:); bcs_hb(1:imax,k,3) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => h2(ip_t:); bcs_ht(1:imax,k,3) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
  ENDDO

! pressure in tmp1, Oy derivative in tmp3
  CALL OPR_POISSON_FXZ(imode_fdm,i2,i3, imax,jmax,kmax,  &
       y,dx,dy,dz, tmp1,tmp3, tmp2,tmp4, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i0, dx, tmp1, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, i0, dz, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient 
! -----------------------------------------------------------------------
#ifdef USE_BLAS
  dummy=-C_1_R
  CALL DAXPY(ilen, dummy, tmp2, 1, h1, 1)
  CALL DAXPY(ilen, dummy, tmp3, 1, h2, 1)
  CALL DAXPY(ilen, dummy, tmp4, 1, h3, 1)

#else
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,isize_field
     h1(ij) = h1(ij) - tmp2(ij)
     h2(ij) = h2(ij) - tmp3(ij)
     h3(ij) = h3(ij) - tmp4(ij)
  ENDDO
!$omp end do
!$omp end parallel

#endif

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
  ibc = 0
  bcs_hb(:,:,1:2) = C_0_R ! default is no-slip (dirichlet)
  bcs_ht(:,:,1:2) = C_0_R
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, h1, bcs_hb(1,1,1),bcs_ht(1,1,1), wrk1d,tmp1,wrk3d)
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, h3, bcs_hb(1,1,2),bcs_ht(1,1,2), wrk1d,tmp1,wrk3d)
  ENDIF

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin 
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     h1(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
     h2(ip_b:ip_b+imax-1) = C_0_R               ! no penetration
     h3(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2); ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     h1(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
     h2(ip_t:ip_t+imax-1) = C_0_R               ! no penetration
     h3(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2); ip_t = ip_t + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_INCOMPRESSIBLE_ALWATER

