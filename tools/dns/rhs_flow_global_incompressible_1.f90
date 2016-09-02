#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2011/10/01 - J.P. Mellado
!#              Adding the buffer zone 
!# 2011/11/28 - C. Ansorge
!#              Combining BLAS and OMP
!#
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
SUBROUTINE  RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1&
     (dte,etime, x,y,z,dx,dy,dz, u,v,w,h1,h2,h3, s, q,hq, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     bcs_hb,bcs_ht,b_ref, vaux, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL
  USE DNS_LOCAL,  ONLY : bcs_flow_jmin, bcs_flow_jmax
  USE DNS_LOCAL,  ONLY : idivergence
  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : buff_type

  IMPLICIT NONE

#include "integers.h"

  TREAL dte,etime
  TREAL, DIMENSION(*)             :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field)   :: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,*) :: s, q,hq
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(*)             :: wrk2d,wrk3d, vaux
  TREAL, DIMENSION(jmax)          :: b_ref
  TREAL, DIMENSION(imax,kmax,*)   :: bcs_hb, bcs_ht

  TARGET tmp2, h2

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip_b, ip_t
  TINTEGER ibc
  TREAL dummy, u_geo, w_geo

  TINTEGER siz, srt, end    !  Variables for OpenMP Partitioning 

  TREAL, DIMENSION(:), POINTER :: p_bcs

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
! Diffusion and convection terms in Ox momentum eqn
! #######################################################################
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, u, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, u, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL OPR_BURGERS_X(i0,i0, imode_fdm, imax,jmax,kmax, &
       g(1), u,u,u, tmp4, i0,i0, i0,i0, tmp1, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction. 
! -----------------------------------------------------------------------

!$omp parallel default( shared ) private( ij, dummy,srt,end,siz )

  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)

  IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
     dummy = coriolis%vector(2)
     DO ij = srt, end 
        h1(ij) = h1(ij) + dummy*( w_geo-w(ij) )  + tmp4(ij) + visc*( tmp6(ij)+tmp5(ij) ) &
             - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) )
     ENDDO
! -----------------------------------------------------------------------
  ELSE
     DO ij = srt,end
        h1(ij) = h1(ij) + tmp4(ij) + visc*( tmp6(ij)+tmp5(ij) ) &
             - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) )
     ENDDO
  ENDIF
!$omp end parallel

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( buoyancy%active(1) ) THEN
     wrk1d(:,1) = C_0_R
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)

!$omp parallel default( shared ) &
#ifdef USE_BLAS
!$omp private( ilen, dummy,srt,end,siz )
#else 
!$omp private( ij,   dummy,srt,end,siz )
#endif 

     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 
     dummy = buoyancy%vector(1)

#ifdef USE_BLAS
     ilen = siz
     CALL DAXPY(ilen, dummy, wrk3d(srt), 1, h1(srt), 1)
#else
     DO ij = srt,end
        h1(ij) = h1(ij) + dummy*wrk3d(ij)
     ENDDO
#endif
!$omp end parallel
     
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
  CALL OPR_BURGERS_Z(i0,i0, imode_fdm, imax,jmax,kmax,&
       g(3), w,w,w, tmp6, i0,i0, i0,i0, tmp3, wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, w, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, w, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction. 
! -----------------------------------------------------------------------

!$omp parallel default( shared ) private( ij, dummy,srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 
  IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
     dummy = coriolis%vector(2)
     DO ij = srt,end
        h3(ij) = h3(ij) + dummy*( u(ij)-u_geo ) + tmp6(ij) + visc*( tmp5(ij)+tmp4(ij) ) &
             - ( v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
     ENDDO
! -----------------------------------------------------------------------
  ELSE
     DO ij = srt,end
        h3(ij) = h3(ij) + tmp6(ij) + visc*( tmp5(ij)+tmp4(ij) ) &
             - ( v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
     ENDDO
  ENDIF
!$omp end parallel

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( buoyancy%active(3) ) THEN
     wrk1d(:,1) = C_0_R
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
     
!$omp parallel default( shared ) &
#ifdef USE_BLAS
!$omp private( ilen, dummy, srt,end,siz)
#else     
!$omp private( ij,   dummy, srt,end,siz )
#endif
     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 
     dummy = buoyancy%vector(3)

#ifdef USE_BLAS
     ilen = siz
     CALL DAXPY(ilen, dummy, wrk3d(srt), 1, h3(srt), 1)
#else
     DO ij = srt,end
        h3(ij) = h3(ij) + dummy*wrk3d(ij)
     ENDDO
#endif

!$omp end parallel
     
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
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, v, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL OPR_BURGERS_Y(i0,i0, imode_fdm, imax,jmax,kmax,&
       g(2), v,v,v, tmp5, i0,i0, i0,i0, tmp2, wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, v, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------

  IF ( buoyancy%active(2) ) THEN
     CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, b_ref)
!$omp parallel default( shared ) &
#ifdef USE_BLAS
!$omp private( ilen, srt,end,siz, dummy )
#else
!$omp private( ij,   srt,end,siz, dummy )
#endif 
     CALL DNS_OMP_PARTITION(isize_field, srt,end,siz)      
     dummy = buoyancy%vector(2)
     DO ij = srt,end
        h2(ij) = h2(ij) + dummy*wrk3d(ij) + tmp5(ij) + visc*( tmp6(ij)+tmp4(ij) ) &
             - ( w(ij)*tmp3(ij) + u(ij)*tmp1(ij) )
     ENDDO
!$omp end parallel 

  ELSE
!$omp parallel default( shared ) &
#ifdef USE_BLAS
!$omp private( ilen, srt,end,siz, dummy )
#else
!$omp private( ij,   srt,end,siz, dummy )
#endif 
     CALL DNS_OMP_PARTITION(isize_field, srt,end,siz) 
     DO ij = srt,end
        h2(ij) = h2(ij) + tmp5(ij) + visc*( tmp6(ij)+tmp4(ij) ) &
             - ( w(ij)*tmp3(ij) + u(ij)*tmp1(ij) )
     ENDDO
!$omp end parallel 
     
  ENDIF

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
  IF ( buff_type .EQ. 1 .OR. buff_type .EQ. 3 ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(&
          vaux(vindex(VA_BUFF_HT)), vaux(vindex(VA_BUFF_HB)), &
          vaux(vindex(VA_BUFF_VI)), vaux(vindex(VA_BUFF_VO)), x,y, q,hq)
  ENDIF

! #######################################################################
! Pressure term
! #######################################################################
  IF ( idivergence .EQ. EQNS_DIVERGENCE ) THEN ! remove residual divergence

!$omp parallel default( shared )&
#ifdef USE_BLAS
!$omp private( ilen, dummy, srt,end,siz )
#else 
!$omp private( ij,   dummy, srt,end,siz )
#endif 

     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
     dummy=C_1_R/dte
     
#ifdef USE_BLAS
     ilen = siz
     CALL DZAXPY(ilen, dummy, v(srt), 1, h2(srt), 1, tmp2(srt), 1)
     CALL DZAXPY(ilen, dummy, u(srt), 1, h1(srt), 1, tmp3(srt), 1)
     CALL DZAXPY(ilen, dummy, w(srt), 1, h3(srt), 1, tmp4(srt), 1)

#else
     DO ij = srt,end 
        tmp2(ij) = h2(ij) + v(ij)*dummy
        tmp3(ij) = h1(ij) + u(ij)*dummy
        tmp4(ij) = h3(ij) + w(ij)*dummy
     ENDDO

#endif
!$omp end parallel

     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp2, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp3, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp4, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

  ELSE
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, h2, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, h1, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, h3, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

  ENDIF

! -----------------------------------------------------------------------
!$omp parallel default( shared ) private( ij,srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
  DO ij = srt,end
     tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
  ENDDO
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
  CALL OPR_POISSON_FXZ(imode_fdm,i2,i3, imax,jmax,kmax, &
       y,dx,dy,dz, tmp1,tmp3, tmp2,tmp4, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

! horizontal derivatives
  CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i0, dx, tmp1, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, i0, dz, tmp1, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient 
! -----------------------------------------------------------------------
!$omp parallel default( shared ) &
#ifdef USE_BLAS
!$omp private( ilen, srt,end,siz,dummy )
#else
!$omp private( ij,   srt,end,siz,dummy )
#endif
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)

#ifdef USE_BLAS
  ilen = siz 
  dummy=-C_1_R
  CALL DAXPY(ilen, dummy, tmp2(srt), 1, h1(srt),1)
  CALL DAXPY(ilen, dummy, tmp3(srt), 1, h2(srt),1)
  CALL DAXPY(ilen, dummy, tmp4(srt), 1, h3(srt),1)
#else
  DO ij = srt,end
     h1(ij) = h1(ij) - tmp2(ij)
     h2(ij) = h2(ij) - tmp3(ij)
     h3(ij) = h3(ij) - tmp4(ij)
  ENDDO
#endif

!$omp end parallel


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
END SUBROUTINE RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1

