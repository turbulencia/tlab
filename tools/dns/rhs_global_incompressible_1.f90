#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2011/11/28 - J.P. Mellado
!#              Created
!# 2012/12/20 - J.P. Mellado
!#              Rearranged to reduce transpositions more
!# 2015/10/09 - J.P. Mellado
!#              Extracting source code into arrays FI_SOURCES_*
!#
!########################################################################
!# DESCRIPTION
!#
!# Derived from RHS_FLOW/SCAL_GLOBAL_INCOMPRESSIBLE_1
!#
!# Momentum equations, nonlinear term in convective form and the 
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# BCs need 2 1st order derivatives in Oy
!# Scalar needed for the buoyancy term
!# It is written such that w transpose is calculated first for the 
!# Oz momentum equation, stored in tmp6 and then used as needed. This
!# saves 2 MPI transpositions. 
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
SUBROUTINE  RHS_GLOBAL_INCOMPRESSIBLE_1&
     (dte,etime, x,y,z,dx,dy,dz, u,v,w,h1,h2,h3, q,hq, s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     bcs_hb,bcs_ht,b_ref, vaux, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_CONSTANTS, ONLY : MAX_NSP
  USE DNS_GLOBAL
  USE DNS_LOCAL,  ONLY : bcs_flow_jmin, bcs_flow_jmax
  USE DNS_LOCAL,  ONLY : bcs_scal_jmin, bcs_scal_jmax
  USE DNS_LOCAL,  ONLY : idivergence
  USE DNS_LOCAL,  ONLY : VA_BUFF_HT, VA_BUFF_HB, VA_BUFF_VO, VA_BUFF_VI, vindex
  USE DNS_LOCAL,  ONLY : buff_type 
  USE DNS_LOCAL,  ONLY : rkm_substep,rkm_endstep,tower_mode 
  USE DNS_TOWER 
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TREAL dte,etime
  TREAL, DIMENSION(*)             :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field)   :: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,*) :: q,hq, s,hs
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(*)             :: wrk2d,wrk3d, vaux
  TREAL, DIMENSION(jmax)          :: b_ref
  TREAL, DIMENSION(imax,kmax,inb_vars) :: bcs_hb, bcs_ht

  TARGET tmp2, h2

! -----------------------------------------------------------------------
  TINTEGER is, ij, k, nxy, ip_b, ip_t
  TINTEGER ibc
  TREAL dummy, diff

  TINTEGER siz, srt, end    !  Variables for OpenMP Partitioning 

  TREAL, DIMENSION(:), POINTER :: p_bcs

#ifdef USE_BLAS
  INTEGER ilen
#endif

! #######################################################################
  nxy   = imax*jmax

#ifdef USE_BLAS
  ilen = isize_field
#endif

! #######################################################################
! Source terms
! #######################################################################
  CALL FI_SOURCES_FLOW(q,s, hq, b_ref, wrk1d,wrk3d)
  CALL FI_SOURCES_SCAL(y,dy, s, hs, tmp1,tmp2, wrk1d,wrk2d,wrk3d)

! #######################################################################
! Ox diffusion and convection terms in Ox momentum eqn
! Initializing tmp5 for the rest of terms
! #######################################################################
  CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, &
       g(1), u,u,u, tmp1, i0,i0, i0,i0, tmp5, wrk2d,wrk3d) ! tmp5 contains u transposed
  
  h1 = h1 + tmp1

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  IF ( kmax_total .GT. 1 ) THEN

  CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax,&
       g(3), w,w,w,    tmp1, i0,i0, i0,i0, tmp6, wrk2d,wrk3d) ! tmp6 contains w transposed
  CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax, &
       g(1), w,u,tmp5, tmp4, i0,i0, i0,i0, tmp2, wrk2d,wrk3d) ! tmp5 contains u transposed
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, w, tmp3, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)  ! tmp2 is used below in BCs

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 
  DO ij = srt,end
     h3(ij) = h3(ij) + tmp1(ij) + tmp4(ij) + visc*tmp3(ij) - v(ij)*tmp2(ij) 
  ENDDO
!$omp end parallel

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
  CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax,&
       g(2), v,v,v,    tmp1, i0,i0, i0,i0, tmp2, wrk2d,wrk3d)
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax,&
       g(3), v,w,tmp6, tmp4, i0,i0, i0,i0, tmp2, wrk2d,wrk3d) ! tmp6 contains w transposed
  CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax,&
       g(1), v,u,tmp5, tmp3, i0,i0, i0,i0, tmp2, wrk2d,wrk3d) ! tmp5 contains u transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field, srt,end,siz) 
  DO ij = srt,end
     h2(ij) = h2(ij) + tmp1(ij) + tmp4(ij) + tmp3(ij)
  ENDDO
!$omp end parallel 

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! The term \nu u'' - u u' has been already added in the beginning
! #######################################################################
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax,&
       g(3), u,w,tmp6, tmp4, i0,i0, i0,i0, tmp2, wrk2d,wrk3d) ! tmp6 contains w transposed
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc, & 
       dy, u, tmp3, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)  ! tmp2 is used below in BCs

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
  DO ij = srt,end
     h1(ij) = h1(ij) + tmp4(ij) + visc*tmp3(ij) - v(ij)*tmp2(ij)
  ENDDO
!$omp end parallel

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
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, i3, &
       tmp1,tmp3, tmp2,tmp4, bcs_hb(1,1,3),bcs_ht(1,1,3), wrk1d,wrk1d(1,5),wrk3d)

! Saving pressure for towers to tmp array 
  IF ( tower_mode .EQ. 1 .AND. rkm_substep .EQ. rkm_endstep ) THEN 
     CALL DNS_TOWER_ACCUMULATE(tmp1,i4,dx,dy,dz,wrk1d) 
  ENDIF

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
! Diffusion and convection terms in scalar eqn
! #######################################################################
  DO is = 1,inb_scal ! Start loop over the N scalars
     
     IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
     ELSE;                                  diff = visc/schmidt(is); ENDIF
        
     CALL OPR_BURGERS_Z(i1,is, imax,jmax,kmax,&
          g(3), s(1,is),w,tmp6, tmp1, i0,i0, i0,i0, tmp2, wrk2d,wrk3d) ! tmp6 contains w transposed
     CALL OPR_BURGERS_X(i1,is, imax,jmax,kmax,&
          g(1), s(1,is),u,tmp5, tmp2, i0,i0, i0,i0, tmp3, wrk2d,wrk3d) ! tmp5 contains u transposed
     CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, s(1,is), tmp3, i0,i0, i0,i0, tmp4, wrk1d,wrk2d,wrk3d)
     
!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
     DO ij = srt,end
        hs(ij,is) = hs(ij,is) + tmp1(ij) + tmp2(ij) + diff*tmp3(ij)  &
             - v(ij)*tmp4(ij)
     ENDDO
!$omp end parallel
     
  ENDDO
  
! #######################################################################
! Boundary conditions
! #######################################################################
  bcs_hb(:,:,1:inb_vars) = C_0_R ! default is no-slip (dirichlet)
  bcs_ht(:,:,1:inb_vars) = C_0_R

! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
  ibc = 0
  IF ( bcs_flow_jmin .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_flow_jmax .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, h1, bcs_hb(1,1,1),bcs_ht(1,1,1), wrk1d,tmp1,wrk3d)
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, h3, bcs_hb(1,1,2),bcs_ht(1,1,2), wrk1d,tmp1,wrk3d)
  ENDIF

  DO is = 1,inb_scal
  ibc = 0
  IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_scal_jmax(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, hs(1,is), bcs_hb(1,1,is+inb_flow),bcs_ht(1,1,is+inb_flow), wrk1d,tmp1,wrk3d)
  ENDIF
  ENDDO

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin 
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     DO is = 1,inb_scal
        hs(ip_b:ip_b+imax-1,is) = bcs_hb(1:imax,k,is+inb_flow) 
     ENDDO
     h1(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,1)
     h2(ip_b:ip_b+imax-1) = C_0_R               ! no penetration
     h3(ip_b:ip_b+imax-1) = bcs_hb(1:imax,k,2); ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     DO is = 1,inb_scal
        hs(ip_t:ip_t+imax-1,is) = bcs_ht(1:imax,k,is+inb_flow)
     ENDDO
     h1(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,1)
     h2(ip_t:ip_t+imax-1) = C_0_R               ! no penetration
     h3(ip_t:ip_t+imax-1) = bcs_ht(1:imax,k,2); ip_t = ip_t + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1

