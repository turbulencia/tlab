#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

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
SUBROUTINE  RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1&
     (dte, u,v,w,h1,h2,h3, q,hq, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_wrk1d, inb_flow
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : visc
  USE DNS_LOCAL,  ONLY : idivergence
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS

  IMPLICIT NONE

#include "integers.h"

  TREAL dte
  TREAL, DIMENSION(isize_field)   :: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,*) :: q,hq
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(*)             :: wrk2d,wrk3d

  TARGET tmp2, h2

! -----------------------------------------------------------------------
  TINTEGER iq, ij, k, nxy, ip_b, ip_t
  TINTEGER ibc, bcs(2,2)
  TREAL dummy

  TINTEGER siz, srt, end    !  Variables for OpenMP Partitioning 

  TREAL, DIMENSION(:), POINTER :: p_bcs

#ifdef USE_ESSL
  INTEGER ilen
#endif

! #######################################################################
  nxy   = imax*jmax

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

#ifdef USE_ESSL
  ilen = isize_field
#endif

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), u, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), u, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, bcs, g(1), u,u,u, tmp4, tmp1, wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
  DO ij = srt,end
     h1(ij) = h1(ij) + tmp4(ij) + visc*( tmp6(ij)+tmp5(ij) ) &
             - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) )
  ENDDO
!$omp end parallel

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  IF ( g(3)%size .GT. 1 ) THEN
  CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax, bcs, g(3), w,w,w, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), w, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), w, tmp4, tmp1, wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij, srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 
  DO ij = srt,end
     h3(ij) = h3(ij) + tmp6(ij) + visc*( tmp5(ij)+tmp4(ij) ) &
          - ( v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO
!$omp end parallel

  ENDIF

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), v, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), v,v,v, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), v, tmp4, tmp1, wrk2d,wrk3d)

!$omp parallel default( shared ) &
#ifdef USE_ESSL
!$omp private( ilen, srt,end,siz)
#else
!$omp private( ij,   srt,end,siz)
#endif 
  CALL DNS_OMP_PARTITION(isize_field, srt,end,siz) 
  DO ij = srt,end
     h2(ij) = h2(ij) + tmp5(ij) + visc*( tmp6(ij)+tmp4(ij) ) &
          - ( w(ij)*tmp3(ij) + u(ij)*tmp1(ij) )
  ENDDO
!$omp end parallel 

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
  IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(q, hq)
  ENDIF

! #######################################################################
! Pressure term
! #######################################################################
  IF ( idivergence .EQ. EQNS_DIVERGENCE ) THEN ! remove residual divergence

!$omp parallel default( shared )&
#ifdef USE_ESSL
!$omp private( ilen, dummy, srt,end,siz )
#else 
!$omp private( ij,   dummy, srt,end,siz )
#endif 

     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
     dummy=C_1_R/dte
     
#ifdef USE_ESSL
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

     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp1, wrk3d,wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp3, tmp2, wrk3d,wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4, tmp3, wrk3d,wrk2d,wrk3d)

  ELSE
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), h2, tmp1, wrk3d,wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), h1, tmp2, wrk3d,wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), h3, tmp3, wrk3d,wrk2d,wrk3d)

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
     p_bcs => h2(ip_b:); BcsFlowJmin%ref(1:imax,k,2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
     p_bcs => h2(ip_t:); BcsFlowJmax%ref(1:imax,k,2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
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
!$omp parallel default( shared ) &
#ifdef USE_ESSL
!$omp private( ilen, srt,end,siz,dummy )
#else
!$omp private( ij,   srt,end,siz,dummy )
#endif
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)

#ifdef USE_ESSL
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
  BcsFlowJmin%ref = C_0_R ! default is no-slip (dirichlet)
  BcsFlowJmax%ref = C_0_R

  DO iq = 1,inb_flow
     ibc = 0
     IF ( BcsFlowJmin%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( BcsFlowJmax%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( ibc .GT. 0 ) THEN
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hq(1,iq), &
             BcsFlowJmin%ref(1,1,iq),BcsFlowJmax%ref(1,1,iq), wrk1d,tmp1,wrk3d)
     ENDIF
  ENDDO

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin 
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     h1(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,1)
     h2(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,2)
     h3(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,3); ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     h1(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,1)
     h2(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,2)
     h3(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,3); ip_t = ip_t + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1

