#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Evolution equations, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!# It is written such that u and w transposes are calculated first for the
!# Ox and Oz momentum equations, stored in tmp5 and tmp6 and then used as needed.
!# This saves 2 MPI transpositions.
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1&
     (dte, u,v,w,h1,h2,h3, q,hq, s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, &
     wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
#ifdef TRACE_ON
  USE DNS_CONSTANTS,ONLY:tfile
#endif
  USE TLAB_VARS, ONLY : imode_eqns
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field, isize_wrk1d
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : rbackground, ribackground
  USE DNS_LOCAL,  ONLY : idivergence
  USE DNS_LOCAL,  ONLY : tower_mode
  USE TIME,       ONLY : rkm_substep,rkm_endstep
  USE DNS_TOWER
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS

  IMPLICIT NONE

#include "integers.h"

  TREAL dte
  TREAL, DIMENSION(isize_field)   :: u,v,w, h1,h2,h3
  TREAL, DIMENSION(isize_field,*) :: q,hq, s,hs
  TREAL, DIMENSION(isize_field)   :: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk3d
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d
  TREAL, DIMENSION(imax,kmax,*)   :: wrk2d

  TARGET h2, hs

! -----------------------------------------------------------------------
  TINTEGER iq, is, ij, k, nxy, ip_b, ip_t
  TINTEGER ibc, bcs(2,2)
  TREAL dummy

  TINTEGER siz, srt, end    !  Variables for OpenMP Partitioning

  TREAL, DIMENSION(:), POINTER :: p_bcs

#ifdef USE_ESSL
  INTEGER ilen
#endif

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'ENTERING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
#endif


! #######################################################################
  nxy = imax*jmax

  bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

#ifdef USE_ESSL
  ilen = isize_field
#endif

! #######################################################################
! Preliminaries for Scalar BC
! (flow BCs initialized below as they are used for pressure in between)
! #######################################################################
! Default is zero
  BcsScalJmin%ref(:,:,:) = C_0_R
  BcsScalJmax%ref(:,:,:) = C_0_R

! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     DO is =1,inb_scal
        IF ( BcsScalJmin%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
             p_bcs => hs(ip_b:,is); BcsScalJmin%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
        IF ( BcsScalJmax%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
             p_bcs => hs(ip_t:,is); BcsScalJmax%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
     ENDDO
     ip_b = ip_b + nxy ! bottom BC address
     ip_t = ip_t + nxy ! top BC address
  ENDDO

! #######################################################################
! Ox diffusion and convection terms in Ox momentum eqn
! Initializing tmp5 for the rest of terms
! #######################################################################
  CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, bcs, g(1), u,u,u, tmp1, tmp5, wrk2d,wrk3d) ! store u transposed in tmp5
  h1 = h1 + tmp1

! #######################################################################
! Oy diffusion and convection terms in Oy momentum eqn
! Initializing tmp4 for the rest of terms
! #######################################################################
  CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), v,v,v, tmp2, tmp4, wrk2d,wrk3d) ! store v transposed in tmp4
  h2 = h2 + tmp2

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
  IF ( g(3)%size .GT. 1 ) THEN

     CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax, bcs, g(1), w,u,tmp5, tmp1, tmp6, wrk2d,wrk3d) ! tmp5 contains u transposed
     CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), w,v,tmp4, tmp2, tmp6, wrk2d,wrk3d) ! tmp4 contains v transposed
     CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax, bcs, g(3), w,w,w,    tmp3, tmp6, wrk2d,wrk3d) ! store w transposed in tmp6

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
     DO ij = srt,end
        h3(ij) = h3(ij) +tmp1(ij) +tmp2(ij) +tmp3(ij)
     ENDDO
!$omp end parallel

  ENDIF

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
  CALL OPR_BURGERS_X(i1,i0, imax,jmax,kmax, bcs, g(1), v,u,tmp5, tmp1, tmp2, wrk2d,wrk3d) ! tmp5 contains u transposed
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax, bcs, g(3), v,w,tmp6, tmp3, tmp2, wrk2d,wrk3d) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field, srt,end,siz)
  DO ij = srt,end
     h2(ij) = h2(ij) + tmp1(ij) +tmp3(ij)
  ENDDO
!$omp end parallel

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! The term \nu u'' - u u' has been already added in the beginning
! #######################################################################
  CALL OPR_BURGERS_Y(i1,i0, imax,jmax,kmax, bcs, g(2), u,v,tmp4, tmp2, tmp1, wrk2d,wrk3d) ! tmp4 contains w transposed
  CALL OPR_BURGERS_Z(i1,i0, imax,jmax,kmax, bcs, g(3), u,w,tmp6, tmp3, tmp1, wrk2d,wrk3d) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
  CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
  DO ij = srt,end
     h1(ij) = h1(ij) +tmp2(ij) +tmp3(ij)
  ENDDO
!$omp end parallel

! #######################################################################
! Diffusion and convection terms in scalar eqns
! #######################################################################
  DO is = 1,inb_scal

     CALL OPR_BURGERS_Y(i1,is, imax,jmax,kmax, bcs, g(2), s(1,is),v,tmp4, tmp2, tmp1, wrk2d,wrk3d) ! Not enough tmp arrays
     hs(:,is) = hs(:,is) +tmp2

     CALL OPR_BURGERS_X(i1,is, imax,jmax,kmax, bcs, g(1), s(1,is),u,tmp5, tmp1, tmp2, wrk2d,wrk3d) ! tmp5 contains u transposed
     CALL OPR_BURGERS_Z(i1,is, imax,jmax,kmax, bcs, g(3), s(1,is),w,tmp6, tmp3, tmp2, wrk2d,wrk3d) ! tmp6 contains w transposed

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
     DO ij = srt,end
        hs(ij,is) = hs(ij,is) +tmp1(ij) +tmp3(ij)
     ENDDO
!$omp end parallel

  ENDDO

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
  IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
     CALL BOUNDARY_BUFFER_RELAX_FLOW(q, hq)
  ENDIF

! #######################################################################
! Pressure term
! #######################################################################
  IF ( idivergence .EQ. EQNS_DIVERGENCE ) THEN ! remove residual divergence

#ifdef USE_ESSL
!$omp parallel default( shared )&
!$omp private( ilen, dummy, srt,end,siz )
#else
!$omp parallel default( shared )&
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

     IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp2)
        CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp3)
        CALL THERMO_ANELASTIC_WEIGHT_INPLACE(imax,jmax,kmax, rbackground, tmp4)
     ENDIF
     CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2,tmp1, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp3,tmp2, wrk3d, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4,tmp3, wrk3d, wrk2d,wrk3d)

  ELSE
     IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
        CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, h2,tmp2)
        CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, h1,tmp3)
        CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, h3,tmp4)
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2,tmp1, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp3,tmp2, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp4,tmp3, wrk3d, wrk2d,wrk3d)
     ELSE
        CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), h2,tmp1, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), h1,tmp2, wrk3d, wrk2d,wrk3d)
        CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), h3,tmp3, wrk3d, wrk2d,wrk3d)
     ENDIF

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

! Adding density in BCs
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     BcsFlowJmin%ref(:,:,2) = BcsFlowJmin%ref(:,:,2) *rbackground(1)
     BcsFlowJmax%ref(:,:,2) = BcsFlowJmax%ref(:,:,2) *rbackground(g(2)%size)
  ENDIF

! pressure in tmp1, Oy derivative in tmp3
  CALL OPR_POISSON_FXZ(.TRUE., imax,jmax,kmax, g, i3, &
       tmp1,tmp3, tmp2,tmp4, BcsFlowJmin%ref(1,1,2),BcsFlowJmax%ref(1,1,2), wrk1d,wrk1d(1,5),wrk3d)

! Saving pressure for towers to tmp array
  IF ( tower_mode .EQ. 1 .AND. rkm_substep .EQ. rkm_endstep ) THEN
     CALL DNS_TOWER_ACCUMULATE(tmp1,i4,wrk1d)
  ENDIF

! horizontal derivatives
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1,tmp2, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp1,tmp4, wrk3d, wrk2d,wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax,jmax,kmax, ribackground, tmp2, h1)
     CALL THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax,jmax,kmax, ribackground, tmp3, h2)
     CALL THERMO_ANELASTIC_WEIGHT_SUBSTRACT(imax,jmax,kmax, ribackground, tmp4, h3)

  ELSE
#ifdef USE_ESSL
!$omp parallel default( shared ) &
!$omp private( ilen, srt,end,siz,dummy )
#else
!$omp parallel default( shared ) &
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
  ENDIF

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
  BcsFlowJmin%ref = C_0_R ! default is no-slip (dirichlet)
  BcsFlowJmax%ref = C_0_R ! Scalar BCs initialized at start of routine

  DO iq = 1,inb_flow
     ibc = 0
     IF ( BcsFlowJmin%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( BcsFlowJmax%type(iq) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( ibc .GT. 0 ) THEN
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hq(1,iq), &
             BcsFlowJmin%ref(1,1,iq),BcsFlowJmax%ref(1,1,iq), wrk1d,tmp1,wrk3d)
     ENDIF
  ENDDO

  DO is = 1,inb_scal
     ibc = 0
     IF ( BcsScalJmin%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( BcsScalJmax%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( ibc .GT. 0 ) THEN
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hs(1,is), &
             BcsScalJmin%ref(1,1,is),BcsScalJmax%ref(1,1,is), wrk1d,tmp1,wrk3d)
     ENDIF

     IF ( BcsScalJmin%type(is) .NE. DNS_SFC_STATIC .OR. &
          BcsScalJmax%type(is) .NE. DNS_SFC_STATIC ) THEN
        CALL BOUNDARY_SURFACE_J(is,bcs,s,hs,tmp1,tmp2,tmp3,wrk1d,wrk2d,wrk3d)
     ENDIF
  ENDDO

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin
! -----------------------------------------------------------------------
  ip_b =                 1
  DO k = 1,kmax
     h1(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,1)
     h2(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,2)
     h3(ip_b:ip_b+imax-1) = BcsFlowJmin%ref(1:imax,k,3)
     DO is = 1,inb_scal
        hs(ip_b:ip_b+imax-1,is) = BcsScalJmin%ref(1:imax,k,is)
     ENDDO
     ip_b = ip_b + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     h1(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,1)
     h2(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,2)
     h3(ip_t:ip_t+imax-1) = BcsFlowJmax%ref(1:imax,k,3)
     DO is = 1,inb_scal
        hs(ip_t:ip_t+imax-1,is) = BcsScalJmax%ref(1:imax,k,is)
     ENDDO
     ip_t = ip_t + nxy
  ENDDO

#ifdef TRACE_ON
  CALL TLAB_WRITE_ASCII(tfile,'LEAVING SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1')
#endif
  RETURN
END SUBROUTINE RHS_GLOBAL_INCOMPRESSIBLE_1
