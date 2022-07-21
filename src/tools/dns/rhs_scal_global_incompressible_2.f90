#include "types.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Scalar equation, nonlinear term in skew-symmetric form and the 
!# diffusion term explicit. 3 2nd order + 3+3 1st order derivatives.
!#
!########################################################################
SUBROUTINE  RHS_SCAL_GLOBAL_INCOMPRESSIBLE_2&
     (is, u,v,w,s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : idiffusion, visc, schmidt
  USE BOUNDARY_BCS,ONLY: BcsScalJmin, BcsScalJmax
  USE TLAB_VARS, ONLY : imode_ibm
  USE DNS_IBM,   ONLY : imode_ibm_scal, ibm_partial

  IMPLICIT NONE

  TINTEGER is
  TREAL, DIMENSION(isize_field) :: u, v, w, s, hs
  TREAL, DIMENSION(isize_field) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk3d
  TREAL, DIMENSION(*)           :: wrk1d
  TREAL, DIMENSION(imax,kmax,2) :: wrk2d

  TARGET hs

  TREAL, DIMENSION(:), POINTER :: p_bcs

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip_b, ip_t, ibc, bcs(2,2)
  TREAL diff

! #######################################################################
  nxy = imax*jmax

  bcs = 0

  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE;                                  diff = visc/schmidt(is); ENDIF


! #######################################################################
! Preliminaries for Scalar BC
! (flow BCs initialized below as they are used for pressure in between)
! #######################################################################
! Default is zero (Dirichlet)
  BcsScalJmin%ref(:,:,is) = C_0_R
  BcsScalJmax%ref(:,:,is) = C_0_R

! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
  ip_b =                 1
  ip_t = imax*(jmax-1) + 1
  DO k = 1,kmax
     IF ( BcsScalJmin%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
        p_bcs => hs(ip_b:); BcsScalJmin%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
     IF ( BcsScalJmax%SfcType(is) .EQ. DNS_SFC_LINEAR ) THEN
        p_bcs => hs(ip_t:); BcsScalJmax%ref(1:imax,k,is) = p_bcs(1:imax); ENDIF
     ip_b = ip_b + nxy ! bottom BC address
     ip_t = ip_t + nxy ! top BC address
  ENDDO

! #######################################################################
! IBM
! #######################################################################
  IF ( imode_ibm == 1 ) THEN
     IF ( imode_ibm_scal == 1 ) THEN ! IBM usage for scalar field
        ! (requirenments: only possible with objects on bottom boundary 
        !  with homogeneous temperature in solid regions)
        ibm_partial = .true.
     ENDIF
  ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
  CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), s, tmp6, tmp3, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), s, tmp5, tmp2, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), s, tmp4, tmp1, wrk2d,wrk3d)

  DO ij = 1,isize_field
     hs(ij) = hs(ij) + diff*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - C_05_R*( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO

! #######################################################################
! Adding the conservative-form terms to complete skew-symmetric formulation
! #######################################################################
  tmp1 = u*s
  tmp2 = v*s
  tmp3 = w*s
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), tmp1, tmp4, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), tmp2, tmp5, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), tmp3, tmp6, wrk3d, wrk2d,wrk3d)
  hs = hs - C_05_R*( tmp6 + tmp5 + tmp4 )

! #######################################################################
! IBM
! #######################################################################
! IBM usage for scalar field, done
  IF ( imode_ibm_scal == 1 ) ibm_partial = .false.

! #######################################################################
! Boundary conditions
! #######################################################################
  ibc = 0
  IF ( BcsScalJmin%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( BcsScalJmax%type(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hs, &
          BcsScalJmin%ref(:,:,is),BcsScalJmax%ref(:,:,is), wrk1d,tmp1,wrk3d)
  ENDIF

  IF ( BcsScalJmin%SfcType(is) .NE. DNS_SFC_STATIC .OR. &
       BcsScalJmax%SfcType(is) .NE. DNS_SFC_STATIC ) THEN
     CALL BOUNDARY_SURFACE_J(is,bcs,s,hs,tmp1,tmp2,tmp3,wrk1d,wrk2d,wrk3d)
  ENDIF
     IF ( imode_ibm == 1 ) CALL IBM_BCS_FIELD(hs(1)) ! set tendency in solid to zero

! -----------------------------------------------------------------------
! Impose bottom at Jmin/max
! -----------------------------------------------------------------------
  ip_b = 1
  ip_t = 1 + imax*(jmax-1)
  DO k = 1,kmax
     hs(ip_b:ip_b+imax-1) = BcsScalJmin%ref(1:imax,k,is); ip_b = ip_b + nxy
     hs(ip_t:ip_t+imax-1) = BcsScalJmax%ref(1:imax,k,is); ip_t = ip_t + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_2
