#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Scalar equation, nonlinear term in convective form and the 
!# diffusion term explicit: 3 2nd order + 3 1st order derivatives.
!# BCs need 1 1st order derivatives in Oy
!#
!# Alberto: Changed to include inietia and differnt difusion. 
!# In this case it will assume that the first scalar is the 
!# enthalpy, the second qt and the third (no equation) ql.
!# File initialed copied from: RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1
!########################################################################
SUBROUTINE  RHS_SCAL_GLOBAL_INCOMPRESSIBLE_AIRWATER&
     (is, dte, dx,dy,dz, u,v,w, s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, wrk1d,wrk2d,wrk3d,h1,h2,h3)
  
#ifdef USE_GLOBAL
  USE OMP_LIB
#endif
  USE THERMO_GLOBAL
  USE DNS_GLOBAL
  USE DNS_LOCAL, ONLY : bcs_scal_jmin, bcs_scal_jmax
    
  IMPLICIT NONE

#include "integers.h"

  TINTEGER is
  TREAL dte
  TREAL, DIMENSION(*)           :: dx, dy, dz
  TREAL, DIMENSION(isize_field) :: u, v, w, hs !,s
  TREAL, DIMENSION(isize_field) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, wrk3d
  TREAL, DIMENSION(*)           :: wrk1d
  TREAL, DIMENSION(imax,kmax,2) :: wrk2d

  TREAL, DIMENSION(isize_field), INTENT(IN) :: h1,h2,h3 ! Alberto rhs of the flow equations

  TREAL, DIMENSION(isize_field,*), TARGET   :: s !Alberto. Now mulidimensinal array (including h,qt,ql)

! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip, ibc
  TREAL diff

!Alberto.Ponter to the qt,ql (2dim array) and enthalpy
  TREAL, DIMENSION(:),  POINTER :: al_qt
  TREAL, DIMENSION(:),  POINTER :: al_ql
  TREAL, DIMENSION(:),  POINTER :: al_h

  TREAL, DIMENSION(:),    POINTER ::al_s !Recieved pointer . Needed in case the mixture is not airwater

!Alberto: Parameters of the equations
  TREAL :: differential_diffusion	   
  TREAL :: dummy, dummy1, dummy2, dummy3

! #######################################################################

  nxy = imax*jmax

!Alberto: check if we have an air water mixture
  IF ( idiffusion .EQ. EQNS_NONE ) THEN; diff = C_0_R
  ELSE
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
        diff = visc/schmidt(1)  !Alberto: changed from schmidt(is) to schmidt(1) to reflect upon Prandtl = Schmidt
        al_h=>s(:,1)  !We point to the scalar array which contains the enthalpy
        al_qt=> s(:,2) ! We point to the scalar array which contains qt
        al_ql=> s(:,3)		! We point to the scalar array which contains ql
        al_s=>s(:,is) !Needed to calculate the laplacian in the part of the routine

        differential_diffusion = (schmidt(1)/schmidt(2) -C_1_R)*diff !Parameter for diff difusion


        IF (is .EQ. 1) THEN !Enthalpy diffusion
!Calculate the laplacian of ql
           CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
                dz, al_ql, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
                dy, al_ql, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
                dx, al_ql, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)


           dummy = ( THERMO_AI(1,1,3) - C_1_R + THERMO_AI(6,1,3) )* differential_diffusion
           dummy1 = stokes*buoyancy%vector(1)*( THERMO_AI(1,1,3) - C_1_R + THERMO_AI(6,1,3) )*settling!Alberto: Parameter multiplying the sedimentation term
           dummy2 = stokes*buoyancy%vector(2)*( THERMO_AI(1,1,3) - C_1_R + THERMO_AI(6,1,3) )*settling !Alberto: Parameter multiplying the sedimentation term
           dummy3 = stokes*buoyancy%vector(3)*( THERMO_AI(1,1,3) - C_1_R + THERMO_AI(6,1,3) )*settling !Alberto: Parameter multiplying the sedimentation term
!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              hs(ij) = hs(ij) + dummy*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) - dummy1*tmp1(ij) - dummy2*tmp2(ij) -dummy3*tmp3(ij) !Alberto. Liqud diff + Sedimentation term 
           ENDDO
!$omp end do
!$omp end parallel

!Start INERTIAL term


!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field
              dummy = stokes*( THERMO_AI(1,1,3) - C_1_R + THERMO_AI(6,1,3) )*(C_1_R - C_2_R*al_ql(ij))
              hs(ij) = hs(ij)  + dummy*h1(ij)*tmp1(ij) + dummy*h2(ij)*tmp2(ij) + dummy*h3(ij)*tmp3(ij) !Add the first inertial term
           ENDDO
!$omp end do
!$omp end parallel

!Second part of the inertial  term (derivative of the convective (nasty and long) term
!ux component
!Velocity derivatives
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, u, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              tmp4(ij) =  al_ql(ij)*( C_1_R - C_2_R*al_ql(ij) )*( u(ij)*tmp1(ij) + v(ij)*tmp2(ij) + w(ij)*tmp3(ij) )
           ENDDO
!$omp end do
!$omp end parallel

!Derivative of the convective term
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp4, tmp5, i0,i0, wrk1d,wrk2d,wrk3d) !Result stored in tmp5


!uy component
!Velocity derivatives
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, v, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              tmp4(ij) =  al_ql(ij)*( C_1_R - C_2_R*al_ql(ij) ) *( u(ij)*tmp1(ij) + v(ij)*tmp2(ij) + w(ij)*tmp3(ij) )
           ENDDO
!$omp end do
!$omp end parallel

!Derivative of the convective term
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp4, tmp6, i0,i0, wrk1d,wrk2d,wrk3d) !Result stored in tmp6

!uz component
!uy component
!Velocity derivatives
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, w, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              tmp4(ij) =  al_ql(ij)*( C_1_R - C_2_R*al_ql(ij) )*( u(ij)*tmp1(ij) + v(ij)*tmp2(ij) + w(ij)*tmp3(ij) )
           ENDDO
!$omp end do
!$omp end parallel
!Derivative of the convective term
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp4, tmp1, i0,i0, wrk1d,wrk2d,wrk3d) !Result stored in tmp1

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              hs(ij) = hs(ij)  + dummy*tmp1(ij) + dummy*tmp5(ij) + dummy*tmp6(ij) !Add the second inertial term
           ENDDO
!$omp end do
!$omp end parallel




        ELSE  !Liquid water diffusion
!Calculate the laplacian of ql (written twice to reuse tmp1...6)
           CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
                dz, al_ql, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
                dy, al_ql, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
                dx, al_ql, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field
              dummy = (C_1_R - al_qt(ij) + al_ql(ij) )* differential_diffusion
              dummy1 = stokes*buoyancy%vector(1)*(C_1_R - al_qt(ij) + al_ql(ij) )*settling
              dummy2 = stokes*buoyancy%vector(2)*(C_1_R - al_qt(ij) + al_ql(ij) )*settling
              dummy3 = stokes*buoyancy%vector(3)*(C_1_R - al_qt(ij) + al_ql(ij) )*settling
              hs(ij) = hs(ij) + dummy*( tmp6(ij)+tmp5(ij)+tmp4(ij)) - dummy1*tmp1(ij) - dummy2*tmp2(ij) - dummy3*tmp3(ij)!Alberto. Liqud diff + Sedimentation term
           ENDDO
!$omp end do
!$omp end parallel

!START INERTIAL term

!Derivative of total water (order ql^2) variation higher
!CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, al_qt, tmp4, i0,i0, wrk1d,wrk2d,wrk3d)
!CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, al_qt, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
!CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, al_qt, tmp6, i0,i0, wrk1d,wrk2d,wrk3d)


           dummy = stokes*( THERMO_AI(1,1,3) - C_1_R + THERMO_AI(6,1,3) )
!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field
!(order ql^2) variation higher
!dummy = (C_1_R - al_qt(ij))* stokes 
!hs(ij) = hs(ij)  + dummy*tmp1(ij)*h1(ij) + dummy*tmp2(ij)*h2(ij) + dummy*tmp3(ij)*h3(ij)  &
!		   - stokes*al_ql(ij)*tmp4(ij)*h1(ij) - stokes*al_ql(ij)*tmp5(ij)*h2(ij) - stokes*al_ql(ij)*tmp6(ij)*h3(ij)
!Approx (nabla qv) * ql sim 0
              dummy = (C_1_R - C_2_R*al_qt(ij))* stokes 
              hs(ij) = hs(ij)  + dummy*tmp1(ij)*h1(ij) + dummy*tmp2(ij)*h2(ij) + dummy*tmp3(ij)*h3(ij)  
           ENDDO
!$omp end do
!$omp end parallel

!Second part of the inertial  term (derivative of the convective (nasty and long) term
!ux component
!Velocity derivatives
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, u, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, u, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, u, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              tmp4(ij) =  al_ql(ij)*(C_1_R-al_qt(ij))*( u(ij)*tmp1(ij) + v(ij)*tmp2(ij) + w(ij)*tmp3(ij) )
           ENDDO
!$omp end do
!$omp end parallel
!Derivative of the convective term
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp4, tmp5, i0,i0, wrk1d,wrk2d,wrk3d) !Result stored in tmp5


!uy component
!Velocity derivatives
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, v, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, v, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, v, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              tmp4(ij) =  al_ql(ij)*(C_1_R-al_qt(ij))*( u(ij)*tmp1(ij) + v(ij)*tmp2(ij) + w(ij)*tmp3(ij) )
           ENDDO
!$omp end do
!$omp end parallel
!Derivative of the convective term
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp4, tmp6, i0,i0, wrk1d,wrk2d,wrk3d) !Result stored in tmp6

!uz component
!uy component
!Velocity derivatives
           CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, w, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, w, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, w, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field			
              tmp4(ij) =  al_ql(ij)*(C_1_R-al_qt(ij))*( u(ij)*tmp1(ij) + v(ij)*tmp2(ij) + w(ij)*tmp3(ij) )
           ENDDO
!$omp end do
!$omp end parallel
!Derivative of the convective term
           CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp4, tmp1, i0,i0, wrk1d,wrk2d,wrk3d) !Result stored in tmp1

!$omp parallel default( shared ) private( ij )
!$omp do
           DO ij = 1,isize_field
              dummy = (C_1_R - al_qt(ij) + al_ql(ij) )* stokes
              hs(ij) = hs(ij)  + dummy*tmp1(ij) + dummy*tmp5(ij) + dummy*tmp6(ij) !Add the second inertial term
           ENDDO
!$omp end do
!$omp end parallel

        ENDIF



     ELSE   !IF not it will have by default one single array in s.
        al_s=>s(:,1)
        diff = visc/schmidt(is)
     ENDIF
  ENDIF

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
! CALL RHS_SCAL_GLOBAL_INCOMPRESSIBLE_MPI&
!    (diff, dx,dy,dz, u,v,w,s,hs, tmp1,tmp2,tmp3,tmp4,tmp5, wrk1d,wrk2d,wrk3d)

  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, al_s, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, al_s, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, al_s, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,isize_field
     hs(ij) = hs(ij) + diff*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO
!$omp end do
!$omp end parallel

! #######################################################################
! Boundary conditions
! #######################################################################
  ibc = 0
  wrk2d(:,:,1:2) = C_0_R ! default is dirichlet
  IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
  IF ( bcs_scal_jmax(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
  IF ( ibc .GT. 0 ) THEN
     CALL BOUNDARY_BCS_NEUMANN_Y(imode_fdm,ibc, imax,jmax,kmax, dy, hs, wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,tmp1,wrk3d)
  ENDIF

! -----------------------------------------------------------------------
! Impose bottom at Jmin 
! -----------------------------------------------------------------------
  ip = 1
  DO k = 1,kmax
     hs(ip:ip+imax-1) = wrk2d(1:imax,k,1); ip = ip + nxy
  ENDDO

! -----------------------------------------------------------------------
! Impose top at Jmax
! -----------------------------------------------------------------------
  ip = 1 + imax*(jmax-1) 
  DO k = 1,kmax
     hs(ip:ip+imax-1) = wrk2d(1:imax,k,2); ip = ip + nxy
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_AIRWATER
