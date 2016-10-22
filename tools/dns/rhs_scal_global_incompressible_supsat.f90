#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2011/07/06 - A de Lozar
!#              Created form rhs_scal_global_incompressible_1
!#
!# 2012/02/13 - A de Lozar
!#              Corrected the liquid intertial flux term (factor q_g term included)
!#
!# 2012/06/04 - A de Lozar
!#              Corrected a q_g in the inertial part of the total water equation
!#
!#              If Stokes<0 Settling >0,  
!#              Added a branch in order to calculate the settling avoiding the calculation 
!#              of inertia (in order to speed the code) 
!#              IMPORTANT: This only works when the settling is in the y direction
!#              From now on Settling includes the Stokes in it. 
!########################################################################
!# DESCRIPTION
!#
!# Rhs of the scalar equation for the supersaturation formalismus
!# There are 3 independent scalar values (h,qt and ql)
!# No termodynamic equilibrium is assumed so that we have a condensation term 
!# (with 1 Damkohler number)
!# All three scalar equations are done at the same time in order to avoid 
! #calcualtion repetitions.
!# (Contrary to the rest of the code where a loop over the different rhs_scal are done)
!# Inertia (preferential condensaion), sedimentation and 
!# differential diffusivity of the liquid are included.(2 Schmidts number are needed)
!#
!# Question: we defined 7 txc but maybe it is possible to do it only with wrk3d ????????
!#
!# IMPORTNANT It is assumed that the Damkohlrt number stores damkohler*qcloud(2/3)
!# IMPORTNANT It is assumed that the stokes number stores stokes*qcloud(-2/3)
!#
!########################################################################
SUBROUTINE  RHS_SCAL_GLOBAL_INCOMPRESSIBLE_SUPSAT&
     (u,v,w, h1,h2,h3, s,hs, tmp1,tmp2,tmp3,tmp4,tmp5,tmp6, tmp7, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif

  USE DNS_GLOBAL
  USE THERMO_GLOBAL
  USE DNS_LOCAL, ONLY : bcs_scal_jmin, bcs_scal_jmax

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(isize_field) :: u, v, w
  TREAL, DIMENSION(isize_field) :: h1, h2, h3  !rhs of the flow equation
  TREAL, DIMENSION(isize_field,*), TARGET :: s !Alberto. Now mulidimensinal array (including h,qt,ql)
  TREAL, DIMENSION(isize_field,*) :: hs !Alberto. Now mulidimensinal array (including h,qt,ql)
  TREAL, DIMENSION(isize_field) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, wrk3d 
  TREAL, DIMENSION(*)           :: wrk1d
  TREAL, DIMENSION(imax,kmax,2) :: wrk2d


! -----------------------------------------------------------------------
  TINTEGER ij, k, nxy, ip, ibc, is, nscal_bound
  TREAL, DIMENSION(2)           :: diff
  TREAL dummy,exp_l,  diff_difusion, cl, Q_latent !Auxialliary variables

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: dx,dy,dz

! -----------------------------------------------------------------------
!Alberto.Ponter to the qt,ql (2dim array) and enthalpy
  TREAL, DIMENSION(:),  POINTER :: al_h
  TREAL, DIMENSION(:),  POINTER :: al_qt
  TREAL, DIMENSION(:),  POINTER :: al_ql

! #######################################################################
! Define pointers
  dx => g(1)%jac(:,1)
  dy => g(2)%jac(:,1)
  dz => g(3)%jac(:,1)

  nxy = imax*jmax

  DO k = 1,2 
     diff(k) = visc/schmidt(k); !factor 1/ReSch of the diffusion 
  ENDDO

  al_h=>s(:,1)  !We point to the scalar array which contains the enthalpy
  al_qt=> s(:,2) ! We point to the scalar array which contains qt
  al_ql=> s(:,3)  ! We point to the scalar array which contains ql

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################

! **************************
!FIRST SCALAR (ENTHALPY)
! **************************

  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, al_h, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, al_h, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, al_h, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,isize_field
     hs(ij,1) = hs(ij,1) + diff(1)*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO
!$omp end do
!$omp end parallel

! **************************
!SECOND SCALAR (TOTAL WATER)
! **************************
  CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
       dz, al_qt, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
       dy, al_qt, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
  CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
       dx, al_qt, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)
!$omp parallel default( shared ) private( ij )
!$omp do
  DO ij = 1,isize_field
     hs(ij,2) = hs(ij,2) + diff(1)*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
          - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
  ENDDO
!$omp end do
!$omp end parallel

! **************************
!THIRD SCALAR (LIQUID WATER)
! **************************

  IF (damkohler(1) .GT. C_0_R) THEN !It is only a convective variable if we allow for no eq
     CALL PARTIAL_ZZ(i1, iunifz, imode_fdm, imax,jmax,kmax, k1bc,&
          dz, al_ql, tmp6, i0,i0, i0,i0, tmp3, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_YY(i1, iunify, imode_fdm, imax,jmax,kmax, j1bc,&
          dy, al_ql, tmp5, i0,i0, i0,i0, tmp2, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_XX(i1, iunifx, imode_fdm, imax,jmax,kmax, i1bc,&
          dx, al_ql, tmp4, i0,i0, i0,i0, tmp1, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,3) = hs(ij,3) + diff(2)*( tmp6(ij)+tmp5(ij)+tmp4(ij) ) &
             - ( w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij) )
     ENDDO
!$omp end do
!$omp end parallel

  ELSE
     IF  (schmidt(1) .NE.  schmidt(2) )   THEN !We need the gradient of ql in tmp1, tmp2, tmp3 for the differential diffusion
        CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, al_ql, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, al_ql, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, al_ql, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)
     ENDIF
  ENDIF

! #######################################################################
! Differential Diffusion total water equation
! First the fluxes are added in tmp4, tmp5, tmp6, tmp7
!We use that the gradient of ql are in tmp1, tmp2, tmp3 (these have to be untouched to be used by the next routine)
! Second the divergence of the flux is done and the terms are added to the rhs
! This is done in two steps (first x_divergenge and then y and z divergecnes) to save the number of arrays needed
! #######################################################################

  IF (schmidt(1) .NE.  schmidt(2) ) THEN

! **************************
!SECOND SCALAR (TOTAL WATER)
! **************************
!CREATE THE X FLUX TERM in tmp4

!$omp parallel default( shared ) private( ij,diff_difusion )
     diff_difusion = diff(2) - diff(1)

!$omp do
     DO ij = 1,isize_field      
        tmp6(ij) = ( C_1_R-al_qt(ij) )*diff_difusion/(C_1_R-al_ql(ij)) !We save some little time putting the multiplying factor in tmp6
        tmp4(ij) = tmp6(ij)*tmp1(ij)
     ENDDO
!$omp end do
!$omp end parallel


!x_Divergence in tmp5
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp4, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs
!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,2) = hs(ij,2) + tmp5(ij)
     ENDDO
!$omp end do
!$omp end parallel

!CREATE THE Y,Z FLUX TERMS in tmp4, tmp5

!$omp parallel default( shared ) private( ij) 
!$omp do
     DO ij = 1,isize_field      
        tmp4(ij) = tmp6(ij)*tmp2(ij)
        tmp6(ij) = tmp6(ij)*tmp3(ij)
     ENDDO
!$omp end do
!$omp end parallel

!x_Divergence in tmp5 and tmp7 
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp4, tmp5, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp6, tmp7, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs
!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,2) = hs(ij,2) + tmp5(ij) + tmp7(ij)
     ENDDO
!$omp end do
!$omp end parallel

  ENDIF

! #######################################################################
! CONDENDSATION EQUATION.
! The condensation term is calculated and then added to the rhs
! From this point tmp7 stores the temperature
! IMPORTNANT It is assumed that the Damkohlrt number stores damkohler*qcloud(1/3)
! #######################################################################
  IF (damkohler(1) .GT. C_0_R) THEN !It supersatturation is allowed

     CALL THERMO_AIRWATER_QSAT(imax,jmax,kmax, s(:,2:3),p_init,al_h, tmp7, tmp4) !Temperature in tmp7 and saturation concentration in tmp4 

!(IN THE FUTURE THIS LOOP CAN BE DONE IN THERMO_AIRWATER_QSAT)
!$omp parallel default( shared ) private( ij, dummy, exp_l )
     dummy = C_3_R*damkohler(1)
     exp_l = C_1_R/C_3_R
!$omp do
     DO ij = 1,isize_field !Condensation term in the rhs
        hs(ij,3) = hs(ij,3) + dummy*( ( al_qt(ij)-al_ql(ij) )/tmp4(ij) - C_1_R )*(al_ql(ij)**exp_l)
     ENDDO
!$omp end do
!$omp end parallel


  ELSE 
     IF ( (schmidt(1) .NE.  schmidt(2) ) .OR. (stokes .GT. C_0_R) .OR. (settling .GT. C_0_R) ) THEN
!Put temperature in tmp7 (mqn_rho, tmp4 and wrk3d are in this case unused)
        CALL THERMO_CALORIC_TEMPERATURE(imax, jmax, kmax, s, tmp4, mean_rho, tmp7, wrk3d) 
     ENDIF
  ENDIF

! #######################################################################
! Differential Diffusion entahlpy water equation
! First the fluxes are added in tmp4, tmp5, tmp6
!We use that the gradient of ql are in tmp1, tmp2, tmp3
!Second the divergence of the flux is done and the terms are added to the rhs
! #######################################################################
  IF (schmidt(1) .NE.  schmidt(2) ) THEN

! **************************
!FIRST SCALAR (ENTHALPY)
! **************************
!CREATE THE FLUX TERM in tmp4, tmp5, tmp6

!$omp parallel default( shared ) private( ij,diff_difusion, cl,Q_latent,dummy )
     diff_difusion = diff(2) - diff(1)
     cl = THERMO_AI(1,1,3) !Specific heats
     Q_latent = THERMO_AI(6,1,3) !latent heat

!$omp do
     DO ij = 1,isize_field
        tmp7(ij) = cl*tmp7(ij) + Q_latent  - al_h(ij)!From now on tmp7 stores hl -h
        dummy = tmp7(ij)*diff_difusion/(C_1_R-al_ql(ij))
        tmp4(ij) = dummy*tmp1(ij)
        tmp5(ij) = dummy*tmp2(ij)
        tmp6(ij) = dummy*tmp3(ij)
     ENDDO
!$omp end do
!$omp end parallel

!Calculate the gradient of the flux in tmp1, tmp2, tmp3
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp4, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp5, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp6, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,1) = hs(ij,1) + tmp1(ij) + tmp2(ij) + tmp3(ij)
     ENDDO
!$omp end do
!$omp end parallel
  ELSE

     IF( (stokes .GT. C_0_R) .OR. (settling .GT. 0) ) THEN
!$omp parallel default( shared ) private( ij,cl,Q_latent )
        cl = THERMO_AI(1,1,3) !Specific heats
        Q_latent = THERMO_AI(6,1,3) !latent heat
!$omp do
        DO ij = 1,isize_field
           tmp7(ij) = cl*tmp7(ij) + Q_latent -al_h(ij) !From now on tmp7 stores hl -h
        ENDDO
!$omp end do
!$omp end parallel
     ENDIF
  ENDIF

  IF(stokes .GT. C_0_R) THEN
! #######################################################################
! INERTIA EQUATION
! First calculate the drift velocity and store it in tmp1,tmp2,tmp3
! IMPORTNANT It is assumed that the stokes number stores stokes*qcloud(-2/3)
! Again we calculate first the flux term and then we do the gradient
! #######################################################################

!NEgative Drift velocity DuDT -g strored in tmp4, tmp5, tmp6
     CALL FI_TRANS_VELOCITY(dx,dy,dz, u,v,w, h1,h2,h3, tmp4,tmp5,tmp6,tmp1,tmp2,tmp3, wrk1d,wrk2d,wrk3d) 


! **************************
!THIRD SCALAR (LIQUID WATER)
! **************************
!CREATE THE FLUX TERM in tmp4, tmp5, tmp6 (we do not store the drift any more because this term  is common to all three equations)
!Except for a q_g term which is taken out later
!$omp parallel default( shared ) private( ij,dummy,exp_l )
     exp_l = C_5_R/C_3_R
!$omp do
     DO ij = 1,isize_field
        dummy = (C_1_R-al_ql(ij))*al_ql(ij)**exp_l
        tmp4(ij) = dummy*tmp4(ij)
        tmp5(ij) = dummy*tmp5(ij)
        tmp6(ij) = dummy*tmp6(ij)
     ENDDO
!$omp end do
!$omp end parallel

     IF (damkohler(1) .GT. C_0_R) THEN !Calculate only the third equation if supersaturation is allowed
!Calculate the gradient of the flux in tmp1, tmp2, tmp3
        CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp4, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp5, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
        CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp6, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
        DO ij = 1,isize_field
           hs(ij,3) = hs(ij,3) + tmp1(ij) + tmp2(ij) + tmp3(ij)
        ENDDO
!$omp end do
!$omp end parallel
     ENDIF

! **************************
!FIRST SCALAR (ENTHALPY)
! **************************
!We have to split the loop due to the lack of tmps. Remember hl -h stored in tmp7
! x component of the flux . q_g is removed from the inertial flux
!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        tmp3(ij) = tmp7(ij)*tmp4(ij)/(C_1_R-al_ql(ij))
     ENDDO
!$omp end do
!$omp end parallel

!Calculate the gradient of the flux in tmp1
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp3, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,1) = hs(ij,1) + tmp1(ij)
     ENDDO
!$omp end do
!$omp end parallel

! y,z component of the flux (hl-h) is not needed any more and can be reuse tmp7.
! q_g is removed from the flux
!$omp parallel default( shared ) private( ij, dummy )
!$omp do
     DO ij = 1,isize_field
        dummy = (C_1_R-al_ql(ij))
        tmp3(ij) = tmp7(ij)*tmp5(ij)/dummy
        tmp7(ij) = tmp7(ij)*tmp6(ij)/dummy
     ENDDO
!$omp end do
!$omp end parallel

!Calculate the gradient of the flux in tmp1, temp2
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp3, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp7, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)


!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,1) = hs(ij,1) + tmp1(ij) + tmp2(ij)
     ENDDO
!$omp end do
!$omp end parallel

! **************************
!SECOND SCALAR (TOTAL WATER)
! **************************
! Now we can use that the drift velocities are not needed any more and do everithing in one loop

!$omp parallel default( shared ) private( ij,dummy )
!$omp do
     DO ij = 1,isize_field
        dummy = (C_1_R - al_qt(ij))/(C_1_R - al_ql(ij))
        tmp4(ij) = dummy*tmp4(ij)
        tmp5(ij) = dummy*tmp5(ij)
        tmp6(ij) = dummy*tmp6(ij)
     ENDDO
!$omp end do
!$omp end parallel

!Calculate the gradient of the flux in tmp1, tmp2, tmp3
     CALL PARTIAL_X(imode_fdm, imax,jmax,kmax, i1bc, dx, tmp4, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp5, tmp2, i0,i0, wrk1d,wrk2d,wrk3d)
     CALL PARTIAL_Z(imode_fdm, imax,jmax,kmax, k1bc, dz, tmp6, tmp3, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,2) = hs(ij,2) + tmp1(ij) + tmp2(ij) + tmp3(ij)
     ENDDO
!$omp end do
!$omp end parallel

! *********************************
!For the case that we want settling velocity but not Inertia. In this case ONLY DIRECTION Y  is considered. Maybe it could be improved to use less TEMP 
!**********************************
  ELSEIF (settling > C_0_R) THEN   

! **************************
!THIRD SCALAR (LIQUID WATER)
! **************************
!CREATE THE FLUX TERM in tmp4, tmp5, tmp6 (we do not store the drift any more because this term  is common to all three equations)
!Except for a q_g term which is taken out later
!$omp parallel default( shared ) private( ij,dummy,exp_l )
     exp_l = C_5_R/C_3_R
     dummy = -settling*buoyancy%vector(2)
!$omp do
     DO ij = 1,isize_field
        tmp4(ij) = dummy*(C_1_R-al_ql(ij))*al_ql(ij)**exp_l      
     ENDDO
!$omp end do
!$omp end parallel

     IF (damkohler(1) .GT. C_0_R) THEN !Calculate only the third equation if supersaturation is allowed
!Calculate the gradient of the flux in tmp2
        CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp4, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
        DO ij = 1,isize_field
           hs(ij,3) = hs(ij,3) + tmp1(ij)
        ENDDO
!$omp end do
!$omp end parallel
     ENDIF

! **************************
!FIRST SCALAR (ENTHALPY)
! **************************
!We have to split the loop due to the lack of tmps. Remember hl -h stored in tmp7
! x component of the flux . q_g is removed from the inertial flux
!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        tmp3(ij) = tmp7(ij)*tmp4(ij)/(C_1_R-al_ql(ij))
     ENDDO
!$omp end do
!$omp end parallel

!Calculate the gradient of the flux in tmp1
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp3, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)
!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,1) = hs(ij,1) + tmp1(ij)
     ENDDO
!$omp end do
!$omp end parallel

! **************************
!SECOND SCALAR (TOTAL WATER)
! **************************
! Now we can use that the drift velocities are not needed any more and do everithing in one loop

!$omp parallel default( shared ) private( ij,dummy )
!$omp do
     DO ij = 1,isize_field
        dummy = (C_1_R - al_qt(ij))/(C_1_R - al_ql(ij))
        tmp4(ij) = dummy*tmp4(ij)
     ENDDO
!$omp end do
!$omp end parallel

!Calculate the gradient of the flux in tmp1, tmp2, tmp3
     CALL PARTIAL_Y(imode_fdm, imax,jmax,kmax, j1bc, dy, tmp4, tmp1, i0,i0, wrk1d,wrk2d,wrk3d)

!Add to the rhs

!$omp parallel default( shared ) private( ij )
!$omp do
     DO ij = 1,isize_field
        hs(ij,2) = hs(ij,2) + tmp1(ij) 
     ENDDO
!$omp end do
!$omp end parallel
  ENDIF
! #######################################################################
! Boundary conditions ALBERTO THINKA ABOUT !!!!! -> LOOP OVER DIMESIONS CHANGING hs->hs(:,cdim)
! #######################################################################

!Check how many boundary conditions have to do
  IF (damkohler(1) .GT. C_0_R) THEN 
     nscal_bound = 3
  ELSE
     nscal_bound = 2
  ENDIF

  DO is =1,nscal_bound
     ibc = 0
     wrk2d(:,:,1:2) = C_0_R ! default is dirichlet
     IF ( bcs_scal_jmin(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 1
     IF ( bcs_scal_jmax(is) .EQ. DNS_BCS_NEUMANN ) ibc = ibc + 2
     IF ( ibc .GT. 0 ) THEN
        CALL BOUNDARY_BCS_NEUMANN_Y(ibc, imax,jmax,kmax, g(2), hs(:,is), wrk2d(1,1,1),wrk2d(1,1,2), wrk1d,tmp1,wrk3d)
     ENDIF

! -----------------------------------------------------------------------
! Impose bottom at Jmin 
! -----------------------------------------------------------------------
     ip = 1
     DO k = 1,kmax
        hs(ip:ip+imax-1,is) = wrk2d(1:imax,k,1); ip = ip + nxy
     ENDDO

! -----------------------------------------------------------------------
! Impose top at Jmax
! -----------------------------------------------------------------------
     ip = 1 + imax*(jmax-1) 
     DO k = 1,kmax
        hs(ip:ip+imax-1,is) = wrk2d(1:imax,k,2); ip = ip + nxy
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RHS_SCAL_GLOBAL_INCOMPRESSIBLE_SUPSAT
