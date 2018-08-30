#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2003/06/12 - J.P. Mellado
!#              Created
!# 2007/08/16 - J.P. Mellado
!#              Case of internal energy formulation added
!# 2007/08/30 - J.P. Mellado
!#              Skewsymmetric formulation added.
!#              Explicit viscous term formulation added.
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# txc   In    3D auxiliar array of size 6 or 9
!#
!########################################################################
SUBROUTINE TIME_SUBSTEP_COMPRESSIBLE(dte, etime, q,hq, s,hs, q_inf,s_inf, txc, wrk1d,wrk2d,wrk3d)

#ifdef USE_OPENMP
  USE OMP_LIB
#endif
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0
  USE BOUNDARY_BUFFER
  USE BOUNDARY_BCS
#ifdef LES
  USE DNS_LOCAL, ONLY : rkm_substep
#endif
#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL, ONLY : iles
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL dte, etime

  TREAL, DIMENSION(isize_field,*) :: q, hq, s, hs, txc
  TREAL, DIMENSION(*)             :: q_inf, s_inf
  TREAL, DIMENSION(*)             :: wrk1d, wrk2d, wrk3d

  TARGET :: q, hq

! -------------------------------------------------------------------
  TREAL rho_ratio, dt_rho_ratio, prefactor
  TREAL M2_max, dummy
  TINTEGER i, is, inb_scal_loc

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w, e, rho, p, T, vis
  TREAL, DIMENSION(:), POINTER :: h0, h1, h2, h3, h4

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING TIME_SUBSTEP')
#endif

! Define pointers
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)

  e   => q(:,4)
  rho => q(:,5)
  p   => q(:,6)
  T   => q(:,7)

  vis => q(:,8)

  h1 => hq(:,1)
  h2 => hq(:,2)
  h3 => hq(:,3)
  h4 => hq(:,4)
  h0 => hq(:,5)
 
! ###################################################################
! Evaluate standard RHS of equations
! global formulation
! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL  .AND. &
       iadvection .EQ. EQNS_SKEWSYMMETRIC .AND. &
       iviscous   .EQ. EQNS_EXPLICIT      .AND. &
       idiffusion .EQ. EQNS_EXPLICIT            ) THEN
     CALL RHS_FLOW_GLOBAL_2(rho,u,v,w,p,e,T,s, h0,h1,h2,h3,h4,hs,&
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
     DO is = 1,inb_scal
        CALL RHS_SCAL_GLOBAL_2(is, rho,u,v,w,s,T, hs, h4,&
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
     ENDDO

  ELSE
! ###################################################################
! Evaluate standard RHS of equations
! split formulation
! ###################################################################
! -------------------------------------------------------------------
! convective terms
! -------------------------------------------------------------------
     IF      ( iadvection .EQ. EQNS_DIVERGENCE    ) THEN
        CALL RHS_FLOW_EULER_DIVERGENCE(rho,u,v,w,p,e, h0,h1,h2,h3,h4,&
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
        DO is = 1,inb_scal
           CALL RHS_SCAL_EULER_DIVERGENCE(rho,u,v,w,s(1,is), hs(1,is),&
                txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
        ENDDO
     ELSE IF ( iadvection .EQ. EQNS_SKEWSYMMETRIC ) THEN
        CALL RHS_FLOW_EULER_SKEWSYMMETRIC(rho,u,v,w,p,e,s, h0,h1,h2,h3,h4,hs,&
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
        DO is = 1,inb_scal
           CALL RHS_SCAL_EULER_SKEWSYMMETRIC(rho,u,v,w,s(1,is), hs(1,is),&
                txc(1,1),txc(1,2),txc(1,3),txc(1,4), wrk2d,wrk3d)
        ENDDO
     ENDIF

! -------------------------------------------------------------------
! viscous terms
! -------------------------------------------------------------------
     IF ( itransport .NE. 1 ) THEN
        CALL IO_WRITE_ASCII(efile,'TIME_SUBSTEP_COMPRESSIBLE. Section requires to allocate array vis.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

     IF      ( iviscous .EQ. EQNS_DIVERGENCE ) THEN
        CALL RHS_FLOW_VISCOUS_DIVERGENCE(vis, u,v,w,p, h1,h2,h3,h4, &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8),txc(1,9), wrk2d,wrk3d)
     ELSE IF ( iviscous .EQ. EQNS_EXPLICIT   ) THEN
        CALL RHS_FLOW_VISCOUS_EXPLICIT(vis, u,v,w,p, h1,h2,h3,h4, &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
     ENDIF

! -------------------------------------------------------------------
! diffusion/conduction terms
! -------------------------------------------------------------------
     IF      ( idiffusion .EQ. EQNS_DIVERGENCE ) THEN
! diffusion transport of enthalpy is accumulated in txc and used then in RHS_FLOW_CONDUCTION
        txc(:,1) = C_0_R; txc(:,2) = C_0_R; txc(:,3) = C_0_R
        DO is = 1,inb_scal
           CALL RHS_SCAL_DIFFUSION_DIVERGENCE(is, vis, s, T, hs, &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
        ENDDO
        CALL RHS_FLOW_CONDUCTION_DIVERGENCE(vis, s, T, h4, &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), wrk2d,wrk3d)
     ELSE IF ( idiffusion .EQ. EQNS_EXPLICIT   ) THEN
        DO is = 1,inb_scal
           CALL RHS_SCAL_DIFFUSION_EXPLICIT(is, vis, s, T, hs, h4, &
                txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6), wrk2d,wrk3d)
        ENDDO
        CALL RHS_FLOW_CONDUCTION_EXPLICIT(vis, s, T, h4, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk2d,wrk3d)
     ENDIF

  ENDIF

! ###################################################################
! Evaluate chemical RHS
! If LES is on, each of these subroutines should have its
! counterpart in LES library
! NOT YET DEVELOPED FOR THE ENERGY FORMULATION !!!
! ###################################################################
#ifdef CHEMISTRY
  IF ( iles .EQ. 0 ) THEN
     IF ( icalc_scal .EQ. 1 ) THEN
        IF      ( imech_type .EQ. CHEM_TYPE_PETERS1991  ) THEN
           CALL CHEM_PETERS1991(i0, rho, T, gama, s, hs, h4, txc)
        ELSE IF ( imech_type .EQ. CHEM_TYPE_PETERS1988  ) THEN
           CALL CHEM_PETERS1988(rho, T, gama, s, hs, h4)
        ELSE IF ( imech_type .EQ. CHEM_TYPE_UNIDECOMP   ) THEN
           CALL CHEM_UNIDECOMP(rho, T, gama, s, hs, h4)
        ELSE IF ( imech_type .EQ. CHEM_TYPE_BSZELDOVICH ) THEN
           CALL CHEM_BSZELDOVICH(rho, T, gama, s, hs, h4)
        ELSE IF ( imech_type .EQ. CHEM_TYPE_ONESTEP     ) THEN
           CALL CHEM_ONESTEP(rho, T, gama, s, hs, h4)
        ELSE IF ( imech_type .EQ. CHEM_TYPE_QUASIBS     ) THEN
           CALL CHEM_QUASIBS(rho, T, gama, s, hs, h4)
        ENDIF
     ENDIF
  ENDIF
#endif

! ###################################################################
! Evaluate LES terms of the RHS of equations
! ###################################################################
#ifdef LES
  IF ( iles .EQ. 1 ) THEN
     CALL LES_RHS(rkm_substep, q,hq, s,hs, txc, vaux, wrk1d,wrk2d,wrk3d)
  ENDIF
#endif

! ###################################################################
! Impose boundary conditions
! Temperature array T is used as auxiliary array because it is no 
! longer used until the fields are updated
! ###################################################################
#define GAMMA_LOC(i) txc(i,6)
#define AUX_LOC(i)   T(i)

  CALL THERMO_GAMMA(imax, jmax, kmax, s, T, GAMMA_LOC(1))

! Maximum Mach for Poinsot & Lele reference pressure BC
  IF ( bcs_euler_drift .EQ. 1 ) THEN
     M2_max = C_0_R
     DO i = 1,isize_field
        dummy = (u(i)*u(i)+v(i)*v(i)+w(i)*w(i))*rho(i)/(GAMMA_LOC(i)*p(i))
        M2_max = MAX(M2_max,dummy)
     ENDDO
#ifdef USE_MPI
     CALL MPI_ALLREDUCE(M2_max, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
     M2_max = dummy
#endif
  ENDIF

  IF ( .NOT. g(2)%periodic ) THEN
     CALL BOUNDARY_BCS_Y(isize_field, M2_max,        rho,u,v,w,p,GAMMA_LOC(1),s, &
          h0,h1,h2,h3,h4,hs,&
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), AUX_LOC(1), wrk2d,wrk3d)
  ENDIF

  IF ( .NOT. g(1)%periodic ) THEN
     CALL BOUNDARY_BCS_X(isize_field, M2_max, etime, rho,u,v,w,p,GAMMA_LOC(1),s, &
          q_inf,s_inf, &
          h0,h1,h2,h3,h4,hs, txc, AUX_LOC(1), wrk2d,wrk3d)
  ENDIF

#undef GAMMA_LOC
#undef AUX_LOC

! ###################################################################
! Impose buffer zone as relaxation terms
! ###################################################################
  IF ( BuffType .EQ. DNS_BUFFER_RELAX .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
     CALL BOUNDARY_BUFFER_RELAXATION_FLOW(q, hq)
     DO is = 1,inb_scal
        CALL BOUNDARY_BUFFER_RELAXATION_SCAL(is, rho,s(1,is), hs(1,is)) 
     ENDDO
  ENDIF

! ###################################################################
! Perform the time stepping 
! ###################################################################
  rho_ratio = C_1_R
  prefactor = (gama0-C_1_R)*mach*mach

  IF ( icalc_flow .EQ. 1 ) THEN
     IF ( icalc_scal .EQ. 1 ) THEN; inb_scal_loc = inb_scal
     ELSE;                          inb_scal_loc = 0;       ENDIF

! -------------------------------------------------------------------
! Total energy formulation
! -------------------------------------------------------------------
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
!$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
!$omp do
        DO i = 1,isize_field
           rho_ratio    = rho(i)
           rho(i)       = rho(i) + dte*h0(i)
           rho_ratio    = rho_ratio/rho(i)
           dt_rho_ratio = dte/rho(i)

           e(i) = rho_ratio*( e(i) + prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i)) ) &
                + dt_rho_ratio*h4(i)

           u(i) = rho_ratio*u(i) + dt_rho_ratio*h1(i)
           v(i) = rho_ratio*v(i) + dt_rho_ratio*h2(i)
           w(i) = rho_ratio*w(i) + dt_rho_ratio*h3(i)

           e(i) = e(i) - prefactor*C_05_R*(u(i)*u(i)+v(i)*v(i)+w(i)*w(i))

           DO is = 1,inb_scal_loc
              s(i,is) = rho_ratio*s(i,is) + dt_rho_ratio*hs(i,is)
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! Internal energy formulation
! -------------------------------------------------------------------
     ELSE IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
!$omp parallel default( shared ) private( i, is, rho_ratio, dt_rho_ratio )
!$omp do
        DO i = 1,isize_field
           rho_ratio    = rho(i)
           rho(i)       = rho(i) + dte*h0(i)
           rho_ratio    = rho_ratio/rho(i)
           dt_rho_ratio = dte/rho(i)

           e(i) = rho_ratio*e(i) + dt_rho_ratio*h4(i)
           u(i) = rho_ratio*u(i) + dt_rho_ratio*h1(i)
           v(i) = rho_ratio*v(i) + dt_rho_ratio*h2(i)
           w(i) = rho_ratio*w(i) + dt_rho_ratio*h3(i)

           DO is = 1,inb_scal_loc
              s(i,is) = rho_ratio*s(i,is) + dt_rho_ratio*hs(i,is)
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel

     ENDIF

  ELSE
     IF ( icalc_scal .EQ. 1 ) THEN
        DO is = 1,inb_scal
!$omp parallel default( shared ) private( i, dt_rho_ratio )
!$omp do
           DO i = 1,isize_field
              dt_rho_ratio = dte/rho(i)
              s(i,is) = rho_ratio*s(i,is) + dt_rho_ratio*hs(i,is)
           ENDDO
!$omp end do
!$omp end parallel
        ENDDO
     ENDIF
  ENDIF

! ###################################################################
! Impose buffer zone as filter
! ###################################################################
  IF ( BuffType .EQ. DNS_BUFFER_FILTER .OR. BuffType .EQ. DNS_BUFFER_BOTH ) THEN
! CALL BOUNDARY_BUFFER_FILTER&
!      (rho,u,v,w,e,s, txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
! OPR_FILTER has changed and this routine needs to be updated
  ENDIF

! ###################################################################
! Calculate other intensive thermodynamic variables
! ###################################################################
  CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
  CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING TIME_SUBSTEP')
#endif

  RETURN
END SUBROUTINE TIME_SUBSTEP_COMPRESSIBLE
