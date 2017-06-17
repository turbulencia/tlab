#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef LES
#include "les_const.h"
#endif

!########################################################################
!# DESCRIPTION
!#
!# Determine the variable time step is a positive cfl is given
!# If negative cfl is prescribed, then constant dtime and this routine
!# calculates CFL and diffustion numbers for log files
!#
!# The reacting case should be reviewed in the new formulation in terms
!# of energy.
!#
!# The diffusion number is fixed in terms of the CFL, which is the input.
!# This depends on the scheme used. From Lele (1992), page 32, we have that, if 
!# the sixth order tridiagonal scheme is used, then the maximum CFL number 
!# for a 4RK is 2.9/1.989, about 1.43. For the (5)4RK from CarpenterKennedy1994 
!# used here we have 3.36/1.989, about 1.69.
!# This holds for periodic case. A safety margin leads to the common value of 1.2.
!#
!# If second order finite different operator is used, then the maximum
!# diffusion number is 2.9/6.857, about 0.42.
!# For the (5)4RK from CarpenterKennedy1994 it is 4.639/6.857 = 0.68
!#
!# If twice the first order finite difference operator is used, then the
!# maximum diffusion number is 2.9/1.989^2, about 0.73.
!# For the (5)4RK from CarpenterKennedy1994 it is 4.639/1.989^2 = 1.17
!#
!# In incompressible mode the arrays rho, p and vis are not used
!#
!########################################################################
SUBROUTINE TIME_COURANT(q,s, wrk3d)

  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_scal, imode_eqns
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : itransport, visc, prandtl, schmidt
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL,     ONLY : cfl, dtime, rkm_mode, logs_data
#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL
#endif
#ifdef CHEMISTRY
  USE CHEM_GLOBAL, ONLY : TGFM
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN)    :: q, s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: wrk3d
  
  TARGET :: q

! -------------------------------------------------------------------
  TINTEGER i,j,k, kdsp, idsp
  TREAL dt_loc
  TREAL cmax,dmax,rmax, dtc,dtd,dtr, cfld,cflr
  TREAL schmidtfactor, viscles, dummy
#ifdef CHEMISTRY
  TREAL mtgfm, zmin, zmax, umin, umax
  TREAL vmin, vmax, wmin, wmax
#endif

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u, v, w, rho, p, vis

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING TIME_COURANT' )
#endif

  u   => q(:,:,:,1)
  v   => q(:,:,:,2)
  w   => q(:,:,:,3)

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     rho => q(:,:,:,5)
     p   => q(:,:,:,6)
     
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,:,:,8)
     
  ENDIF

#ifdef USE_MPI
  idsp = ims_offset_i 
  kdsp = ims_offset_k 
#else 
  idsp = 0
  kdsp = 0
#endif

! diffusion number set equal to a factor of the given CFL number
! A factor 1/4 is used. See header of the routine. 
  cfld = C_025_R*cfl
  cflr = C_05_R *cfl ! this one for the reaction I'm not sure

! So that the minimum non-zero determines dt at the end
  dtc = C_BIG_R 
  dtd = C_BIG_R 
  dtr = C_BIG_R 

! ###################################################################
! CFL number condition
! ###################################################################
  cmax = C_0_R

! -------------------------------------------------------------------
! Incompressible
! Calculate global maximum of u/dx + v/dy + w/dz 
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     IF ( g(3)%size .GT. 1 ) THEN
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = ABS(u(i,j,k)) *g(1)%jac(i+idsp,3) &
                        + ABS(v(i,j,k)) *g(2)%jac(j,3)      &
                        + ABS(w(i,j,k)) *g(3)%jac(k+kdsp,3)
        ENDDO; ENDDO; ENDDO
     ELSE
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = ABS(u(i,j,k)) *g(1)%jac(i+idsp,3) &
                        + ABS(v(i,j,k)) *g(2)%jac(j,3)
        ENDDO; ENDDO; ENDDO
     ENDIF
     
! -------------------------------------------------------------------
! Compressible
! Calculate global maximum of (u+c)/dx + (v+c)/dy + (w+c)/dz 
! -------------------------------------------------------------------
  ELSE
     wrk3d = SQRT(gama0*p/rho) ! sound speed; positiveness of p and rho checked in routine DNS_CONTROL
     IF ( g(3)%size .GT. 1 ) THEN
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = (ABS(u(i,j,k))+wrk3d(i,j,k)) *g(1)%jac(i+idsp,3) &
                        + (ABS(v(i,j,k))+wrk3d(i,j,k)) *g(2)%jac(j,3)      &
                        + (ABS(w(i,j,k))+wrk3d(i,j,k)) *g(3)%jac(k+kdsp,3)
        ENDDO; ENDDO; ENDDO
     ELSE
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = (ABS(u(i,j,k))+wrk3d(i,j,k)) *g(1)%jac(i+idsp,3) &
                        + (ABS(v(i,j,k))+wrk3d(i,j,k)) *g(2)%jac(j,3)
        ENDDO; ENDDO; ENDDO
     ENDIF

  ENDIF

  cmax = MAXVAL(wrk3d)
  
#ifdef USE_MPI
  CALL MPI_ALLREDUCE(cmax, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  cmax = dummy
#endif

  IF ( cmax .GT. C_0_R ) dtc = cfl/cmax ! Set time step for the given CFL number
     
! ###################################################################
! Diffusion number condition
! ###################################################################
  dmax = C_0_R

! maximum diffusivities
  schmidtfactor = C_1_R
  dummy         = C_1_R/prandtl
  schmidtfactor = MAX(schmidtfactor, dummy)
  dummy         = C_1_R/MINVAL(schmidt(1:inb_scal))
  schmidtfactor = MAX(schmidtfactor, dummy)

#ifdef LES
  IF ( iles .EQ. 1 ) THEN
     IF ( iles_type_regu .EQ. LES_REGU_NONE ) THEN
        viscles = C_0_R
     ELSE IF ( iles_type_regu.EQ.LES_REGU_SMGDYN .OR. iles_type_regu.EQ.LES_REGU_SMGDYNRMS ) THEN
        viscles = MAX(smg_udiff,smg_ediff)
        IF ( icalc_scal .EQ. 1 ) viscles = MAX(viscles,smg_zdiff)
     ELSE IF ( iles_type_regu .EQ. LES_REGU_SMGSTA ) THEN
        viscles = C_1_R
        dummy = C_1_R/smg_prandtl
        viscles = MAX(viscles, dummy)
        dummy = C_1_R/smg_schmidt
        viscles = MAX(viscles, dummy)
        viscles = viscles * smg_udiff
     ENDIF
     IF ( iviscous .EQ. EQNS_NONE ) schmidtfactor = C_0_R
  ELSE
#endif
     viscles = C_0_R
#ifdef LES
  ENDIF
#endif

! -------------------------------------------------------------------
! Incompressible
! Calculate global maximum of \nu*(1/dx^2 + 1/dy^2 + 1/dz^2)
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     IF ( g(3)%size .GT. 1 ) THEN
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = g(1)%jac(i+idsp,4) + g(2)%jac(j,4) + g(3)%jac(k+kdsp,4)
        ENDDO; ENDDO; ENDDO
     ELSE
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = g(1)%jac(i+idsp,4) + g(2)%jac(j,4)
        ENDDO; ENDDO; ENDDO
     ENDIF
     dmax = (schmidtfactor*visc+viscles) *MAXVAL(wrk3d)

! -------------------------------------------------------------------
! Compressible
! Calculate global maximum of \mu/rho*(1/dx^2 + 1/dy^2 + 1/dz^2)
! -------------------------------------------------------------------
  ELSE
     IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        IF ( g(3)%size .GT. 1 ) THEN
           DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
              wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) + g(3)%jac(k+kdsp,4) )&
                            *( schmidtfactor *visc *vis(i,j,k) /rho(i,j,k) +viscles )
           ENDDO; ENDDO; ENDDO
        ELSE
           DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
              wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) )&
                            *( schmidtfactor *visc *vis(i,j,k) /rho(i,j,k) +viscles )
           ENDDO; ENDDO; ENDDO
        ENDIF
        
     ELSE ! constant dynamic viscosity
        IF ( g(3)%size .GT. 1 ) THEN
           DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
              wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) + g(3)%jac(k+kdsp,4) )&
                            *( schmidtfactor *visc /rho(i,j,k) +viscles )
           ENDDO; ENDDO; ENDDO
        ELSE
           DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
              wrk3d(i,j,k) = ( g(1)%jac(i+idsp,4) + g(2)%jac(j,4) )&
                            *( schmidtfactor *visc /rho(i,j,k) +viscles )
           ENDDO; ENDDO; ENDDO
        ENDIF

     ENDIF
     
     dmax = MAXVAL(wrk3d)

  ENDIF

#ifdef USE_MPI
  CALL MPI_ALLREDUCE(dmax, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  dmax = dummy
#endif

  IF ( dmax .GT. C_0_R ) dtd = cfld/dmax ! Set time step for the given diffusion number

! ###################################################################
! Reacting time step control
! ###################################################################
  rmax = C_0_R

#ifdef CHEMISTRY
! -------------------------------------------------------------------
! Infinitely fast
! -------------------------------------------------------------------
  IF ( ireactive .EQ. CHEM_INFINITE ) THEN
     CALL MINMAX(imax,jmax,kmax, u, umin,umax)
     CALL MINMAX(imax,jmax,kmax, v, vmin,vmax)
     CALL MINMAX(imax,jmax,kmax, w, wmin,wmax)
     CALL MINMAX(imax,jmax,kmax, s(1,1,1,inb_scal), zmin,zmax)

     umax = umax-umin
     vmax = vmax-vmin
     wmax = wmax-wmin

     umax = MAX(umax, vmax)
     umax = MAX(umax, wmax)

     mtgfm = TGFM*(zmax-zmin)/umax

     IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        DO k=1, kmax
           DO ij=1, imax*jmax
              index = ij-1
              wrk2d(ij,1) = C_1_R/(g(1)%jac(MOD(index,imax)+1+idsp,1)**3)
              wrk2d(ij,1) = wrk2d(ij,1) + C_1_R/(g(2)%jac((ij-1)/imax+1,1)**3)
              wrk2d(ij,1) = (mtgfm*vis(ij,1,k)/rho(ij,1,k))*wrk2d(ij,1)
           ENDDO
           
           IF ( g(3)%size .GT. 1 ) THEN
              DO ij=1, imax*jmax
                 wrk2d(ij,1) = wrk2d(ij,1) + (mtgfm*vis(ij,1,k)/rho(ij,1,k))/(g(3)%jac(k+kdsp,1)**3)
              ENDDO
           ENDIF

           bmax = MAXVAL(wrk2d(1:nxy,1))
           IF (bmax .ge. rmax) rmax = bmax
        ENDDO

     ELSE
        DO k=1, kmax
           DO ij=1, imax*jmax
              index = ij-1
              wrk2d(ij,1) = C_1_R/(g(1)%jac(MOD(index,imax)+1+idsp,1)**3)
              wrk2d(ij,1) = wrk2d(ij,1) + C_1_R/(g(2)%jac((ij-1)/imax+1,1)**3)
              wrk2d(ij,1) = (mtgfm/rho(ij,1,k))*wrk2d(ij,1)
           ENDDO
           
           IF ( g(3)%size .GT. 1 ) THEN
              DO ij=1, imax*jmax
                 wrk2d(ij,1) = wrk2d(ij,1) + (mtgfm/rho(ij,1,k))/(g(3)%jac(k+kdsp,1)**3)
              ENDDO
           ENDIF

           bmax = MAXVAL(wrk2d(1:nxy,1))
           IF (bmax .ge. rmax) rmax = bmax
        ENDDO

     ENDIF

#ifdef USE_MPI
     CALL MPI_ALLREDUCE(rmax, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
     rmax = dummy
#endif

     IF ( rmax .GT. C_0_R ) dtr = SQRT(cflr/rmax)

! -------------------------------------------------------------------
! Finite rate
! -------------------------------------------------------------------
  ELSE IF ( ireactive .EQ. CHEM_FINITE ) THEN
! I obtain TGFM/(\rho*c^2) (not TGFM/(\rho*u^2)) as eigenvalue if I 
! assume that omega depends only on \rho as \rho*S, S constant.     
     DO k=1, kmax
        DO ij=1, imax*jmax
           wrk2d(ij,1) = rho(ij,1,k)*(u(ij,1,k)*u(ij,1,k)&
                +v(ij,1,k)*v(ij,1,k))
        ENDDO

        IF ( g(3)%size .GT. 1 ) THEN
           DO ij=1, imax*jmax
              wrk2d(ij,1) = wrk2d(ij,1) + rho(ij,1,k)*w(ij,1,k)*w(ij,1,k)
           ENDDO
        ENDIF

        bmax = MAXVAL(wrk2d(1:nxy,1))
        IF (bmax .ge. rmax) rmax = bmax
     ENDDO
     rmax = TGFM/rmax

#ifdef USE_MPI
     CALL MPI_ALLREDUCE(rmax, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
     rmax = dummy
#endif

     IF ( rmax .GT. C_0_R ) dtr = cflr/rmax

  ENDIF

#endif

! ###################################################################
! Final operations
! ###################################################################
  IF ( cfl .GT. C_0_R ) THEN
     IF ( rkm_mode .EQ. RKM_EXP3 .OR. rkm_mode .EQ. RKM_EXP4 ) THEN 
        dt_loc = MIN(dtc, dtd)
     ELSE 
        dt_loc = dtc 
     ENDIF
     dt_loc = MIN(dt_loc, dtr)
#ifdef USE_MPI
     CALL MPI_ALLREDUCE(dt_loc, dummy, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ims_err)
     dt_loc = dummy
#endif
     
     dtime = dt_loc

  ENDIF

! Real CFL and diffusion numbers being used
  cmax = dtime * cmax
  dmax = dtime * dmax
#ifdef CHEMISTRY
  IF      ( ireactive .EQ. CHEM_INFINITE ) THEN; rmax = (dtime**2) *rmax
  ELSE IF ( ireactive .EQ. CHEM_FINITE   ) THEN; rmax =  dtime     *rmax; ENDIF
#endif

  logs_data(2) = cmax
  logs_data(3) = dmax
  logs_data(4) = rmax

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING TIME_COURANT' )
#endif

  RETURN

END SUBROUTINE TIME_COURANT
