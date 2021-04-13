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
!# If the extension by Lamballais et al is used, then the maximum
!# diffusion number is 2.9/pi^2, about 0.29.
!# For the (5)4RK from CarpenterKennedy1994 it is 4.639/pi^2 = 0.47.
!#
!# If twice the first order finite difference operator is used, then the
!# maximum diffusion number is 2.9/1.989^2, about 0.73.
!# For the (5)4RK from CarpenterKennedy1994 it is 4.639/1.989^2 = 1.17
!#
!# In incompressible mode the arrays rho, p and vis are not used
!#
!########################################################################
SUBROUTINE TIME_COURANT(q,s, wrk3d)

#ifdef TRACE_ON
  USE DNS_CONSTANTS, ONLY : tfile
#endif
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_scal, imode_eqns
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : inb_flow_array, inb_scal_array
  USE DNS_GLOBAL,    ONLY : itransport, visc, prandtl, schmidt
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL,     ONLY : cfla, cfld, cflr, dtime, rkm_mode, logs_data
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

  TREAL, DIMENSION(imax,jmax,kmax,inb_flow_array), INTENT(IN) :: q
  TREAL, DIMENSION(imax,jmax,kmax,inb_scal_array), INTENT(IN) :: s
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: wrk3d

  TARGET :: q

! -------------------------------------------------------------------
  TINTEGER i,j,k, kdsp, idsp, ipmax
  TREAL dt_loc
  TREAL pmax(3), dtc,dtd,dtr
  TREAL schmidtfactor, viscles, dummy
#ifdef CHEMISTRY
  TREAL mtgfm, zmin, zmax, umin, umax
  TREAL vmin, vmax, wmin, wmax
#endif
#ifdef USE_MPI
  TREAL pmax_aux(3)
#endif

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: rho, p, vis

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING TIME_COURANT' )
#endif

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

! So that the minimum non-zero determines dt at the end
  dtc = C_BIG_R
  dtd = C_BIG_R
  dtr = C_BIG_R

! Initialize counter of time restrictions
  ipmax = 0

! ###################################################################
! CFL number condition
! ###################################################################
  ipmax = ipmax +1
  pmax(1) = C_0_R

! -------------------------------------------------------------------
! Incompressible: Calculate global maximum of u/dx + v/dy + w/dz
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     IF ( g(3)%size .GT. 1 ) THEN
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = ABS(q(i,j,k,1)) *g(1)%jac(i+idsp,3) &
                        + ABS(q(i,j,k,2)) *g(2)%jac(j,3)      &
                        + ABS(q(i,j,k,3)) *g(3)%jac(k+kdsp,3)
        ENDDO; ENDDO; ENDDO
     ELSE
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = ABS(q(i,j,k,1)) *g(1)%jac(i+idsp,3) &
                        + ABS(q(i,j,k,2)) *g(2)%jac(j,3)
        ENDDO; ENDDO; ENDDO
     ENDIF

! -------------------------------------------------------------------
! Compressible: Calculate global maximum of (u+c)/dx + (v+c)/dy + (w+c)/dz
! -------------------------------------------------------------------
  ELSE
     wrk3d = SQRT(gama0*p/rho) ! sound speed; positiveness of p and rho checked in routine DNS_CONTROL
     IF ( g(3)%size .GT. 1 ) THEN
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = (ABS(q(i,j,k,1))+wrk3d(i,j,k)) *g(1)%jac(i+idsp,3) &
                        + (ABS(q(i,j,k,2))+wrk3d(i,j,k)) *g(2)%jac(j,3)      &
                        + (ABS(q(i,j,k,3))+wrk3d(i,j,k)) *g(3)%jac(k+kdsp,3)
        ENDDO; ENDDO; ENDDO
     ELSE
        DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
           wrk3d(i,j,k) = (ABS(q(i,j,k,1))+wrk3d(i,j,k)) *g(1)%jac(i+idsp,3) &
                        + (ABS(q(i,j,k,2))+wrk3d(i,j,k)) *g(2)%jac(j,3)
        ENDDO; ENDDO; ENDDO
     ENDIF

  ENDIF

  pmax(1) = MAXVAL(wrk3d)

! ###################################################################
! Diffusion number condition
! ###################################################################
  ipmax = ipmax +1
  pmax(2) = C_0_R

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
! Incompressible: Calculate global maximum of \nu*(1/dx^2 + 1/dy^2 + 1/dz^2)
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
     pmax(2) = (schmidtfactor*visc+viscles) *MAXVAL(wrk3d)

! -------------------------------------------------------------------
! Compressible: Calculate global maximum of \mu/rho*(1/dx^2 + 1/dy^2 + 1/dz^2)
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

     pmax(2) = MAXVAL(wrk3d)

  ENDIF

#ifdef CHEMISTRY
! ###################################################################
! Reacting time step control
! ###################################################################
  ipmax = ipmax +1
  pmax(ipmax) = C_0_R

! -------------------------------------------------------------------
! Infinitely fast
! -------------------------------------------------------------------
  IF ( ireactive .EQ. CHEM_INFINITE ) THEN
     CALL MINMAX(imax,jmax,kmax, q(1,1,1,1), umin,umax)
     CALL MINMAX(imax,jmax,kmax, q(1,1,1,2), vmin,vmax)
     CALL MINMAX(imax,jmax,kmax, q(1,1,1,3), wmin,wmax)
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
           IF (bmax .ge. pmax(ipmax)) pmax(ipmax) = bmax
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
           IF (bmax .ge. pmax(ipmax)) pmax(ipmax) = bmax
        ENDDO

     ENDIF

! -------------------------------------------------------------------
! Finite rate
! -------------------------------------------------------------------
  ELSE IF ( ireactive .EQ. CHEM_FINITE ) THEN
! I obtain TGFM/(\rho*c^2) (not TGFM/(\rho*u^2)) as eigenvalue if I
! assume that omega depends only on \rho as \rho*S, S constant.
     IF ( g(3)%size .GT. 1 ) THEN
        wrk3d = rho *( q(:,:,:,1) *q(:,:,:,1) +q(:,:,:,2) *q(:,:,:,2) +q(:,:,:,3) *q(:,:,:,3) )
     ELSE
        wrk3d = rho *( q(:,:,:,1) *q(:,:,:,1) +q(:,:,:,2) *q(:,:,:,2) )
     ENDIF

     bmax = MAXVAL(wrk3d)
     IF (bmax .ge. pmax(ipmax)) pmax(ipmax) = bmax
     pmax(ipmax) = TGFM/pmax(ipmax)

  ENDIF

#endif

! ###################################################################
! Final operations
! ###################################################################
#ifdef USE_MPI
  CALL MPI_ALLREDUCE(pmax, pmax_aux, ipmax, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
  pmax(1:ipmax) = pmax_aux(1:ipmax)
#endif

  IF ( pmax(1) .GT. C_0_R ) dtc = cfla /pmax(1) ! Set time step for the given CFL number
  IF ( pmax(2) .GT. C_0_R ) dtd = cfld/pmax(2) ! Set time step for the given diffusion number
#ifdef CHEMISTRY
  IF ( pmax(ipmax) .GT. C_0_R ) dtr = cflr/pmax(ipmax)
#endif

! -------------------------------------------------------------------
  IF ( cfla .GT. C_0_R ) THEN
     IF ( rkm_mode .EQ. RKM_EXP3 .OR. rkm_mode .EQ. RKM_EXP4 ) THEN
        dt_loc = MIN(dtc, dtd)
     ELSE
        dt_loc = dtc
     ENDIF
     dt_loc = MIN(dt_loc, dtr)

     dtime = dt_loc

  ENDIF

! Real CFL and diffusion numbers being used, for the logfile
  logs_data(2) = dtime * pmax(1)
  logs_data(3) = dtime * pmax(2)
#ifdef CHEMISTRY
  IF      ( ireactive .EQ. CHEM_INFINITE ) THEN; logs_data(4) = (dtime**2) *pmax(3)
  ELSE IF ( ireactive .EQ. CHEM_FINITE   ) THEN; logs_data(4) =  dtime     *pmax(3); ENDIF
#endif

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING TIME_COURANT' )
#endif

  RETURN

END SUBROUTINE TIME_COURANT
