#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef LES
#include "les_const.h"
#endif

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2008/02/09 - J.P. Mellado
!#              Cleaned
!#
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
!# ARGUMENTS 
!#
!########################################################################

SUBROUTINE TIME_COURANT(q,s, wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : gama0
  USE DNS_LOCAL
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

  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(IN) :: q, s
  TREAL, DIMENSION(imax,jmax,kmax)               :: wrk3d
  TREAL, DIMENSION(imax,jmax)                    :: wrk2d
  
  TARGET :: q

! -------------------------------------------------------------------
  TINTEGER k, ij, index, nxy, kdsp, idsp
  TREAL dt_loc
  TREAL cmax, bmax, dmax, dtc, cfld, dtd, rmax, dtr, cflr
  TREAL vscfct, vsctmp, vscles
#ifdef CHEMISTRY
  TREAL mtgfm, zmin, zmax, umin, umax
  TREAL vmin, vmax, wmin, wmax
#endif
#ifdef USE_MPI
  TREAL dummy
#endif

! Pointers to existing allocated space
  TREAL, DIMENSION(:,:,:), POINTER :: u, v, w, rho, p, vis
  TREAL, DIMENSION(:),     POINTER :: dx, dy, dz

! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING TIME_COURANT' )
#endif

! Define pointers
  dx => g(1)%jac(:,1)
  dy => g(2)%jac(:,1)
  dz => g(3)%jac(:,1)

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

  nxy = imax*jmax

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
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE   .OR. &
       imode_eqns .EQ. DNS_EQNS_ANELASTIC             ) THEN
     DO k = 1,kmax
        DO ij = 1,imax*jmax
           index = ij-1
           wrk2d(ij,1) = ABS(u(ij,1,k))/dx(MOD(index,imax)+1+idsp)
           wrk2d(ij,1) = wrk2d(ij,1) + ABS(v(ij,1,k))/dy((ij-1)/imax+1)
        ENDDO
        IF ( kmax_total .GT. 1 ) THEN
           DO ij = 1,imax*jmax
              wrk2d(ij,1) = wrk2d(ij,1) + ABS(w(ij,1,k))/dz(k+kdsp)
           ENDDO
        ENDIF
        
        bmax = MAXVAL(wrk2d(1:nxy,1))
        cmax = MAX(cmax,bmax)
        
     ENDDO
  
! -------------------------------------------------------------------
! Compressible
! Calculate global maximum of (u+c)/dx + (v+c)/dy + (w+c)/dz 
! -------------------------------------------------------------------
  ELSE
! sound speed; positiveness of pressure and density is check in routine DNS_CONTROL
     wrk3d = SQRT(gama0*p/rho)

     DO k = 1,kmax
        DO ij = 1,imax*jmax
           index = ij-1
           wrk2d(ij,1) = (ABS(u(ij,1,k))+wrk3d(ij,1,k))/dx(MOD(index,imax)+1+idsp)
           wrk2d(ij,1) = wrk2d(ij,1) + &
                (ABS(v(ij,1,k))+wrk3d(ij,1,k))/dy((ij-1)/imax+1)
        ENDDO
        IF ( kmax_total .GT. 1 ) THEN
           DO ij = 1,imax*jmax
              wrk2d(ij,1) = wrk2d(ij,1) + (ABS(w(ij,1,k))+wrk3d(ij,1,k))/dz(k+kdsp)
           ENDDO
        ENDIF
        
        bmax = MAXVAL(wrk2d(1:nxy,1))
        cmax = MAX(cmax,bmax)
        
     ENDDO
     
  ENDIF

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
  vscfct = C_1_R
  vsctmp = C_1_R/prandtl
  vscfct = MAX(vscfct, vsctmp)
  vsctmp = C_1_R/MINVAL(schmidt(1:inb_scal))
  vscfct = MAX(vscfct, vsctmp)

#ifdef LES
  IF ( iles .EQ. 1 ) THEN
     IF ( iles_type_regu .EQ. LES_REGU_NONE ) THEN
        vscles = C_0_R
     ELSE IF ( iles_type_regu.EQ.LES_REGU_SMGDYN .OR. iles_type_regu.EQ.LES_REGU_SMGDYNRMS ) THEN
        vscles = MAX(smg_udiff,smg_ediff)
        IF ( icalc_scal .EQ. 1 ) vscles = MAX(vscles,smg_zdiff)
     ELSE IF ( iles_type_regu .EQ. LES_REGU_SMGSTA ) THEN
        vscles = C_1_R
        vsctmp = C_1_R/smg_prandtl
        vscles = MAX(vscles, vsctmp)
        vsctmp = C_1_R/smg_schmidt
        vscles = MAX(vscles, vsctmp)
        vscles = vscles * smg_udiff
     ENDIF
     IF ( iviscous .EQ. EQNS_NONE ) vscfct = C_0_R
  ELSE
#endif
     vscles = C_0_R
#ifdef LES
  ENDIF
#endif

! -------------------------------------------------------------------
! Incompressible
! Calculate global maximum of \nu*(1/dx^2 + 1/dy^2 + 1/dz^2)
! -------------------------------------------------------------------
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE   .OR. &
       imode_eqns .EQ. DNS_EQNS_ANELASTIC      ) THEN
     DO k = 1,kmax
        DO ij = 1,imax*jmax
           index = ij-1
           wrk2d(ij,1) = C_1_R/(dx(MOD(index,imax)+1+idsp)**2) + C_1_R/(dy((ij-1)/imax+1)**2)
           wrk2d(ij,1) = (vscfct*visc+vscles)*wrk2d(ij,1)
        ENDDO
        IF ( kmax_total .GT. 1 ) THEN
           DO ij=1, imax*jmax
              wrk2d(ij,1) = wrk2d(ij,1) + (vscfct*visc+vscles)/(dz(k+kdsp)**2)
           ENDDO
        ENDIF

        bmax = MAXVAL(wrk2d(1:nxy,1))
        dmax = MAX(dmax,bmax)

     ENDDO

! -------------------------------------------------------------------
! Compressible
! Calculate global maximum of \mu/rho*(1/dx^2 + 1/dy^2 + 1/dz^2)
! -------------------------------------------------------------------
  ELSE
     IF ( itransport .EQ. EQNS_TRANS_POWERLAW ) THEN
        DO k = 1,kmax
           DO ij = 1,imax*jmax
              index = ij-1
              wrk2d(ij,1) = C_1_R/(dx(MOD(index,imax)+1+idsp)**2) + C_1_R/(dy((ij-1)/imax+1)**2)
              wrk2d(ij,1) = (vscfct*visc*vis(ij,1,k)/rho(ij,1,k)+vscles)*wrk2d(ij,1)
           ENDDO
           IF ( kmax_total .GT. 1 ) THEN
              DO ij=1, imax*jmax
                 wrk2d(ij,1) = wrk2d(ij,1) + &
                   (vscfct*visc*vis(ij,1,k)/rho(ij,1,k)+vscles)/(dz(k+kdsp)**2)
              ENDDO
           ENDIF

           bmax = MAXVAL(wrk2d(1:nxy,1))
           dmax = MAX(dmax,bmax)

        ENDDO
        
     ELSE ! constant dynamic viscosity
        DO k = 1,kmax
           DO ij = 1,imax*jmax
              index = ij-1
              wrk2d(ij,1) = C_1_R/(dx(MOD(index,imax)+1+idsp)**2) + C_1_R/(dy((ij-1)/imax+1)**2)
              wrk2d(ij,1) = (vscfct*visc/rho(ij,1,k)+vscles)*wrk2d(ij,1)
           ENDDO
           IF ( kmax_total .GT. 1 ) THEN
              DO ij=1, imax*jmax
                 wrk2d(ij,1) = wrk2d(ij,1) + &
                      (vscfct*visc/rho(ij,1,k)+vscles)/(dz(k+kdsp)**2)
              ENDDO
           ENDIF

           bmax = MAXVAL(wrk2d(1:nxy,1))
           dmax = MAX(dmax,bmax)

        ENDDO

     ENDIF
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
              wrk2d(ij,1) = C_1_R/(dx(MOD(index,imax)+1+idsp)**3)
              wrk2d(ij,1) = wrk2d(ij,1) + C_1_R/(dy((ij-1)/imax+1)**3)
              wrk2d(ij,1) = (mtgfm*vis(ij,1,k)/rho(ij,1,k))*wrk2d(ij,1)
           ENDDO
           
           IF ( kmax_total .GT. 1 ) THEN
              DO ij=1, imax*jmax
                 wrk2d(ij,1) = wrk2d(ij,1) + (mtgfm*vis(ij,1,k)/rho(ij,1,k))/(dz(k+kdsp)**3)
              ENDDO
           ENDIF

           bmax = MAXVAL(wrk2d(1:nxy,1))
           IF (bmax .ge. rmax) rmax = bmax
        ENDDO

     ELSE
        DO k=1, kmax
           DO ij=1, imax*jmax
              index = ij-1
              wrk2d(ij,1) = C_1_R/(dx(MOD(index,imax)+1+idsp)**3)
              wrk2d(ij,1) = wrk2d(ij,1) + C_1_R/(dy((ij-1)/imax+1)**3)
              wrk2d(ij,1) = (mtgfm/rho(ij,1,k))*wrk2d(ij,1)
           ENDDO
           
           IF ( kmax_total .GT. 1 ) THEN
              DO ij=1, imax*jmax
                 wrk2d(ij,1) = wrk2d(ij,1) + (mtgfm/rho(ij,1,k))/(dz(k+kdsp)**3)
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

        IF ( kmax_total .GT. 1 ) THEN
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
