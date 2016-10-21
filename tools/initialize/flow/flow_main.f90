#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INIFLOW"

!########################################################################
!# Tool/Library INIT/FLOW
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!# 2007/03/17 - J.P. Mellado
!#              Adding shear with gravity
!# 2007/04/30 - J.P. Mellado
!#              Splitting calculation into homentropic fluctuations
!#              and mean profiles, to allow easier introduction of
!#              buoyancy cases
!# 2007/10/30 - J.P. Mellado
!#              Jet added.
!#
!########################################################################
!# DESCRIPTION
!#
!# Create an initial condition that approximately satisfies the steady
!# inviscid Navier-Stokes equations:
!# 1) Mean profiles
!# 2) Add perturbation
!# The velocity fluctuation in first calculated to reduce aux arrays.
!#
!########################################################################
PROGRAM INIFLOW

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE FLOW_LOCAL
#ifdef USE_CGLOC
  USE CG_GLOBAL, ONLY : cg_unif, cg_ord
#endif
#ifdef CHEMISTRY
  USE CHEM_GLOBAL
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_FFTW
#include "fftw3.f"
#endif

! -------------------------------------------------------------------
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: q, s, txc
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk1d,wrk2d,wrk3d
#ifdef USE_CGLOC
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: ci,cj,ck, ipos,jpos,kpos
#endif

  TARGET q, wrk3d
  TREAL, DIMENSION(:),   POINTER :: u, v, w, p, rho

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
  CHARACTER*64 str, line
  CHARACTER*32 inifile

  TREAL r0
  TINTEGER iread_flow, iread_scal, isize_wrk3d, ierr

!########################################################################
!########################################################################
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  CALL FLOW_READ_LOCAL(inifile)
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(inifile)
#endif

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Definitions
! -------------------------------------------------------------------
  itime = 0; rtime = C_0_R
  
  inb_wrk2d   = MAX(inb_wrk2d,6)

  r0   = C_0_R

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------      
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

  iread_flow = 1
  iread_scal = 0

  IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL   )THEN
     iread_scal = 1
  ENDIF

  IF ( flag_u .EQ. 0 ) THEN; inb_txc = 2
  ELSE;                      inb_txc = 8; ENDIF

  isize_wrk3d = isize_txc_field
!  isize_wrk3d = MAX(isize_wrk1d*300,isize_wrk3d)
#include "dns_alloc_arrays.h"

#ifdef USE_CGLOC
  IF ( flag_u .NE. 0 )THEN
     ALLOCATE(ci(isize_field*2))
     ALLOCATE(cj(isize_field*2))
     ALLOCATE(ck(isize_field*2))
     ALLOCATE(ipos(imax*4))
     ALLOCATE(jpos(jmax*4))
     ALLOCATE(kpos(kmax*4))
  ENDIF
#endif

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! -------------------------------------------------------------------
! Initialize Poisson Solver
! -------------------------------------------------------------------
  IF ( flag_u .NE. 0 ) THEN
     IF ( ifourier .EQ. 1 .AND. g(1)%periodic .AND. g(3)%periodic ) THEN ! Doubly periodic in xOz 
        CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)

     ELSE
#ifdef USE_CGLOC
        IF ( .NOT. g(1)%uniform .NOT. .OR. g(2)%uniform ) THEN
           CALL IO_WRITE_ASCII(lfile, 'Initializing conjugate gradient, non-uniform grid, second-order.')
           cg_unif = 1; cg_ord = 2
           CALL CGBC2(cg_unif, imode_fdm, imax,jmax,kmax,kmax_total, &
                i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
        ELSE
           CALL IO_WRITE_ASCII(lfile, 'Initializing conjugate gradient, uniform grid, fourth-order.')
           cg_unif = 0; cg_ord = 4
           CALL CGBC4(cg_unif, imax,jmax,kmax,kmax_total, &
                i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
        ENDIF
#else
        CALL IO_WRITE_ASCII(efile, 'INIFLOW: CG routines needed.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
#endif
     ENDIF

  ENDIF
  
! ###################################################################
  q(:,:)  = C_0_R

  u => q(:,1)
  v => q(:,2)
  w => q(:,3)
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     p   => q(:,4)
     rho => q(:,5)
  ELSE
     rho => wrk3d ! array not used
  ENDIF

! ###################################################################
! Pressure and density mean fields
! Array txc1 used in this section for the temperature
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 1')
#endif

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     CALL PRESSURE_MEAN(p,txc(1,1),s, wrk1d,wrk2d,wrk3d)
     
#ifdef CHEMISTRY
     IF ( ireactive .EQ. CHEM_NONE ) THEN
#endif
        CALL DENSITY_MEAN(rho,p,txc(1,1),s, txc(1,2), wrk1d,wrk2d,wrk3d)
        
#ifdef CHEMISTRY
     ELSE
        IF ( icalc_scal .EQ. 1 ) THEN
           CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, &
                inb_scal,inb_scal, isize_wrk3d, s(1,inb_scal), wrk3d)
           IF ( ireactive .EQ. CHEM_FINITE .AND. ichem_config .NE. CHEM_PREMIXED .AND.&
                flag_mixture .EQ. 2 ) THEN
! Initialize density from flame
              CALL CHEM_READ_TEXT(flame_ini_file, i11, i11, i2, isize_field, r0, s(1,inb_scal), rho, isize_wrk3d, wrk3d)
           ELSE
              CALL THERMO_BURKESCHUMANN(rho, s(1,inb_scal))
           ENDIF
        ENDIF
        CALL IO_WRITE_ASCII(efile, 'INIFLOW: Chemistry part to be checked')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
#endif

  ENDIF

! ###################################################################
! Velocity mean fields
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 2')
#endif
  CALL VELOCITY_MEAN(rho,u,v,w, wrk1d,wrk3d)

! ###################################################################
! Velocity perturbation fields
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 3')
#endif
  IF ( flag_u .NE. 0 ) THEN
     txc(:,1:3) = C_0_R

     IF      ( flag_u .EQ. 1 ) THEN
        CALL VELOCITY_DISCRETE(i1, txc(1,1),txc(1,2),txc(1,3))

     ELSE IF ( flag_u .GT. 1 ) THEN
        CALL VELOCITY_BROADBAND(flag_u, txc(1,1),txc(1,2),txc(1,3), &
             txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk1d,wrk2d,wrk3d)
        
     ENDIF

     IF ( norm_ini_u .GE. C_0_R ) THEN
        CALL NORMALIZE(imax,jmax,kmax, txc(1,1),txc(1,2),txc(1,3), norm_ini_u)
     ENDIF
     
     q(1:isize_field,1:3) =  q(1:isize_field,1:3) + txc(1:isize_field,1:3)

  ENDIF

! ###################################################################
! Pressure and density fluctuation fields
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 4')
#endif

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     IF ( flag_u .NE. 0 ) THEN
        CALL PRESSURE_FLUCTUATION(u,v,w,rho,p,txc(1,1), &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
     ENDIF
     
     IF ( imixture .GT. 0 ) THEN
        CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)
     ENDIF
     
     IF ( flag_t .EQ. 4 .OR. flag_t .EQ. 5 ) THEN
        CALL DENSITY_FLUCTUATION(flag_t, s,p,rho, txc(1,1),txc(1,2), wrk2d,wrk3d)
     ENDIF

  ENDIF

! ###################################################################
! Output file
! Array s should contain the species fields at this point.
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 5')
#endif
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     CALL THERMO_THERMAL_TEMPERATURE(imax,jmax,kmax, s, p, rho, txc(1,1))
     CALL THERMO_CALORIC_ENERGY(imax,jmax,kmax, s, txc(1,1), p)
  ENDIF
  CALL DNS_WRITE_FIELDS('flow.ics', i2, imax,jmax,kmax, inb_flow, isize_wrk3d, q, wrk3d)
  
  CALL DNS_END(0)

  STOP
END PROGRAM INIFLOW
