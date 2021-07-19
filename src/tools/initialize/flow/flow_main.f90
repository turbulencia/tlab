#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INIFLOW"

PROGRAM INIFLOW

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE TLAB_ARRAYS
  USE THERMO_GLOBAL, ONLY : imixture
  USE FLOW_LOCAL
#ifdef USE_CGLOC
  USE CG_GLOBAL, ONLY : cg_unif, cg_ord
#endif
#ifdef CHEMISTRY
  USE CHEM_GLOBAL
#endif

  IMPLICIT NONE

  ! -------------------------------------------------------------------
  ! Additional local arrays
#ifdef USE_CGLOC
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE         :: ci,cj,ck, ipos,jpos,kpos
#endif

  TREAL, DIMENSION(:),   POINTER :: e, rho, p, T

  !########################################################################
  CALL DNS_START()

  CALL DNS_READ_GLOBAL(ifile)
  CALL FLOW_READ_LOCAL(ifile)
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(ifile)
#endif

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  inb_wrk2d=MAX(inb_wrk2d,3)
  isize_wrk3d = isize_txc_field

  IF ( flag_u .EQ. 0 ) THEN; inb_txc = 2
  ELSE;                      inb_txc = 8
  ENDIF

  CALL TLAB_ALLOCATE(C_FILE_LOC)

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

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

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Initializing flow fiels.')

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
    e   => q(:,4)
    rho => q(:,5)
    p   => q(:,6)
    T   => q(:,7)
  ENDIF

  IF ( flag_u .NE. 0 ) THEN ! Initialize Poisson Solver
    IF ( ifourier .EQ. 1 .AND. g(1)%periodic .AND. g(3)%periodic ) THEN
      CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)

    ELSE
#ifdef USE_CGLOC
      IF ( .NOT. g(1)%uniform .NOT. .OR. g(2)%uniform ) THEN
        CALL IO_WRITE_ASCII(lfile, 'Initializing conjugate gradient, non-uniform grid, second-order.')
        cg_unif = 1; cg_ord = 2
        ! to be rewritten in terms of grid derived type
        ! CALL CGBC2(cg_unif, imode_fdm, imax,jmax,kmax,g(3)%size, &
        !      i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
      ELSE
        CALL IO_WRITE_ASCII(lfile, 'Initializing conjugate gradient, uniform grid, fourth-order.')
        cg_unif = 0; cg_ord = 4
        ! CALL CGBC4(cg_unif, imax,jmax,kmax,g(3)%size, &
        !      i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
      ENDIF
#else
      CALL IO_WRITE_ASCII(efile, 'INIFLOW: CG routines needed.')
      CALL DNS_STOP(DNS_ERROR_OPTION)
#endif
    ENDIF

  ENDIF

  itime = 0; rtime = C_0_R
  q = C_0_R

  ! ###################################################################
  ! Pressure and density mean fields
  ! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 1')
#endif

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
    CALL PRESSURE_MEAN(p,T,s, wrk1d,wrk2d,wrk3d)

#ifdef CHEMISTRY
    IF ( ireactive .EQ. CHEM_NONE ) THEN
#endif
      CALL DENSITY_MEAN(rho,p,T,s, txc, wrk1d,wrk2d,wrk3d)

#ifdef CHEMISTRY
    ELSE
      IF ( icalc_scal .EQ. 1 ) THEN
        CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal,inb_scal, isize_wrk3d, s(1,inb_scal), wrk3d)
        IF ( ireactive .EQ. CHEM_FINITE .AND. ichem_config .NE. CHEM_PREMIXED .AND. flag_mixture .EQ. 2 ) THEN ! Initialize density from flame
          r0   = C_0_R
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
! Velocity
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 2')
#endif

  CALL VELOCITY_MEAN( q(1,1),q(1,2),q(1,3), wrk1d,wrk3d )

  SELECT CASE( flag_u )
  CASE( 1 )
    CALL VELOCITY_DISCRETE(txc(1,1),txc(1,2),txc(1,3), wrk1d, wrk2d)
    q(1:isize_field,1:3) =  q(1:isize_field,1:3) + txc(1:isize_field,1:3)

  CASE( 2,3,4 )
    CALL VELOCITY_BROADBAND(txc(1,1),txc(1,2),txc(1,3), txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk1d,wrk2d,wrk3d)
    q(1:isize_field,1:3) =  q(1:isize_field,1:3) + txc(1:isize_field,1:3)

  END SELECT

! ###################################################################
! Pressure and density fluctuation fields
! ###################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'INIFLOW: Section 3')
#endif

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
    IF ( flag_u .NE. 0 ) THEN
      CALL PRESSURE_FLUCTUATION(q(1,1),q(1,2),q(1,3),rho,p,txc(1,1), &
      txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
    ENDIF

    IF ( imixture .GT. 0 ) THEN
      CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal,i0, isize_wrk3d, s, wrk3d)
    ENDIF

    IF ( flag_t .EQ. 4 .OR. flag_t .EQ. 5 ) THEN
      CALL DENSITY_FLUCTUATION(flag_t, s,p,rho, txc(1,1),txc(1,2), wrk2d,wrk3d)
    ENDIF

    ! Calculate specfic energy. Array s should contain the species fields at this point.
    CALL THERMO_THERMAL_TEMPERATURE(imax,jmax,kmax, s, p, rho, txc(1,1))
    CALL THERMO_CALORIC_ENERGY(imax,jmax,kmax, s, txc(1,1), e)

  ENDIF

  ! ###################################################################
  CALL DNS_WRITE_FIELDS('flow.ics', i2, imax,jmax,kmax, inb_flow, isize_wrk3d, q, wrk3d)

  CALL DNS_STOP(0)
END PROGRAM INIFLOW
