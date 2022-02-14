#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INIFLOW"

PROGRAM INIFLOW

  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE TLAB_ARRAYS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_PROCS
#endif
  USE THERMO_VARS, ONLY : imixture
  USE FLOW_LOCAL
  USE IO_FIELDS
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
  TREAL, DIMENSION(:), ALLOCATABLE, SAVE :: ci,cj,ck, ipos,jpos,kpos
#endif

  TREAL, DIMENSION(:), POINTER :: e, rho, p, T

  !########################################################################
  CALL TLAB_START()

  CALL DNS_READ_GLOBAL(ifile)
  CALL FLOW_READ_LOCAL(ifile)
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(ifile)
#endif

#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

  inb_wrk2d = MAX(inb_wrk2d,3)
  isize_wrk3d = isize_txc_field

  IF ( flag_u == 0 ) THEN; inb_txc = 2
  ELSE;                    inb_txc = 8
  END IF

  CALL TLAB_ALLOCATE(C_FILE_LOC)

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  IF ( flag_u /= 0 ) THEN ! Initialize Poisson Solver
     IF ( ifourier == 1 .AND. g(1)%periodic .AND. g(3)%periodic ) THEN
        CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)

     ELSE
#ifdef USE_CGLOC
        ALLOCATE(ci(isize_field*2))
        ALLOCATE(cj(isize_field*2))
        ALLOCATE(ck(isize_field*2))
        ALLOCATE(ipos(imax*4))
        ALLOCATE(jpos(jmax*4))
        ALLOCATE(kpos(kmax*4))

        IF ( .NOT. g(1)%uniform .OR. .NOT. g(2)%uniform ) THEN
           CALL TLAB_WRITE_ASCII(lfile, 'Initializing conjugate gradient, non-uniform grid, second-order.')
           cg_unif = 1; cg_ord = 2
           ! to be rewritten in terms of grid derived type
           ! CALL CGBC2(cg_unif, imode_fdm, imax,jmax,kmax,g(3)%size, &
           !      i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
        ELSE
           CALL TLAB_WRITE_ASCII(lfile, 'Initializing conjugate gradient, uniform grid, fourth-order.')
           cg_unif = 0; cg_ord = 4
           ! CALL CGBC4(cg_unif, imax,jmax,kmax,g(3)%size, &
           !      i1bc,j1bc,k1bc, scalex,scaley,scalez, dx,dy,dz, ipos,jpos,kpos,ci,cj,ck, wrk2d)
        END IF
#else
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. CG routines needed.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
#endif
     END IF

  END IF

  ! ###################################################################
  itime = 0; rtime = C_0_R
  q = C_0_R

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Initializing velocity.')

  CALL VELOCITY_MEAN( q(1,1),q(1,2),q(1,3), wrk1d )

  SELECT CASE( flag_u )
  CASE( 1 )
     CALL VELOCITY_DISCRETE(txc(1,1),txc(1,2),txc(1,3), wrk1d, wrk2d)
     q(1:isize_field,1:3) =  q(1:isize_field,1:3) + txc(1:isize_field,1:3)

  CASE( 2,3,4 )
     CALL VELOCITY_BROADBAND(txc(1,1),txc(1,2),txc(1,3), txc(1,4),txc(1,5),txc(1,6),txc(1,7),txc(1,8), wrk1d,wrk2d,wrk3d)
     q(1:isize_field,1:3) =  q(1:isize_field,1:3) + txc(1:isize_field,1:3)

  END SELECT

  ! ###################################################################
  IF ( imode_eqns == DNS_EQNS_TOTAL .OR. imode_eqns == DNS_EQNS_INTERNAL ) THEN
     CALL TLAB_WRITE_ASCII(lfile,'Initializing pressure and density.')
     e   => q(:,4)
     rho => q(:,5)
     p   => q(:,6)
     T   => q(:,7)

     CALL PRESSURE_MEAN(p,T,s, wrk1d)
     CALL DENSITY_MEAN(rho,p,T,s, txc, wrk1d,wrk2d,wrk3d)

     IF ( flag_u /= 0 ) THEN
        CALL PRESSURE_FLUCTUATION(q(1,1),q(1,2),q(1,3),rho,p,txc(1,1), &
             txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)
     END IF

     IF ( imixture > 0 ) THEN
        CALL IO_READ_FIELDS(TRIM(ADJUSTL(tag_scal))//'ics', IO_SCAL, imax,jmax,kmax, inb_scal,i0, s, wrk3d)
     END IF

     IF ( flag_t == 4 .OR. flag_t == 5 ) THEN
        CALL DENSITY_FLUCTUATION(flag_t, s,p,rho, txc(1,1),txc(1,2), wrk2d,wrk3d)
     END IF

     ! Calculate specfic energy. Array s should contain the species fields at this point.
     CALL THERMO_THERMAL_TEMPERATURE(imax,jmax,kmax, s, p, rho, txc(1,1))
     CALL THERMO_CALORIC_ENERGY(imax,jmax,kmax, s, txc(1,1), e)

  END IF

  ! ###################################################################
  CALL IO_WRITE_FIELDS(TRIM(ADJUSTL(tag_flow))//'ics', IO_FLOW, imax,jmax,kmax, inb_flow, q, wrk3d)

  CALL TLAB_STOP(0)
END PROGRAM INIFLOW
