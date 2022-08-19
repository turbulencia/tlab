#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "INISCAL"

PROGRAM INISCAL

  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE TLAB_ARRAYS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_PROCS
#endif
  USE THERMO_VARS, ONLY : imixture
  USE IO_FIELDS
  USE SCAL_LOCAL

  IMPLICIT NONE

! -------------------------------------------------------------------
  TINTEGER is, inb_scal_loc

! ###################################################################
  CALL TLAB_START()

  CALL IO_READ_GLOBAL(ifile)
  CALL SCAL_READ_LOCAL(ifile)
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(ifile)
#endif

#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

  isize_wrk3d = isize_field

  IF ( flag_s .EQ. 1 .OR. flag_s .EQ. 3 .OR. radiation%type .NE. EQNS_NONE ) THEN; inb_txc = 1
  ELSE;                                                                            inb_txc = 0
  ENDIF

  CALL TLAB_ALLOCATE(C_FILE_LOC)

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Initializing scalar fiels.')

  CALL FI_PROFILES_INITIALIZE(wrk1d)

  itime = 0; rtime = C_0_R
  s = C_0_R

#ifdef CHEMISTRY
  IF ( ireactive .EQ. CHEM_NONE ) THEN
#endif

  inb_scal_loc = inb_scal
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
    IF ( damkohler(1) .GT. C_0_R .AND. flag_mixture .EQ. 1 ) THEN
      inb_scal_loc = inb_scal - 1
    ENDIF
  ENDIF

  DO is = 1,inb_scal_loc
    CALL SCAL_MEAN(is, s(1,is), wrk1d,wrk2d,wrk3d)

    SELECT CASE( flag_s )
    CASE( 1,2,3 )
      CALL SCAL_FLUCTUATION_VOLUME( is, s(1,is), txc, wrk1d, wrk2d, wrk3d )

    CASE( 4,5,6,7,8,9 )
      CALL SCAL_FLUCTUATION_PLANE(is, s(1,is), wrk2d)

    END SELECT

  ENDDO

  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN ! Initial liquid in equilibrium; overwrite previous values
    IF ( damkohler(3) .GT. C_0_R .AND. flag_mixture .EQ. 1 ) THEN
      CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2), s(1,1), epbackground,pbackground)
    ENDIF
  ENDIF

  ! ###################################################################
  IF ( radiation%type .NE. EQNS_NONE ) THEN         ! Initial radiation effect as an accumulation during a certain interval of time

    IF ( ABS(radiation%parameters(1)) .GT. C_0_R ) THEN
      radiation%parameters(3) = radiation%parameters(3) /radiation%parameters(1) *norm_ini_radiation
    ENDIF
    radiation%parameters(1) = norm_ini_radiation
    IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Calculate q_l
      CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2), s(1,1), epbackground,pbackground)
    ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
      CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))
    ENDIF
    DO is = 1,inb_scal
      IF ( radiation%active(is) ) THEN
         CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(is)), txc, wrk1d,wrk3d)
         s(1:isize_field,is) = s(1:isize_field,is) + txc(1:isize_field,1)
      ENDIF
    ENDDO

  ENDIF

! ###################################################################
  CALL IO_WRITE_FIELDS('scal.ics', IO_SCAL, imax,jmax,kmax, inb_scal, s, wrk3d)

  CALL TLAB_STOP(0)
END PROGRAM INISCAL
