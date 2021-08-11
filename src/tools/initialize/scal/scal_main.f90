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
  USE TLAB_MPI_PROCS
#endif
  USE THERMO_VARS, ONLY : imixture
  USE SCAL_LOCAL

  IMPLICIT NONE

! -------------------------------------------------------------------
  TINTEGER is, inb_scal_loc

! ###################################################################
  CALL TLAB_START()

  CALL DNS_READ_GLOBAL(ifile)
  CALL SCAL_READ_LOCAL(ifile)
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(ifile)
#endif

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
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

#ifdef USE_MPI
  CALL SCAL_MPIO_AUX ! Needed for options 4, 6, 8
#else
  io_aux(1)%offset = 52 ! header size in bytes
#endif

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
#ifdef CHEMISTRY
  ELSE
    is = inb_scal

    IF      ( flag_mixture .EQ. 0 ) THEN            ! pasive scalar field
      CALL SCAL_MEAN(is, s(1,is), wrk1d,wrk2d,wrk3d)
    ELSE IF ( flag_mixture .EQ. 2 ) THEN
      CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, i1,i1, isize_wrk3d, s(1,is), wrk3d)
    ENDIF

    IF      ( ireactive .EQ. CHEM_FINITE ) THEN     ! species mass fractions
      CALL SCREACT_FINITE(x, s, isize_wrk3d, wrk3d)
    ELSE IF ( ireactive .EQ. CHEM_INFINITE .AND. inb_scal .GT. 1 ) THEN
      CALL SCREACT_INFINITE(x, s, isize_wrk3d, wrk3d)
    ENDIF
    CALL TLAB_WRITE_ASCII(efile, 'INISCAL. Chemistry part to be checked')
    CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)

  ENDIF
#endif

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
  CALL DNS_WRITE_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, isize_wrk3d, s, wrk3d)

  CALL TLAB_STOP(0)
END PROGRAM INISCAL

! ###################################################################
! ###################################################################
#ifdef USE_MPI

SUBROUTINE SCAL_MPIO_AUX()

  USE TLAB_VARS, ONLY : imax,kmax
  USE TLAB_VARS, ONLY : io_aux
  USE TLAB_MPI_VARS

  IMPLICIT NONE

#include "mpif.h"

! -----------------------------------------------------------------------
  TINTEGER                :: ndims, id
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset

! #######################################################################
  io_aux(:)%active = .FALSE. ! defaults
  io_aux(:)%offset = 0

! ###################################################################
! Subarray information to read plane data
! ###################################################################
  id = 1

  io_aux(id)%active = .TRUE.
  io_aux(id)%communicator = MPI_COMM_WORLD
  io_aux(:)%offset  = 52 ! size of header in bytes

  ndims = 3
  sizes(1)  =imax *ims_npro_i; sizes(2)   = 1; sizes(3)   = kmax *ims_npro_k
  locsize(1)=imax;             locsize(2) = 1; locsize(3) = kmax
  offset(1) =ims_offset_i;     offset(2)  = 0; offset(3)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
       MPI_ORDER_FORTRAN, MPI_REAL8, io_aux(id)%subarray, ims_err)
  CALL MPI_Type_commit(io_aux(id)%subarray, ims_err)

  RETURN
END SUBROUTINE SCAL_MPIO_AUX

#endif
