#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library INIT/SCAL
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!# 2007/05/09 - J.P. Mellado
!#              Adding multispecies
!# 2007/08/17 - J.P. Mellado
!#              Adding plane perturbation
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM INISCAL

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE SCAL_LOCAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
  TREAL, DIMENSION(:),   ALLOCATABLE :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(:,:), ALLOCATABLE :: s
  TREAL, DIMENSION(:),   ALLOCATABLE :: wrk1d, wrk2d, wrk3d, txc

  TREAL dummy
  TINTEGER isize_wrk3d, is, ierr

  CHARACTER*64 str, line
  CHARACTER*32 inifile

! ###################################################################
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  CALL SCAL_READ_LOCAL(inifile)
#ifdef CHEMISTRY
  CALL CHEM_READ_GLOBAL(inifile)
#endif

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  CALL IO_WRITE_ASCII(lfile,'Initializing scalar fiels.')

  itime = 0; rtime = C_0_R
  
  isize_wrk3d = isize_field
  isize_wrk3d = MAX(isize_wrk1d*300, isize_wrk3d)

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------      
  ALLOCATE(x(imax_total))
  ALLOCATE(y(jmax_total))
  ALLOCATE(z(kmax_total))
  ALLOCATE(dx(imax_total*inb_grid))
  ALLOCATE(dy(jmax_total*inb_grid))
  ALLOCATE(dz(kmax_total*inb_grid))

  WRITE(str,*) inb_scal_array; line = 'Allocating array scal. Size '//TRIM(ADJUSTL(str))//'x'
  WRITE(str,*) isize_field;    line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
  CALL IO_WRITE_ASCII(lfile,line)
  ALLOCATE(s(isize_field,inb_scal_array),stat=ierr)
  IF ( ierr .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, 'INISCAL. Not enough memory for scal.')
     CALL DNS_STOP(DNS_ERROR_ALLOC)
  ENDIF
  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk3d(isize_wrk3d))

  IF ( flag_s .EQ. 2 .OR. flag_s .EQ. 3 .OR. &
       iradiation .NE. EQNS_NONE ) THEN
     WRITE(str,*) i1;          line = 'Allocating array txc. Size '//TRIM(ADJUSTL(str))//'x'
     WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
     CALL IO_WRITE_ASCII(lfile,line)
     ALLOCATE(txc(isize_field),stat=ierr)
     IF ( ierr .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile, 'INISCAL. Not enough memory for txc.')
        CALL DNS_STOP(DNS_ERROR_ALLOC)
     ENDIF
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN; ALLOCATE(wrk2d(isize_wrk2d*5))
  ELSE;                                        ALLOCATE(wrk2d(isize_wrk2d  ));  ENDIF

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! ###################################################################
  s = C_0_R

! ###################################################################
! Non-reacting case
! ###################################################################
  IF ( ireactive .EQ. CHEM_NONE ) THEN

! -------------------------------------------------------------------
! Mean
! -------------------------------------------------------------------
     DO is = 1,inb_scal
        !Start the supersaturation simulation with a concentration in equilibrium
        IF ( (is .EQ. 3) .AND. (imixture .EQ. MIXT_TYPE_SUPSAT) .AND. (flag_mixture .EQ. 1) ) THEN 
           dummy = p_init/C_1_R ! MRATIO = 1
           CALL  THERMO_AIRWATER_PHAL(imax,jmax,kmax, s(1,2), dummy, s(1,1))
        ELSE
           CALL SCAL_MEAN(is, x,y,z, s(1,is), wrk1d,wrk2d,wrk3d)
        END IF
     END DO

! -------------------------------------------------------------------
! Fluctuation field
! -------------------------------------------------------------------
      DO is = 1,inb_scal
        IF ((is .EQ. 3) .AND. (imixture .EQ. MIXT_TYPE_SUPSAT) .AND. (flag_mixture .EQ. 1)) THEN !In Eq concentration do nothing
        ELSE
           IF      ( flag_s .EQ. 1 ) THEN
              CALL SCAL_VOLUME_DISCRETE(is, x,y,z, s(1,is))
           ELSE IF ( flag_s .EQ. 2 ) THEN
              CALL SCAL_VOLUME_BROADBAND(is, y,dx,dz, s(1,is), txc, wrk3d)
           ELSE IF ( flag_s .GE. 4 ) THEN ! rest of the cases
              CALL SCAL_PLANE(flag_s, is, x,y,z,dx,dz, s(1,is), wrk2d)
           ENDIF
           
        ENDIF
     ENDDO

! ###################################################################
! Reacting case
! ###################################################################
  ELSE
#ifdef CHEMISTRY
     is = inb_scal

! pasive scalar field
     IF      ( flag_mixture .EQ. 0 ) THEN
        CALL SCAL_MEAN(is, x,y,z, s(1,is), wrk1d,wrk2d,wrk3d)
     ELSE IF ( flag_mixture .EQ. 2 ) THEN
        CALL DNS_READ_FIELDS('scal.ics', i1, imax,jmax,kmax, i1,i1, isize_wrk3d, s(1,is), wrk3d)
     ENDIF

! species mass fractions
     IF      ( ireactive .EQ. CHEM_FINITE ) THEN
        CALL SCREACT_FINITE(x, s, isize_wrk3d, wrk3d)
     ELSE IF ( ireactive .EQ. CHEM_INFINITE .AND. inb_scal .GT. 1 ) THEN
        CALL SCREACT_INFINITE(x, s, isize_wrk3d, wrk3d)
     ENDIF
#else
     CALL IO_WRITE_ASCII(efile, 'INISCAL. Chemistry part to be checked')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
#endif

  ENDIF

! ------------------------------------------------------------------
! Add Radiation component after the fluctuation field
! ------------------------------------------------------------------
  IF ( iradiation .NE. EQNS_NONE ) THEN
     
! An initial effect of radiation is imposed as an accumulation during a certain interval of time
     rad_param(1) = norm_ini_radiation
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))
     ENDIF
     CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(1,inb_scal_array), txc, wrk1d,wrk3d)
     s(1:isize_field,irad_scalar) = s(1:isize_field,irad_scalar) + txc(1:isize_field)

  ENDIF

! ###################################################################
! Output file
! ###################################################################
  CALL DNS_WRITE_FIELDS('scal.ics', i1, imax,jmax,kmax, inb_scal, isize_wrk3d, s, wrk3d)

  CALL DNS_END(0)

  STOP
END PROGRAM INISCAL
