!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/05/18 - J.P. Mellado
!#              Cleaned. Removing routines WRT_FLOWRES, WRT_SCALRES.
!# 2010/11/01 - J.P. Mellado
!#              Split into flow and scal routines
!#
!########################################################################
!# DESCRIPTION
!#
!# Write log file.
!#
!# Check values of density and pressure to make sure they are positive;
!# Parameters ?_bound_min are used as threshold.
!#
!# Only in compressible mode
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

! ###################################################################
! Flow
! ###################################################################
SUBROUTINE DNS_CONTROL_FLOW(q,s, wrk3d)

  USE DNS_GLOBAL,ONLY : icalc_flow, icalc_scal
  USE DNS_GLOBAL,ONLY : isize_field, imax,jmax,kmax, inb_flow, inb_scal
  USE DNS_GLOBAL,ONLY : itime
  USE DNS_CONSTANTS,ONLY : efile, tag_flow, tag_scal
  USE DNS_LOCAL, ONLY : ilimit_flow, p_bound_min,p_bound_max, r_bound_min,r_bound_max
  USE DNS_LOCAL, ONLY : logs_data

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(isize_field,*), TARGET :: q
  TREAL, DIMENSION(*)                     :: s, wrk3d

! -------------------------------------------------------------------
  TREAL r_min_loc, r_max_loc, p_min_loc, p_max_loc
  TINTEGER ij, ierr_loc
  CHARACTER*32 fname

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: rho, p

! ###################################################################
  ierr_loc = 0

! Define pointers
  rho => q(:,5)
  p   => q(:,6)

! ###################################################################
  CALL MINMAX(imax,jmax,kmax, rho, r_min_loc,r_max_loc)
  CALL MINMAX(imax,jmax,kmax, p,   p_min_loc,p_max_loc)

! -------------------------------------------------------------------
! Check density
! -------------------------------------------------------------------
  IF ( r_min_loc .LE. r_bound_min .OR. r_max_loc .GE. r_bound_max ) THEN
     IF ( ilimit_flow .EQ. 1 ) THEN
        DO ij = 1,isize_field
           rho(ij) = MIN(MAX(rho(ij),r_bound_min), r_bound_max)
        ENDDO
        
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_CONTROL_FLOW. Density out of bounds.')
        ierr_loc = 1
           
     ENDIF
  ENDIF
     
! -------------------------------------------------------------------
! Check pressure
! -------------------------------------------------------------------
  IF ( p_min_loc .LE. p_bound_min .OR. p_max_loc .GE. p_bound_max ) THEN
     IF ( ilimit_flow .EQ. 1 ) THEN
        DO ij = 1,isize_field
           p(ij) = MIN(MAX(p(ij),p_bound_min), p_bound_max)
        ENDDO
           
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_CONTROL_FLOW. Pressure out of bounds.')
        ierr_loc = 2
        
     ENDIF
  ENDIF

! -------------------------------------------------------------------
! Pass data to global vairables for logfiles
! -------------------------------------------------------------------
  logs_data(1) = ierr_loc  ! Status
  logs_data(5) = p_min_loc
  logs_data(6) = p_max_loc
  logs_data(7) = r_min_loc
  logs_data(8) = r_max_loc
  
! -------------------------------------------------------------------
! Error control
! -------------------------------------------------------------------
  IF ( ierr_loc .GT. 0 ) THEN
     IF ( icalc_flow .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i2, imax,jmax,kmax, inb_flow, isize_field, q, wrk3d)
     ENDIF
     IF ( icalc_scal .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        CALL DNS_WRITE_FIELDS(fname, i1, imax,jmax,kmax, inb_scal, isize_field, s, wrk3d)
     ENDIF
     
     CALL DNS_LOGS(i2)
     
     IF      ( ierr_loc .EQ. 1 ) THEN; CALL DNS_STOP(DNS_ERROR_NEGDENS)
     ELSE IF ( ierr_loc .EQ. 2 ) THEN; CALL DNS_STOP(DNS_ERROR_NEGPRESS); ENDIF
        
  ENDIF
     
  RETURN
END SUBROUTINE DNS_CONTROL_FLOW

! ###################################################################
! Check scalar
! ###################################################################
SUBROUTINE DNS_CONTROL_SCAL(s)

  USE DNS_GLOBAL,ONLY : imode_eqns, isize_field, icalc_scal
  USE DNS_GLOBAL,ONLY : inb_scal
  USE DNS_LOCAL, ONLY : ilimit_scal
  USE DNS_LOCAL, ONLY : z_bound_min, z_bound_max

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*) :: s

! -------------------------------------------------------------------
  TINTEGER ij, is

! ###################################################################
  IF ( ilimit_scal .EQ. 1 .AND. icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal
        DO ij = 1,isize_field
           s(ij,is) = MIN(MAX(s(ij,is),z_bound_min), z_bound_max)
        ENDDO
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE DNS_CONTROL_SCAL
