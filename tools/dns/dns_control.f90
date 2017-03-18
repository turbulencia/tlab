#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# DESCRIPTION
!#
!# Check values of density and pressure to make sure they are positive;
!# Parameters ?_bound_min are used as threshold.
!#
!########################################################################
SUBROUTINE DNS_CONTROL(flag_dilatation, q,s, txc, wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,ONLY : imode_eqns, icalc_scal, inb_scal
  USE DNS_GLOBAL,ONLY : isize_field, imax,jmax,kmax, itime
  USE DNS_LOCAL, ONLY : ilimit_flow, p_bound_min,p_bound_max, r_bound_min,r_bound_max, d_bound_max
  USE DNS_LOCAL, ONLY : ilimit_scal, s_bound_min,s_bound_max
  USE DNS_LOCAL, ONLY : logs_data

  IMPLICIT NONE

  TINTEGER,                        INTENT(IN)    :: flag_dilatation
  TREAL, DIMENSION(isize_field,*), INTENT(INOUT) :: q,s
  TREAL, DIMENSION(isize_field,2), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d

! -------------------------------------------------------------------
  TREAL r_min_loc,r_max_loc, p_min_loc,p_max_loc
  TINTEGER ij, is

! Pointers to existing allocated space
  TARGET q
  TREAL, DIMENSION(:), POINTER :: rho, p

! ###################################################################
! Scalars
! ###################################################################
  IF ( ilimit_scal .EQ. 1 .AND. icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal
        DO ij = 1,isize_field
           s(ij,is) = MIN(MAX(s(ij,is),s_bound_min(is)), s_bound_max(is))
        ENDDO
     ENDDO
  ENDIF
  
! ###################################################################
! Incompressible flow
! ###################################################################
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     IF ( flag_dilatation .EQ. 0 ) THEN
        CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), wrk2d,wrk3d)
        
        CALL MINMAX(imax,jmax,kmax, txc(1,1), logs_data(11),logs_data(10))
        logs_data(10)=-logs_data(10); logs_data(11)=-logs_data(11)
        
        IF ( MAX(ABS(logs_data(10)),ABS(logs_data(11))) .GT. d_bound_max ) THEN
           CALL IO_WRITE_ASCII(efile, 'DNS_CONTROL. Dilatation out of bounds.')
           logs_data(1) = DNS_ERROR_DILATATION
        ENDIF
     ENDIF
     
  ELSE
! ###################################################################
! Compressible flow
! ###################################################################
! Define pointers
     rho => q(:,5)
     p   => q(:,6)

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
           CALL IO_WRITE_ASCII(efile, 'DNS_CONTROL. Density out of bounds.')
           logs_data(1) = DNS_ERROR_NEGDENS
           
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
           CALL IO_WRITE_ASCII(efile, 'DNS_CONTROL. Pressure out of bounds.')
           logs_data(1) = DNS_ERROR_NEGPRESS
        
        ENDIF
     ENDIF
     
! -------------------------------------------------------------------
! Pass data to global vairables for logfiles
! -------------------------------------------------------------------
     logs_data(5) = p_min_loc
     logs_data(6) = p_max_loc
     logs_data(7) = r_min_loc
     logs_data(8) = r_max_loc
           
  ENDIF
  
  RETURN
END SUBROUTINE DNS_CONTROL
