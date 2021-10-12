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

  USE TLAB_CONSTANTS, ONLY : efile, lfile
  USE TLAB_PROCS
  USE TLAB_VARS,ONLY : imode_eqns, imode_ibm, icalc_scal, inb_scal
  USE TLAB_VARS,ONLY : isize_field, imax,jmax,kmax
  USE TLAB_VARS,ONLY : rbackground
  USE DNS_LOCAL, ONLY : ilimit_flow, p_bound_min,p_bound_max, r_bound_min,r_bound_max, d_bound_max
  USE DNS_LOCAL, ONLY : ilimit_scal, s_bound_min,s_bound_max
  USE DNS_LOCAL, ONLY : logs_data
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_offset_i, ims_offset_k
#endif

  IMPLICIT NONE
  
#include "integers.h"

  TINTEGER,                        INTENT(IN)    :: flag_dilatation
  TREAL, DIMENSION(isize_field,*), INTENT(INOUT) :: q,s
  TREAL, DIMENSION(isize_field,5), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: wrk3d

! -------------------------------------------------------------------
  TREAL r_min_loc,r_max_loc, p_min_loc,p_max_loc, dummy
  TINTEGER ij, is, idummy(3)
  CHARACTER*128 line
  CHARACTER*32 str

! Pointers to existing allocated space
  TARGET q, wrk3d
  TREAL, DIMENSION(:),     POINTER :: rho, p
  TREAL, DIMENSION(:,:,:), POINTER :: loc_max

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
        IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, q(1,1),txc(1,3))
           CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, q(1,2),txc(1,4))
           CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, q(1,3),txc(1,5))
           CALL FI_INVARIANT_P(imax,jmax,kmax, txc(1,3),txc(1,4),txc(1,5), txc(1,1),txc(1,2), wrk2d,wrk3d)
        ELSE
           CALL FI_INVARIANT_P(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), wrk2d,wrk3d)
        ENDIF

        IF ( imode_ibm == 1 ) CALL IBM_BCS_FLOW(txc(1,1),i1) ! IBM - zeros in solid

        CALL MINMAX(imax,jmax,kmax, txc(1,1), logs_data(11),logs_data(10))
        logs_data(10)=-logs_data(10); logs_data(11)=-logs_data(11)

        IF ( MAX(ABS(logs_data(10)),ABS(logs_data(11))) .GT. d_bound_max ) THEN
           CALL TLAB_WRITE_ASCII(efile, 'DNS_CONTROL. Dilatation out of bounds.')
           logs_data(1) = DNS_ERROR_DILATATION

! Locating the points where the maximum dilatation occurs
           wrk3d = -txc(:,1)
           loc_max(1:imax,1:jmax,1:kmax) => wrk3d(1:imax*jmax*kmax)

           dummy = MAXVAL(wrk3d)
           IF ( ABS(dummy) .GT. d_bound_max ) THEN
              idummy = MAXLOC(loc_max)
              WRITE(str,1000) dummy; line = 'Maximum dilatation '//TRIM(ADJUSTL(str))
#ifdef USE_MPI
              idummy(1) = idummy(1) +ims_offset_i
              idummy(3) = idummy(3) +ims_offset_k
#endif
              WRITE(str,*) idummy(1); line = TRIM(ADJUSTL(line))//' at grid node '//TRIM(ADJUSTL(str))
              WRITE(str,*) idummy(2); line = TRIM(ADJUSTL(line))//':'//TRIM(ADJUSTL(str))
              WRITE(str,*) idummy(3); line = TRIM(ADJUSTL(line))//':'//TRIM(ADJUSTL(str))//'.'
              CALL TLAB_WRITE_ASCII(lfile, line, .TRUE.)
           ENDIF

           dummy = MINVAL(wrk3d)
           IF ( ABS(dummy) .GT. d_bound_max ) THEN
              idummy = MINLOC(loc_max)
              WRITE(str,1000) dummy; line = 'Minimum dilatation '//TRIM(ADJUSTL(str))
#ifdef USE_MPI
              idummy(1) = idummy(1) +ims_offset_i
              idummy(3) = idummy(3) +ims_offset_k
#endif
              WRITE(str,*) idummy(1); line = TRIM(ADJUSTL(line))//' at grid node '//TRIM(ADJUSTL(str))
              WRITE(str,*) idummy(2); line = TRIM(ADJUSTL(line))//':'//TRIM(ADJUSTL(str))
              WRITE(str,*) idummy(3); line = TRIM(ADJUSTL(line))//':'//TRIM(ADJUSTL(str))//'.'
              CALL TLAB_WRITE_ASCII(lfile, line, .TRUE.)
           ENDIF

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
           CALL TLAB_WRITE_ASCII(efile, 'DNS_CONTROL. Density out of bounds.')
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
           CALL TLAB_WRITE_ASCII(efile, 'DNS_CONTROL. Pressure out of bounds.')
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

1000 FORMAT(G_FORMAT_R)

END SUBROUTINE DNS_CONTROL
