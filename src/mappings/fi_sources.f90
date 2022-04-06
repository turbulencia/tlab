#include "types.h"
#include "dns_const.h"

!########################################################################
!#
!# Sources from the evolution equations.
!#
!########################################################################
SUBROUTINE FI_SOURCES_FLOW(q,s, hq, tmp1, wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field, isize_wrk1d
  USE TLAB_VARS, ONLY : buoyancy, coriolis, subsidence, random
  USE TLAB_VARS, ONLY : bbackground, pbackground, rbackground, epbackground

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*), INTENT(IN   ) :: q,s
  TREAL, DIMENSION(isize_field,*), INTENT(  OUT) :: hq
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1
  TREAL, DIMENSION(isize_wrk1d,*), INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: wrk3d

  ! -----------------------------------------------------------------------
  TINTEGER ij, iq
  TREAL dummy, u_geo, w_geo, dtr1, dtr3

  TINTEGER siz, srt, END    !  Variables for OpenMP Partitioning

#ifdef USE_BLAS
  INTEGER ILEN
#endif

  ! #######################################################################
#ifdef USE_BLAS
  ILEN = isize_field
#endif

  ! -----------------------------------------------------------------------
  ! Coriolis. Remember that coriolis%vector already contains the Rossby #.
  ! -----------------------------------------------------------------------
  SELECT CASE ( coriolis%type )
  CASE( EQNS_EXPLICIT )
     DO ij = 1, isize_field
       hq(ij,1) = hq(ij,1) + coriolis%vector(3)*q(ij,2) - coriolis%vector(2)*q(ij,3)
       hq(ij,2) = hq(ij,2) + coriolis%vector(1)*q(ij,3) - coriolis%vector(3)*q(ij,1)
       hq(ij,3) = hq(ij,3) + coriolis%vector(2)*q(ij,1) - coriolis%vector(1)*q(ij,2)
     END DO

  CASE( EQNS_COR_NORMALIZED ) ! So far, rotation only in the Oy direction.
     u_geo = COS(coriolis%parameters(1)) *coriolis%parameters(2)
     w_geo =-SIN(coriolis%parameters(1)) *coriolis%parameters(2)

!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
     CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)

     dummy = coriolis%vector(2)
     dtr3=C_0_R; dtr1=C_0_R
     DO ij = srt,END
        hq(ij,1) = hq(ij,1) + dummy*( w_geo-q(ij,3) )
        hq(ij,3) = hq(ij,3) + dummy*( q(ij,1)-u_geo )
     END DO
!$omp end parallel

  END SELECT

  ! -----------------------------------------------------------------------
  DO iq = 1,3

     ! -----------------------------------------------------------------------
     ! Buoyancy. Remember that buoyancy%vector already contains the Froude #.
     ! -----------------------------------------------------------------------
     IF ( buoyancy%active(iq) ) THEN

        IF ( buoyancy%type == EQNS_EXPLICIT ) THEN
           CALL THERMO_ANELASTIC_BUOYANCY(imax,jmax,kmax, s, epbackground,pbackground,rbackground, wrk3d)

        ELSE
           IF ( iq == 2 ) THEN
              CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, bbackground)
           ELSE
              wrk1d(:,1) = C_0_R
              CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
           END IF
        END IF

#ifdef USE_BLAS
!$omp parallel default( shared ) &
!$omp private( ilen, dummy, srt,end,siz)
#else
!$omp parallel default( shared ) &
!$omp private( ij,   dummy, srt,end,siz )
#endif
        CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)

        dummy = buoyancy%vector(iq)
#ifdef USE_BLAS
        ILEN = siz
        CALL DAXPY(ILEN, dummy, wrk3d(srt), 1, hq(srt,iq), 1)
#else
        DO ij = srt,END
           hq(ij,iq) = hq(ij,iq) + dummy* wrk3d(ij)
        END DO
#endif
!$omp end parallel

     END IF

     ! -----------------------------------------------------------------------
     ! Subsidence
     ! -----------------------------------------------------------------------
     IF ( subsidence%active(iq) ) THEN
        CALL FI_SUBSIDENCE(subsidence, imax,jmax,kmax, q(1,iq), tmp1, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)

        DO ij = srt,END
           hq(ij,iq) = hq(ij,iq) + tmp1(ij)
        END DO
!$omp end parallel

     END IF

     ! -----------------------------------------------------------------------
     ! Random forcing
     ! -----------------------------------------------------------------------
     IF ( random%active(iq) ) THEN
        CALL FI_RANDOM(random,imax,jmax,kmax,hq(1,iq),tmp1)
     ENDIF
  END DO

  RETURN
END SUBROUTINE FI_SOURCES_FLOW

! #######################################################################
! #######################################################################
SUBROUTINE FI_SOURCES_SCAL(s, hs, tmp1,tmp2, wrk1d,wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : imax,jmax,kmax, inb_scal, isize_field, isize_wrk1d
  USE TLAB_VARS, ONLY : imode_eqns
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : radiation, transport, chemistry, subsidence
  USE TLAB_VARS, ONLY : rbackground, ribackground

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: s
  TREAL, DIMENSION(isize_field,*), INTENT(OUT)   :: hs
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(isize_wrk1d,*), INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: wrk3d

  ! -----------------------------------------------------------------------
  TINTEGER ij, is, flag_grad

  TINTEGER siz, srt, END    !  Variables for OpenMP Partitioning

#ifdef USE_BLAS
  INTEGER ILEN
#endif

  ! #######################################################################
#ifdef USE_BLAS
  ILEN = isize_field
#endif

  DO is = 1,inb_scal ! Start loop over the N scalars

     ! -----------------------------------------------------------------------
     ! Radiation
     ! -----------------------------------------------------------------------
     IF ( radiation%active(is) ) THEN
        IF ( imode_eqns == DNS_EQNS_ANELASTIC ) THEN
           CALL THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax,jmax,kmax, rbackground, s(1,radiation%scalar(is)), tmp2)
           CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), tmp2,                      tmp1, wrk1d,wrk3d)
           CALL THERMO_ANELASTIC_WEIGHT_ADD(imax,jmax,kmax, ribackground, tmp1, hs(1,is))

        ELSE
           CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(is)), tmp1, wrk1d,wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
           CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)

           DO ij = srt,END
              hs(ij,is) = hs(ij,is) + tmp1(ij)
           END DO
!$omp end parallel

        END IF

     END IF

     ! -----------------------------------------------------------------------
     ! Transport, such as settling
     ! array tmp2 should not be used inside the loop on is
     ! -----------------------------------------------------------------------
     IF ( transport%active(is) ) THEN
        IF ( is == 1 ) THEN; flag_grad = 1;
        ELSE;                  flag_grad = 0
        END IF
        CALL FI_TRANSPORT(transport, flag_grad, imax,jmax,kmax, is, s,tmp1, tmp2, wrk2d,wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)

        DO ij = srt,END
           hs(ij,is) = hs(ij,is) + tmp1(ij)
        END DO
!$omp end parallel

     END IF

     ! -----------------------------------------------------------------------
     ! Chemistry
     ! -----------------------------------------------------------------------
     IF ( chemistry%active(is) ) THEN
        CALL FI_CHEM(chemistry, imax,jmax,kmax, is, s, tmp1)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)

        DO ij = srt,END
           hs(ij,is) = hs(ij,is) + tmp1(ij)
        END DO
!$omp end parallel

     END IF

     ! -----------------------------------------------------------------------
     ! Subsidence
     ! -----------------------------------------------------------------------
     IF ( subsidence%active(is) ) THEN
        CALL FI_SUBSIDENCE(subsidence, imax,jmax,kmax, s(1,is), tmp1, wrk1d,wrk2d,wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,END,siz)

        DO ij = srt,END
           hs(ij,is) = hs(ij,is) + tmp1(ij)
        END DO
!$omp end parallel

     END IF

  END DO

  RETURN
END SUBROUTINE FI_SOURCES_SCAL
