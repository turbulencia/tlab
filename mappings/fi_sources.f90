#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2015/10/09 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Sources from the evolution equations.
!#
!########################################################################
SUBROUTINE FI_SOURCES_FLOW(q,s, hq, wrk1d,wrk3d)

  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, isize_field, isize_wrk1d
  USE DNS_GLOBAL,    ONLY : buoyancy, coriolis
  USE DNS_GLOBAL,    ONLY : bbackground, pbackground, rbackground, epbackground
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: q,s
  TREAL, DIMENSION(isize_field,*), INTENT(OUT)   :: hq
  TREAL, DIMENSION(isize_wrk1d,*), INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: wrk3d

! -----------------------------------------------------------------------
  TINTEGER ij, iq
  TREAL dummy, u_geo, w_geo

  TINTEGER siz, srt, end    !  Variables for OpenMP Partitioning 

#ifdef USE_BLAS
  INTEGER ilen
#endif

! #######################################################################
#ifdef USE_BLAS
  ilen = isize_field
#endif

! -----------------------------------------------------------------------
! Coriolis. So far, rotation only in the Oy direction. 
! -----------------------------------------------------------------------
  IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
     u_geo = COS(coriolis%parameters(1))
     w_geo =-SIN(coriolis%parameters(1))

!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 

     dummy = coriolis%vector(2)
     DO ij = srt,end
        hq(ij,1) = hq(ij,1) + dummy*( w_geo-q(ij,3) )
        hq(ij,3) = hq(ij,3) + dummy*( q(ij,1)-u_geo ) 
     ENDDO
  ENDIF

! -----------------------------------------------------------------------
! Buoyancy. Remember that buoyancy%vector contains the Froude # already.
! -----------------------------------------------------------------------
  DO iq = 1,3
     IF ( buoyancy%active(iq) ) THEN

        IF ( buoyancy%type .EQ. EQNS_EXPLICIT ) THEN
           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
              CALL THERMO_AIRWATER_BUOYANCY(imax,jmax,kmax, s(1,2),s(1,1), epbackground,pbackground,rbackground, wrk3d)
           ENDIF
        ELSE
           IF ( iq .EQ. 2 ) THEN
              CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, bbackground)
           ELSE
              wrk1d(:,1) = C_0_R
              CALL FI_BUOYANCY(buoyancy, imax,jmax,kmax, s, wrk3d, wrk1d)
           ENDIF
        ENDIF

!$omp parallel default( shared ) &
#ifdef USE_BLAS
!$omp private( ilen, dummy, srt,end,siz)
#else     
!$omp private( ij,   dummy, srt,end,siz )
#endif
        CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 
        
        dummy = buoyancy%vector(iq)
#ifdef USE_BLAS
        ilen = siz
        CALL DAXPY(ilen, dummy, wrk3d(srt), 1, hq(srt,iq), 1)
#else
        DO ij = srt,end
           hq(ij,iq) = hq(ij,iq) + dummy* wrk3d(ij)
        ENDDO
#endif
!$omp end parallel
        
     ENDIF

  ENDDO

  RETURN
END SUBROUTINE FI_SOURCES_FLOW

! #######################################################################
! #######################################################################
SUBROUTINE FI_SOURCES_SCAL(s, hs, tmp1,tmp2, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, inb_scal, isize_field, isize_wrk1d
  USE DNS_GLOBAL, ONLY : g
  USE DNS_GLOBAL, ONLY : radiation, transport, chemistry

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: s
  TREAL, DIMENSION(isize_field,*), INTENT(OUT)   :: hs
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2
  TREAL, DIMENSION(isize_wrk1d,*), INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: wrk2d,wrk3d
  
! -----------------------------------------------------------------------
  TINTEGER ij, is, flag_grad

  TINTEGER siz, srt, end    !  Variables for OpenMP Partitioning 
  
#ifdef USE_BLAS
  INTEGER ilen
#endif

! #######################################################################
#ifdef USE_BLAS
  ilen = isize_field
#endif

  DO is = 1,inb_scal ! Start loop over the N scalars  

! -----------------------------------------------------------------------
! Radiation
! -----------------------------------------------------------------------
     IF ( radiation%active(is) ) THEN
        CALL OPR_RADIATION(radiation, imax,jmax,kmax, g(2), s(1,radiation%scalar(is)), tmp1, wrk1d,wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
        
        DO ij = srt,end
           hs(ij,is) = hs(ij,is) + tmp1(ij)
        ENDDO
!$omp end parallel
        
     ENDIF
     
! -----------------------------------------------------------------------
! Transport, such as settling 
! array tmp2 should not be used inside the loop on is
! -----------------------------------------------------------------------
     IF ( transport%active(is) ) THEN
        IF ( is .EQ. 1 ) THEN; flag_grad = 1;
        ELSE;                  flag_grad = 0; ENDIF
        CALL FI_TRANS_FLUX(transport, flag_grad, imax,jmax,kmax, is, s,tmp1, tmp2, wrk2d,wrk3d)
        
!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
        
        DO ij = srt,end
           hs(ij,is) = hs(ij,is) + tmp1(ij)
        ENDDO
!$omp end parallel
        
     ENDIF

! -----------------------------------------------------------------------
! Chemistry
! -----------------------------------------------------------------------
     IF ( chemistry%active(is) ) THEN
        CALL FI_CHEM(chemistry, imax,jmax,kmax, is, s, tmp1)
        
!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)

        DO ij = srt,end
           hs(ij,is) = hs(ij,is) + tmp1(ij)
        ENDDO
!$omp end parallel

     ENDIF
     
  ENDDO
  
  RETURN
END SUBROUTINE FI_SOURCES_SCAL
