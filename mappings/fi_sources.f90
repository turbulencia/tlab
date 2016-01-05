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
SUBROUTINE FI_SOURCES_FLOW(q,s, hq, b_ref, wrk1d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, isize_field, isize_wrk1d
  USE DNS_GLOBAL, ONLY : ibodyforce, body_param, body_vector
  USE DNS_GLOBAL, ONLY : icoriolis, rotn_param, rotn_vector

  IMPLICIT NONE

  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: q,s
  TREAL, DIMENSION(isize_field,*), INTENT(OUT)   :: hq
  TREAL, DIMENSION(jmax),          INTENT(IN)    :: b_ref
  TREAL, DIMENSION(isize_wrk1d,*), INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk3d

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
  IF ( icoriolis .EQ. EQNS_COR_NORMALIZED ) THEN
     u_geo = COS(rotn_param(1))
     w_geo =-SIN(rotn_param(1))

!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
     CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 

     dummy = rotn_vector(2)
     DO ij = srt,end
        hq(ij,1) = hq(ij,1) + dummy*( w_geo-q(ij,3) )
        hq(ij,3) = hq(ij,3) + dummy*( q(ij,1)-u_geo ) 
     ENDDO
  ENDIF

! -----------------------------------------------------------------------
! Buoyancy. Remember that body_vector contains the Froude # already.
! -----------------------------------------------------------------------
  IF ( ibodyforce .NE. EQNS_NONE ) THEN
     DO iq = 1,3
        IF ( ABS(body_vector(iq)) .GT. C_0_R ) THEN
           
           IF ( iq .EQ. 2 ) THEN
              CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s, wrk3d, b_ref)
           ELSE
              wrk1d(:,1) = C_0_R
              CALL FI_BUOYANCY(ibodyforce, imax,jmax,kmax, body_param, s, wrk3d, wrk1d)
           ENDIF
           
!$omp parallel default( shared ) &
#ifdef USE_BLAS
!$omp private( ilen, dummy, srt,end,siz)
#else     
!$omp private( ij,   dummy, srt,end,siz )
#endif
           CALL DNS_OMP_PARTITION(isize_field,srt,end,siz) 
           
           dummy = body_vector(iq)
#ifdef USE_BLAS
           ilen = siz
           CALL DAXPY(ilen, dummy, wrk3d(srt), 1, hq(srt,iq), 1)
#else
           DO ij = srt,end
              hq(ij,iq) = hq(ij,iq) + dummy*wrk3d(ij)
           ENDDO
#endif
!$omp end parallel

        ENDIF

     ENDDO
  ENDIF

  RETURN
END SUBROUTINE FI_SOURCES_FLOW

! #######################################################################
! #######################################################################
SUBROUTINE FI_SOURCES_SCAL(y,dy, s, hs, tmp1,tmp2,tmp3,tmp4, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, inb_scal, inb_scal_array, isize_field, isize_wrk1d
  USE DNS_GLOBAL, ONLY : scaley, ycoor_i
  USE DNS_GLOBAL, ONLY : iradiation, rad_param, irad_scalar
  USE DNS_GLOBAL, ONLY : itransport, trans_param, settling
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

  TREAL, DIMENSION(*),             INTENT(IN)    :: y, dy
  TREAL, DIMENSION(isize_field,*), INTENT(IN)    :: s
  TREAL, DIMENSION(isize_field,*), INTENT(OUT)   :: hs
  TREAL, DIMENSION(isize_field),   INTENT(INOUT) :: tmp1,tmp2,tmp3,tmp4
  TREAL, DIMENSION(isize_wrk1d,*), INTENT(INOUT) :: wrk1d
  TREAL, DIMENSION(*),             INTENT(INOUT) :: wrk2d,wrk3d
  
! -----------------------------------------------------------------------
  TINTEGER ij, is, i,j,k
  TREAL xi, ycenter, yrel, dummy

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
     IF ( irad_scalar .EQ. is ) THEN
        CALL OPR_RADIATION(iradiation, imax,jmax,kmax, dy, rad_param, s(:,inb_scal_array), tmp4, wrk1d,wrk3d)
        
!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
        
        DO ij = srt,end
           hs(ij,is) = hs(ij,is) + tmp4(ij)
        ENDDO
!$omp end parallel
        
     ENDIF
     
! -----------------------------------------------------------------------
! Relaxation
! -----------------------------------------------------------------------
     IF ( imixture .EQ. MIXT_TYPE_BILAIRWATERSTRAT .AND. is .EQ. 2 ) THEN
        
        ycenter=y(1)+scaley*ycoor_i(3)+0.005 !postition of s3 plus constant
        DO i=1,imax
           DO k=1,kmax
              DO j=1,jmax
                 yrel=y(j)-ycenter
                 xi=yrel/0.2 !thicknes constant
                 wrk3d(i+(j-1)*imax+(k-1)*imax*jmax) = (1+TANH(xi))/2 !strength constant
              ENDDO
           ENDDO
        ENDDO
        
!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
        CALL DNS_OMP_PARTITION(isize_field,srt,end,siz)
        
        dummy = -C_1_R/0.02
        DO ij = srt,end
           hs(ij,is) = hs(ij,is) + dummy*wrk3d(ij)*s(ij,is)
        ENDDO
!$omp end parallel
        
     ENDIF

! -----------------------------------------------------------------------
! Transport, such as settling 
! array tmp1 should not be used inside the loop on is
! -----------------------------------------------------------------------
     IF ( itransport .GT. 0 ) THEN  
        CALL FI_TRANS_FLUX(itransport, imax,jmax,kmax, is, inb_scal_array, trans_param, settling, &
             dy, s,tmp1, tmp4, wrk1d,wrk2d,wrk3d)
        
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
