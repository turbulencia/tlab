#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

SUBROUTINE STATS_TEMPORAL_LAGRANGIAN(x,y,z,dx,dy,dz, q,s,hq, l_q,l_hq,l_txc,l_tags, txc, mean, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*),                 INTENT(IN)    :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(isize_field,    *), INTENT(IN)    :: q,s
  TREAL, DIMENSION(isize_field,    *), INTENT(INOUT) :: hq ! Used as aux array
  TREAL, DIMENSION(isize_txc_field,*), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),                 INTENT(INOUT) :: wrk1d,wrk2d,wrk3d, mean

  TREAL, DIMENSION(isize_particle,inb_particle),     INTENT(IN)    :: l_q, l_hq
  TREAL, DIMENSION(isize_particle,inb_particle_txc), INTENT(INOUT) :: l_txc
  INTEGER(8), DIMENSION(*)                                         :: l_tags

! -------------------------------------------------------------------
  TINTEGER is
  CHARACTER*32 fname

! ###################################################################
! Particle calculations
! ###################################################################
! Lagrange Liquid and Liquid without diffusion
  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN

     l_txc = C_1_R; ! We want density
     CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, txc(1,5))
     
     hq(:,1) = hq(:,1) + 0.00000001
     DO is = inb_scal_array+2,inb_scal_particle+inb_scal_array+1
        mean_i(is) = C_1_R; delta_i(is) = C_0_R; ycoor_i(is) = ycoor_i(1); schmidt(is) = schmidt(1)
        l_txc(:,1)=l_q(:,3+is-inb_scal_array-1) !!! DO WE WANT l_txc(:,is) ???
        CALL PARTICLE_TO_FIELD(l_q,l_txc,x,y,z,wrk1d,wrk2d,wrk3d, hq(1,2))   
        hq(:,2) = hq(:,2) /hq(:,1)
        CALL AVG_SCAL_XZ(is, q,s, hq(1,2), &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4), txc(1,5),txc(1,6),mean, wrk1d,wrk2d,wrk3d)
     ENDDO
  ENDIF
  
! Save particle pathlines for particle_pdf
  IF ( icalc_particle_pdf .EQ. 1) THEN
     number_of_bins = particle_pdf_max/particle_pdf_interval
     WRITE(fname,*) itime; fname = "particle_pdf."//TRIM(ADJUSTL(fname))
     CALL PARTICLE_PDF(fname,s,wrk1d,wrk2d,wrk3d,x,y,z,l_txc,l_tags,l_hq,l_q)
  END IF
  
! Save particle residence times
  IF ( ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN
     WRITE(fname,*) itime; fname = "residence_pdf."//TRIM(ADJUSTL(fname))
     CALL RESIDENCE_PDF(fname,l_hq,l_q)
  END IF
  
  RETURN
END SUBROUTINE STATS_TEMPORAL_LAGRANGIAN
