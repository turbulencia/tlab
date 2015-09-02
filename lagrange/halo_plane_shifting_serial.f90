#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
!########################################################################
!# Tool/Library
!#
!########################################################################
!# DESCRIPTION
!#
!# shifting of data to ghost plane
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE  halo_plane_shifting_serial(q,txc,halo_field_1, halo_field_2, halo_field_3)

USE DNS_GLOBAL, ONLY: imax,jmax,kmax, imax_total, kmax_total, inb_particle
USE LAGRANGE_GLOBAL

IMPLICIT NONE
#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax,*),INTENT(in)       :: q
  TREAL, DIMENSION(imax,jmax,kmax,inb_lag_aux_field)  :: txc ! Rest of the arrays to be interpolated
  TREAL, DIMENSION(2,jmax,kmax,inb_lag_total_interp)          :: halo_field_1
  TREAL, DIMENSION(imax,jmax,2,inb_lag_total_interp)          :: halo_field_2
  TREAL, DIMENSION(2,jmax,2,inb_lag_total_interp)             :: halo_field_3
  TINTEGER i, iloc 


  DO i=1,3 
    halo_field_1(1,1:jmax,1:kmax,i)=q(imax_total,1:jmax,1:kmax,i)  
    halo_field_1(2,1:jmax,1:kmax,i)=q(1,1:jmax,1:kmax,i)  

    halo_field_2(1:imax,1:jmax,1,i)=q(1:imax,1:jmax,kmax_total,i)
    halo_field_2(1:imax,1:jmax,2,i)=q(1:imax,1:jmax,1,i)
    
    halo_field_3(1,1:jmax,1,i)=q(imax_total,1:jmax,kmax_total,i)
    halo_field_3(1,1:jmax,2,i)=q(imax_total,1:jmax,1,i)
    halo_field_3(2,1:jmax,1,i)=q(1,1:jmax,kmax_total,i)
    halo_field_3(2,1:jmax,2,i)=q(1,1:jmax,1,i)
  END DO

  DO i=4,inb_lag_total_interp
    iloc=i - 3
    halo_field_1(1,1:jmax,1:kmax,i)=txc(imax_total,1:jmax,1:kmax,iloc)  
    halo_field_1(2,1:jmax,1:kmax,i)=txc(1,1:jmax,1:kmax,iloc)  

    halo_field_2(1:imax,1:jmax,1,i)=txc(1:imax,1:jmax,kmax_total,iloc)
    halo_field_2(1:imax,1:jmax,2,i)=txc(1:imax,1:jmax,1,iloc)
    
    halo_field_3(1,1:jmax,1,i)=txc(imax_total,1:jmax,kmax_total,iloc)
    halo_field_3(1,1:jmax,2,i)=txc(imax_total,1:jmax,1,iloc)
    halo_field_3(2,1:jmax,1,i)=txc(1,1:jmax,kmax_total,iloc)
    halo_field_3(2,1:jmax,2,i)=txc(1,1:jmax,1,iloc)
  END DO


 
  RETURN
END SUBROUTINE halo_plane_shifting_serial
