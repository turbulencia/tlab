#include "types.h"
#include "dns_const.h"
#include "avgij_map.h"

!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2008/01/08 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Accumulating information for running statistics in spatial mode
!#
!########################################################################
SUBROUTINE DNS_SAVE(icount_stat, q,h,z1, txc, vaux, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER icount_stat
  TREAL, DIMENSION(imax*jmax*kmax,*) :: q, h, z1, txc
  TREAL, DIMENSION(*)                :: wrk1d, wrk2d, wrk3d, vaux

  TARGET q

! -------------------------------------------------------------------
  TINTEGER ie, is, bcs(2,1)

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w, T, rho, p, vis

! ###################################################################
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)
  rho => q(:,5)
  p   => q(:,6)
  T   => q(:,7)
  vis => q(:,8)

  bcs = 0
  
! ###################################################################
! Save running averages
! ###################################################################
  IF ( frunstat .EQ. 1 ) THEN
     nstatavg_points = nstatavg_points + kmax_total

     CALL DNS_SAVE_AVGIJ(rho,u,v,w,p,vis,T, &
          txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), &
          h(1,1), h(1,2), h(1,3), h(1,4), h(1,5), & 
          vaux(vindex(VA_MEAN_WRK)), wrk1d,wrk2d,wrk3d)

     IF ( icalc_scal .EQ. 1 ) THEN
        CALL DNS_SAVE_SCBDGIJ(rho,u,v,w,p,z1,vis, &
             txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6),txc(1,7), &
             h(1,1), h(1,2), h(1,3), h(1,4), h(1,5), & 
             vaux(vindex(VA_MEAN_WRK)+MA_MOMENTUM_SIZE*nstatavg*jmax), wrk1d,wrk2d,wrk3d)
     ENDIF

  ENDIF

! ###################################################################
! Save data into plane and line arrays
! ###################################################################
  CALL DNS_SAVE_T(icount_stat, vaux(vindex(VA_TIMES)))

  IF ( frunline .EQ. 1 ) THEN
     CALL DNS_SAVE_IJ(icount_stat, rho, u, v, w, p, z1, vaux(vindex(VA_LINE_SPA_WRK)) )
  ENDIF

  IF ( frunplane .EQ. 1 ) THEN
     CALL DNS_SAVE_I(icount_stat, rho, u, v, w, p, z1, T, vaux(vindex(VA_PLANE_SPA_WRK)))

     IF ( nstatplnextra .GT. 0 ) THEN
        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), rho, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL DNS_SAVE_EXTRA_I(icount_stat, i1, txc, vaux(vindex(VA_PLANE_SPA_WRK)))

        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), u, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL DNS_SAVE_EXTRA_I(icount_stat, i2, txc, vaux(vindex(VA_PLANE_SPA_WRK)))

        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), v, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL DNS_SAVE_EXTRA_I(icount_stat, i3, txc, vaux(vindex(VA_PLANE_SPA_WRK)))

        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), w, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL DNS_SAVE_EXTRA_I(icount_stat, i4, txc, vaux(vindex(VA_PLANE_SPA_WRK)))

        CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), p, txc(1,1), wrk3d, wrk2d,wrk3d)
        CALL DNS_SAVE_EXTRA_I(icount_stat, i5, txc, vaux(vindex(VA_PLANE_SPA_WRK)))

        DO is = 1,inb_scal
           ie = is + 5
           CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), z1(1,is), txc(1,1), wrk3d, wrk2d,wrk3d)
           CALL DNS_SAVE_EXTRA_I(icount_stat, ie, txc, vaux(vindex(VA_PLANE_SPA_WRK)))
        ENDDO
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE DNS_SAVE
