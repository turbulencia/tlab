SUBROUTINE DNS_SAVE_AVGKIN(rho, u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, mean1d, wrk2d)

! ##############################################
! # Running statistics for kinetic energy before
! # filtering (Favre averages)
! #
! # 10/12/2000 Juan Pedro Mellado
! ##############################################
  
#include "types.h"
#include "avgij_map.h"

  USE DNS_GLOBAL

  IMPLICIT NONE

  TREAL, DIMENSION(imax,jmax,kmax) :: rho, u, v, w
  TREAL, DIMENSION(imax,jmax,kmax) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7

  TREAL mean1d(nstatavg,jmax,*)
  TREAL wrk2d(isize_wrk2d,*)

  TINTEGER j
  TINTEGER NNstat

  NNstat = jmax*nstatavg

! #################################
! # Temporary array storage
! #
! # tmp1 = rho
! # tmp2 = rho*u*u 
! # tmp3 = rho*v*v
! # tmp4 = rho*w*w
! # tmp5 = rho*u
! # tmp6 = rho*v
! # tmp7 = rho*w
! #################################

  CALL REDUCE( imax, jmax, kmax, rho, nstatavg, statavg, tmp1 )
  CALL REDUCE( imax, jmax, kmax, u,   nstatavg, statavg, tmp2 )
  CALL REDUCE( imax, jmax, kmax, v,   nstatavg, statavg, tmp3 )
  CALL REDUCE( imax, jmax, kmax, w,   nstatavg, statavg, tmp4 )

  DO j = 1,NNstat*kmax
     tmp5(j,1,1) = tmp1(j,1,1)*tmp2(j,1,1)
     tmp6(j,1,1) = tmp1(j,1,1)*tmp3(j,1,1)
     tmp7(j,1,1) = tmp1(j,1,1)*tmp4(j,1,1)
     tmp2(j,1,1) = tmp5(j,1,1)*tmp2(j,1,1)
     tmp3(j,1,1) = tmp6(j,1,1)*tmp3(j,1,1)
     tmp4(j,1,1) = tmp7(j,1,1)*tmp4(j,1,1)
  ENDDO

  CALL SUM1V1D_V( NNstat, kmax, tmp5, wrk2d(1,1), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp6, wrk2d(1,2), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp7, wrk2d(1,3), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp2, wrk2d(1,4), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp3, wrk2d(1,5), wrk2d(1,11) )
  CALL SUM1V1D_V( NNstat, kmax, tmp4, wrk2d(1,6), wrk2d(1,11) )

  DO j = 1,NNstat
     MA_FLT_RU(j) = MA_FLT_RU(j) + wrk2d(j,1)
     MA_FLT_RV(j) = MA_FLT_RV(j) + wrk2d(j,2)
     MA_FLT_RW(j) = MA_FLT_RW(j) + wrk2d(j,3)
     MA_FLT_RUU(j) = MA_FLT_RUU(j) + wrk2d(j,4)
     MA_FLT_RVV(j) = MA_FLT_RVV(j) + wrk2d(j,5)
     MA_FLT_RWW(j) = MA_FLT_RWW(j) + wrk2d(j,6)
  ENDDO

  RETURN
END SUBROUTINE DNS_SAVE_AVGKIN
