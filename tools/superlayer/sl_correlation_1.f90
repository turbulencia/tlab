#include "types.h"

!########################################################################
!# Tool/Library SUPERLAYER
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE SL_CORRELATION_1(ilog, y, dx, dy, dz, u, v, w, z1, corr, &
     strain, vorticity, gradient, tmp1, tmp2, wrk1d, wrk2d, wrk3d)

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TINTEGER ilog
  TREAL y(jmax)
  TREAL dx(imax)
  TREAL dy(jmax)
  TREAL dz(kmax_total)
  TREAL u(imax,jmax,kmax)
  TREAL v(imax,jmax,kmax)
  TREAL w(imax,jmax,kmax)
  TREAL z1(imax,jmax,kmax)
  TREAL strain(imax,jmax,kmax)
  TREAL vorticity(imax,jmax,kmax)
  TREAL gradient(imax,jmax,kmax)

  TREAL tmp1(imax,jmax,kmax)
  TREAL tmp2(imax,jmax,kmax)

  TREAL corr(jmax,*)

  TREAL wrk1d(jmax,*)
  TREAL wrk2d(imax,kmax,*)
  TREAL wrk3d(imax,jmax,kmax)

! -------------------------------------------------------------------
  TINTEGER j
  TREAL mean_1, mean_2, var_1, var_2, delta_w
  TREAL AVG1V2D, COV2V2D

  CHARACTER*32 fname
  CHARACTER*250 line1
  CHARACTER*550 line2

! ###################################################################
  IF ( delta_u .EQ. C_0_R ) THEN
     delta_w = C_1_R
  ELSE
     DO j = 1,jmax
        wrk1d(j,1) = AVG1V2D(imax, jmax, kmax, j, i1, u)
     ENDDO
     CALL PARTIAL_Y(imode_fdm, i1, jmax, i1, j1bc, dy, wrk1d(1,1), &
          wrk1d(1,2), i0, i0, wrk1d(1,3), wrk2d, wrk3d)
     delta_w = delta_u/MAXVAL(ABS(wrk1d(1:jmax,2)))
  ENDIF

! ###################################################################
! Define fields
! ###################################################################
  CALL FI_STRAIN(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
       dx, dy, dz, u, v, w, strain, tmp1, tmp2, wrk1d, wrk2d, wrk3d)
  CALL FI_VORTICITY(imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
       dx, dy, dz, u, v, w, vorticity, tmp1, tmp2, wrk1d, wrk2d, wrk3d)
  CALL FI_GRADIENT(imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, dx,dy,dz, z1,gradient, tmp1, wrk1d,wrk2d,wrk3d)

  IF ( ilog .EQ. 1 ) THEN
     strain    = log(strain)
     vorticity = log(vorticity)
     gradient  = log(gradient)
  ENDIF

! ###################################################################
! Compute plane correlations
! ###################################################################
  DO j = 1,jmax

     corr(j,1) = COV2V2D(imax, jmax, kmax, j, vorticity, strain)
     mean_1 = AVG1V2D(imax, jmax, kmax, j, i1, vorticity)
     var_1  = AVG1V2D(imax, jmax, kmax, j, i2, vorticity)-mean_1*mean_1
     mean_2 = AVG1V2D(imax, jmax, kmax, j, i1, strain)
     var_2  = AVG1V2D(imax, jmax, kmax, j, i2, strain)-mean_2*mean_2
     IF ( var_1 .GT. C_0_R .AND. var_2 .GT. C_0_R ) THEN
        corr(j,1) = (corr(j,1)-mean_1*mean_2)/SQRT(var_1*var_2)
     ELSE
        corr(j,1) = C_2_R
     ENDIF

     corr(j,2) = COV2V2D(imax, jmax, kmax, j, vorticity, gradient)
     mean_1 = AVG1V2D(imax, jmax, kmax, j, i1, vorticity)
     var_1  = AVG1V2D(imax, jmax, kmax, j, i2, vorticity)-mean_1*mean_1
     mean_2 = AVG1V2D(imax, jmax, kmax, j, i1, gradient)
     var_2  = AVG1V2D(imax, jmax, kmax, j, i2, gradient)-mean_2*mean_2
     IF ( var_1 .GT. C_0_R .AND. var_2 .GT. C_0_R ) THEN
        corr(j,2) = (corr(j,2)-mean_1*mean_2)/SQRT(var_1*var_2)
     ELSE
        corr(j,2) = C_2_R
     ENDIF

     corr(j,3) = COV2V2D(imax, jmax, kmax, j, gradient, strain)
     mean_1 = AVG1V2D(imax, jmax, kmax, j, i1, gradient)
     var_1  = AVG1V2D(imax, jmax, kmax, j, i2, gradient)-mean_1*mean_1
     mean_2 = AVG1V2D(imax, jmax, kmax, j, i1, strain)
     var_2  = AVG1V2D(imax, jmax, kmax, j, i2, strain)-mean_2*mean_2
     IF ( var_1 .GT. C_0_R .AND. var_2 .GT. C_0_R ) THEN
        corr(j,3) = (corr(j,3)-mean_1*mean_2)/SQRT(var_1*var_2)
     ELSE
        corr(j,3) = C_2_R
     ENDIF

  ENDDO

! ###################################################################
! TkStat Output
! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     WRITE(fname,*) itime; fname='slc'//TRIM(ADJUSTL(fname))

     OPEN(UNIT=23, FILE=fname, STATUS='unknown')
     WRITE(23, '(A8,E14.7E3)') 'RTIME = ', rtime
     WRITE(23, '(A8)') 'IMAX = 1'
     WRITE(23, '(A7,I8)') 'JMAX = ', jmax

! -------------------------------------------------------------------
! Header
! -------------------------------------------------------------------
! Independent variables
     line2='I J Y SW'

! Main fields
     line1 = 'W-S W-G S-G'
     WRITE(i23,1010) 'GROUP = MainFields '//TRIM(ADJUSTL(line1))
     line2 = TRIM(ADJUSTL(line2))//' '//TRIM(ADJUSTL(line1))

     WRITE(i23,1010) TRIM(ADJUSTL(line2))

1010 FORMAT(A)

! -------------------------------------------------------------------
! Body
! -------------------------------------------------------------------
     DO j = 1,jmax
        WRITE(23,1020) 1, j, y(j), (y(j)-g(2)%scale *qbg(1)%ymean-y(1))/delta_w,&
             corr(j,1), corr(j,2), corr(j,3)
     ENDDO

1020 FORMAT(I3,1X,I3,2(1X,E12.5E3),3(1X,E12.5E3))

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE SL_CORRELATION_1
