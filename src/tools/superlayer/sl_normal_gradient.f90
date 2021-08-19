#include "types.h"
#include "dns_const.h"

SUBROUTINE SL_NORMAL_GRADIENT(isl, nmax, istep, kstep, ibuffer_npy, &
     u, v, w, z1, a, sl, profiles, txc, wrk1d, wrk2d, wrk3d)
  
  USE TLAB_VARS

  IMPLICIT NONE

#include "integers.h"

#define L_NFIELDS_MAX 1

  TINTEGER isl, nmax, istep, kstep, ibuffer_npy
  TREAL u(*), v(*), w(*), z1(*), a(*), sl(*)
  TREAL profiles(L_NFIELDS_MAX,nmax,imax/istep,kmax/kstep)
  TREAL txc(imax*jmax*kmax,*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TREAL vmin, vmax
  TINTEGER ij, i, k, n, ifield, nfield, jmin_loc, jmax_loc
  CHARACTER*32 fname

! ###################################################################
  nfield = L_NFIELDS_MAX
  jmin_loc = MAX(1,2*ibuffer_npy)
  jmax_loc = MIN(jmax,jmax - 2*ibuffer_npy +1)

! Calculate scalar gradient field sqrt(G_iG_i), put it in array z1
  CALL FI_GRADIENT(imax,jmax,kmax, z1,a, txc(1,1), wrk2d,wrk3d)
  DO ij = 1,imax*jmax*kmax
     a(ij) = SQRT(a(ij))
  ENDDO

! Calculate boundaries, upper or lower depending on flag isl
  CALL MINMAX(imax,jmax,kmax, a, vmin,vmax)
  vmin = vmin + C_1EM2_R*(vmax-vmin)
  IF ( isl .EQ. 1 ) THEN
     CALL SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, vmin, g(2)%nodes, a, txc(1,1), sl, wrk2d)
  ELSE IF ( isl .EQ. 2 ) THEN
     CALL SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, vmin, g(2)%nodes, a, txc(1,1), sl, wrk2d)
  ENDIF

! -------------------------------------------------------------------
! Normal analysis
! -------------------------------------------------------------------
! Calculate gradient of conditioning field; normal stored in u,v,w
  CALL OPR_PARTIAL_X(OPR_P1, imax,jmax,kmax, bcs, g(1), a, u, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Y(OPR_P1, imax,jmax,kmax, bcs, g(2), a, v, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,jmax,kmax, bcs, g(3), a, w, wrk3d, wrk2d,wrk3d)

! -------------------------------------------------------------------
! TkStat file
! -------------------------------------------------------------------
  WRITE(fname,*) itime; fname = 'slg'//TRIM(ADJUSTL(fname))

  OPEN(unit=21,file=fname)
  WRITE(21,'(A8,E14.7E3)') 'RTIME = ', rtime
  WRITE(21,*) 'I J N G'

  DO k = 1,kmax/kstep
     DO i = 1,imax/istep
        DO n = 1,nmax
           WRITE(21,1020) i, k, M_REAL(n-1-nmax/2)*(g(1)%nodes(2)-g(1)%nodes(1)), &
                (profiles(ifield,n,i,k),ifield=1,nfield)
        ENDDO
     ENDDO
  ENDDO
1020 FORMAT(I3,1X,I3,1X,E10.3E3,L_NFIELDS_MAX(1X,E10.3E3))

  CLOSE(21)

  RETURN
END SUBROUTINE SL_NORMAL_GRADIENT
