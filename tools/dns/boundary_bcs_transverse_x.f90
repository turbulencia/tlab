!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/11/23 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate transverse terms at Ox_min (array tmin) and at Ox_max (array tmax)
!# according to Lodato et al, JCP 227 (2008), 5105-5143
!# The sign is the opposite to that paper
!#
!########################################################################
!# ARGUMENTS 
!#
!# tmin    In    Transverse term at OxMin
!# tmax    In    Transverse term at OxMax
!#
!########################################################################
#include "types.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

SUBROUTINE BOUNDARY_BCS_TRANSVERSE_X(dx, dy, dz, u, v, w, p, r, gamma, z1, &
     tmin, mmin, tmax, mmax, tmp1, ddy, ddz, wrk1d, wrk2d, wrk3d)

  USE DNS_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif
 
  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(*)                                      :: dx, dy, dz
  TREAL, DIMENSION(imax,jmax,kmax)                         :: u, v, w, p, r, gamma
#ifdef USE_MPI
  TREAL, DIMENSION(ims_bcs_imax,               jmax,kmax)  :: tmp1, ddy, ddz
#else
  TREAL, DIMENSION(2*(inb_flow+inb_scal_array),jmax,kmax)  :: tmp1, ddy, ddz
#endif
  TREAL, DIMENSION(imax,jmax,kmax,*)                       :: z1
  TREAL, DIMENSION(jmax,kmax,*)                            :: tmin, tmax, mmin, mmax

  TREAL, DIMENSION(*) :: wrk1d, wrk2d, wrk3d

! -----------------------------------------------------------------------
  TINTEGER ip, j, k, is
  TREAL c
#ifdef USE_MPI
  TINTEGER imode_fdm_loc
#endif

! #######################################################################
! -------------------------------------------------------------------
! Arrange data
! -------------------------------------------------------------------
  ip = 0

! BCs at x_min
  DO k = 1,kmax
     DO j = 1,jmax
        tmp1(ip+1,j,k) = u(1,j,k)
        tmp1(ip+2,j,k) = v(1,j,k)
        tmp1(ip+3,j,k) = w(1,j,k)
        tmp1(ip+4,j,k) = p(1,j,k)
        tmp1(ip+5,j,k) = r(1,j,k)
        DO is = 1,inb_scal_array
           tmp1(ip+5+is,j,k) = z1(1,j,k,is)
        ENDDO
     ENDDO
  ENDDO
  ip = ip + inb_flow + inb_scal_array

! BCs at x_max
  DO k = 1,kmax
     DO j = 1,jmax
        tmp1(ip+1,j,k) = u(imax,j,k)
        tmp1(ip+2,j,k) = v(imax,j,k)
        tmp1(ip+3,j,k) = w(imax,j,k)
        tmp1(ip+4,j,k) = p(imax,j,k)
        tmp1(ip+5,j,k) = r(imax,j,k)
        DO is = 1,inb_scal_array
           tmp1(ip+5+is,j,k) = z1(imax,j,k,is)
        ENDDO
     ENDDO
  ENDDO
  ip = ip + inb_flow + inb_scal_array
  
! -------------------------------------------------------------------
! Construct t1-t5
! -------------------------------------------------------------------
#ifdef USE_MPI
  CALL PARTIAL_Y(imode_fdm, ims_bcs_imax, jmax, kmax, j1bc,&
       dy, tmp1, ddy, i0, i0, wrk1d, wrk2d, wrk3d)
  imode_fdm_loc = imode_fdm + (DNS_MPI_K_NRBCX-1)*100
  CALL PARTIAL_Z(imode_fdm_loc, ims_bcs_imax, jmax, kmax, k1bc,&
       dz, tmp1, ddz, i0, i0, wrk1d, wrk2d, wrk3d)
#else
  CALL PARTIAL_Y(imode_fdm, ip, jmax, kmax, j1bc,&
       dy, tmp1, ddy, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, ip, jmax, kmax, k1bc,&
       dz, tmp1, ddz, i0, i0, wrk1d, wrk2d, wrk3d)
#endif
  ip = 0

! BCs at x_min
  DO k = 1,kmax
     DO j = 1,jmax
        tmin(j,k,1) = r(1,j,k)*ddy(ip+2,j,k) + v(1,j,k)*ddy(ip+5,j,k) &
                    + r(1,j,k)*ddz(ip+3,j,k) + w(1,j,k)*ddz(ip+5,j,k)
        tmin(j,k,2) = v(1,j,k)*ddy(ip+1,j,k) + w(1,j,k)*ddz(ip+1,j,k)
        tmin(j,k,3) = v(1,j,k)*ddy(ip+2,j,k) + w(1,j,k)*ddz(ip+2,j,k) &
                    + ddy(ip+4,j,k)/r(1,j,k) - buoyancy%vector(2)
        tmin(j,k,4) = v(1,j,k)*ddy(ip+3,j,k) + w(1,j,k)*ddz(ip+3,j,k) &
                    + ddz(ip+4,j,k)/r(1,j,k) - buoyancy%vector(3)
        tmin(j,k,5) = v(1,j,k)*ddy(ip+4,j,k) + w(1,j,k)*ddz(ip+4,j,k) &
                    + gamma(1,j,k)*p(1,j,k)*( ddy(ip+2,j,k) + ddz(ip+3,j,k) )
        DO is = 1,inb_scal_array
           tmin(j,k,5+is) = v(1,j,k)*ddy(ip+5+is,j,k) + w(1,j,k)*ddz(ip+5+is,j,k)
        ENDDO
! constructing M1-M5
        c = SQRT(gamma(1,j,k)*p(1,j,k)/r(1,j,k))
        mmin(j,k,1) = (v(1,j,k)-c)*( ddy(ip+4,j,k)     - ddy(ip+2,j,k)*r(1,j,k)*c )
        mmin(j,k,2) =  v(1,j,k)   *( ddy(ip+1,j,k) )
        mmin(j,k,3) =  v(1,j,k)   *( ddy(ip+5,j,k)*c*c - ddy(ip+4,j,k) )
        mmin(j,k,4) =  v(1,j,k)   *( ddy(ip+3,j,k) )
        mmin(j,k,5) = (v(1,j,k)+c)*( ddy(ip+4,j,k)     + ddy(ip+2,j,k)*r(1,j,k)*c )
        DO is = 1,inb_scal_array
           mmin(j,k,5+is) = v(1,j,k)*ddy(ip+5+is,j,k)
        ENDDO
     ENDDO
  ENDDO
  ip = ip + inb_flow + inb_scal_array
     
! BCs at x_max
  DO k = 1,kmax
     DO j = 1,jmax
        tmax(j,k,1) = r(imax,j,k)*ddy(ip+2,j,k) + v(imax,j,k)*ddy(ip+5,j,k) &
                    + r(imax,j,k)*ddz(ip+3,j,k) + w(imax,j,k)*ddz(ip+5,j,k)
        tmax(j,k,2) = v(imax,j,k)*ddy(ip+1,j,k) + w(imax,j,k)*ddz(ip+1,j,k)
        tmax(j,k,3) = v(imax,j,k)*ddy(ip+2,j,k) + w(imax,j,k)*ddz(ip+2,j,k) &
                    + ddy(ip+4,j,k)/r(imax,j,k) - buoyancy%vector(2)
        tmax(j,k,4) = v(imax,j,k)*ddy(ip+3,j,k) + w(imax,j,k)*ddz(ip+3,j,k) &
                    + ddz(ip+4,j,k)/r(imax,j,k) - buoyancy%vector(3)
        tmax(j,k,5) = v(imax,j,k)*ddy(ip+4,j,k) + w(imax,j,k)*ddz(ip+4,j,k) &
                    + gamma(imax,j,k)*p(imax,j,k)*( ddy(ip+2,j,k) + ddz(ip+3,j,k) )
        DO is = 1,inb_scal_array
           tmax(j,k,5+is) = v(imax,j,k)*ddy(ip+5+is,j,k) + w(imax,j,k)*ddz(ip+5+is,j,k)
        ENDDO
! constructing M1-M5
        c = SQRT(gamma(imax,j,k)*p(imax,j,k)/r(imax,j,k))
        mmax(j,k,1) = (v(imax,j,k)-c)*( ddy(ip+4,j,k)     - ddy(ip+2,j,k)*r(imax,j,k)*c)
        mmax(j,k,2) =  v(imax,j,k)   *( ddy(ip+1,j,k) )
        mmax(j,k,3) =  v(imax,j,k)   *( ddy(ip+5,j,k)*c*c - ddy(ip+4,j,k) )
        mmax(j,k,4) =  v(imax,j,k)   *( ddy(ip+3,j,k) )
        mmax(j,k,5) = (v(imax,j,k)+c)*( ddy(ip+4,j,k)     + ddy(ip+2,j,k)*r(imax,j,k)*c)
        DO is = 1,inb_scal_array
           mmax(j,k,5+is) = v(imax,j,k)*ddy(ip+5+is,j,k)
        ENDDO
     ENDDO
  ENDDO
  ip = ip + inb_flow + inb_scal_array 
  
! -------------------------------------------------------------------
! Change sign 
! -------------------------------------------------------------------
  DO is = 1, inb_flow + inb_scal_array
     DO k = 1,kmax
        DO j = 1,jmax
           tmin(j,k,is) =-tmin(j,k,is)
           tmax(j,k,is) =-tmax(j,k,is)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE BOUNDARY_BCS_TRANSVERSE_X
