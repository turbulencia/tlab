#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

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
SUBROUTINE BOUNDARY_BCS_TRANSVERSE_X(u, v, w, p, r, gamma, z1, &
     tmin, mmin, tmax, mmax, tmp1, ddy, ddz, wrk2d, wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE TLAB_VARS,    ONLY : g
  USE TLAB_VARS,    ONLY : imax,jmax,kmax, inb_flow, inb_scal_array
  USE TLAB_VARS,    ONLY : buoyancy
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax)                         :: u, v, w, p, r, gamma
#ifdef USE_MPI
  TREAL, DIMENSION(ims_bcs_imax,               jmax,kmax)  :: tmp1, ddy, ddz
#else
  TREAL, DIMENSION(2*(inb_flow+inb_scal_array),jmax,kmax)  :: tmp1, ddy, ddz
#endif
  TREAL, DIMENSION(imax,jmax,kmax,*)                       :: z1
  TREAL, DIMENSION(jmax,kmax,*)                            :: tmin, tmax, mmin, mmax

  TREAL, DIMENSION(*) :: wrk2d, wrk3d

! -----------------------------------------------------------------------
  TINTEGER ip, j, k, is, bcs(2,2)
  TREAL c

! #######################################################################
  bcs = 0

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
  CALL OPR_PARTIAL_Y(OPR_P1,     ims_bcs_imax,jmax,kmax, bcs, g(2), tmp1, ddy, wrk3d, wrk2d,wrk3d)
! Needs to be checked
  CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_BCS_TRANSVERSE_X. To be checked')
  CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
!  imode_fdm_loc = imode_fdm + (DNS_MPI_K_NRBCX-1)*100
  CALL OPR_PARTIAL_Z(OPR_P1_BCS, ims_bcs_imax,jmax,kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d,wrk3d)
#else
  CALL OPR_PARTIAL_Y(OPR_P1, ip,jmax,kmax, bcs, g(2), tmp1, ddy, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, ip,jmax,kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d,wrk3d)
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
