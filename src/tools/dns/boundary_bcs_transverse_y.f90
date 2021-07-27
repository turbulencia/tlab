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
!# Calculate transverse terms at Oy_min (array tmin) and at Oy_max (array tmax)
!# according to Lodato et al, JCP 227 (2008), 5105-5143.
!# The sign is the opposite to that paper
!#
!########################################################################
!# ARGUMENTS
!#
!# tmin    In    Transverse term at OyMin
!# tmin    In    Transverse term at OyMin
!#
!########################################################################
SUBROUTINE BOUNDARY_BCS_TRANSVERSE_Y(u,v,w,p,r,gamma,z1, &
     tmin,lmin,tmax,lmax, tmp1,ddx,ddz, wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : g
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow, inb_scal_array
  USE DNS_GLOBAL,    ONLY : buoyancy
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"

  TREAL, DIMENSION(imax,jmax,kmax)                        :: u, v, w, p, r, gamma
#ifdef USE_MPI
  TREAL, DIMENSION(imax,ims_bcs_jmax,              kmax)  :: tmp1, ddx, ddz
#else
  TREAL, DIMENSION(imax,2*(inb_flow+inb_scal_array),kmax) :: tmp1, ddx, ddz
#endif
  TREAL, DIMENSION(imax,jmax,kmax,*)                      :: z1
  TREAL, DIMENSION(imax,kmax,*)                           :: tmin, lmin, tmax, lmax

  TREAL, DIMENSION(*)              :: wrk2d, wrk3d

! -----------------------------------------------------------------------
  TINTEGER ip, i, k, is, bcs(2,2)
  TREAL c

! #######################################################################
  bcs = 0

! -------------------------------------------------------------------
! Arrange data
! -------------------------------------------------------------------
  ip = 0

! BCs at y_min
  DO k = 1,kmax; DO i = 1,imax
     tmp1(i,ip+1,k) = u(i,1,k)
     tmp1(i,ip+2,k) = v(i,1,k)
     tmp1(i,ip+3,k) = w(i,1,k)
     tmp1(i,ip+4,k) = p(i,1,k)
     tmp1(i,ip+5,k) = r(i,1,k)
     DO is = 1,inb_scal_array
        tmp1(i,ip+5+is,k) = z1(i,1,k,is)
     ENDDO
  ENDDO; ENDDO
  ip = ip + inb_flow + inb_scal_array

! BCs at y_max
  DO k = 1,kmax; DO i = 1,imax
     tmp1(i,ip+1,k) = u(i,jmax,k)
     tmp1(i,ip+2,k) = v(i,jmax,k)
     tmp1(i,ip+3,k) = w(i,jmax,k)
     tmp1(i,ip+4,k) = p(i,jmax,k)
     tmp1(i,ip+5,k) = r(i,jmax,k)
     DO is = 1,inb_scal_array
        tmp1(i,ip+5+is,k) = z1(i,jmax,k,is)
     ENDDO
  ENDDO; ENDDO
  ip = ip + inb_flow + inb_scal_array

! -------------------------------------------------------------------
! Construct t1-t5
! -------------------------------------------------------------------
#ifdef USE_MPI
  CALL OPR_PARTIAL_X(OPR_P1,     imax,ims_bcs_jmax,kmax, bcs, g(1), tmp1, ddx, wrk3d, wrk2d,wrk3d)
! Needs to be checked
  CALL TLAB_WRITE_ASCII(efile,'BOUNDARY_BCS_TRANSVERSE_Y. To be checked')
  CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
!  imode_fdm_loc = imode_fdm + (DNS_MPI_K_NRBCY-1)*100
  CALL OPR_PARTIAL_Z(OPR_P1_BCS, imax,ims_bcs_jmax,kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d,wrk3d)
#else
  CALL OPR_PARTIAL_X(OPR_P1, imax,ip,kmax, bcs, g(1), tmp1, ddx, wrk3d, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1, imax,ip,kmax, bcs, g(3), tmp1, ddz, wrk3d, wrk2d,wrk3d)
#endif
  ip = 0

! BCs at y_min
  DO k = 1,kmax; DO i = 1,imax
     tmin(i,k,1) = r(i,1,k)*ddx(i,ip+1,k) + u(i,1,k)*ddx(i,ip+5,k) &
                 + r(i,1,k)*ddz(i,ip+3,k) + w(i,1,k)*ddz(i,ip+5,k)
     tmin(i,k,2) = u(i,1,k)*ddx(i,ip+1,k) + w(i,1,k)*ddz(i,ip+1,k) &
                 + ddx(i,ip+4,k)/r(i,1,k) - buoyancy%vector(1)
     tmin(i,k,3) = u(i,1,k)*ddx(i,ip+2,k) + w(i,1,k)*ddz(i,ip+2,k)
     tmin(i,k,4) = u(i,1,k)*ddx(i,ip+3,k) + w(i,1,k)*ddz(i,ip+3,k) &
                 + ddz(i,ip+4,k)/r(i,1,k) - buoyancy%vector(3)
     tmin(i,k,5) = u(i,1,k)*ddx(i,ip+4,k) + w(i,1,k)*ddz(i,ip+4,k) &
                 + gamma(i,1,k)*p(i,1,k)*( ddx(i,ip+1,k) + ddz(i,ip+3,k) )
     DO is = 1,inb_scal_array
        tmin(i,k,5+is) = u(i,1,k)*ddx(i,ip+5+is,k) + w(i,1,k)*ddz(i,ip+5+is,k)
     ENDDO
! constructing L1-L5
     c = SQRT(gamma(i,1,k)*p(i,1,k)/r(i,1,k))
     lmin(i,k,1) = (u(i,1,k)-c)*( ddx(i,ip+4,k)     - ddx(i,ip+1,k)*r(i,1,k)*c )
     lmin(i,k,2) =  u(i,1,k)   *( ddx(i,ip+5,k)*c*c - ddx(i,ip+4,k) )
     lmin(i,k,3) =  u(i,1,k)   *( ddx(i,ip+2,k) )
     lmin(i,k,4) =  u(i,1,k)   *( ddx(i,ip+3,k) )
     lmin(i,k,5) = (u(i,1,k)+c)*( ddx(i,ip+4,k)     + ddx(i,ip+1,k)*r(i,1,k)*c )
     DO is = 1,inb_scal_array
        lmin(i,k,5+is) = u(i,1,k)*ddx(i,ip+5+is,k)
     ENDDO
  ENDDO; ENDDO
  ip = ip + inb_flow + inb_scal_array

! BCs at y_max
  DO k = 1,kmax; DO i = 1,imax
     tmax(i,k,1) = r(i,jmax,k)*ddx(i,ip+1,k) + u(i,jmax,k)*ddx(i,ip+5,k) &
                 + r(i,jmax,k)*ddz(i,ip+3,k) + w(i,jmax,k)*ddz(i,ip+5,k)
     tmax(i,k,2) = u(i,jmax,k)*ddx(i,ip+1,k) + w(i,jmax,k)*ddz(i,ip+1,k) &
                 + ddx(i,ip+4,k)/r(i,jmax,k) - buoyancy%vector(1)
     tmax(i,k,3) = u(i,jmax,k)*ddx(i,ip+2,k) + w(i,jmax,k)*ddz(i,ip+2,k)
     tmax(i,k,4) = u(i,jmax,k)*ddx(i,ip+3,k) + w(i,jmax,k)*ddz(i,ip+3,k) &
                 + ddz(i,ip+4,k)/r(i,jmax,k) - buoyancy%vector(3)
     tmax(i,k,5) = u(i,jmax,k)*ddx(i,ip+4,k) + w(i,jmax,k)*ddz(i,ip+4,k) &
                 + gamma(i,jmax,k)*p(i,jmax,k)*( ddx(i,ip+1,k) + ddz(i,ip+3,k) )
     DO is = 1,inb_scal_array
        tmax(i,k,5+is) = u(i,jmax,k)*ddx(i,ip+5+is,k) + w(i,jmax,k)*ddz(i,ip+5+is,k)
     ENDDO
! constructing L1-L5
     c = SQRT(gamma(i,jmax,k)*p(i,jmax,k)/r(i,jmax,k))
     lmax(i,k,1) = (u(i,jmax,k)-c)*( ddx(i,ip+4,k)     - ddx(i,ip+1,k)*r(i,jmax,k)*c)
     lmax(i,k,2) =  u(i,jmax,k)   *( ddx(i,ip+5,k)*c*c - ddx(i,ip+4,k) )
     lmax(i,k,3) =  u(i,jmax,k)   *( ddx(i,ip+2,k) )
     lmax(i,k,4) =  u(i,jmax,k)   *( ddx(i,ip+3,k) )
     lmax(i,k,5) = (u(i,jmax,k)+c)*( ddx(i,ip+4,k)     + ddx(i,ip+1,k)*r(i,jmax,k)*c)
     DO is = 1,inb_scal_array
        lmax(i,k,5+is) = u(i,jmax,k)*ddx(i,ip+5+is,k)
     ENDDO
  ENDDO; ENDDO
  ip = ip + inb_flow + inb_scal_array

! -------------------------------------------------------------------
! Change sign
! -------------------------------------------------------------------
  DO is = 1,inb_flow + inb_scal_array
     DO k = 1,kmax; DO i = 1,imax
        tmin(i,k,is) =-tmin(i,k,is)
        tmax(i,k,is) =-tmax(i,k,is)
     ENDDO; ENDDO
  ENDDO

  RETURN
END SUBROUTINE BOUNDARY_BCS_TRANSVERSE_Y
