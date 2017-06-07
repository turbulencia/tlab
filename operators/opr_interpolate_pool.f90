#include "types.h"
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
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2011/08/01 - J.P. Mellado
!#              Parallelized in Oz using MPI
!# 2012/11/01 - J.P. Mellado
!#              Parallelized in Ox using MPI
!#
!########################################################################
!# DESCRIPTION
!#
!# All three routines called the same kernel routine INTERPOLATE_1D, 
!# calls in turn the routines from the library spline
!#
!########################################################################

! #######################################################################
! Interpolation in X
! #######################################################################
SUBROUTINE OPR_INTERPOLATE_X(nx,ny,nz, nx_dst, periodic, scalex, &
     x_org,x_dst, u_org,u_dst, u_tmp1,u_tmp2, isize_wrk3d, wrk3d)

  USE DNS_GLOBAL, ONLY : isize_txc_field
#ifdef USE_MPI
  USE DNS_MPI 
#endif

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER nx,ny,nz, nx_dst, isize_wrk3d
  TREAL scalex
  TREAL, DIMENSION(*)                       :: x_org, x_dst
  TREAL, DIMENSION(nx    *ny*nz),    TARGET :: u_org
  TREAL, DIMENSION(nx_dst*ny*nz),    TARGET :: u_dst
  TREAL, DIMENSION(isize_txc_field), TARGET :: u_tmp1, u_tmp2
  TREAL, DIMENSION(isize_wrk3d)             :: wrk3d

! -----------------------------------------------------------------------
  TINTEGER nyz, nx_total, nx_total_dst
#ifdef USE_MPI
  TINTEGER id
#endif

  TREAL, DIMENSION(:), POINTER :: p_a, p_b
 
! #######################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
!     id = DNS_MPI_I_PARTIAL
     id = DNS_MPI_I_AUX1
     CALL DNS_MPI_TRPF_I(u_org, u_tmp1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

     p_a => u_tmp1
     p_b => u_tmp2

     nyz          = ims_size_i(id)
     nx_total     = nx    *ims_npro_i
     nx_total_dst = nx_dst*ims_npro_i

  ELSE
#endif
     p_a => u_org
     p_b => u_dst

     nyz          = ny*nz
     nx_total     = nx
     nx_total_dst = nx_dst

#ifdef USE_MPI
  ENDIF
#endif

! -----------------------------------------------------------------------
  CALL INTERPOLATE_1D(periodic, nx_total,nyz, nx_total_dst, scalex, x_org,x_dst, p_a,p_b, isize_wrk3d, wrk3d)
 
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_i .GT. 1 ) THEN
     id = DNS_MPI_I_AUX2
     CALL DNS_MPI_TRPB_I(u_tmp2, u_dst, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
  ENDIF
#endif
  NULLIFY(p_a,p_b)

  RETURN
END SUBROUTINE OPR_INTERPOLATE_X

! ###################################################################
! Interpolation in Oz direction
! ###################################################################
SUBROUTINE OPR_INTERPOLATE_Z(nx,ny,nz, nz_dst, periodic, scalez, &
     z_org,z_dst, u_org,u_dst, u_tmp1,u_tmp2, isize_wrk3d, wrk3d)

  USE DNS_GLOBAL, ONLY : isize_txc_field
#ifdef USE_MPI
  USE DNS_MPI 
#endif

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER nx,ny,nz, nz_dst, isize_wrk3d
  TREAL scalez
  TREAL, DIMENSION(*)                       :: z_org, z_dst
  TREAL, DIMENSION(nx*ny*nz    ),    TARGET :: u_org
  TREAL, DIMENSION(nx*ny*nz_dst),    TARGET :: u_dst
  TREAL, DIMENSION(isize_txc_field), TARGET :: u_tmp1, u_tmp2
  TREAL, DIMENSION(isize_wrk3d)             :: wrk3d

! -----------------------------------------------------------------------
  TINTEGER nxy, nz_total, nz_total_dst
#ifdef USE_MPI
  TINTEGER id
#endif

  TREAL, DIMENSION(:), POINTER :: p_a, p_b
 
! #######################################################################
! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     id = DNS_MPI_K_AUX1
     CALL DNS_MPI_TRPF_K(u_org, u_tmp2, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

     p_a => u_tmp2
     p_b => u_tmp1

     nxy          = ims_size_k(id)
     nz_total     = nz    *ims_npro_k
     nz_total_dst = nz_dst*ims_npro_k

  ELSE
#endif
     p_a => u_org
     p_b => u_dst

     nxy          = nx*ny
     nz_total     = nz
     nz_total_dst = nz_dst
#ifdef USE_MPI
  ENDIF
#endif

! -------------------------------------------------------------------
! Make z direction the first one 
! -------------------------------------------------------------------
#ifdef USE_ESSL
     CALL DGETMO(p_a, nxy, nxy, nz_total, u_tmp1, nz_total)
#else
     CALL DNS_TRANSPOSE(p_a, nxy, nz_total, nxy, u_tmp1, nz_total)
#endif

! -----------------------------------------------------------------------
  CALL INTERPOLATE_1D(periodic, nz_total,nxy, nz_total_dst, scalez, z_org,z_dst, u_tmp1,u_tmp2, isize_wrk3d, wrk3d)

! -------------------------------------------------------------------
! Put arrays back in the right order
! -------------------------------------------------------------------
#ifdef USE_ESSL
     CALL DGETMO(u_tmp2, nz_total_dst, nz_total_dst, nxy, p_b, nxy)
#else
     CALL DNS_TRANSPOSE(u_tmp2, nz_total_dst, nxy, nz_total_dst, p_b, nxy)
#endif

! -------------------------------------------------------------------
! Transposition
! -------------------------------------------------------------------
#ifdef USE_MPI         
  IF ( ims_npro_k .GT. 1 ) THEN
     id = DNS_MPI_K_AUX2
     CALL DNS_MPI_TRPB_K(u_tmp1, u_dst, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif
  NULLIFY(p_a,p_b)

  RETURN
END SUBROUTINE OPR_INTERPOLATE_Z

! #######################################################################
! Interpolation in Y
! #######################################################################
SUBROUTINE OPR_INTERPOLATE_Y(nx,ny,nz, ny_dst, periodic, scaley, &
     y_org,y_dst, u_org,u_dst, u_tmp1,u_tmp2, isize_wrk3d, wrk3d)

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER nx,ny,nz, ny_dst, isize_wrk3d
  TREAL scaley
  TREAL y_org(ny+1)
  TREAL y_dst(ny_dst)
  TREAL, DIMENSION(nx,ny,    nz) :: u_org, u_tmp1
  TREAL, DIMENSION(nx,ny_dst,nz) :: u_dst, u_tmp2
  TREAL, DIMENSION(isize_wrk3d)  :: wrk3d

! -----------------------------------------------------------------------
  TINTEGER ikmax, nyz

! #######################################################################
! -------------------------------------------------------------------
! Make y direction the first one (x direction the last one)
! -------------------------------------------------------------------
  nyz = ny*nz
#ifdef USE_ESSL
  CALL DGETMO(u_org, nx, nx, nyz,        u_tmp1, nyz)
#else
  CALL DNS_TRANSPOSE(u_org, nx, nyz,        nx, u_tmp1, nyz)
#endif

! -----------------------------------------------------------------------
  ikmax = nx*nz
  CALL INTERPOLATE_1D(periodic, ny,ikmax, ny_dst, scaley, y_org,y_dst, u_tmp1,u_tmp2, isize_wrk3d, wrk3d)

! -------------------------------------------------------------------
! Put arrays back in the right order
! -------------------------------------------------------------------
  nyz = ny_dst*nz
#ifdef USE_ESSL
  CALL DGETMO(u_tmp2, nyz, nyz,        nx, u_dst, nx)
#else
  CALL DNS_TRANSPOSE(u_tmp2, nyz, nx, nyz,        u_dst, nx)
#endif

  RETURN
END SUBROUTINE OPR_INTERPOLATE_Y

! #######################################################################
! #######################################################################

! #######################################################################
! Interpolation in 1D
! #######################################################################
SUBROUTINE INTERPOLATE_1D(periodic, imax,kmax, imax_dst, scalex, x_org,x_dst, u_org,u_dst, isize_wrk, wrk)

  USE DNS_CONSTANTS, ONLY : efile

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER imax,kmax, imax_dst, isize_wrk
  TREAL scalex
  TREAL x_org(imax+1)
  TREAL x_dst(imax_dst)
  TREAL u_org(imax,*)
  TREAL u_dst(imax_dst,*)

  TREAL wrk(isize_wrk)

! -----------------------------------------------------------------------
  TINTEGER i,k
  TINTEGER iopt, kx, nest, lwrk, nx 
  TINTEGER ip1, ip2, ip3, ip4, ip5, ip6, ier, imax1

  TREAL xb, xe, s, fp

  CHARACTER*64 line

! #######################################################################
  iopt = 0
  xb = x_org(1); xe = x_org(imax)                  ! just used in the non-periodic case
  imax1 = imax+1; x_org(imax1) = x_org(1) + scalex ! the periodic case
  kx = 5                                           ! order of the splines interpolation
  s = C_0_R
  nest = imax1 + 2*kx                              ! the periodic case requires maximum sizes
  lwrk = imax1*(kx+1)+nest*(8+5*kx)

! Working array relative possitions
  ip1 = 1           ! w
  ip2 = ip1 + imax1 ! t
  ip3 = ip2 + nest  ! c
  ip4 = ip3 + nest  ! wrk
  ip5 = ip4 + lwrk  ! iwkr
  ip6 = ip5 + nest  ! returned value

  IF ( isize_wrk .LT. ip6 ) THEN
     CALL IO_WRITE_ASCII(efile, 'INTERPOLATE_1D. Temporary Array not large enough')
     CALL DNS_STOP(DNS_ERROR_CURFIT)
  ENDIF

! #######################################################################
  DO i=1,imax
     wrk(ip1+i-1) = C_1_R
  ENDDO

  DO k = 1,kmax
     IF ( periodic ) THEN
        CALL percur(iopt, imax1, x_org, u_org(1,k), wrk(ip1),        kx, s, &
             nest, nx, wrk(ip2), wrk(ip3), fp, wrk(ip4), lwrk, wrk(ip5), ier)

        IF ( ier .NE. 0 .AND. ier .NE. -1 ) THEN
           WRITE(line, *) 'INTERPOLATE_1D. Percur error code = ', ier
           CALL IO_WRITE_ASCII(efile, line)
           CALL DNS_STOP(DNS_ERROR_CURFIT)
        ENDIF
     ELSE
        CALL curfit(iopt, imax,  x_org, u_org(1,k), wrk(ip1), xb,xe, kx, s, &
             nest, nx, wrk(ip2), wrk(ip3), fp, wrk(ip4), lwrk, wrk(ip5), ier)

        IF ( ier .NE. 0 .AND. ier .NE. -1 ) THEN
           WRITE(line, *) 'INTERPOLATE_1D. Curfit error code = ', ier
           CALL IO_WRITE_ASCII(efile, line)
           CALL DNS_STOP(DNS_ERROR_CURFIT)
        ENDIF
     ENDIF


     CALL splev(wrk(ip2),nx,wrk(ip3),kx,x_dst,u_dst(1,k), imax_dst, ier)

  ENDDO

  RETURN
END SUBROUTINE INTERPOLATE_1D
