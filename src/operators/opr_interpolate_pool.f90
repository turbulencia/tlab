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

  USE TLAB_VARS, ONLY : isize_txc_field
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
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
     id = TLAB_MPI_I_AUX1
     u_tmp2(1:nx*ny*nz) = u_org(1:nx*ny*nz) ! Need additional space for transposition
     CALL TLAB_MPI_TRPF_I(u_tmp2, u_tmp1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

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
     id = TLAB_MPI_I_AUX2
     CALL TLAB_MPI_TRPB_I(u_tmp2, u_tmp1, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     u_dst(1:nx_dst*ny*nz) = u_tmp1(1:nx_dst*ny*nz)
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

  USE TLAB_VARS, ONLY : isize_txc_field
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_PROCS
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
     id = TLAB_MPI_K_AUX1
     CALL TLAB_MPI_TRPF_K(u_org, u_tmp2, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

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
     id = TLAB_MPI_K_AUX2
     CALL TLAB_MPI_TRPB_K(u_tmp1, u_dst, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
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

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_PROCS

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER imax,kmax, imax_dst, isize_wrk, k
  TREAL scalex
  TREAL x_org(imax+1)
  TREAL x_dst(imax_dst)
  TREAL u_org(imax,*)
  TREAL u_dst(imax_dst,*)

  TREAL wrk(isize_wrk), rdum

  TINTEGER, dimension(2) :: CSpline_BCType
  TREAL,    dimension(2) :: CSpline_BCVal

! #######################################################################

  IF ( isize_wrk .LT. 12*imax +1 ) THEN
     CALL TLAB_WRITE_ASCII(efile, 'INTERPOLATE_1D. Temporary Array not large enough')
     CALL TLAB_STOP(DNS_ERROR_CURFIT)
  ENDIF
  !------------------------------------------! the periodic case
     IF ( periodic ) THEN                    !
        Cspline_BCType(1) = CS_BCS_PERIODIC  !
        Cspline_BCType(2) = CSpline_BCType(1)!
        CSpline_BCVal(:) = C_0_R             !
        x_org(imax+1) = x_org(1) + scalex    ! extend x_org for periodic point
        DO k = 1,kmax-1                      !
           rdum=u_org(imax+1,k)              ! use u_org(imax+1,k) to extend array and
           u_org(imax+1,k) = u_org(1,k)      ! avoid the copy of the whole line
           CALL CUBIC_SPLINE(CSpline_BCType,CSpline_BCVal,&
                imax+1,imax_dst,x_org,u_org(1,k),x_dst,u_dst(1,k),&
                wrk(imax+2))                 !
           u_org(imax+1,k) = rdum            ! set u_org back to stored value rdum
        ENDDO                                !
        wrk(1:imax) = u_org(1:imax,kmax)     ! cannot avoid the copy for the last line
        wrk(imax+1) = u_org(1,kmax)          ! as u_org(imax+1,kmax) is out of bounds
        CALL CUBIC_SPLINE(CSpline_BCType,CSpline_BCVal,&
             imax+1,imax_dst,x_org,wrk,x_dst,u_dst(1,kmax),&
             wrk(imax+2))
     !---------------------------------------! the aperiodic case
     ELSE
        CSpline_BCType(1) = CS_BCS_NATURAL; CSpline_BCVal(1) = C_0_R
        CSpline_BCType(2) = CS_BCS_NATURAL; CSpline_BCVal(2) = C_0_R
        DO k=1,kmax
           CALL CUBIC_SPLINE(CSpline_BCType,CSpline_BCVal,&
                imax,imax_dst,x_org,u_org(1,k),x_dst,u_dst(1,k),&
                wrk)
        ENDDO
     ENDIF

  RETURN
END SUBROUTINE INTERPOLATE_1D
