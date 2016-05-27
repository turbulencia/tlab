#include "types.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculates turbulent dissipation per unit volume \rho \epsilon = \tau'_{ij}u'_{ij}
!# It assumes constant visocsity
!#
!########################################################################
!# ARGUMENTS 
!#
!# ifluc    In     Flag to consider total or fluctuation field:
!#                 0 => tau_ji * u_i,j
!#                 1 => tau'_ij * u'_i,j 
!#
!########################################################################
SUBROUTINE FI_DISSIPATION(ifluc, imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
     area, visc, dx, dy, dz, u, v, w, eps, &
     tmp1, tmp2, tmp3, mean2d, wrk1d, wrk2d, wrk3d)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER imode_fdm, ifluc
  TINTEGER imax, jmax, kmax
  TINTEGER i1bc, j1bc, k1bc
  TREAL area, visc

  TREAL, DIMENSION(*)              :: dx, dy, dz
  TREAL, DIMENSION(imax,jmax,kmax) :: u, v, w, eps
  TREAL, DIMENSION(imax,jmax,kmax) :: tmp1, tmp2, tmp3, wrk3d
  TREAL, DIMENSION(*)              :: wrk1d, wrk2d

  TREAL mean2d(jmax,5)

! -------------------------------------------------------------------
  TREAL c23, dil
  TREAL tau11, tau22, tau33, tau12, tau13, tau23
  TREAL upy, vpy, wpy
  TREAL AVG_IK
  TINTEGER i, j, k

#define rU(j)     mean2d(j,1)
#define rV(j)     mean2d(j,1)
#define rW(j)     mean2d(j,1)

#define rU_y(j)   mean2d(j,2)
#define rV_y(j)   mean2d(j,2)
#define rW_y(j)   mean2d(j,2)

#define Tau_xx(j) mean2d(j,3)
#define Tau_yy(j) mean2d(j,4)
#define Tau_zz(j) mean2d(j,5)
#define Tau_xy(j) mean2d(j,3)
#define Tau_xz(j) mean2d(j,3)
#define Tau_yz(j) mean2d(j,3)

! ###################################################################
  c23 = C_2_R/C_3_R

! -------------------------------------------------------------------
! Diagonal terms
! -------------------------------------------------------------------
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, w, tmp3, i0, i0, wrk1d, wrk2d, wrk3d)

! calculate mean profiles
  IF ( ifluc .EQ. 1 ) THEN
     DO j=1,jmax
        DO k=1,kmax
           DO i=1,imax
!              vs = vis(i,j,k)
              dil = (tmp1(i,j,k)+tmp2(i,j,k)+tmp3(i,j,k))*c23                  
!              wrk3d(i,1,k) = vs*(C_2_R*tmp1(i,j,k)-dil)
!              wrk3d(i,2,k) = vs*(C_2_R*tmp2(i,j,k)-dil)
!              wrk3d(i,3,k) = vs*(C_2_R*tmp3(i,j,k)-dil)
              wrk3d(i,1,k) = C_2_R*tmp1(i,j,k)-dil
              wrk3d(i,2,k) = C_2_R*tmp2(i,j,k)-dil
              wrk3d(i,3,k) = C_2_R*tmp3(i,j,k)-dil
           ENDDO
        ENDDO
        Tau_xx(j) = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        Tau_yy(j) = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)
        Tau_zz(j) = AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)
     ENDDO

     DO j=1,jmax
        rV(j) = AVG_IK(imax, jmax, kmax, j, v, dx, dz, area)
     ENDDO
     CALL PARTIAL_Y(imode_fdm, i1, jmax, i1, j1bc, dy, rV(1), &
          rV_y(1), i0, i0, wrk1d, wrk2d, wrk3d)
  ELSE
     DO j = 1, jmax
        Tau_xx(j) = C_0_R
        Tau_yy(j) = C_0_R
        Tau_zz(j) = C_0_R

        rV_y(j) = C_0_R
     ENDDO
  ENDIF

! calculate dissipation
  DO k=1, kmax
     DO j=1, jmax
        DO i=1, imax
!           vs = vis(i,j,k)
           dil = (tmp1(i,j,k)+tmp2(i,j,k)+tmp3(i,j,k))*c23
!           tau11 = vs*(C_2_R*tmp1(i,j,k)-dil)-Tau_xx(j)
!           tau22 = vs*(C_2_R*tmp2(i,j,k)-dil)-Tau_yy(j)
!           tau33 = vs*(C_2_R*tmp3(i,j,k)-dil)-Tau_zz(j)
           tau11 = (C_2_R*tmp1(i,j,k)-dil)-Tau_xx(j)
           tau22 = (C_2_R*tmp2(i,j,k)-dil)-Tau_yy(j)
           tau33 = (C_2_R*tmp3(i,j,k)-dil)-Tau_zz(j)
           vpy = tmp2(i,j,k)-rV_y(j)

           eps(i,j,k) = tau11*tmp1(i,j,k)+tau22*vpy+tau33*tmp3(i,j,k)
        ENDDO
     ENDDO
  ENDDO

! -------------------------------------------------------------------
! 12 contribution
! -------------------------------------------------------------------
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, v, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)

! calculate mean profiles
  IF ( ifluc .EQ. 1 ) THEN
     DO j=1,jmax
        DO k=1,kmax
           DO i=1,imax
!              wrk3d(i,4,k) = vis(i,j,k)*(tmp1(i,j,k)+tmp2(i,j,k))
              wrk3d(i,4,k) = tmp1(i,j,k)+tmp2(i,j,k)
           ENDDO
        ENDDO
        Tau_xy(j) = AVG_IK(imax, jmax, kmax, i4, wrk3d, dx, dz, area)

     ENDDO

     DO j=1,jmax
        rU(j) = AVG_IK(imax, jmax, kmax, j, u, dx, dz, area)
     ENDDO
     CALL PARTIAL_Y(imode_fdm, i1, jmax, i1, j1bc, dy, rU(1), &
          rU_y(1), i0, i0, wrk1d, wrk2d, wrk3d)
  ELSE
     DO j = 1,jmax
        Tau_xy(j) = C_0_R
        rU_y(j) = C_0_R
     ENDDO
  ENDIF

! calculate dissipation
  DO k=1, kmax
     DO j=1, jmax
        DO i=1, imax
!           tau12 = vis(i,j,k)*(tmp1(i,j,k)+tmp2(i,j,k))-Tau_xy(j)
           tau12 = (tmp1(i,j,k)+tmp2(i,j,k))-Tau_xy(j)
           upy = tmp1(i,j,k)-rU_y(j)

           eps(i,j,k) = eps(i,j,k) + tau12*(upy+tmp2(i,j,k))
        ENDDO
     ENDDO
  ENDDO

! -------------------------------------------------------------------
! 13 term
! -------------------------------------------------------------------
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, u, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_X(imode_fdm, imax, jmax, kmax, i1bc,&
       dx, w, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)

! calculate mean profiles
  IF ( ifluc .EQ. 1 ) THEN
     DO j=1,jmax
        DO k=1,kmax
           DO i=1,imax
!              wrk3d(i,5,k) = vis(i,j,k)*(tmp1(i,j,k)+tmp2(i,j,k))
              wrk3d(i,5,k) = (tmp1(i,j,k)+tmp2(i,j,k))
           ENDDO
        ENDDO

        Tau_xz(j) = AVG_IK(imax, jmax, kmax, i5, wrk3d, dx, dz, area)
     ENDDO
  ELSE
     DO j = 1,jmax
        Tau_xz(j) = C_0_R
     ENDDO
  ENDIF

! calculate dissipation
  DO k=1, kmax
     DO j=1, jmax
        DO i=1, imax
!           tau13 = vis(i,j,k)*(tmp1(i,j,k)+tmp2(i,j,k))-Tau_xz(j)
           tau13 = (tmp1(i,j,k)+tmp2(i,j,k))-Tau_xz(j)
           eps(i,j,k) = eps(i,j,k) + tau13*(tmp1(i,j,k)+tmp2(i,j,k))
        ENDDO
     ENDDO
  ENDDO

! -------------------------------------------------------------------
! 23 term
! -------------------------------------------------------------------      
  CALL PARTIAL_Z(imode_fdm, imax, jmax, kmax, k1bc,&
       dz, v, tmp1, i0, i0, wrk1d, wrk2d, wrk3d)
  CALL PARTIAL_Y(imode_fdm, imax, jmax, kmax, j1bc,&
       dy, w, tmp2, i0, i0, wrk1d, wrk2d, wrk3d)

! calculate mean profiles
  IF ( ifluc .EQ. 1 ) THEN
     DO j=1,jmax
        DO k=1,kmax
           DO i=1,imax
!              wrk3d(i,6,k) = vis(i,j,k)*(tmp1(i,j,k)+tmp2(i,j,k))
              wrk3d(i,6,k) = (tmp1(i,j,k)+tmp2(i,j,k))
           ENDDO
        ENDDO

        Tau_yz(j) = AVG_IK(imax, jmax, kmax, i6, wrk3d, dx, dz, area) 
     ENDDO

     DO j=1,jmax
        rW(j) = AVG_IK(imax, jmax, kmax, j, w, dx, dz, area)
     ENDDO
     CALL PARTIAL_Y(imode_fdm, i1, jmax, i1, j1bc, dy, rW(1), &
          rW_y(1), i0, i0, wrk1d, wrk2d, wrk3d)
  ELSE
     DO j = 1,jmax
        Tau_yz(j) = C_0_R
        rW_y(j) = C_0_R
     ENDDO

  ENDIF

! calculate dissipation
  DO k=1, kmax
     DO j=1, jmax
        DO i=1, imax               
!           tau23 = vis(i,j,k)*(tmp1(i,j,k)+tmp2(i,j,k))-Tau_yz(j)
           tau23 = (tmp1(i,j,k)+tmp2(i,j,k))-Tau_yz(j)
           wpy = tmp2(i,j,k)-rW_y(j)

           eps(i,j,k) = eps(i,j,k) + tau23*(tmp1(i,j,k)+wpy)
        ENDDO
     ENDDO
  ENDDO

! -------------------------------------------------------------------
! Final calculation
! -------------------------------------------------------------------      
  DO k=1, kmax
     DO j=1, jmax
        DO i=1, imax
           eps(i,j,k) = visc*eps(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FI_DISSIPATION
