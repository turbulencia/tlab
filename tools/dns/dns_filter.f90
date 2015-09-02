!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"

SUBROUTINE DNS_FILTER(flag_save, y,dx,dy,dz, q,s, txc, vaux, wrk1d,wrk2d,wrk3d)

  USE DNS_GLOBAL
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  LOGICAL flag_save

  TREAL, DIMENSION(*)             :: y, dx, dy, dz
  TREAL, DIMENSION(isize_field,*) :: q, s, txc
  TREAL, DIMENSION(*)             :: wrk2d, wrk3d, vaux
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d

  TARGET q

! -----------------------------------------------------------------------
  TINTEGER ij, id, is, ibc_x(4), ibc_y(4), ibc_z(4)

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: u, v, w, e, rho, p, T, vis

! #######################################################################
! Define pointers
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)

  e   => q(:,4)
  rho => q(:,5)
  p   => q(:,6)
  T   => q(:,7)

  vis => q(:,8)

! BCs for the filters (see routine FILTER)
  ibc_x(1) = ifilt_x; ibc_x(2) = i1bc; ibc_x(3) = 0; ibc_x(4) = 0
  ibc_y(1) = ifilt_y; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0 
  ibc_z(1) = ifilt_z; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0 

! #######################################################################
! Domain filter
! #######################################################################
! -----------------------------------------------------------------------
! Flow
! -----------------------------------------------------------------------
! Statistics
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. frunstat .EQ. 1 ) THEN
     CALL DNS_SAVE_AVGKIN(rho, u, v, w, txc(1,1), txc(1,2), txc(1,3), txc(1,4), &
          txc(1,5), txc(1,6), txc(1,7), vaux(vindex(VA_MEAN_WRK)), wrk2d)     
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. ffltdmp .EQ. 1 .AND. flag_save ) THEN
     CALL FI_DISSIPATION(i1, imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
          area, visc, dx, dy, dz, rho, u, v, w, txc(1,1), &
          txc(1,2), txc(1,3), txc(1,4), wrk1d, wrk1d(1,6), wrk2d, wrk3d)         
     CALL DNS_AVG_KIN(i1, itime, rtime, imax, jmax, kmax, y, &
          dx, dz, area, rho, u, v, w, txc, vaux(vindex(VA_MEAN_WRK)), wrk3d)
  ENDIF
  
! Filter     
  DO ij = 1,isize_field
     txc(ij,1) = rho(ij)*u(ij)
     txc(ij,2) = rho(ij)*v(ij)
     txc(ij,3) = rho(ij)*w(ij)
     txc(ij,4) = rho(ij)*e(ij)
  ENDDO
     
  IF ( ifilt_domain .EQ. DNS_FILTER_6E        .OR. &
       ifilt_domain .EQ. DNS_FILTER_4E        .OR. &
       ifilt_domain .EQ. DNS_FILTER_COMPACT ) THEN
     id = 1
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, rho,     &
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,1),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,2),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,3),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,wrk3d)
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,4),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,wrk3d)
     
  ELSE  IF ( ifilt_domain .EQ. DNS_FILTER_ADM ) THEN ! 2 txc array needed
     id = 1
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, rho,     &
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc(1,5))
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,1),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc(1,5))
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,2),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc(1,5))
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,3),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc(1,5))
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,4),&
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc(1,5))
     
  ENDIF
     
  DO ij = 1,isize_field
     u(ij) = txc(ij,1)/rho(ij)
     v(ij) = txc(ij,2)/rho(ij)
     w(ij) = txc(ij,3)/rho(ij)
     e(ij) = txc(ij,4)/rho(ij)
  ENDDO
     
! Statistics
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. ffltdmp .EQ. 1 .AND. flag_save ) THEN
     CALL FI_DISSIPATION(i1, imode_fdm, imax, jmax, kmax, i1bc, j1bc, k1bc, &
          area, visc, dx, dy, dz, rho, u, v, w, txc(1,1), &
          txc(1,2), txc(1,3), txc(1,4), wrk1d, wrk1d(1,6), wrk2d, wrk3d)         
     CALL DNS_AVG_KIN(i2, itime, rtime, imax, jmax, kmax, y, &
          dx, dz, area, rho, u, v, w, txc, vaux(vindex(VA_MEAN_WRK)), wrk3d)
  ENDIF
  
! recalculation of p and T
  CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
  CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
  IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)
  
! -------------------------------------------------------------------
! Scalar
! -------------------------------------------------------------------
  IF ( icalc_scal .EQ. 1 .AND. ifilt_scalar .EQ. 1 ) THEN
     DO is = 1,inb_scal
        DO ij = 1,isize_field
           txc(ij,1) = rho(ij)*s(ij,is)
        ENDDO
        
        IF ( ifilt_domain .EQ. DNS_FILTER_4E        .OR. &
             ifilt_domain .EQ. DNS_FILTER_6E        .OR. &
             ifilt_domain .EQ. DNS_FILTER_COMPACT ) THEN
           id = 1
           CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,1), &
                vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,wrk3d)
           
        ELSE IF (ifilt_domain .EQ. DNS_FILTER_ADM ) THEN ! 2 txc array needed
           id = 1
           CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, txc(1,1), &
                vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc(1,5))
                      
        ENDIF
        
        DO ij = 1,isize_field
           s(ij,is) = txc(ij,1)/rho(ij)
        ENDDO
        
        IF ( ilimit_scal .EQ. 1 ) THEN
           DO ij = 1,isize_field
              s(ij,is) = MIN(MAX(s(ij,is),z_bound_min), z_bound_max)
           ENDDO
        ENDIF
        
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE DNS_FILTER
     
! #######################################################################
! #######################################################################
SUBROUTINE DNS_AVG_KIN(iflag, itime, rtime, imax,jmax,kmax, y, dx, dz, &
     area, rho, u, v, w, eps, mean2d, wrk3d)

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER iflag, itime
  TINTEGER imax, jmax, kmax
  TREAL rtime, area

  TREAL, DIMENSION(*)              :: y, dx, dz
  TREAL, DIMENSION(imax,jmax,kmax) :: rho, u, v, w, eps, wrk3d

  TREAL mean2d(jmax, 8)

  TINTEGER i, j, k
  TREAL AVG_IK
  CHARACTER*32 fname
#ifdef USE_MPI
  INTEGER ims_err, ims_pro
#endif

#define rR(j)     mean2d(j,1)
#define fU(j)     mean2d(j,2)
#define fV(j)     mean2d(j,3)
#define fW(j)     mean2d(j,4)
#define Kin_0(j)  mean2d(j,5)
#define Kin_1(j)  mean2d(j,6)
#define Eps_0(j)  mean2d(j,7)
#define Eps_1(j)  mean2d(j,8)

  DO j=1,jmax

     rR(j) = AVG_IK(imax, jmax, kmax, j, rho, dx, dz, area)

     DO k=1, kmax
        DO i=1, imax
           wrk3d(i,1,k) = rho(i,j,k)*u(i,j,k)
           wrk3d(i,2,k) = rho(i,j,k)*v(i,j,k)
           wrk3d(i,3,k) = rho(i,j,k)*w(i,j,k)
        ENDDO
     ENDDO

     fU(j) = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)/rR(j)
     fV(j) = AVG_IK(imax, jmax, kmax, i2, wrk3d, dx, dz, area)/rR(j)
     fW(j) = AVG_IK(imax, jmax, kmax, i3, wrk3d, dx, dz, area)/rR(j)

     DO k=1, kmax
        DO i=1, imax
           wrk3d(i,1,k) = ((u(i,j,k)-fU(j))**2 + (v(i,j,k)-fV(j))**2 + (w(i,j,k)-fW(j))**2)*C_05_R
        ENDDO
     ENDDO

     IF ( iflag .EQ. 1 ) THEN
        Kin_0(j) = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        Eps_0(j) = AVG_IK(imax, jmax, kmax, j, eps, dx, dz, area)
     ELSE
        Kin_1(j) = AVG_IK(imax, jmax, kmax, i1, wrk3d, dx, dz, area)
        Eps_1(j) = AVG_IK(imax, jmax, kmax, j, eps, dx, dz, area)
     ENDIF

  ENDDO

  IF ( iflag .EQ. 2 ) THEN

#ifdef USE_MPI
     CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)

     IF ( ims_pro .EQ. 0 ) THEN
#endif
        WRITE(fname,*) itime; fname='kin'//TRIM(ADJUSTL(fname))

        OPEN(UNIT=21, FILE=fname)

        WRITE(21, '(A8,E14.7E3)') 'RTIME = ', rtime
        WRITE(21, '(A7,I8)') 'IMAX = ', i1
        WRITE(21, '(A7,I8)') 'JMAX = ', jmax

        WRITE(21, '(A)') 'I J Y Eps_before Eps_after Kin_before Kin_after'

        DO j=1, jmax
           WRITE(21,1000) 1, j, y(j), Eps_0(j), Eps_1(j), Kin_0(j), Kin_1(j)
        ENDDO

        CLOSE(21)

1000    FORMAT(I5,(1X,I5),5(1X,E12.5E3))

#ifdef USE_MPI
     ENDIF
#endif

  ENDIF

  RETURN
END SUBROUTINE DNS_AVG_KIN
