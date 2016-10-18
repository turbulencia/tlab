#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!# 2016/05/27 - J.P. Mellado
!#              Cleaning
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
SUBROUTINE DNS_FILTER(flag_save, q,s, txc, vaux, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : lfile
  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL

  IMPLICIT NONE

#include "integers.h"

  LOGICAL flag_save
  TREAL, DIMENSION(isize_field,*) :: q, s
  TREAL, DIMENSION(isize_field,*) :: txc
  TREAL, DIMENSION(*)             :: wrk2d, wrk3d, vaux
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d

  TARGET q

! -----------------------------------------------------------------------
  TINTEGER id, iq,is,ij, ibc_x(4), ibc_y(4), ibc_z(4)
  CHARACTER*250 line

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: e, rho, p, T, vis
  TREAL, DIMENSION(:), POINTER :: y, dx,dy,dz

! #######################################################################
! Define pointers
                   dx => g(1)%jac(:,1)
  y => g(2)%nodes; dy => g(2)%jac(:,1)
                   dz => g(3)%jac(:,1)

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     e   => q(:,4)
     rho => q(:,5)
     p   => q(:,6)
     T   => q(:,7)
     
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,8)

  ENDIF
  
! BCs for the filters (see routine FILTER)
  ibc_x(1) = ifilt_x; ibc_x(2) = i1bc; ibc_x(3) = 0; ibc_x(4) = 0
  ibc_y(1) = ifilt_y; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0 
  ibc_z(1) = ifilt_z; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0 

! #######################################################################
! Domain filter
! #######################################################################
  WRITE(line,*) itime; line = 'Filtering fields at It'//TRIM(ADJUSTL(line))//'.'
  CALL IO_WRITE_ASCII(lfile,line)

! Statistics
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. frunstat .EQ. 1 ) THEN
     CALL DNS_SAVE_AVGKIN(q(1,5), q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), txc(1,3), txc(1,4), &
          txc(1,5), txc(1,6), txc(1,7), vaux(vindex(VA_MEAN_WRK)), wrk2d)     
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. ffltdmp .EQ. 1 .AND. flag_save ) THEN
     CALL FI_DISSIPATION(i1, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          area, visc, dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1), &
          txc(1,2),txc(1,3),txc(1,4), wrk1d, wrk1d(1,6), wrk2d, wrk3d)         
     CALL DNS_AVG_KIN(i1, itime, rtime, imax,jmax,kmax, y, dx,dz, area, q, txc, vaux(vindex(VA_MEAN_WRK)), wrk3d)
  ENDIF
  
! filtering
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! contruct fields per unit volume
     DO iq = 1,inb_flow-1
        q(:,iq) = q(:,iq) *q(:,inb_flow)
     ENDDO
     DO is = 1,inb_scal
        s(:,is) = s(:,is) *q(:,inb_flow)
     ENDDO
  ENDIF
  
  id = 1
  DO iq = 1,inb_flow
     CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, q(1,iq),     &
          vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc)
  ENDDO
  IF ( icalc_scal .EQ. 1 .AND. ifilt_scalar .EQ. 1 ) THEN
     DO is = 1,inb_scal
        CALL OPR_FILTER(ifilt_domain, imax,jmax,kmax, ibc_x,ibc_y,ibc_z, id, s(1,is),     &
             vaux(vindex(VA_FLT_CX)),vaux(vindex(VA_FLT_CY)),vaux(vindex(VA_FLT_CZ)), wrk1d,wrk2d,txc)
     ENDDO
  ENDIF
  
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! re-contruct fields per unit mass
     DO iq = 1,inb_flow-1
        q(:,iq) = q(:,iq) /q(:,inb_flow)
     ENDDO
     DO is = 1,inb_scal
        s(:,is) = s(:,is) /q(:,inb_flow)
     ENDDO
  ENDIF

  IF ( icalc_scal .EQ. 1 .AND. ifilt_scalar .EQ. 1 ) THEN
     IF ( ilimit_scal .EQ. 1 ) THEN
        DO ij = 1,isize_field
           s(ij,is) = MIN(MAX(s(ij,is),s_bound_min(is)), s_bound_max(is))
        ENDDO
     ENDIF
  ENDIF

! statistics
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. ffltdmp .EQ. 1 .AND. flag_save ) THEN
     CALL FI_DISSIPATION(i1, imode_fdm, imax,jmax,kmax, i1bc,j1bc,k1bc, &
          area, visc, dx,dy,dz, q(1,1),q(1,2),q(1,3), txc(1,1), &
          txc(1,2),txc(1,3),txc(1,4), wrk1d, wrk1d(1,6), wrk2d, wrk3d)         
     CALL DNS_AVG_KIN(i2, itime, rtime, imax,jmax,kmax, y, dx,dz, area, q, txc, vaux(vindex(VA_MEAN_WRK)), wrk3d)
  ENDIF
  
! recalculation of diagnostic variables
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN 
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))        
     ENDIF

  ELSE
     CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
     CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) CALL THERMO_VISCOSITY(imax,jmax,kmax, T, vis)

  ENDIF

  RETURN
END SUBROUTINE DNS_FILTER
     
! #######################################################################
! #######################################################################
#define rR(j)     mean2d(j,1)
#define fU(j)     mean2d(j,2)
#define fV(j)     mean2d(j,3)
#define fW(j)     mean2d(j,4)
#define Kin_0(j)  mean2d(j,5)
#define Kin_1(j)  mean2d(j,6)
#define Eps_0(j)  mean2d(j,7)
#define Eps_1(j)  mean2d(j,8)

SUBROUTINE DNS_AVG_KIN(iflag, itime, rtime, nx,ny,nz, y, dx,dz, area, q, eps, mean2d, wrk3d)

  USE DNS_GLOBAL, ONLY : imode_eqns

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER iflag, itime, nx,ny,nz
  TREAL rtime, area
  TREAL, DIMENSION(*),          INTENT(IN)    :: y, dx, dz
  TREAL, DIMENSION(nx,ny,nz,*), INTENT(IN)    :: q
  TREAL, DIMENSION(nx,ny,nz),   INTENT(IN)    :: eps
  TREAL, DIMENSION(nx,ny,nz),   INTENT(INOUT) :: wrk3d
  TREAL, DIMENSION(ny,8),       INTENT(INOUT) :: mean2d

! -----------------------------------------------------------------------
  TINTEGER i, j, k
  TREAL AVG_IK
  CHARACTER*32 fname
#ifdef USE_MPI
  INTEGER ims_err, ims_pro
#endif

! #######################################################################
  DO j = 1,ny

     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        rR(j) = AVG_IK(nx, ny, nz, j, q(1,1,1,5), dx, dz, area)
        DO k = 1,nz; DO i = 1,nx
           wrk3d(i,1,k) = q(i,j,k,5) *q(i,j,k,1)
           wrk3d(i,2,k) = q(i,j,k,5) *q(i,j,k,2)
           wrk3d(i,3,k) = q(i,j,k,5) *q(i,j,k,3)
        ENDDO; ENDDO

     ELSE
        rR(j) = C_1_R
        DO k=1, nz; DO i=1, nx
           wrk3d(i,1,k) = q(i,j,k,1)
           wrk3d(i,2,k) = q(i,j,k,2)
           wrk3d(i,3,k) = q(i,j,k,3)
        ENDDO; ENDDO
     
     ENDIF

     fU(j) = AVG_IK(nx,ny,nz, i1, wrk3d, dx,dz, area) /rR(j)
     fV(j) = AVG_IK(nx,ny,nz, i2, wrk3d, dx,dz, area) /rR(j)
     fW(j) = AVG_IK(nx,ny,nz, i3, wrk3d, dx,dz, area) /rR(j)

     DO k = 1,nz; DO i = 1,nx
        wrk3d(i,1,k) = rR(j) *((q(i,j,k,1)-fU(j))**2 + (q(i,j,k,2)-fV(j))**2 + (q(i,j,k,3)-fW(j))**2) *C_05_R
     ENDDO; ENDDO

     IF ( iflag .EQ. 1 ) THEN
        Kin_0(j) = AVG_IK(nx,ny,nz, i1, wrk3d, dx,dz, area)
        Eps_0(j) = AVG_IK(nx,ny,nz, j,  eps,   dx,dz, area)
     ELSE
        Kin_1(j) = AVG_IK(nx,ny,nz, i1, wrk3d, dx,dz, area)
        Eps_1(j) = AVG_IK(nx,ny,nz, j,  eps,   dx,dz, area)
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
        WRITE(21, '(A7,I8)') 'JMAX = ', ny

        WRITE(21, '(A)') 'I J Y Eps_before Eps_after Kin_before Kin_after'

        DO j=1, ny
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
