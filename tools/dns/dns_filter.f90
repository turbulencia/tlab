#include "types.h"
#include "dns_const.h"

!########################################################################
!########################################################################
SUBROUTINE DNS_FILTER(flag_save, q,s, txc, vaux, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : lfile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow,inb_scal,inb_scal_array, isize_field, isize_wrk1d
  USE DNS_GLOBAL,    ONLY : imode_eqns, imode_sim, itransport
  USE DNS_GLOBAL,    ONLY : itime, rtime
  USE DNS_GLOBAL,    ONLY : epbackground, pbackground, damkohler
  USE DNS_GLOBAL,    ONLY : FilterDomain
  USE THERMO_GLOBAL, ONLY : imixture
  USE DNS_LOCAL,     ONLY : nitera_stats_spa, ffltdmp
  USE DNS_LOCAL,     ONLY : vindex, VA_MEAN_WRK
  USE DNS_LOCAL,     ONLY : ilimit_scal, s_bound_min, s_bound_max  

  IMPLICIT NONE

#include "integers.h"

  LOGICAL flag_save
  TREAL, DIMENSION(isize_field,*) :: q, s
  TREAL, DIMENSION(isize_field,*) :: txc
  TREAL, DIMENSION(*)             :: wrk2d, wrk3d, vaux
  TREAL, DIMENSION(isize_wrk1d,*) :: wrk1d

  TARGET q

! -----------------------------------------------------------------------
  TINTEGER iq,is,ij
  CHARACTER*250 line

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: e, rho, p, T, vis

! #######################################################################
! Define pointers
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     e   => q(:,4)
     rho => q(:,5)
     p   => q(:,6)
     T   => q(:,7)
     
     IF ( itransport .EQ. EQNS_TRANS_SUTHERLAND .OR. itransport .EQ. EQNS_TRANS_POWERLAW ) vis => q(:,8)

  ENDIF
  
! Might be better to filter the pressure instead of the energy in compressible flows
  ! IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
  !    iq_loc = (/ 5,1,2,3,6 /) ! Filtered variables: rho, u,v,w, p
  ! ELSE
  !    iq_loc = (/ 1,2,3 /)
  ! ENDIF

! #######################################################################
! Domain filter
! #######################################################################
  WRITE(line,*) itime; line = 'Filtering fields at It'//TRIM(ADJUSTL(line))//'.'
  CALL IO_WRITE_ASCII(lfile,line)

! Statistics
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL .AND. nitera_stats_spa .GT. 0 ) THEN
     CALL DNS_SAVE_AVGKIN(q(1,5), q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2), txc(1,3), txc(1,4), &
          txc(1,5), txc(1,6), txc(1,7), vaux(vindex(VA_MEAN_WRK)), wrk2d)     
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. ffltdmp .EQ. 1 .AND. flag_save ) THEN
     CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)         
     CALL DNS_AVG_KIN(i1, itime, rtime, imax,jmax,kmax,q, txc, vaux(vindex(VA_MEAN_WRK)), wrk3d)
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
  
  DO iq = 1,inb_flow
     CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, q(1,iq), wrk1d,wrk2d,txc)
  ENDDO
  DO is = 1,inb_scal
     CALL OPR_FILTER(imax,jmax,kmax, FilterDomain, s(1,is), wrk1d,wrk2d,txc)
  ENDDO
  
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN ! re-contruct fields per unit mass
     DO iq = 1,inb_flow-1
        q(:,iq) = q(:,iq) /q(:,inb_flow)
     ENDDO
     DO is = 1,inb_scal
        s(:,is) = s(:,is) /q(:,inb_flow)
     ENDDO
  ENDIF

  IF ( ilimit_scal .EQ. 1 ) THEN
     DO ij = 1,isize_field
        s(ij,is) = MIN(MAX(s(ij,is),s_bound_min(is)), s_bound_max(is))
     ENDDO
  ENDIF

! statistics
  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. ffltdmp .EQ. 1 .AND. flag_save ) THEN
     CALL FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5), wrk1d,wrk2d,wrk3d)         
     CALL DNS_AVG_KIN(i2, itime, rtime, imax,jmax,kmax, q, txc, vaux(vindex(VA_MEAN_WRK)), wrk3d)
  ENDIF
  
! #######################################################################
! recalculation of diagnostic variables
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN
        CALL THERMO_AIRWATER_PH(imax,jmax,kmax, s(1,2), s(1,1), epbackground,pbackground)
        
     ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR                        ) THEN 
        CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(1,inb_scal_array))
        
     ENDIF

  ELSE
     ! IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     !    CALL THERMO_AIRWATER_RP(imax,jmax,kmax, s, p, rho, T, wrk3d)
     ! ELSE
     !    CALL THERMO_THERMAL_TEMPERATURE(imax,jmax,kmax, s, p, rho, T)
     ! ENDIF
     ! CALL THERMO_CALORIC_ENERGY(imax,jmax,kmax, s, T, e)

! This recalculation of T and p is made to make sure that the same numbers are
! obtained in statistics postprocessing as in the simulation; avg* files
! can then be compared with diff command.
!     IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     CALL THERMO_CALORIC_TEMPERATURE(imax,jmax,kmax, s, e, rho, T, wrk3d)
     CALL THERMO_THERMAL_PRESSURE(imax,jmax,kmax, s, rho, T, p)
!     ENDIF
     
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

SUBROUTINE DNS_AVG_KIN(iflag, itime, rtime, nx,ny,nz, q, eps, mean2d, wrk3d)

  USE DNS_GLOBAL, ONLY : imode_eqns
  USE DNS_GLOBAL, ONLY : g, area

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER iflag, itime, nx,ny,nz
  TREAL rtime
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
        rR(j) = AVG_IK(nx, ny, nz, j, q(1,1,1,5), g(1)%jac, g(3)%jac, area)
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

     fU(j) = AVG_IK(nx,ny,nz, i1, wrk3d, g(1)%jac,g(3)%jac, area) /rR(j)
     fV(j) = AVG_IK(nx,ny,nz, i2, wrk3d, g(1)%jac,g(3)%jac, area) /rR(j)
     fW(j) = AVG_IK(nx,ny,nz, i3, wrk3d, g(1)%jac,g(3)%jac, area) /rR(j)

     DO k = 1,nz; DO i = 1,nx
        wrk3d(i,1,k) = rR(j) *((q(i,j,k,1)-fU(j))**2 + (q(i,j,k,2)-fV(j))**2 + (q(i,j,k,3)-fW(j))**2) *C_05_R
     ENDDO; ENDDO

     IF ( iflag .EQ. 1 ) THEN
        Kin_0(j) = AVG_IK(nx,ny,nz, i1, wrk3d, g(1)%jac,g(3)%jac, area)
        Eps_0(j) = AVG_IK(nx,ny,nz, j,  eps,   g(1)%jac,g(3)%jac, area)
     ELSE
        Kin_1(j) = AVG_IK(nx,ny,nz, i1, wrk3d, g(1)%jac,g(3)%jac, area)
        Eps_1(j) = AVG_IK(nx,ny,nz, j,  eps,   g(1)%jac,g(3)%jac, area)
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
           WRITE(21,1000) 1, j, g(2)%nodes(j), Eps_0(j), Eps_1(j), Kin_0(j), Kin_1(j)
        ENDDO

        CLOSE(21)

1000    FORMAT(I5,(1X,I5),5(1X,E12.5E3))

#ifdef USE_MPI
     ENDIF
#endif

  ENDIF

  RETURN
END SUBROUTINE DNS_AVG_KIN
