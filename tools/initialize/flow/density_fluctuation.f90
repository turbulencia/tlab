#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/10/08 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Setting up a perturbation of the thermodynamic fields by a 
!# displacement of the reference center plane.
!#
!# Array s enters with the scalar total field, including fluctuations.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE DENSITY_FLUCTUATION(code, x,y,z,dx,dz, s, p, rho, T, h, disp, wrk3d)

  USE DNS_GLOBAL
  USE THERMO_GLOBAL, ONLY : imixture
  USE FLOW_LOCAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER code

  TREAL, DIMENSION(*)                :: x, y, z, dx, dz
  TREAL, DIMENSION(imax,jmax,kmax)   :: T, h, rho, p, wrk3d
  TREAL, DIMENSION(imax,jmax,kmax,*) :: s
  TREAL, DIMENSION(imax,kmax)        :: disp

! -------------------------------------------------------------------
  TINTEGER i, j, k, inx2d, inx3d, inz3d
  TINTEGER idummy, idsp, iprof_loc
  TREAL wx, wz, wxloc, wzloc, dummy, ycenter, mean_loc, delta_loc
  TREAL AVG_IK, FLOW_SHEAR_TEMPORAL
  TREAL xcenter, amplify

! ###################################################################
  disp = C_0_R

! ###################################################################
! Center plane displacement
! ###################################################################
! -------------------------------------------------------------------
! Broadband case
! -------------------------------------------------------------------
  IF ( code .EQ. 4 ) THEN
     idummy = jmax_total; jmax_total = 1
     CALL DNS_READ_FIELDS('scal.rand', i1, imax,i1,kmax, i1,i0, isize_field, disp, wrk3d)
     jmax_total = idummy
! remove mean
     dummy = AVG_IK(imax, i1, kmax, i1, disp, dx, dz, area)
     DO k = 1,kmax
        DO i = 1,imax
           disp(i,k) = disp(i,k)-dummy
        ENDDO
     ENDDO
  ENDIF

! -------------------------------------------------------------------
! Discrete case
! -------------------------------------------------------------------
  IF ( code .EQ. 5 ) THEN
     wx = C_2_R*C_PI_R/scalex
     wz = C_2_R*C_PI_R/scalez

! 1D perturbation along X
     DO inx2d = 1,nx2d
        wxloc = M_REAL(inx2d)*wx
        DO k = 1,kmax
           DO i = 1,imax
              disp(i,k) = disp(i,k) + &
                   A2d(inx2d)*COS(wxloc*x(i)+Phix2d(inx2d))
           ENDDO
        ENDDO
     ENDDO

! 2D perturbation along X and Z
     IF (kmax_total .GT. 1) THEN
#ifdef USE_MPI
        idsp = ims_offset_k 
#else
        idsp = 0
#endif
        DO inx3d = 1,nx3d
           DO inz3d = 1,nz3d
              wxloc = M_REAL(inx3d)*wx
              wzloc = M_REAL(inz3d)*wz
              DO k = 1,kmax
                 DO i = 1,imax
                    disp(i,k) = disp(i,k) + &
                         A3d(inx3d)*COS(wxloc*x(i)+Phix3d(inx3d))* &
                         COS(wzloc*z(idsp+k)+Phiz3d(inz3d))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDIF

  ENDIF

! -------------------------------------------------------------------
! Modulation
! -------------------------------------------------------------------
  IF ( frc_delta .GT. C_0_R ) THEN
     DO k = 1,kmax
        DO i = 1,imax
           xcenter   = x(i) - scalex*C_05_R - x(1)
           amplify   = EXP(-(C_05_R*xcenter/frc_delta)**2)
           disp(i,k) = disp(i,k)*amplify
        ENDDO
     ENDDO
  ENDIF

! ###################################################################
! Perturbation in the thermodynamic fields
! ###################################################################
  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
! -------------------------------------------------------------------
! temperature
! -------------------------------------------------------------------
     IF ( iprof_rho .EQ. PROFILE_NONE ) THEN

! temperature/mixture profile is given
        IF ( iprof_tem .GT. 0 ) THEN
           DO k = 1,kmax
              DO i = 1,imax
                 delta_loc = delta_tem + (prof_tem(2)-prof_tem(1))*disp(i,k)*scaley
                 mean_loc  = mean_tem  + C_05_R*(prof_tem(2)+prof_tem(1))*disp(i,k)*scaley
                 ycenter   = y(1) + scaley*ycoor_tem + disp(i,k)
                 DO j = 1,jmax
                    T(i,j,k) =  FLOW_SHEAR_TEMPORAL&
                         (iprof_tem, thick_tem, delta_loc, mean_loc, ycenter, prof_tem, y(j))
                 ENDDO
              ENDDO
           ENDDO

           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
              CALL THERMO_AIRWATER_PT(imax, jmax, kmax, s, p, T)
           ENDIF

! enthalpy/mixture profile is given
        ELSE IF ( iprof_tem .LT. 0 ) THEN
           DO k = 1,kmax
              DO i = 1,imax
                 delta_loc = delta_tem + (prof_tem(2)-prof_tem(1))*disp(i,k)*scaley
                 mean_loc  = mean_tem  + C_05_R*(prof_tem(2)+prof_tem(1))*disp(i,k)*scaley
                 ycenter   = y(1) + scaley*ycoor_tem + disp(i,k)
                 iprof_loc =-iprof_tem
                 DO j = 1,jmax
                    h(i,j,k) =  FLOW_SHEAR_TEMPORAL&
                         (iprof_loc, thick_tem, delta_loc, mean_loc, ycenter, prof_tem, y(j))
                 ENDDO
              ENDDO
           ENDDO

           IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
!              CALL THERMO_AIRWATER_PH(imax, jmax, kmax, s, p, h, T, wrk3d)
              CALL THERMO_AIRWATER_PH2(imax,jmax,kmax, s, p, h, T)
           ENDIF

        ENDIF

! compute perturbation in density
        CALL THERMO_THERMAL_DENSITY(imax,jmax,kmax, s, p, T, rho)

     ELSE

! -------------------------------------------------------------------
! density
! -------------------------------------------------------------------
! to be developed
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE DENSITY_FLUCTUATION
