#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2007/05/08 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Evaluate the integral \int_ycenter^y dx/H(x), where H(x) is the scale 
!# height in the system, where ycenter is some given value
!#
!########################################################################
SUBROUTINE FI_HYDROSTATIC(imax,jmax,kmax, ycenter, y, integral)

  IMPLICIT NONE

  TINTEGER imax, jmax, kmax
  TREAL ycenter
  TREAL y(jmax)
  TREAL integral(imax,jmax,kmax)

! -------------------------------------------------------------------
  TINTEGER i, j, k, jcenter
  TREAL FI_HYDROSTATIC_SCALEHEIGHT_INV
  EXTERNAL FI_HYDROSTATIC_SCALEHEIGHT_INV

  TREAL ABSERR,EPSABS,EPSREL,RESULT
  TINTEGER IER,NEVAL
  TREAL WORK, params(10)
  TINTEGER IWORK,KEY,LAST,LENW,LIMIT
  DIMENSION IWORK(100),WORK(400)

! ###################################################################
  DO j = 1,jmax
     IF ( y(j  ) .LE. ycenter .AND. &
          y(j+1) .GT. ycenter ) THEN
        jcenter = j
        EXIT
     ENDIF
  ENDDO

! Convergence constants
  EPSABS = C_1EM10_R
  EPSREL = C_1EM10_R
  KEY = 6
  LIMIT = 100
  LENW = LIMIT*4

  params(1) = y(1)
  params(2) = y(1)
  params(3) = y(jmax)
  params(4) = C_0_R
  params(5) = C_0_R
  params(6) = C_0_R
  params(7) = C_0_R

  DO j = jcenter+1,jmax
     CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,ycenter,y(j),EPSABS,&
          EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
     DO k = 1,kmax
        DO i = 1,imax
           integral(i,j,k) = RESULT
        ENDDO
     ENDDO
  ENDDO

  DO j = jcenter,1,-1
     CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(j),ycenter,EPSABS,&
          EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
     DO k = 1,kmax
        DO i = 1,imax
           integral(i,j,k) = -RESULT
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FI_HYDROSTATIC
      
!########################################################################
!# Modification from FI_HYDROSTATIC to allow passing a function through
!# parameters; nvar variables enter through wrk1d
!#
!########################################################################
SUBROUTINE FI_HYDROSTATIC_STEPS(imax,jmax,kmax, nvar, ycenter, y, integral, wrk1d)
  
  IMPLICIT NONE

  TINTEGER imax, jmax, kmax, nvar
  TREAL ycenter
  TREAL y(jmax)
  TREAL integral(imax,jmax,kmax)
  TREAL wrk1d(jmax,nvar)

! -------------------------------------------------------------------
  TINTEGER i, j, k, iv, jcenter
  TREAL FI_HYDROSTATIC_SCALEHEIGHT_INV
  EXTERNAL FI_HYDROSTATIC_SCALEHEIGHT_INV

  TREAL ABSERR,EPSABS,EPSREL,RESULT
  TINTEGER IER,NEVAL
  TREAL WORK, params(10)
  TINTEGER IWORK,KEY,LAST,LENW,LIMIT
  DIMENSION IWORK(100),WORK(400)

! ###################################################################
  DO iv = 1,10
     params(iv) = C_0_R
  ENDDO

  DO j = 1,jmax
     IF ( y(j  ) .LE. ycenter .AND. &
          y(j+1) .GT. ycenter ) THEN
        jcenter = j
        EXIT
     ENDIF
  ENDDO

! Convergence constants
  EPSABS = C_1EM10_R
  EPSREL = C_1EM10_R
  KEY = 6
  LIMIT = 100
  LENW = LIMIT*4

  params(1) = y(1)

! from center, upwards
  params(2) = y(jcenter)
  params(3) = y(jcenter+1)
  DO iv = 1,nvar
     params((iv-1)*2+4) = wrk1d(jcenter  ,iv)
     params((iv-1)*2+5) = wrk1d(jcenter+1,iv)
  ENDDO
  CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,ycenter,y(jcenter+1),EPSABS,&
       EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
  DO k = 1,kmax
     DO i = 1,imax
        integral(i,jcenter+1,k) = RESULT
     ENDDO
  ENDDO
  DO j = jcenter+2,jmax
     params(2) = y(j-1)
     params(3) = y(j)
     DO iv = 1,nvar
        params((iv-1)*2+4) = wrk1d(j-1,iv)
        params((iv-1)*2+5) = wrk1d(j  ,iv)
     ENDDO
     CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(j-1),y(j),EPSABS,&
          EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
     DO k = 1,kmax
        DO i = 1,imax
           integral(i,j,k) = integral(i,j-1,k) + RESULT
        ENDDO
     ENDDO
  ENDDO

! from center, downwards
  params(2) = y(jcenter)
  params(3) = y(jcenter+1)
  DO iv = 1,nvar
     params((iv-1)*2+4) = wrk1d(jcenter  ,iv)
     params((iv-1)*2+5) = wrk1d(jcenter+1,iv)
  ENDDO
  CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(jcenter),ycenter,EPSABS,&
       EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
  DO k = 1,kmax
     DO i = 1,imax
        integral(i,jcenter,k) = -RESULT
     ENDDO
  ENDDO
  DO j = jcenter-1,1,-1
     params(2) = y(j)
     params(3) = y(j+1)
     DO iv = 1,nvar
        params((iv-1)*2+4) = wrk1d(j  ,iv)
        params((iv-1)*2+5) = wrk1d(j+1,iv)
     ENDDO
     CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(j),y(j+1),EPSABS,&
          EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
     DO k = 1,kmax
        DO i = 1,imax
           integral(i,j,k) = integral(i,j+1,k)-RESULT
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE FI_HYDROSTATIC_STEPS
       
!########################################################################
! Compute hydrostatic equilibrium from profiles (h,q_t).
!########################################################################
SUBROUTINE FI_HYDROSTATIC_AIRWATER_H(nmax, y, s,h,T,p, wrk1d)

  USE DNS_GLOBAL, ONLY : p_init, ycoor_tem, scaley

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nmax
  TREAL, DIMENSION(nmax)   :: y
  TREAL, DIMENSION(nmax)   :: h, T, p
  TREAL, DIMENSION(nmax,*) :: s, wrk1d

! -------------------------------------------------------------------
  TINTEGER iter, niter
  TREAL ycenter

! ###################################################################
  niter = 10

! ###################################################################
! initialize iteration
! ###################################################################
! pressure from hydrostatic equilibrium
  ycenter = y(1) + ycoor_tem*scaley      
  wrk1d(:,1) = p_init
  CALL FI_HYDROSTATIC_STEPS(i1, nmax, i1, i1, ycenter, y, p, wrk1d)
  p(:) = p_init*EXP(p(:))

! ###################################################################
! iteration
! ###################################################################
  DO iter = 1,niter
! pressure from hydrostatic equilibrium
     ycenter = y(1) + ycoor_tem*scaley      
     wrk1d(:,1) = p(:)
     CALL FI_HYDROSTATIC_STEPS(i1, nmax, i1, i1, ycenter, y, p, wrk1d)
     p(:) = p_init*EXP(p(:))

  ENDDO

! ###################################################################
! compute equilibrium values of q_l and T
! ###################################################################
  CALL THERMO_AIRWATER_PH2(i1, nmax, i1, s, p, h, T)

  RETURN
END SUBROUTINE FI_HYDROSTATIC_AIRWATER_H
      
!########################################################################
!Compute hydrostatic equilibrium from profiles (T,q_t).
!########################################################################
SUBROUTINE FI_HYDROSTATIC_AIRWATER_T(y, dy, z1, T, p, rho, wrk1d, wrk2d, wrk3d)

  USE DNS_GLOBAL, ONLY : p_init, ycoor_tem, scaley
  USE DNS_GLOBAL, ONLY : imode_fdm, jmax, j1bc
  USE DNS_GLOBAL, ONLY : buoyancy

  USE THERMO_GLOBAL, ONLY : dsmooth, WGHT_INV

  IMPLICIT NONE

#include "integers.h"

  TREAL y(jmax), dy(jmax)
  TREAL z1(jmax,*), T(*), p(*), rho(*)
  TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! -------------------------------------------------------------------
  TINTEGER ij, iter, niter
  TREAL ycenter, qsat, dsmooth_loc

! ###################################################################
  niter = 3

! ###################################################################
! initialize iteration from condition q_l = 0
! ###################################################################
  DO ij = 1,jmax
     z1(ij,2) = C_0_R
  ENDDO

! ###################################################################
! iteration
! ###################################################################
  DO iter = 1,niter
! pressure 
     ycenter = y(1) + ycoor_tem*scaley
     CALL FI_HYDROSTATIC_STEPS(i1, jmax, i1, i1, ycenter, y, p, z1(1,2))
     DO ij = 1,jmax
        p(ij) = p_init*EXP(p(ij))
     ENDDO

! density profile from pressure gradient
     CALL PARTIAL_Y(imode_fdm, i1, jmax, i1, j1bc, dy, p, rho, i0, i0, wrk1d, wrk2d, wrk3d)
     DO ij = 1,jmax
        rho(ij) = rho(ij)/buoyancy%vector(2)
     ENDDO

! equilibrium from state (rho,T)
     CALL THERMO_POLYNOMIAL_PSAT(i1, jmax, i1, T, z1(1,2))

     DO ij = 1,jmax
        qsat = z1(ij,2)/(rho(ij)*T(ij)*WGHT_INV(1))

        IF ( qsat .GE. z1(ij,1) ) THEN
           z1(ij,2) = 0
        ELSE
           z1(ij,2) = z1(ij,1)-qsat
        ENDIF
        IF ( dsmooth .GT. C_0_R ) THEN
           dsmooth_loc = dsmooth*qsat
           z1(ij,2) = C_05_R*dsmooth_loc&
                *LOG(EXP(C_2_R*(z1(ij,1)-qsat)/dsmooth_loc)+C_1_R)
        ENDIF

     ENDDO

  ENDDO

  RETURN
END SUBROUTINE FI_HYDROSTATIC_AIRWATER_T

!########################################################################
!# Computes the INVERSE of the scale height of the system at position Y
!########################################################################
!# ARGUMENTS 
!# p(1)   In   Value of the grid point y(1)
!#
!# In case of AIRWATER mixture:
!# p(2)   In   Value of grid point y_l   
!# p(3)   In   Value of grid point y_r
!# p(4)   In   Value of function_1 at point y_l
!# p(5)   In   Value of function_1 at point y_r
!# p(6)   In   Value of function_2 at point y_l
!# p(7)   In   Value of function_2 at point y_r
!#
!########################################################################
FUNCTION FI_HYDROSTATIC_SCALEHEIGHT_INV(y,p)

  USE DNS_CONSTANTS, ONLY : MAX_NSP
  USE DNS_GLOBAL, ONLY : imode_flow
  USE DNS_GLOBAL, ONLY : iprof_tem, mean_tem, delta_tem, thick_tem, ycoor_tem, prof_tem
  USE DNS_GLOBAL, ONLY : iprof_i, mean_i, delta_i, thick_i, ycoor_i, prof_i
  USE DNS_GLOBAL, ONLY : inb_scal, scaley
  USE DNS_GLOBAL, ONLY : buoyancy
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TREAL y, p(*), FI_HYDROSTATIC_SCALEHEIGHT_INV

! -------------------------------------------------------------------
  TREAL t_loc, y_i_loc(MAX_NSP+1), ycenter, FLOW_SHEAR_TEMPORAL
  TREAL r1, p_loc, h_loc !, dummy
  TINTEGER is, iprof_loc
  EXTERNAL FLOW_SHEAR_TEMPORAL

! ###################################################################
  r1 = C_1_R

  IF ( imode_flow .EQ. DNS_FLOW_SHEAR ) THEN
     ycenter = p(1) + scaley*ycoor_tem
     t_loc = FLOW_SHEAR_TEMPORAL&
          (iprof_tem, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y)
     DO is = 1,inb_scal
        ycenter = p(1) + scaley*ycoor_i(is)
        y_i_loc(is) =  FLOW_SHEAR_TEMPORAL&
             (iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i, y)
     ENDDO
  ENDIF

! AIRWATER case
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. iprof_tem .GT. 0 ) THEN
     y_i_loc(2) = p(4)+(p(5)-p(4))/(p(3)-p(2))*(y-p(2))
  ENDIF
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. iprof_tem .LT. 0 ) THEN
     ycenter = p(1) + scaley*ycoor_tem
     iprof_loc =-iprof_tem
     h_loc = FLOW_SHEAR_TEMPORAL&
          (iprof_loc, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y)
     p_loc = p(4)+(p(5)-p(4))/(p(3)-p(2))*(y-p(2))
!     CALL THERMO_AIRWATER_PH(i1, i1, i1, y_i_loc, p_loc, h_loc)
     ! CALL THERMO_AIRWATER_TEMPERATURE(i1, i1, i1, y_i_loc, h_loc, t_loc)
     CALL THERMO_AIRWATER_PH2(i1, i1, i1, y_i_loc, p_loc, h_loc, t_loc)
  ENDIF

! Setting the pressure entry to 1, THERMO_RHO give the thermodynamic part
! of the inverse of the scale height W/(R_0 T)
  r1 = C_1_R
  CALL THERMO_THERMAL_DENSITY(i1, i1, i1, y_i_loc, r1, t_loc, FI_HYDROSTATIC_SCALEHEIGHT_INV)

! adding the volumetric force part
  FI_HYDROSTATIC_SCALEHEIGHT_INV = buoyancy%vector(2)*FI_HYDROSTATIC_SCALEHEIGHT_INV

  RETURN
END FUNCTION FI_HYDROSTATIC_SCALEHEIGHT_INV
