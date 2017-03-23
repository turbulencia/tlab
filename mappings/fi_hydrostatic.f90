#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2007/05/08 - J.P. Mellado
!#              Created
!# 2017/03/23 - J.P. Mellado
!#              Rewritten in terms of FDM operators
!#
!########################################################################
!# DESCRIPTION
!#
!# Evaluate the integral \int_ycenter^y dx/H(x), where H(x) is the scale 
!# height in the system, where ycenter is some given value
!#
!########################################################################

!########################################################################
! Compute hydrostatic equilibrium from profiles s=(h,q_t).
!########################################################################
SUBROUTINE FI_HYDROSTATIC_H(g, s, e, T,p, wrk1d)

  USE DNS_TYPES, ONLY : grid_structure

  USE DNS_GLOBAL, ONLY : imode_eqns
  USE DNS_GLOBAL, ONLY : p_init, ycoor_p, p_scale_height, damkohler, buoyancy
  USE THERMO_GLOBAL, ONLY : imixture

  IMPLICIT NONE

#include "integers.h"

  TYPE(grid_structure),       INTENT(IN)    :: g
  TREAL, DIMENSION(g%size),   INTENT(IN)    :: e
  TREAL, DIMENSION(g%size),   INTENT(OUT)   :: T,p
  TREAL, DIMENSION(g%size,*), INTENT(INOUT) :: s      ! We calculate equilibrium composition
  TREAL, DIMENSION(g%size,*), INTENT(INOUT) :: wrk1d

! -------------------------------------------------------------------
  TINTEGER iter, niter, ibc, j, jcenter
  TREAL dummy, ycenter
  
! ###################################################################
! Get the center  
  ycenter = g%nodes(1) + ycoor_p *g%scale
  DO j = 1,g%size
     IF ( g%nodes(j  ) .LE. ycenter .AND. &
          g%nodes(j+1) .GT. ycenter ) THEN
        jcenter = j
        EXIT
     ENDIF
  ENDDO

! Setting the pressure entry to 1 to get 1/RT
  wrk1d(:,6) = C_1_R
  
! Prepare the pentadiagonal system
  ibc = 1     ! Boundary condition at the bottom for integral calulation 
  CALL INT_C1N6_LHS(g%size, ibc, wrk1d(1,1),wrk1d(1,2),wrk1d(1,3),wrk1d(1,4),wrk1d(1,5))
  CALL PENTADFS(g%size-1,        wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5))

  niter = 10

  p(:) = p_init        ! initialize iteration
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Get the composition, if necessary
     s(:,3) = C_0_R
  ENDIF
  DO iter = 1,niter    ! iterate
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
! Setting the pressure entry to 1 to get 1/RT.
        CALL THERMO_ANELASTIC_DENSITY(i1,g%size,i1, s, e,wrk1d(1,6), wrk1d(1,7))
        dummy = C_1_R / SIGN(p_scale_height,buoyancy%vector(2))
     ELSE
        CALL THERMO_AIRWATER_PH_RE(i1,g%size, i1, s(1,2), p, s(1,1), T)
! Setting the pressure entry to 1 to get 1/RT.
        CALL THERMO_THERMAL_DENSITY(i1,g%size,i1, s(1,2),wrk1d(1,6),T, wrk1d(1,7))
        dummy = buoyancy%vector(2)
     ENDIF
     wrk1d(:,7) = dummy *wrk1d(:,7)

! Calculate integral
     CALL INT_C1N6_RHS(g%size,i1, ibc, g%jac, wrk1d(1,7), p)
     CALL PENTADSS(g%size-1,i1, wrk1d(2,1),wrk1d(2,2),wrk1d(2,3),wrk1d(2,4),wrk1d(2,5), p(2))
     p(1) = C_0_R

! Calculate pressure and normalize s.t. p=p_init at y=ycoor_p
     p(:) = EXP(p(:))
     IF ( ABS(ycenter-g%nodes(jcenter)) .EQ. C_0_R ) THEN
        dummy = p(jcenter)
     ELSE
        dummy = p(jcenter) + (p(jcenter+1)      -p(jcenter)      ) &
                           / (g%nodes(jcenter+1)-g%nodes(jcenter)) *(ycenter-g%nodes(jcenter))
     ENDIF
     dummy = p_init /dummy     
     p(:)  = dummy *p(:)

     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN ! Get the composition, if necessary
        IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
           CALL THERMO_AIRWATER_PH(i1,g%size,i1, s(1,2),s(1,1), e,p)
        ELSE
           CALL THERMO_AIRWATER_PH_RE(i1,g%size,i1, s(1,2), p, s(1,1), T)
        ENDIF        
     ENDIF
     
  ENDDO

! compute equilibrium values of T
  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     CALL THERMO_ANELASTIC_TEMPERATURE(i1,g%size,i1, s, e, T)
  ENDIF
  
  RETURN
END SUBROUTINE FI_HYDROSTATIC_H

! ! ###################################################################
! ! TO BE REMOVED
! ! ###################################################################

! SUBROUTINE FI_HYDROSTATIC(imax,jmax,kmax, ycenter, y, integral)

!   IMPLICIT NONE

!   TINTEGER imax, jmax, kmax
!   TREAL ycenter
!   TREAL y(jmax)
!   TREAL integral(imax,jmax,kmax)

! ! -------------------------------------------------------------------
!   TINTEGER i, j, k, jcenter
!   TREAL FI_HYDROSTATIC_SCALEHEIGHT_INV
!   EXTERNAL FI_HYDROSTATIC_SCALEHEIGHT_INV

!   TREAL ABSERR,EPSABS,EPSREL,RESULT
!   TINTEGER IER,NEVAL
!   TREAL WORK, params(10)
!   TINTEGER IWORK,KEY,LAST,LENW,LIMIT
!   DIMENSION IWORK(100),WORK(400)

! ! ###################################################################
!   DO j = 1,jmax
!      IF ( y(j  ) .LE. ycenter .AND. &
!           y(j+1) .GT. ycenter ) THEN
!         jcenter = j
!         EXIT
!      ENDIF
!   ENDDO

! ! Convergence constants
!   EPSABS = C_1EM10_R
!   EPSREL = C_1EM10_R
!   KEY = 6
!   LIMIT = 100
!   LENW = LIMIT*4

!   params(1) = y(1)
!   params(2) = y(1)
!   params(3) = y(jmax)
!   params(4) = C_0_R
!   params(5) = C_0_R
!   params(6) = C_0_R
!   params(7) = C_0_R

!   DO j = jcenter+1,jmax
!      CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,ycenter,y(j),EPSABS,&
!           EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
!      DO k = 1,kmax
!         DO i = 1,imax
!            integral(i,j,k) = RESULT
!         ENDDO
!      ENDDO
!   ENDDO

!   DO j = jcenter,1,-1
!      CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(j),ycenter,EPSABS,&
!           EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
!      DO k = 1,kmax
!         DO i = 1,imax
!            integral(i,j,k) = -RESULT
!         ENDDO
!      ENDDO
!   ENDDO

!   RETURN
! END SUBROUTINE FI_HYDROSTATIC
      
! !########################################################################
! !# Modification from FI_HYDROSTATIC to allow passing a function through
! !# parameters; nvar variables enter through wrk1d
! !#
! !########################################################################
! SUBROUTINE FI_HYDROSTATIC_STEPS(imax,jmax,kmax, nvar, ycenter, y, integral, wrk1d)
  
!   IMPLICIT NONE

!   TINTEGER imax, jmax, kmax, nvar
!   TREAL ycenter
!   TREAL y(jmax)
!   TREAL integral(imax,jmax,kmax)
!   TREAL wrk1d(jmax,nvar)

! ! -------------------------------------------------------------------
!   TINTEGER i, j, k, iv, jcenter
!   TREAL FI_HYDROSTATIC_SCALEHEIGHT_INV
!   EXTERNAL FI_HYDROSTATIC_SCALEHEIGHT_INV

!   TREAL ABSERR,EPSABS,EPSREL,RESULT
!   TINTEGER IER,NEVAL
!   TREAL WORK, params(10)
!   TINTEGER IWORK,KEY,LAST,LENW,LIMIT
!   DIMENSION IWORK(100),WORK(400)

! ! ###################################################################
!   DO iv = 1,10
!      params(iv) = C_0_R
!   ENDDO

!   DO j = 1,jmax
!      IF ( y(j  ) .LE. ycenter .AND. &
!           y(j+1) .GT. ycenter ) THEN
!         jcenter = j
!         EXIT
!      ENDIF
!   ENDDO

! ! Convergence constants
!   EPSABS = C_1EM10_R
!   EPSREL = C_1EM10_R
!   KEY = 6
!   LIMIT = 100
!   LENW = LIMIT*4

!   params(1) = y(1)

! ! from center, upwards
!   params(2) = y(jcenter)
!   params(3) = y(jcenter+1)
!   DO iv = 1,nvar
!      params((iv-1)*2+4) = wrk1d(jcenter  ,iv)
!      params((iv-1)*2+5) = wrk1d(jcenter+1,iv)
!   ENDDO
!   CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,ycenter,y(jcenter+1),EPSABS,&
!        EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
!   DO k = 1,kmax
!      DO i = 1,imax
!         integral(i,jcenter+1,k) = RESULT
!      ENDDO
!   ENDDO
!   DO j = jcenter+2,jmax
!      params(2) = y(j-1)
!      params(3) = y(j)
!      DO iv = 1,nvar
!         params((iv-1)*2+4) = wrk1d(j-1,iv)
!         params((iv-1)*2+5) = wrk1d(j  ,iv)
!      ENDDO
!      CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(j-1),y(j),EPSABS,&
!           EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
!      DO k = 1,kmax
!         DO i = 1,imax
!            integral(i,j,k) = integral(i,j-1,k) + RESULT
!         ENDDO
!      ENDDO
!   ENDDO

! ! from center, downwards
!   params(2) = y(jcenter)
!   params(3) = y(jcenter+1)
!   DO iv = 1,nvar
!      params((iv-1)*2+4) = wrk1d(jcenter  ,iv)
!      params((iv-1)*2+5) = wrk1d(jcenter+1,iv)
!   ENDDO
!   CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(jcenter),ycenter,EPSABS,&
!        EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
!   DO k = 1,kmax
!      DO i = 1,imax
!         integral(i,jcenter,k) = -RESULT
!      ENDDO
!   ENDDO
!   DO j = jcenter-1,1,-1
!      params(2) = y(j)
!      params(3) = y(j+1)
!      DO iv = 1,nvar
!         params((iv-1)*2+4) = wrk1d(j  ,iv)
!         params((iv-1)*2+5) = wrk1d(j+1,iv)
!      ENDDO
!      CALL QUADAG(FI_HYDROSTATIC_SCALEHEIGHT_INV,params,y(j),y(j+1),EPSABS,&
!           EPSREL,KEY,RESULT,ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
!      DO k = 1,kmax
!         DO i = 1,imax
!            integral(i,j,k) = integral(i,j+1,k)-RESULT
!         ENDDO
!      ENDDO
!   ENDDO

!   RETURN
! END SUBROUTINE FI_HYDROSTATIC_STEPS
       
! !########################################################################
! ! Compute hydrostatic equilibrium from profiles (h,q_t).
! !########################################################################
! SUBROUTINE FI_HYDROSTATIC_H_OLD(nmax, y, s, e, T,p, wrk1d)

!   USE DNS_GLOBAL, ONLY : imode_eqns
!   USE DNS_GLOBAL, ONLY : p_init, scaley, damkohler, ycoor_p !ycoor_tem, ycoor_i
!   USE THERMO_GLOBAL, ONLY : imixture

!   IMPLICIT NONE

! #include "integers.h"

!   TINTEGER nmax
!   TREAL, DIMENSION(nmax)   :: y
!   TREAL, DIMENSION(nmax)   :: T, p, e
!   TREAL, DIMENSION(nmax,*) :: s, wrk1d

! ! -------------------------------------------------------------------
!   TINTEGER iter, niter
!   TREAL ycenter

! ! ###################################################################
!   niter = 10

!   ycenter = y(1) + ycoor_p*scaley      
!   IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN 
! !     ycenter = y(1) + ycoor_i(1)*scaley      
!   ELSE
! !     ycenter = y(1) + ycoor_tem*scaley
!      wrk1d(:,2) = e(:)
!   ENDIF

!   p(:) = p_init        ! initialize iteration
!   DO iter = 1,niter    ! iterate
!      wrk1d(:,1) = p(:) ! pressure from hydrostatic equilibrium
!      CALL FI_HYDROSTATIC_STEPS(i1,nmax,i1, i2, ycenter, y, p, wrk1d)
!      p(:) = p_init*EXP(p(:))

!   ENDDO

! ! compute equilibrium values of q_l and T
!   IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
!      IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R )  THEN
!         CALL THERMO_AIRWATER_PH(i1,nmax,i1, s(:,2),s(:,1), e,p)
!      ENDIF
!      CALL THERMO_ANELASTIC_TEMPERATURE(i1,nmax,i1, s, e, T)
!   ELSE
!      IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R )  THEN
!         CALL THERMO_AIRWATER_PH_RE(i1,nmax,i1, s(:,2), p, s(:,1), T)
!      ENDIF
!   ENDIF
  
!   RETURN
! END SUBROUTINE FI_HYDROSTATIC_H_OLD

! !########################################################################
! ! Compute hydrostatic equilibrium from profiles (T,q_t).
! !########################################################################
! SUBROUTINE FI_HYDROSTATIC_AIRWATER_T(y, dy, s, T, p, rho, wrk1d, wrk2d, wrk3d)

!   USE DNS_GLOBAL, ONLY : p_init, ycoor_tem, scaley
!   USE DNS_GLOBAL, ONLY : imode_fdm, jmax, j1bc
!   USE DNS_GLOBAL, ONLY : buoyancy

!   USE THERMO_GLOBAL, ONLY : dsmooth, WGHT_INV

!   IMPLICIT NONE

! #include "integers.h"

!   TREAL y(jmax), dy(jmax)
!   TREAL s(jmax,*), T(*), p(*), rho(*)
!   TREAL wrk1d(*), wrk2d(*), wrk3d(*)

! ! -------------------------------------------------------------------
!   TINTEGER ij, iter, niter
!   TREAL ycenter, qsat, dsmooth_loc

! ! ###################################################################
!   niter = 3

!   s(:,2) = C_0_R ! initialize iteration from condition q_l = 0

! ! ###################################################################
! ! iteration
! ! ###################################################################
!   DO iter = 1,niter
! ! pressure 
!      ycenter = y(1) + ycoor_tem*scaley
!      CALL FI_HYDROSTATIC_STEPS(i1, jmax, i1, i1, ycenter, y, p, s(1,2))
!      DO ij = 1,jmax
!         p(ij) = p_init*EXP(p(ij))
!      ENDDO

! ! density profile from pressure gradient
!      CALL PARTIAL_Y(imode_fdm, i1, jmax, i1, j1bc, dy, p, rho, i0, i0, wrk1d, wrk2d, wrk3d)
!      DO ij = 1,jmax
!         rho(ij) = rho(ij)/buoyancy%vector(2)
!      ENDDO

! ! equilibrium from state (rho,T)
!      CALL THERMO_POLYNOMIAL_PSAT(i1, jmax, i1, T, s(1,2))

!      DO ij = 1,jmax
!         qsat = s(ij,2)/(rho(ij)*T(ij)*WGHT_INV(1))

!         IF ( qsat .GE. s(ij,1) ) THEN
!            s(ij,2) = 0
!         ELSE
!            s(ij,2) = s(ij,1)-qsat
!         ENDIF
!         IF ( dsmooth .GT. C_0_R ) THEN
!            dsmooth_loc = dsmooth*qsat
!            s(ij,2) = C_05_R*dsmooth_loc&
!                 *LOG(EXP(C_2_R*(s(ij,1)-qsat)/dsmooth_loc)+C_1_R)
!         ENDIF

!      ENDDO

!   ENDDO

!   RETURN
! END SUBROUTINE FI_HYDROSTATIC_AIRWATER_T

! !########################################################################
! !# Computes the INVERSE of the scale height of the system at position Y
! !########################################################################
! !# ARGUMENTS 
! !# p(1)   In   Value of the grid point y(1)
! !#
! !# In case of AIRWATER mixture:
! !# p(2)   In   Value of grid point y_l   
! !# p(3)   In   Value of grid point y_r
! !# p(4)   In   Value of function_1 at point y_l
! !# p(5)   In   Value of function_1 at point y_r
! !# p(6)   In   Value of function_2 at point y_l
! !# p(7)   In   Value of function_2 at point y_r
! !#
! !########################################################################
! FUNCTION FI_HYDROSTATIC_SCALEHEIGHT_INV(y,p)

!   USE DNS_CONSTANTS, ONLY : MAX_NSP
!   USE DNS_GLOBAL, ONLY : imode_eqns
!   USE DNS_GLOBAL, ONLY : iprof_tem, mean_tem, delta_tem, thick_tem, ycoor_tem, prof_tem
!   USE DNS_GLOBAL, ONLY : p_scale_height
!   USE DNS_GLOBAL, ONLY : iprof_i, mean_i, delta_i, thick_i, ycoor_i, prof_i
!   USE DNS_GLOBAL, ONLY : inb_scal, scaley, damkohler
!   USE DNS_GLOBAL, ONLY : buoyancy
!   USE THERMO_GLOBAL, ONLY : imixture

!   IMPLICIT NONE

! #include "integers.h"

!   TREAL y, p(*), FI_HYDROSTATIC_SCALEHEIGHT_INV

! ! -------------------------------------------------------------------
!   TREAL t_loc, y_i_loc(MAX_NSP+1), ycenter, FLOW_SHEAR_TEMPORAL
!   TREAL r1, e_loc, p_loc, h_loc
!   TINTEGER is, iprof_loc
!   EXTERNAL FLOW_SHEAR_TEMPORAL

! ! ###################################################################
!   r1 = C_1_R

!   p_loc = p(4)+(p(5)-p(4))/(p(3)-p(2))*(y-p(2))
!   e_loc = p(6)+(p(7)-p(6))/(p(3)-p(2))*(y-p(2))

!   DO is = 1,inb_scal
!      ycenter = p(1) + scaley*ycoor_i(is)
!      y_i_loc(is) =  FLOW_SHEAR_TEMPORAL&
!           (iprof_i(is), thick_i(is), delta_i(is), mean_i(is), ycenter, prof_i, y)
!   ENDDO
  
!   IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN 
!      IF      ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R )  THEN
!         CALL THERMO_AIRWATER_PH(i1,i1,i1, y_i_loc(2),y_i_loc(1), e_loc,p_loc)
!      ENDIF        
! ! Setting the pressure entry to 1, THERMO_RHO gives the thermodynamic part
! ! of the inverse of the scale height W/(R_0 T)
!      CALL THERMO_ANELASTIC_DENSITY(i1,i1,i1, y_i_loc(1), e_loc,r1, FI_HYDROSTATIC_SCALEHEIGHT_INV)
!      FI_HYDROSTATIC_SCALEHEIGHT_INV = SIGN(FI_HYDROSTATIC_SCALEHEIGHT_INV,buoyancy%vector(2)) /p_scale_height

!   ELSE     
!      ycenter = p(1) + scaley*ycoor_tem
!      t_loc = FLOW_SHEAR_TEMPORAL&
!           (iprof_tem, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y)

!      IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. iprof_tem .GT. 0 ) THEN
!         y_i_loc(2) = p_loc ! p(4)+(p(5)-p(4))/(p(3)-p(2))*(y-p(2))
!      ENDIF
!      IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. iprof_tem .LT. 0 ) THEN
!         ycenter = p(1) + scaley*ycoor_tem
!         iprof_loc =-iprof_tem
!         h_loc = FLOW_SHEAR_TEMPORAL&
!              (iprof_loc, thick_tem, delta_tem, mean_tem, ycenter, prof_tem, y)
!         CALL THERMO_AIRWATER_PH_RE(i1, i1, i1, y_i_loc, p_loc, h_loc, t_loc)
!      ENDIF

! ! Setting the pressure entry to 1, THERMO_RHO gives the thermodynamic part
! ! of the inverse of the scale height W/(R_0 T)
!      CALL THERMO_THERMAL_DENSITY(i1, i1, i1, y_i_loc, r1, t_loc, FI_HYDROSTATIC_SCALEHEIGHT_INV)

!      FI_HYDROSTATIC_SCALEHEIGHT_INV = buoyancy%vector(2)*FI_HYDROSTATIC_SCALEHEIGHT_INV
     
!   ENDIF

!   RETURN
! END FUNCTION FI_HYDROSTATIC_SCALEHEIGHT_INV

