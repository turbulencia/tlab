#include "types.h"
#include "dns_const.h"

PROGRAM STATE
  
  USE DNS_GLOBAL
  USE THERMO_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL p, ps, t, qs, qv, qt, ql, r, e, h, z1(2), dummy, dqldqt, ep, theta, theta_e
  TREAL heat1, heat2, cp1, cp2, alpha, as, bs
  TREAL r1, h1, s(3)
  TREAL Cd, Cdv, Cvl
  TINTEGER iopt

! ###################################################################
  CALL DNS_INITIALIZE

  imixture = MIXT_TYPE_AIRWATER
  CALL THERMO_INITIALIZE
  MRATIO = C_1_R
  IF ( gama0 .GT. C_0_R ) GRATIO = (gama0-C_1_R)/gama0
  ep = C_0_R
  dsmooth = C_0_R
     
  Cd = THERMO_AI(1,1,2)
  Cdv= THERMO_AI(1,1,1) - THERMO_AI(1,1,2)
  Cvl= THERMO_AI(1,1,3) - THERMO_AI(1,1,1)

  WRITE(*,*) 'Case p-t (1) or d-e (2) or p-h (3)?'
  READ(*,*) iopt

  IF ( iopt .EQ. 1 ) THEN
     WRITE(*,*) 'temperature (K) ?'
     READ(*,*) t
     t=t/TREF !t = (t+273.15)/TREF
     WRITE(*,*) 'pressure (bar) ?'
     READ(*,*) p

  ELSE IF ( iopt .EQ. 2 ) THEN
     WRITE(*,*) 'density ?'
     READ(*,*) r
     WRITE(*,*) 'energy ?'
     READ(*,*) e

  ELSE IF ( iopt .EQ. 3 ) THEN
     WRITE(*,*) 'enthalpy (kJ/kg)?'
     READ(*,*) h
     WRITE(*,*) 'pressure (bar) ?'
     READ(*,*) p

  ENDIF

  WRITE(*,*) 'water specific humidity (g/kg) ?'
  READ(*,*) qt
  qt = qt*C_1EM3_R

! ###################################################################
  IF ( iopt .EQ. 1 ) THEN
     CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, ps)
     qs = C_1_R/(MRATIO*p/ps-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qs = qs/(C_1_R+qs)
     IF ( qt .GT. qs ) THEN
        qv = qs*(1-qt)
        ql = qt-qv
     ELSE
        qv = qt
        ql = C_0_R
     ENDIF
     z1(1) = qt
     z1(2) = ql
     CALL THERMO_CALORIC_ENTHALPY(i1, i1, i1, z1, t, h)
     CALL THERMO_CALORIC_ENERGY(i1, i1, i1, z1, t, e)
     CALL THERMO_THERMAL_DENSITY(i1, i1, i1, z1, p, t, r)

     s(1) = h; s(2:3) = z1(1:2)
     CALL THERMO_ANELASTIC_THETA(i1,i1,i1, s, ep,p, theta)
     CALL THERMO_ANELASTIC_THETA_E(i1,i1,i1, s, ep,p, theta_e)

  ELSE IF ( iopt .EQ. 2 ) THEN
     z1(1) = qt
     CALL THERMO_CALORIC_TEMPERATURE(i1, i1, i1, z1, e, r, T, dqldqt)
     ql = z1(2)
     qv = qt - ql
     qs = qv ! initial condition for next routine
!     CALL THERMO_CALORIC_QSAT(i1, i1, i1, e, r, qs, qs)
     CALL THERMO_THERMAL_PRESSURE(i1, i1, i1, z1, r, t, p)
     CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, ps)
     qs = C_1_R/(MRATIO*p/ps-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qs = qs/(C_1_R+qs)
     CALL THERMO_CALORIC_ENTHALPY(i1, i1, i1, z1, t, h)

     s(1) = h; s(2:3) = z1(1:2)
     CALL THERMO_ANELASTIC_THETA(i1,i1,i1, s, ep,p, theta)
     CALL THERMO_ANELASTIC_THETA_E(i1,i1,i1, s, ep,p, theta_e)

  ELSE IF ( iopt .EQ. 3 ) THEN
     h = h/TREF/1.007
     z1(1) = qt
     CALL THERMO_AIRWATER_PH(i1,i1,i1, z1,h, ep,p)
     s(1) = h; s(2:3) = z1(1:2)
     CALL THERMO_ANELASTIC_TEMPERATURE(i1,i1,i1, s, ep, T)
     ! CALL THERMO_AIRWATER_PH_RE(i1,i1,i1, z1, p, h, T)
     ql = z1(2)
     qv = qt - ql

     CALL THERMO_POLYNOMIAL_PSAT(i1,i1,i1, T, ps)
     qs = C_1_R/(MRATIO*p/ps-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qs = qs/(C_1_R+qs)
     CALL THERMO_THERMAL_DENSITY(i1,i1,i1, z1,p,T, r)
     CALL THERMO_CALORIC_ENERGY(i1,i1,i1, z1, T, e)
     CALL THERMO_ANELASTIC_THETA(i1,i1,i1, s, ep,p, theta)
     CALL THERMO_ANELASTIC_THETA_E(i1,i1,i1, s, ep,p, theta_e)

! check
     CALL THERMO_ANELASTIC_DENSITY(i1,i1,i1, s, ep,p, r1)
!     r2 = p/(T*(1- qt +WGHT_INV(1)/WGHT_INV(2)*qv ) )
     CALL THERMO_CALORIC_ENTHALPY(i1,i1,i1, z1,T,h1)

  ENDIF

  WRITE(*,1000) 'Saturation specific humidity ......:', qs
  WRITE(*,1000) 'Vapor specific humidity ...........:', qv
  WRITE(*,1000) 'Liquid specific humidity ..........:', ql
  WRITE(*,1000) 'Density ...........................:', r
  WRITE(*,1000) 'Pressure (bar) ....................:', p
  WRITE(*,1000) 'Saturation pressure (bar) .........:', ps
  WRITE(*,1000) 'Temperature (K) ...................:', t*TREF !- 273.15
  WRITE(*,1000) 'Specific heat capacity ............:', Cd + qt *Cdv + ql *Cvl
  WRITE(*,1000) 'Specific energy ...................:', e
  WRITE(*,1000) 'Specific enthalpy .................:', h
  WRITE(*,1000) 'Reference latent heat (kJ/kg) .....:', -THERMO_AI(6,1,3) *1.007 *TREF 
  WRITE(*,1000) 'Latent heat (kJ/kg) ...............:', (-THERMO_AI(6,1,3)-t*(THERMO_AI(1,1,3)-THERMO_AI(1,1,1)) ) *1.007 *TREF 
  WRITE(*,1000) 'Liquid-water potential T (K) ......:', theta*TREF
  WRITE(*,1000) 'Equivalent potential T (K) ........:', theta_e*TREF
  IF ( iopt .EQ. 3 ) THEN
     WRITE(*,1000) 'Density ...........................:', r1
     WRITE(*,1000) 'Specific enthalpy .................:', h1
  ENDIF
    
! ###################################################################
  WRITE(*,*) ' '
  WRITE(*,*) 'Calculate reversal linear coefficients (1-yes/0-no) ?'
  READ(*,*) iopt

  IF ( iopt .EQ. 1 .AND. ql .GT. C_0_R ) THEN
     heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
          (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t
     heat2 =  heat1*( C_1_R + qv/(C_1_R-qt) ) -&
          (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t

     cp1 = (C_1_R-qt)*THERMO_AI(1,1,2) + qv*THERMO_AI(1,1,1) + ql*THERMO_AI(1,1,3)
     dummy = (heat1**2)*qv/((t**2)*cp1*WGHT_INV(1)*GRATIO)
     cp2 = cp1*( C_1_R + dummy*( C_1_R + qv/(C_1_R-qt)*WGHT_INV(1)/WGHT_INV(2) ) )

     alpha = C_1_R + heat1*qv/((C_1_R-qt)*GRATIO*WGHT_INV(2)*t)

     as =-alpha/cp2/t
     bs = heat2*as + C_1_R/(C_1_R-qt)
     WRITE(*,*) 'Enthalpy coefficient ..........:', as
     WRITE(*,*) 'Water fraction coefficient ....:', bs

  ELSE IF ( iopt .EQ. 1 .AND. ql .EQ. C_0_R ) THEN
     cp1 = THERMO_AI(1,1,2) + qt*(THERMO_AI(1,1,1)-THERMO_AI(1,1,2))

     as =-C_1_R/cp1/t
     bs = (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))/cp1&
          - ( WGHT_INV(1)-WGHT_INV(2) )/( WGHT_INV(2) + qt*(WGHT_INV(1)-WGHT_INV(2)) )
     WRITE(*,*) 'Enthalpy coefficient ..........:', as
     WRITE(*,*) 'Water fraction coefficient ....:', bs
  ENDIF

  STOP

1000 FORMAT(A,G_FORMAT_R)

END PROGRAM STATE
