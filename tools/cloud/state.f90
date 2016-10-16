PROGRAM STATE
  
#include "types.h"
#include "dns_const.h"

  USE DNS_GLOBAL
  USE THERMO_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL p, t, qs, qv, qt, ql, r, e, h, z1(2), dummy, dqldqt
  TREAL heat1, heat2, cp1, cp2, alpha, as, bs, t_eq, l, cp_ref
  TREAL r1,r2,r3, h2
  TINTEGER iopt

! ###################################################################
  CALL DNS_INITIALIZE

  imixture = MIXT_TYPE_AIRWATER
  CALL THERMO_INITIALIZE
  MRATIO = C_1_R
  dsmooth = C_0_R

  WRITE(*,*) 'Case p-t (1) or d-e (2) or p-h (3)?'
  READ(*,*) iopt

  IF ( iopt .EQ. 1 ) THEN
     WRITE(*,*) 'temperature (C) ?'
     READ(*,*) t
     t = (t+273.15)/TREF
     WRITE(*,*) 'pressure (bar) ?'
     READ(*,*) p

  ELSE IF ( iopt .EQ. 2 ) THEN
     WRITE(*,*) 'density ?'
     READ(*,*) r
     WRITE(*,*) 'energy ?'
     READ(*,*) e

  ELSE IF ( iopt .EQ. 3 ) THEN
     WRITE(*,*) 'enthapy ?'
     READ(*,*) h
     WRITE(*,*) 'pressure (bar) ?'
     READ(*,*) p

  ENDIF

  WRITE(*,*) 'water specific humidity (g/kg) ?'
  READ(*,*) qt
  qt = qt*C_1EM3_R

! ###################################################################
  IF ( iopt .EQ. 1 ) THEN
     CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, dummy)
     dummy = C_1_R/(MRATIO*p/dummy-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qs = dummy/(C_1_R+dummy)
     IF ( qt .GT. qs ) THEN
        qv = dummy*(1-qt)
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

  ELSE IF ( iopt .EQ. 2 ) THEN
     z1(1) = qt
     CALL THERMO_CALORIC_TEMPERATURE(i1, i1, i1, z1, e, r, T, dqldqt)
     ql = z1(2)
     qv = qt - ql
     qs = qv ! initial condition for next routine
!     CALL THERMO_AIRWATER_QSAT_RE(i1, i1, i1, e, r, qs, dummy)
     CALL THERMO_THERMAL_PRESSURE(i1, i1, i1, z1, r, t, p)
     CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, dummy)
     dummy = C_1_R/(MRATIO*p/dummy-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qs = dummy/(C_1_R+dummy)
     CALL THERMO_CALORIC_ENTHALPY(i1, i1, i1, z1, t, h)

  ELSE IF ( iopt .EQ. 3 ) THEN
     h = h/TREF/1.007
     z1(1) = qt
     CALL THERMO_AIRWATER_PH2(i1, i1, i1, z1, p, h, T)

     ql = z1(2)
     qv = qt - ql
     qs = qv ! initial condition for next routine
  
     CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, T, dummy)
     WRITE(*,*) 'Saturation pressure................:', dummy
     dummy = C_1_R/(MRATIO*p/dummy-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qs = dummy/(C_1_R+dummy)
 
     r1 = p/(T*(1- qt +WGHT_INV(1)/WGHT_INV(2)*qv ) )
     CALL THERMO_THERMAL_DENSITY(i1,i1,i1,z1,p,T, r2)
     CALL THERMO_CALORIC_ENTHALPY(i1,i1,i1,z1,T,h2)
     CALL THERMO_AIRWATER_DENSITY(i1,i1,i1, z1,p,h2,r3)

     WRITE(*,*) 'Latent heat..............................', (THERMO_AI(6,1,3))*1.007*TREF 
  ENDIF

  WRITE(*,*) 'Saturation specific humidity ......:', qs
  WRITE(*,*) 'Vapor specific humidity ...........:', qv
  WRITE(*,*) 'Liquid specific humidity ..........:', ql
  WRITE(*,*) 'Density ...........................:', r
  WRITE(*,*) 'Pressure ..........................:', p
  WRITE(*,*) 'Temperature .......................:', t*TREF - 273.15 ! 273.16
  WRITE(*,*) 'Specific energy ...................:', e
  WRITE(*,*) 'Specific enthalpy .................:', h
  IF ( iopt .EQ. 3 ) THEN
     WRITE(*,*) 'Density ...........................:', r1,r2,r3
     WRITE(*,*) 'Specific enthalpy .................:', h*1.007*TREF,  h2*1.007*TREF
  ENDIF
  cp_ref = (1-qt)*THERMO_AI(1,1,2) + qt*THERMO_AI(1,1,3)
  l      = THERMO_AI(6,1,1)-THERMO_AI(6,1,3)
  t_eq   = t*(C_1_R/p)**((1-qt)*GRATIO*WGHT_INV(2)/cp_ref)
  t_eq   = t_eq * EXP(qv*l/cp_ref/t) 

  WRITE(*,*) 'Equivalent potential temperature ..:', t_eq*TREF

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
END PROGRAM STATE
