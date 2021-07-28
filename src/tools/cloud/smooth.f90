PROGRAM SMOOTH

#include "types.h"
#include "dns_const.h"

  USE TLAB_VARS
  USE TLAB_PROCS
  USE THERMO_VARS

  IMPLICIT NONE

#include "integers.h"

  TREAL qt_min, qt_max, qt_del, qt, qs, dqldqt
  TREAL z1(2), e, rho, p, T, h, ep, s(3)
  TINTEGER opt

! ###################################################################
  CALL TLAB_START

  imixture = MIXT_TYPE_AIRWATER
  CALL THERMO_INITIALIZE
  MRATIO = C_1_R
  IF ( gama0 .GT. C_0_R ) GRATIO = (gama0-C_1_R)/gama0
  ep = C_0_R

  WRITE(*,*) 'Case d-e (1) or d-p (2) or p-h (3) ?'
  READ(*,*) opt

  WRITE(*,*) 'Minimum qt ?'
  READ(*,*) qt_min
  WRITE(*,*) 'Maximum qt ?'
  READ(*,*) qt_max
  WRITE(*,*) 'Increment qt ?'
  READ(*,*) qt_del

  IF ( opt .EQ. 1 ) THEN
     WRITE(*,*) 'Density value ?'
     READ(*,*) rho
     WRITE(*,*) 'Energy value ?'
     READ(*,*) e
  ELSE IF ( opt .EQ. 2 ) THEN
     WRITE(*,*) 'Density value ?'
     READ(*,*) rho
     WRITE(*,*) 'Pressure value ?'
     READ(*,*) p
  ELSE IF ( opt .EQ. 3 ) THEN
     WRITE(*,*) 'enthalpy?'
     READ(*,*) h
     WRITE(*,*) 'pressure?'
     READ(*,*) p
  ENDIF

  WRITE(*,*) 'Smoothing factor ?'
  READ(*,*) dsmooth

! ###################################################################
  OPEN(21,file='vapor.dat')
  WRITE(21,*) '# qt, ql, qv, qs(T), r, T, p, e, h'

  qt = qt_min
  DO WHILE ( qt .LE. qt_max )

     z1(1) = qt
     IF ( opt .EQ. 1 ) THEN
        CALL THERMO_CALORIC_TEMPERATURE(i1, i1, i1, z1, e, rho, T, dqldqt)
        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, T, qs)
        qs = qs/(rho*T*WGHT_INV(1))
        CALL THERMO_THERMAL_PRESSURE(i1, i1, i1, z1, rho, T, p)
        CALL THERMO_CALORIC_ENTHALPY(i1, i1, i1, z1, T, h)

     ELSE IF ( opt .EQ. 2 ) THEN
        CALL THERMO_AIRWATER_RP(i1, i1, i1, z1, p, rho, T, dqldqt)
        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, T, qs)
        qs = qs/(rho*T*WGHT_INV(1))
        CALL THERMO_CALORIC_ENERGY(i1, i1, i1, z1, T, e)
        CALL THERMO_CALORIC_ENTHALPY(i1, i1, i1, z1, T, h)

     ELSE IF ( opt .EQ. 3 ) THEN
        CALL THERMO_AIRWATER_PH(i1,i1,i1, z1,h, ep,p)
        s(1) = h; s(2:3) = z1(1:2)
        CALL THERMO_ANELASTIC_TEMPERATURE(i1,i1,i1, s, ep, T)
!        CALL THERMO_AIRWATER_PH_RE(i1,i1,i1, z1, p, h, T)
        CALL THERMO_POLYNOMIAL_PSAT(i1,i1,i1, T, qs)
        qs = C_1_R/(MRATIO*p/qs-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qs = qs/(C_1_R+qs)
        CALL THERMO_THERMAL_DENSITY(i1,i1,i1, z1, p, T, rho)
        CALL THERMO_CALORIC_ENERGY(i1,i1,i1, z1, T, e)

     ENDIF
     WRITE(21,1000) qt, z1(2), qt-z1(2), qs, rho, T, p, e, h

     qt = qt+qt_del
  ENDDO

  CLOSE(21)

  STOP

1000 FORMAT(9(G_FORMAT_R))

END PROGRAM SMOOTH
