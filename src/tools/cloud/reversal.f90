PROGRAM REVERSAL

#include "types.h"
#include "dns_const.h"

  USE DNS_GLOBAL
  USE THERMO_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL qt_1, qt_2, h_1, h_2, qt, h, x, qsat, dqldqt
  TREAL qvqd, ba_ratio, heat1, heat2, alpha, dummy
  TREAL z1(2), e, rho, p, T, ep, s(3)
  TREAL t_1, t_2, c0, c1, c2
  TREAL r_1, r_2, r_max, r_old, x_c, x_max
  TINTEGER n, nmax, iopt, iup

! ###################################################################
  CALL DNS_START

  imixture = MIXT_TYPE_AIRWATER
  CALL THERMO_INITIALIZE
  MRATIO = C_1_R
  dsmooth = C_0_R

  WRITE(*,*) '1 - Density profile from nondimensional state'
  WRITE(*,*) '2 - Density profile from dimensional state'
  WRITE(*,*) '3 - Coefficient as (enthalpy) table'
  WRITE(*,*) '4 - Coefficient bs (water content) table'
  WRITE(*,*) '5 - Coefficients ratio table'
  WRITE(*,*) '6 - Correction temperature table'
  WRITE(*,*) '7 - Correction temperature table (negative branch)'

  READ(*,*) iopt

  IF ( iopt .EQ. 1 ) THEN
     WRITE(*,*) 'Initial enthalpy, h_1 ?'
     READ(*,*) h_1
     WRITE(*,*) 'Final enthalpy, h_2 ?'
     READ(*,*) h_2
     WRITE(*,*) 'Initial water, qt_1 ?'
     READ(*,*) qt_1
     WRITE(*,*) 'Final water, qt_2 ?'
     READ(*,*) qt_2

     WRITE(*,*) 'pressure (bar) ?'
     READ(*,*) p

  ELSE IF ( iopt .EQ. 2 ) THEN
     WRITE(*,*) 'Lower temperature (C) ?'
     READ(*,*) t_1
     t_1 = (t_1+273.15)/TREF
     WRITE(*,*) 'Upper temperature (C) ?'
     READ(*,*) t_2
     t_2 = (t_2+273.15)/TREF
     WRITE(*,*) 'Lower water content (g/kg) ?'
     READ(*,*) qt_1
     qt_1 = qt_1*C_1EM3_R
     WRITE(*,*) 'Upper water content (g/kg) ?'
     READ(*,*) qt_2
     qt_2 = qt_2*C_1EM3_R

     WRITE(*,*) 'pressure (bar) ?'
     READ(*,*) p

! calculate nondimensional limits
     z1(1) = qt_1
     CALL THERMO_AIRWATER_PT(i1, i1, i1, z1, p, t_1)
     CALL THERMO_CALORIC_ENTHALPY(i1, i1, i1, z1, t_1, h_1)
     CALL THERMO_THERMAL_DENSITY(i1, i1, i1, z1, p, t_1, r_1)

     z1(1) = qt_2
     CALL THERMO_AIRWATER_PT(i1, i1, i1, z1, p, t_2)
     CALL THERMO_CALORIC_ENTHALPY(i1, i1, i1, z1, t_2, h_2)
     CALL THERMO_THERMAL_DENSITY(i1, i1, i1, z1, p, t_2, r_2)

  ELSE IF ( iopt .EQ. 3 .OR. iopt .EQ. 4 .OR. iopt .EQ. 5 .OR. iopt .EQ. 6 .OR. iopt .EQ. 7 ) THEN
     WRITE(*,*) 'pressure (bar) ?'
     READ(*,*) p

     WRITE(*,*) 'Minimum temperature (C) ?'
     READ(*,*) t_1
     t_1 = (t_1+273.15)/TREF
     WRITE(*,*) 'Maximum temperature (C) ?'
     READ(*,*) t_2
     t_2 = (t_2+273.15)/TREF

     IF ( iopt .EQ. 3 ) THEN
        WRITE(*,*) 'Coefficient -as ?'
        READ(*,*) ba_ratio
     ELSE IF ( iopt .EQ. 4 ) THEN
        WRITE(*,*) 'Coefficient -bs ?'
        READ(*,*) ba_ratio
     ELSE IF ( iopt .EQ. 5 ) THEN
        WRITE(*,*) 'Coefficients ratio ?'
        READ(*,*) ba_ratio
     ELSE
        WRITE(*,*) 'Correction temperature difference (K) ?'
        READ(*,*) ba_ratio
        ba_ratio=ba_ratio/TREF
     ENDIF

  ENDIF

  WRITE(*,*) 'number of points ?'
  READ(*,*) nmax

  OPEN(21,file='reversal.dat')
! ###################################################################
  IF ( iopt .EQ. 1 .OR. iopt .EQ. 2 ) THEN
     WRITE(21,*) '# x, qt, h, ql, qv, qsat(T), r, T, p, e'
     r_max = C_0_R
     r_old = r_1
     x_c   =-C_1_R
     iup   = 0
     DO n = 1,nmax
        x = M_REAL(n-1)/M_REAL(nmax-1)

        qt = qt_1 + x*(qt_2-qt_1)
        h  = h_1  + x*(h_2-h_1)
        ep = C_0_R

        z1(1) = qt
        CALL THERMO_AIRWATER_PH(i1, i1, i1, z1, h, ep,p)
        s(1) = h; s(2:3) = z1(1:2)
        CALL THERMO_ANELASTIC_TEMPERATURE(i1, i1, i1, s, ep,T)
        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, T, qsat)
        qsat = C_1_R/(MRATIO*p/qsat-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qsat/(C_1_R+qsat)
        CALL THERMO_THERMAL_DENSITY(i1, i1, i1, z1, p, T, rho)
        CALL THERMO_CALORIC_ENERGY(i1, i1, i1, z1, T, e)

        WRITE(21,1010) x, qt, h, z1(2), qt-z1(2), qsat, rho, T, p, e

        IF ( (rho-r_old) .GT. 0 .AND. iup .EQ. 0 ) iup=1
        IF ( rho .LT. r_1 .AND. iup .EQ. 1 .AND. x_c .LT. C_0_R ) x_c=x
        IF ( rho .GT. r_max ) THEN
           r_max = rho
           x_max = x
        ENDIF

     ENDDO

     WRITE(*,*) 'Maximum density  r_max .................:',  r_max
     WRITE(*,*) 'Maximum density diference r_max-r_1 ....:',  r_max-r_1
     WRITE(*,*) 'Ratio (r_max-r_1)/(r_1-r_2) ............:', (r_max-r_1)/(r_1-r_2)
     WRITE(*,*) 'Maximum mixing fraction ................:',  x_max
     WRITE(*,*) 'Cross-over mixing fraction .............:',  x_c

! ###################################################################
  ELSE IF ( iopt .EQ. 3 ) THEN
     WRITE(21,*) '# T (C), T/298K, qt (g/kg)'

     DO n = 1,nmax
        t = t_1 + (t_2-t_1)*M_REAL(n-1)/M_REAL(nmax-1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t
        heat2 =  heat1*( C_1_R + qvqd ) -&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t

        alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

        dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))*&
             qvqd*(C_1_R+qvqd*WGHT_INV(1)/WGHT_INV(2))
        dummy = dummy +&
             THERMO_AI(1,1,2) + qvqd*(THERMO_AI(1,1,1)-THERMO_AI(1,1,3)) - THERMO_AI(1,1,3)

        dummy = (alpha/(ba_ratio*t) - THERMO_AI(1,1,3))/dummy

        qt = C_1_R-dummy

        IF ( qt .LT. qsat ) EXIT

        WRITE(21,*) t*TREF-273.15, t, qt*1000

     ENDDO

! look for crossing point with saturation curve using bracketing
     dummy = (t_2-t_1)/M_REAL(nmax-1)
     t_2 = t
     t_1 = t - dummy
     DO WHILE ( (t_2-t_1)/t_1 .GT. C_1EM6_R )

        t = C_05_R*(t_2+t_1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t
        heat2 =  heat1*( C_1_R + qvqd ) -&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t

        alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

        dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))*&
             qvqd*(C_1_R+qvqd*WGHT_INV(1)/WGHT_INV(2))
        dummy = dummy +&
             THERMO_AI(1,1,2) + qvqd*(THERMO_AI(1,1,1)-THERMO_AI(1,1,3)) - THERMO_AI(1,1,3)

        dummy = (alpha/(ba_ratio*t) - THERMO_AI(1,1,3))/dummy

        qt = C_1_R-dummy

        IF ( qt .GT. qsat ) t_1 = t
        IF ( qt .LE. qsat ) t_2 = t

     ENDDO

     WRITE(21,*) t*TREF-273.15, t, qt*1000

! ###################################################################
  ELSE IF ( iopt .EQ. 4 ) THEN
     WRITE(21,*) '# T (C), T/298K, qt (g/kg)'

     DO n = 1,nmax
        t = t_1 + (t_2-t_1)*M_REAL(n-1)/M_REAL(nmax-1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t
        heat2 =  heat1*( C_1_R + qvqd ) -&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t

        alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

        dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))*&
             qvqd*(C_1_R+qvqd*WGHT_INV(1)/WGHT_INV(2))
        dummy = dummy +&
             THERMO_AI(1,1,2) + qvqd*(THERMO_AI(1,1,1)-THERMO_AI(1,1,3)) - THERMO_AI(1,1,3)

        c2 = ba_ratio*dummy
        c1 =-( dummy*(C_1_R+ba_ratio) + (dummy+THERMO_AI(1,1,3))*ba_ratio - alpha*heat2/t )
        c0 = (C_1_R+ba_ratio)*(dummy+THERMO_AI(1,1,3)) - alpha*heat2/t

        qt = (-c1 + SQRT(c1*c1-C_4_R*c0*c2)) / (C_2_R*c2)

        IF ( (-c1 - SQRT(c1*c1-C_4_R*c0*c2)) / (C_2_R*c2) .LT. C_1_R ) THEN
           PRINT*,'error, second root also less than 1'
        ENDIF

        IF ( qt .LT. qsat ) EXIT

        WRITE(21,*) t*TREF-273.15, t, qt*1000

     ENDDO

! look for crossing point with saturation curve using bracketing
     dummy = (t_2-t_1)/M_REAL(nmax-1)
     t_2 = t
     t_1 = t - dummy
     DO WHILE ( (t_2-t_1)/t_1 .GT. C_1EM6_R )

        t = C_05_R*(t_2+t_1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t
        heat2 =  heat1*( C_1_R + qvqd ) -&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t

        alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

        dummy = heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))*&
             qvqd*(C_1_R+qvqd*WGHT_INV(1)/WGHT_INV(2))
        dummy = dummy +&
             THERMO_AI(1,1,2) + qvqd*(THERMO_AI(1,1,1)-THERMO_AI(1,1,3)) - THERMO_AI(1,1,3)

        c2 = ba_ratio*dummy
        c1 =-( dummy*(C_1_R+ba_ratio) + (dummy+THERMO_AI(1,1,3))*ba_ratio - alpha*heat2/t )
        c0 = (C_1_R+ba_ratio)*(dummy+THERMO_AI(1,1,3)) - alpha*heat2/t

        qt = (-c1 + SQRT(c1*c1-C_4_R*c0*c2)) / (C_2_R*c2)

        IF ( (-c1 - SQRT(c1*c1-C_4_R*c0*c2)) / (C_2_R*c2) .LT. C_1_R ) THEN
           PRINT*,'error, second root also less than 1'
        ENDIF

        IF ( qt .GT. qsat ) t_1 = t
        IF ( qt .LE. qsat ) t_2 = t

     ENDDO

     WRITE(21,*) t*TREF-273.15, t, qt*1000

! ###################################################################
  ELSE IF ( iopt .EQ. 5 ) THEN
     WRITE(21,*) '# T (C), T/298K, qt (g/kg)'

     DO n = 1,nmax
        t = t_1 + (t_2-t_1)*M_REAL(n-1)/M_REAL(nmax-1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t
        heat2 =  heat1*( C_1_R + qvqd ) -&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t

        alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

        dummy = (heat2-ba_ratio)/t*alpha
        dummy = dummy -&
             heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))*&
             qvqd*(C_1_R+qvqd*WGHT_INV(1)/WGHT_INV(2))
        dummy = dummy - &
             THERMO_AI(1,1,2) - qvqd*(THERMO_AI(1,1,1)-THERMO_AI(1,1,3))
        dummy = dummy/THERMO_AI(1,1,3)

        qt = dummy/(C_1_R+dummy)

        IF ( qt .LT. qsat ) EXIT

        WRITE(21,*) t*TREF-273.15, t, qt*1000

     ENDDO

! look for crossing point with saturation curve using bracketing
     dummy = (t_2-t_1)/M_REAL(nmax-1)
     t_2 = t
     t_1 = t - dummy
     DO WHILE ( (t_2-t_1)/t_1 .GT. C_1EM6_R )

        t = C_05_R*(t_2+t_1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t
        heat2 =  heat1*( C_1_R + qvqd ) -&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t

        alpha = C_1_R + heat1*qvqd/(GRATIO*WGHT_INV(2)*t)

        dummy = (heat2-ba_ratio)/t*alpha
        dummy = dummy -&
             heat1*heat1/(GRATIO*(t**2)*WGHT_INV(1))*&
             qvqd*(C_1_R+qvqd*WGHT_INV(1)/WGHT_INV(2))
        dummy = dummy - &
             THERMO_AI(1,1,2) - qvqd*(THERMO_AI(1,1,1)-THERMO_AI(1,1,3))
        dummy = dummy/THERMO_AI(1,1,3)

        qt = dummy/(C_1_R+dummy)

        IF ( qt .GT. qsat ) t_1 = t
        IF ( qt .LE. qsat ) t_2 = t

     ENDDO

     WRITE(21,*) t*TREF-273.15, t, qt*1000

! ###################################################################
  ELSE IF ( iopt .EQ. 6 ) THEN
     WRITE(21,*) '# T (C), T/298K, qt (g/kg)'

     DO n = 1,nmax
        t = t_1 + (t_2-t_1)*M_REAL(n-1)/M_REAL(nmax-1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t

        qt = (qsat*heat1 + ba_ratio)/(heat1-(THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t)

        IF ( qt .GT. qsat ) WRITE(21,*) t*TREF-273.15, t, qt*1000

     ENDDO

! ###################################################################
  ELSE IF ( iopt .EQ. 7 ) THEN
     WRITE(21,*) '# T (C), T/298K, qt (g/kg)'

     DO n = 1,nmax
        t = t_1 + (t_2-t_1)*M_REAL(n-1)/M_REAL(nmax-1)

        CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t, qvqd)
        qvqd = C_1_R/(MRATIO*p/qvqd-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        qsat = qvqd/(C_1_R+qvqd)

        heat1 = THERMO_AI(6,1,1) - THERMO_AI(6,1,3) +&
             (THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t

        qt = -ba_ratio/((THERMO_AI(1,1,1)-THERMO_AI(1,1,2))*t)

        IF ( qt .LT. qsat ) WRITE(21,*) t*TREF-273.15, t, qt*1000

     ENDDO
  ENDIF

  CLOSE(21)

  CALL DNS_STOP

  STOP
1010 FORMAT(10(1X,G_FORMAT_R))
END PROGRAM REVERSAL
