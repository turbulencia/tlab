PROGRAM SATURATION
  
#include "types.h"
#include "dns_const.h"

  USE DNS_GLOBAL
  USE THERMO_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL t_min, t_max, t_del, t, psat, qsat, dummy, t_loc, p, l
  TINTEGER iopt, ipsat

! ###################################################################
  CALL DNS_INITIALIZE

  imixture = MIXT_TYPE_AIRWATER
  CALL THERMO_INITIALIZE
  MRATIO = C_1_R

  WRITE(*,*) '1 - Saturation pressure as a function of T'
  WRITE(*,*) '2 - Saturation specific humidity as a function of T-p'
  READ(*,*) iopt

  WRITE(*,*) 'Minimum T (cetigrade) ?'
  READ(*,*) t_min
  WRITE(*,*) 'Maximum T (centigrade) ?'
  READ(*,*) t_max
  WRITE(*,*) 'Increment T ?'
  READ(*,*) t_del

  IF ( iopt .EQ. 2 ) THEN
     WRITE(*,*) 'Pressure (bar) ?'
     READ(*,*) p
  ENDIF

! ###################################################################
  OPEN(21,file='vapor.dat')
  IF ( iopt .EQ. 1 ) THEN
     WRITE(21,*) '# T (C), T (K), psat (bar), L-ps (J/kg), L-cp (J/kg)'
  ELSE IF ( iopt .EQ. 2 ) THEN
     WRITE(21,*) '# T (C), T (K), qsat (g/kg)'
  ENDIF

  t = t_min
  DO WHILE ( t .LE. t_max ) 

     t_loc = (t+273.15)/TREF
     CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t_loc, psat)
     dummy = C_1_R/(MRATIO*p/psat-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qsat = dummy/(C_1_R+dummy)
     IF ( iopt .EQ. 1 ) THEN
        l = C_0_R
        DO ipsat = NPSAT,2,-1
           l = l*t_loc + THERMO_PSAT(ipsat)*M_REAL(ipsat-1)
        ENDDO        
        WRITE(21,*) t, t_loc*TREF, psat, l*t_loc**2/psat*WGHT_INV(1)*(RGAS*TREF/WREF), &
             ((THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t_loc+THERMO_AI(6,1,1)-THERMO_AI(6,1,3))*(RGAS*TREF/WREF)/GRATIO
     ELSE IF ( iopt .EQ. 2 ) THEN
        WRITE(21,*) t, t_loc, qsat*1.d3
     ENDIF

     t = t+t_del
  ENDDO

  CLOSE(21)

  CALL DNS_STOP

  STOP
END PROGRAM SATURATION
