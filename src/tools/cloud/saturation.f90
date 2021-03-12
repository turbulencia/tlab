PROGRAM SATURATION
  
#include "types.h"
#include "dns_const.h"

  USE DNS_GLOBAL
  USE THERMO_GLOBAL

  IMPLICIT NONE

#include "integers.h"

  TREAL t_min, t_max, t_del, t, psat, qsat, dummy, t_loc, p, dpsat!, dpsat2
  TINTEGER iopt!, ipsat

! ###################################################################
  CALL DNS_INITIALIZE

  imixture = MIXT_TYPE_AIRWATER
  CALL THERMO_INITIALIZE
  MRATIO = C_1_R
  IF ( gama0 .GT. C_0_R ) GRATIO = (gama0-C_1_R)/gama0

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
     CALL THERMO_POLYNOMIAL_PSAT(i1, i1, i1, t_loc,  psat)
     CALL THERMO_POLYNOMIAL_DPSAT(i1, i1, i1, t_loc, dpsat)
     dummy = C_1_R/(MRATIO*p/psat-C_1_R)*WGHT_INV(2)/WGHT_INV(1)
     qsat = dummy/(C_1_R+dummy)
     IF ( iopt .EQ. 1 ) THEN
        ! dpsat2 = C_0_R
        ! DO ipsat = NPSAT,2,-1
        !    dpsat2 = dpsat2 *t_loc + THERMO_PSAT(ipsat)*M_REAL(ipsat-1)
        ! ENDDO
        ! PRINT*,dpsat-dpsat2
        WRITE(21,1000) t, t_loc*TREF, psat, dpsat*t_loc**2/psat*WGHT_INV(1)*(RGAS*TREF/WREF), &
             ((THERMO_AI(1,1,1)-THERMO_AI(1,1,3))*t_loc+THERMO_AI(6,1,1)-THERMO_AI(6,1,3))*(RGAS*TREF/WREF)/GRATIO
     ELSE IF ( iopt .EQ. 2 ) THEN
        WRITE(21,2000) t, t_loc, qsat*1.d3
     ENDIF

     t = t+t_del
  ENDDO

  CLOSE(21)

  STOP
  
1000 FORMAT(5(1X,G_FORMAT_R))
2000 FORMAT(3(1X,G_FORMAT_R))

END PROGRAM SATURATION
