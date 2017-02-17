#include "types.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!# logs_data01 State (>0 if error)
!#
!# logs_data02 Maximum CFL number
!# logs_data03 Maximum diffusion number
!# logs_data04 Maximum source number
!#
!# logs_data05 Minimum pressure
!# logs_data06 Maximum pressure
!# logs_data07 Minimum density
!# logs_data08 Maximum density
!# logs_data09 NEWTONRAPHSON_ERROR
!#
!# logs_data10 Minimum dilatation
!# logs_data11 Maximum dilatation
!#
!########################################################################
SUBROUTINE DNS_LOGS(iflag)

  USE DNS_CONSTANTS, ONLY : ofile
  USE DNS_GLOBAL,    ONLY : imode_eqns
  USE DNS_GLOBAL,    ONLY : itime, rtime, visc
  USE DNS_GLOBAL,    ONLY : damkohler
  USE DNS_LOCAL,     ONLY : logs_data, dtime

  USE THERMO_GLOBAL, ONLY : imixture
  USE THERMO_GLOBAL, ONLY : NEWTONRAPHSON_ERROR  

  IMPLICIT NONE

  TINTEGER iflag

! -----------------------------------------------------------------------
  TINTEGER ip
  CHARACTER*256 line1
  CHARACTER*256 line2

! #######################################################################
! Headers
! #######################################################################
  IF ( iflag .EQ. 1 ) THEN
     line1 = '#'; ip = 1
     line1 = line1(1:ip)//' '//' Itn.'; ip = ip + 1 + 7
     line1 = line1(1:ip)//' '//' time'; ip = ip + 1 + 13
     line1 = line1(1:ip)//' '//' dt';   ip = ip + 1 + 10
     line1 = line1(1:ip)//' '//' CFL#'; ip = ip + 1 + 10
     line1 = line1(1:ip)//' '//' D#';   ip = ip + 1 + 10
     line1 = line1(1:ip)//' '//' visc'; ip = ip + 1 + 10
     
#ifdef CHEMISTRY
     IF ( ireactive .NE. CHEM_NONE ) THEN
        line1 = line1(1:ip)//' '//' R#'; ip = ip + 1 + 10
     ENDIF
#endif
     
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. &
          imode_eqns .EQ. DNS_EQNS_ANELASTIC     )THEN
        line1 = line1(1:ip)//' '//' DilMin'; ip = ip + 1 + 13
        line1 = line1(1:ip)//' '//' DilMax'; ip = ip + 1 + 13
     ELSE
        line1 = line1(1:ip)//' '//' PMin'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' PMax'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' RMin'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' RMax'; ip = ip + 1 + 10
     ENDIF

     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN
        line1 = line1(1:ip)//' '//' NewtonRs'; ip = ip + 1 + 10
     ENDIF
     
     line1 = line1(1:ip-1)//'#'
     CALL IO_WRITE_ASCII(ofile, REPEAT('#',LEN_TRIM(line1)))
     CALL IO_WRITE_ASCII(ofile, TRIM(ADJUSTL(line1)))
     CALL IO_WRITE_ASCII(ofile, REPEAT('#',LEN_TRIM(line1)))
     
  ENDIF

! #######################################################################
! Construct line with data
! #######################################################################
  IF ( iflag .EQ. 2 ) THEN
     WRITE(line1,100) INT(logs_data(1)), itime, rtime, dtime, (logs_data(ip),ip=2,3), visc
100  FORMAT((1X,I1),(1X,I7),(1X,E13.6),4(1X,E10.3))
     
#ifdef CHEMISTRY
     IF ( ireactive .NE. CHEM_NONE ) THEN
        WRITE(line2,101) logs_data(4)
101     FORMAT(1(1X,E10.3))
        line1 = TRIM(line1)//TRIM(line2)

     ENDIF
#endif
     
     IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. &
          imode_eqns .EQ. DNS_EQNS_ANELASTIC     )THEN
        WRITE(line2,200) logs_data(10), logs_data(11)
200     FORMAT(2(1X,E13.6))
        line1 = TRIM(line1)//TRIM(line2)
        
     ELSE
        WRITE(line2,300) (logs_data(ip),ip=5,8)
300     FORMAT(4(1X,E10.3))
        line1 = TRIM(line1)//TRIM(line2)
        
     ENDIF
  
     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .AND. damkohler(3) .LE. C_0_R ) THEN
        WRITE(line2,400) NEWTONRAPHSON_ERROR
400     FORMAT(1(1X,E10.3))
        line1 = TRIM(line1)//TRIM(line2)
     ENDIF

! Output
     CALL IO_WRITE_ASCII(ofile, TRIM(ADJUSTL(line1)))

  ENDIF
     
  RETURN
END SUBROUTINE DNS_LOGS
