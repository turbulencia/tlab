#include "types.h"
#include "dns_error.h"

SUBROUTINE RAND_READ_LOCAL(inifile)

  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE RAND_LOCAL

  IMPLICIT NONE

#include "integers.h"

  CHARACTER*(*) inifile

! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*32 bakfile

  TINTEGER :: idummy
  TREAL    :: rdummy(6)

! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL IO_WRITE_ASCII(lfile, 'Reading local input data')

! ###################################################################
  CALL IO_WRITE_ASCII(bakfile,'#')
  CALL IO_WRITE_ASCII(bakfile,'#[Broadband]')
  CALL IO_WRITE_ASCII(bakfile,'#Type=<None/Physical/Phase>')
  CALL IO_WRITE_ASCII(bakfile,'#Seed=<random seed>')
  CALL IO_WRITE_ASCII(bakfile,'#Spectrum=<none/uniform/quartic/quadratic/gaussian>')
  CALL IO_WRITE_ASCII(bakfile,'#f0=<frequencies>')
  CALL IO_WRITE_ASCII(bakfile,'#Distribution=<uniform/gaussian>')
  CALL IO_WRITE_ASCII(bakfile,'#Covariance=<Rxx,Ryy,Rzz,Rxy,Rxz,Ryz>')

  CALL SCANINICHAR(bakfile, inifile, 'Broadband', 'Type', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'     ) THEN; flag_type = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'physical' ) THEN; flag_type = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'phase'    ) THEN; flag_type = 2
  ELSE
     CALL IO_WRITE_ASCII(efile, 'RAND_READ_LOCAL. Broadband: Type unknown.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINIINT(bakfile, inifile, 'Broadband', 'Seed', '0', seed)

! Spectral density
  CALL SCANINICHAR(bakfile, inifile, 'Broadband', 'Spectrum', 'quartic', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; ispectrum = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'uniform'   ) THEN; ispectrum = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'quartic'   ) THEN; ispectrum = 3
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'quadratic' ) THEN; ispectrum = 4
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'gaussian'  ) THEN; ispectrum = 6; ENDIF

  spc_param = C_0_R
  CALL SCANINICHAR(bakfile, inifile, 'Broadband', 'f0',    '1.0', sRes)
  idummy = i3
  CALL LIST_REAL(sRes, idummy, spc_param)

  CALL SCANINIREAL(bakfile, inifile, 'Broadband', 'Sigma', '1.0', spc_param(4))

! distribution function
  CALL SCANINICHAR(bakfile, inifile, 'Broadband', 'Distribution', 'uniform', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'uniform'  ) THEN; ipdf = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'gaussian' ) THEN; ipdf = 2
  ELSE
     CALL IO_WRITE_ASCII(efile, 'RAND_READ_LOCAL. Broadband: Distribution type unknown.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! covariance matrix: Rxx,Ryy,Rzz, Rxy,Rxz,Ryz
  ucov(1:3) = C_1_R ! diagonal terms
  ucov(4:6) = C_0_R ! off-diagonal terms
  CALL SCANINICHAR(bakfile, inifile, 'Broadband', 'Covariance', '-1', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .NE. '-1' ) THEN
     idummy = i6
     CALL LIST_REAL(sRes, idummy, rdummy)
     IF ( idummy .EQ. 6 ) THEN; ucov(1:6) = rdummy(1:6)
     ELSE
        CALL IO_WRITE_ASCII(efile, 'RAND_READ_LOCAL. Broadband: Incorrect number of variances.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF
  CALL IO_WRITE_ASCII(bakfile,'Covariance matrix:')
  WRITE(sRes,'(3E11.4)') ucov(1), ucov(4), ucov(5); CALL IO_WRITE_ASCII(bakfile,sRes) 
  WRITE(sRes,'(3E11.4)') ucov(4), ucov(2), ucov(6); CALL IO_WRITE_ASCII(bakfile,sRes) 
  WRITE(sRes,'(3E11.4)') ucov(5), ucov(6), ucov(3); CALL IO_WRITE_ASCII(bakfile,sRes) 

! impose a symmetric field
  CALL SCANINICHAR(bakfile, inifile, 'Broadband', 'Symmetric', 'no', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'yes' ) THEN; isymmetric = 1
  ELSE;                                       isymmetric = 0; ENDIF

  RETURN
END SUBROUTINE RAND_READ_LOCAL
