#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

SUBROUTINE SCAL_READ_LOCAL(inifile)

  USE DNS_CONSTANTS, ONLY : efile, lfile, wfile, MAX_NSP
  USE DNS_GLOBAL,    ONLY : inb_scal
  USE DNS_GLOBAL,    ONLY : sbg
  USE DNS_TYPES,     ONLY : MAX_MODES
  USE SCAL_LOCAL

  IMPLICIT NONE

  CHARACTER*(*) inifile

! -------------------------------------------------------------------
  TREAL dummy(MAX_NSP)
  TINTEGER idummy
  CHARACTER*512 sRes
  CHARACTER*32 bakfile

! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL IO_WRITE_ASCII(lfile, 'Reading local input data')

! ###################################################################
! Scalar fluctuation field
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile,'#')
  CALL IO_WRITE_ASCII(bakfile,'#[IniFields]')
  CALL IO_WRITE_ASCII(bakfile,'#Scalar=<option>')
  CALL IO_WRITE_ASCII(bakfile,'#ThickIniS=<values>')
  CALL IO_WRITE_ASCII(bakfile,'#ProfileIniS=<values>')
  CALL IO_WRITE_ASCII(bakfile,'#YCoorIniS=<values>')
  CALL IO_WRITE_ASCII(bakfile,'#NormalizeS=<values>')
  CALL IO_WRITE_ASCII(bakfile,'#Mixture=<string>')

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Scalar', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'           ) THEN; flag_s = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'layerbroadband' ) THEN; flag_s = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'layerdiscrete'  ) THEN; flag_s = 2
!  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'both'           ) THEN; flag_s = 3
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planebroadband' ) THEN; flag_s = 4
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planediscrete'  ) THEN; flag_s = 5
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'deltabroadband' ) THEN; flag_s = 6
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'deltadiscrete'  ) THEN; flag_s = 7
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'fluxbroadband'  ) THEN; flag_s = 8
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'fluxdiscrete'   ) THEN; flag_s = 9; ENDIF

! Geometry and scaling of perturbation
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ProfileIniS', 'GaussianSurface', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'                  ) THEN; Sini(:)%type = PROFILE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussian'              ) THEN; Sini(:)%type = PROFILE_GAUSSIAN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussiansurface'       ) THEN; Sini(:)%type = PROFILE_GAUSSIAN_SURFACE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussiansymmetric'     ) THEN; Sini(:)%type = PROFILE_GAUSSIAN_SYM
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussianantisymmetric' ) THEN; Sini(:)%type = PROFILE_GAUSSIAN_ANTISYM
  ELSE
    CALL IO_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Wrong ProfileIni parameter.')
    CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIniS', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) &    ! backwards compatilibity
       CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIni', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN; Sini(:)%thick = sbg(:)%thick;
  ELSE
     dummy = C_0_R; idummy = MAX_NSP
     CALL LIST_REAL(sRes, idummy, dummy)
     Sini(:)%thick = dummy(:)
     IF ( idummy .NE. inb_scal ) THEN         ! Consistency check
        IF ( idummy .EQ. 1 ) THEN
           Sini(2:)%thick = Sini(1)%thick
           CALL IO_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using ThickIniS(1) for all scalars.')
        ELSE
           CALL IO_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. ThickIniS size does not match number of scalars.')
           CALL DNS_STOP(DNS_ERROR_OPTION)
        ENDIF
     ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIniS', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) &    ! backwards compatilibity
       CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIni', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN; Sini(:)%ymean = sbg(:)%ymean;
  ELSE
     dummy = C_0_R; idummy = MAX_NSP
     CALL LIST_REAL(sRes, idummy, dummy)
     Sini(:)%ymean = dummy(:)
     IF ( idummy .NE. inb_scal ) THEN         ! Consistency check
        IF ( idummy .EQ. 1 ) THEN
           Sini(2:)%ymean = Sini(1)%ymean
           CALL IO_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using YCoorIniS(1) for all scalars.')
        ELSE
           CALL IO_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. YCoorIniS size does not match number of scalars.')
           CALL DNS_STOP(DNS_ERROR_OPTION)
        ENDIF
     ENDIF
  ENDIF

  DO idummy = 1,inb_scal; Sini(idummy)%parameters(:) = C_0_R; ENDDO

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'NormalizeS', '1.0', sRes)
  norm_ini_s(:) = C_0_R; idummy = MAX_NSP
  CALL LIST_REAL(sRes, idummy, norm_ini_s)
  IF ( idummy .NE. inb_scal ) THEN            ! Consistency check
     IF ( idummy .EQ. 1 ) THEN
        norm_ini_s(2:) = norm_ini_s(1)
        CALL IO_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using NormalizeS(1) for all scalars.')
     ELSE
        CALL IO_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. NormalizeS size does not match number of scalars.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

  CALL SCANINIREAL(bakfile,inifile,'IniFields', 'NormalizeR', '0.0', norm_ini_radiation) ! Radiation field

! Additional parameters
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Mixture', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'        ) THEN; flag_mixture = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'equilibrium' ) THEN; flag_mixture = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'loadfields'  ) THEN; flag_mixture = 2
  ENDIF

! ###################################################################
! Discrete Forcing
! ###################################################################
CALL IO_WRITE_ASCII(bakfile, '#')
CALL IO_WRITE_ASCII(bakfile, '#[Discrete]')
CALL IO_WRITE_ASCII(bakfile, '#Type=<Varicose/Sinuous/Gaussian/Step>')
CALL IO_WRITE_ASCII(bakfile, '#Aplitude=<value>')
CALL IO_WRITE_ASCII(bakfile, '#ModeX=<value>')
CALL IO_WRITE_ASCII(bakfile, '#ModeZ=<value>')
CALL IO_WRITE_ASCII(bakfile, '#PhaseX=<value>')
CALL IO_WRITE_ASCII(bakfile, '#PhaseZ=<value>')
CALL IO_WRITE_ASCII(bakfile, '#Broadening=<value>')

CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'Amplitude', 'void', sRes)
IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) &        ! backwards compatilibity
CALL SCANINICHAR(bakfile, inifile, 'Discrete', '2DAmpl', '0.0', sRes)
fp%amplitude(:)=C_0_R; fp%size = MAX_MODES
CALL LIST_REAL(sRes, fp%size, fp%amplitude)

CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'ModeX', 'void', sRes)
IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN     ! Default
  DO idummy = 1,fp%size; fp%modex(idummy) = idummy; ENDDO
  ELSE
    idummy = MAX_MODES
    CALL LIST_INTEGER(sRes, idummy, fp%modex)
    IF ( idummy .NE. fp%size ) THEN
      CALL IO_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeX.')
      CALL DNS_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'ModeZ', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN     ! Default
    fp%modez = 0
  ELSE
    idummy = MAX_MODES
    CALL LIST_INTEGER(sRes, idummy, fp%modez)
    IF ( idummy .NE. fp%size ) THEN
      CALL IO_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeZ.')
      CALL DNS_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'PhaseX', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) &        ! backwards compatilibity
  CALL SCANINICHAR(bakfile, inifile, 'Discrete', '2DPhi', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN     ! Default
    fp%phasex = 0
  ELSE
    idummy = MAX_MODES
    CALL LIST_REAL(sRes, idummy, fp%phasex)
    IF ( idummy .NE. fp%size ) THEN
      CALL IO_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseX.')
      CALL DNS_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'PhaseZ', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN     ! Default
    fp%phasez = 0
  ELSE
    idummy = MAX_MODES
    CALL LIST_REAL(sRes, idummy, fp%phasez)
    IF ( idummy .NE. fp%size ) THEN
      CALL IO_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseZ.')
      CALL DNS_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'Type', 'none', sRes) ! Modulation type
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; fp%type = 0
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'gaussian' ) THEN; fp%type = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'step' )     THEN; fp%type = 2
  ELSE
    CALL IO_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Error in Discrete.Type.')
    CALL DNS_STOP(DNS_ERROR_INFDISCR)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', fp%parameters(1))
  CALL SCANINIREAL(bakfile, inifile, 'Discrete', 'ThickStep',  '-1.0', fp%parameters(2))

  RETURN
END SUBROUTINE SCAL_READ_LOCAL
