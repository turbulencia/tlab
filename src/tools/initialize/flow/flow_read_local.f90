#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE FLOW_READ_LOCAL(inifile)

  USE DNS_CONSTANTS, ONLY : efile, lfile, wfile
  USE TLAB_VARS,    ONLY : qbg
  USE TLAB_PROCS
  USE TLAB_TYPES,     ONLY : MAX_MODES
  USE FLOW_LOCAL

  IMPLICIT NONE

  TINTEGER :: bcs_flow_jmin, bcs_flow_jmax ! Boundary conditions

  CHARACTER*(*) inifile

  ! -------------------------------------------------------------------
  TREAL dummy(1)
  TINTEGER idummy
  CHARACTER*512 sRes
  CHARACTER*32 bakfile

  ! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL TLAB_WRITE_ASCII(lfile, 'Reading local input data')

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile,'#')
  CALL TLAB_WRITE_ASCII(bakfile,'#[IniFields]')
  CALL TLAB_WRITE_ASCII(bakfile,'#Velocity=<option>')
  CALL TLAB_WRITE_ASCII(bakfile,'#Temperature=<option>')
  CALL TLAB_WRITE_ASCII(bakfile,'#ForceDilatation=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile,'#ProfileIniK=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,'#ThickIniK=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,'#YCoorIniK=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,'#NormalizeK=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,'#NormalizeP=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,'#Mixture=<string>')

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Velocity', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'              ) THEN; flag_u = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'velocitydiscrete'  ) THEN; flag_u = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'velocitybroadband' ) THEN; flag_u = 2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'vorticitybroadband') THEN; flag_u = 3
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'potentialbroadband') THEN; flag_u = 4
  ELSE
    CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Velocity forcing type unknown')
    CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ForceDilatation', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'no' ) THEN; flag_dilatation=0
  ELSE;                                      flag_dilatation=1; ENDIF

  Kini = qbg(1) ! default geometry and scaling of perturbation

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ProfileIniK', 'GaussianSurface', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'             ) THEN; Kini%type = PROFILE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussian'         ) THEN; Kini%type = PROFILE_GAUSSIAN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussianvaricose' ) THEN; Kini%type = PROFILE_GAUSSIAN_ANTISYM
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussiansinuous'  ) THEN; Kini%type = PROFILE_GAUSSIAN_SYM
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussiansurface'  ) THEN; Kini%type = PROFILE_GAUSSIAN_SURFACE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'parabolicsurface' ) THEN; Kini%type = PROFILE_PARABOLIC_SURFACE
  ELSE
    CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Wrong ProfileIni parameter.')
    CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIniK', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) & ! backwards compatilibity
    CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIni', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .NE. 'void' ) THEN
    dummy(1) = C_1_R; idummy = 1
    CALL LIST_REAL(sRes, idummy, dummy)
    Kini%thick = dummy(1)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIniK', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) & ! backwards compatilibity
    CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIni', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .NE. 'void' ) THEN
    dummy(1) = C_1_R; idummy = 1
    CALL LIST_REAL(sRes, idummy, dummy)
    Kini%ymean = dummy(1)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeK', '-1.0', norm_ini_u)
  CALL SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeP', '-1.0', norm_ini_p)

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Temperature', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'           ) THEN; flag_t = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planebroadband' ) THEN; flag_t = 4
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planediscrete'  ) THEN; flag_t = 5
  ELSE
    CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Temperature forcing type unknown')
    CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  ! Additional parameters
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Mixture', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'        ) THEN; flag_mixture = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'equilibrium' ) THEN; flag_mixture = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'loadfields'  ) THEN; flag_mixture = 2
  ENDIF

  ! Boundary conditions
  flag_wall = 0
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmin', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_jmin = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_jmin = DNS_BCS_DIRICHLET; flag_wall=flag_wall+1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_jmin = DNS_BCS_NEUMANN
  ELSE
    CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. BoundaryConditions.VelocityJmin.')
    CALL TLAB_STOP(DNS_ERROR_IBC)
  ENDIF
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmax', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_jmax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_jmax = DNS_BCS_DIRICHLET; flag_wall=flag_wall+2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_jmax = DNS_BCS_NEUMANN
  ELSE
    CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. BoundaryConditions.VelocityJmax.')
    CALL TLAB_STOP(DNS_ERROR_IBC)
  ENDIF

  ! ###################################################################
  ! Discrete Forcing
  ! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Discrete]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Aplitude=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ModeX=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ModeZ=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#PhaseX=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#PhaseZ=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Type=<Gaussian>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Broadening=<value>')

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'Amplitude', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) & ! backwards compatilibity
  CALL SCANINICHAR(bakfile, inifile, 'Discrete', '2DAmpl', '0.0', sRes)
  fp%amplitude(:)=C_0_R; fp%size = MAX_MODES
  CALL LIST_REAL(sRes, fp%size, fp%amplitude)

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'ModeX', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN ! Default
    DO idummy = 1,fp%size; fp%modex(idummy) = idummy; ENDDO
  ELSE
    idummy = MAX_MODES
    CALL LIST_INTEGER(sRes, idummy, fp%modex)
    IF ( idummy .NE. fp%size ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeX.')
      CALL TLAB_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'ModeZ', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN ! Default
    fp%modez = 0
  ELSE
    idummy = MAX_MODES
    CALL LIST_INTEGER(sRes, idummy, fp%modez)
    IF ( idummy .NE. fp%size ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeZ.')
      CALL TLAB_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'PhaseX', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) & ! backwards compatilibity
  CALL SCANINICHAR(bakfile, inifile, 'Discrete', '2DPhi', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN ! Default
    fp%phasex = 0
  ELSE
    idummy = MAX_MODES
    CALL LIST_REAL(sRes, idummy, fp%phasex)
    IF ( idummy .NE. fp%size ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseX.')
      CALL TLAB_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'PhaseZ', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN ! Default
    fp%phasez = 0
  ELSE
    idummy = MAX_MODES
    CALL LIST_REAL(sRes, idummy, fp%phasez)
    IF ( idummy .NE. fp%size ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseZ.')
      CALL TLAB_STOP(DNS_ERROR_INFDISCR)
    ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'Type', 'None', sRes) ! Modulation type
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; fp%type = 0
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'gaussian' ) THEN; fp%type = 1
  ELSE
    CALL TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Error in Discrete.Type.')
    CALL TLAB_STOP(DNS_ERROR_INFDISCR)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', fp%parameters(1))

  RETURN
END SUBROUTINE FLOW_READ_LOCAL
