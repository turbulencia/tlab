#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

SUBROUTINE FLOW_READ_LOCAL(inifile)

  USE DNS_CONSTANTS, ONLY : efile, lfile, wfile
  USE DNS_GLOBAL,    ONLY : ycoor_u, thick_u
  USE FLOW_LOCAL

  IMPLICIT NONE

! Boundary conditions
  TINTEGER :: bcs_flow_jmin, bcs_flow_jmax

#include "integers.h"

  CHARACTER*(*) inifile

! -------------------------------------------------------------------
  TINTEGER idummy
  CHARACTER*512 sRes
  CHARACTER*32 bakfile

! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL IO_WRITE_ASCII(lfile, 'Reading local input data')

! ###################################################################
  CALL IO_WRITE_ASCII(bakfile,'#')
  CALL IO_WRITE_ASCII(bakfile,'#[IniFields]')
  CALL IO_WRITE_ASCII(bakfile,'#Velocity=<option>')
  CALL IO_WRITE_ASCII(bakfile,'#Temperature=<option>')
  CALL IO_WRITE_ASCII(bakfile,'#ForceDilatation=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile,'#ThickIniK=<value>')
  CALL IO_WRITE_ASCII(bakfile,'#YCoorIniK=<value>')
  CALL IO_WRITE_ASCII(bakfile,'#NormalizeK=<value>')
  CALL IO_WRITE_ASCII(bakfile,'#NormalizeP=<value>')
  CALL IO_WRITE_ASCII(bakfile,'#Mixture=<string>')

! velocity
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Velocity', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'              ) THEN; flag_u = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'velocitydiscrete'  ) THEN; flag_u = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'velocitybroadband' ) THEN; flag_u = 2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'vorticitybroadband') THEN; flag_u = 3
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'potentialbroadband') THEN; flag_u = 4
  ELSE 
     CALL IO_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Velocity forcing type unknown')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ForceDilatation', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .eq. 'no' ) THEN; flag_dilatation=0
  ELSE;                                      flag_dilatation=1; ENDIF

! Geometry and scaling of perturbation
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIniK', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) & ! backwards compatilibity
       CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIni', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN; thick_ini = thick_u;
  ELSE
     thick_ini = C_1_R; idummy = 1
     CALL LIST_REAL(sRes, idummy, thick_ini)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIniK', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) & ! backwards compatilibity
       CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIni', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'void' ) THEN; ycoor_ini = ycoor_u;
  ELSE
     ycoor_ini = C_1_R; idummy = 1
     CALL LIST_REAL(sRes, idummy, ycoor_ini)
  ENDIF
  
  CALL SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeK', '1.0', norm_ini_u)
  CALL SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeP', '1.0', norm_ini_p)

! Temperature
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Temperature', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'           ) THEN; flag_t = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planebroadband' ) THEN; flag_t = 4
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planediscrete'  ) THEN; flag_t = 5
  ELSE 
     CALL IO_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Temperature forcing type unknown')
     CALL DNS_STOP(DNS_ERROR_OPTION)
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
     CALL IO_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. BoundaryConditions.VelocityJmin.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF
  CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmax', 'freeslip', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'     ) THEN; bcs_flow_jmax = DNS_BCS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; bcs_flow_jmax = DNS_BCS_DIRICHLET; flag_wall=flag_wall+2
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; bcs_flow_jmax = DNS_BCS_NEUMANN
  ELSE
     CALL IO_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. BoundaryConditions.VelocityJmax.')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF

! ###################################################################
! Discrete Forcing
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Discrete]')
  CALL IO_WRITE_ASCII(bakfile, '#Type=<Varicose/Sinuous/Gaussian>')
  CALL IO_WRITE_ASCII(bakfile, '#2DAmpl=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#3DAmpl=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#2DPhi=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#3DXPhi=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#3DZPhi=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Broadening=<value>')

  CALL SCANINICHAR(bakfile,inifile,'Discrete','2DPhi', '0.0',sRes)
  Phix2D(:)=C_0_R; nx2d = MAX_FRC_FREC
  CALL LIST_REAL(sRes, nx2d, Phix2D)
  CALL SCANINICHAR(bakfile,inifile,'Discrete','2DAmpl','0.0',sRes)
  A2D(:)=C_0_R; nx2d = MAX_FRC_FREC ! The amplitude sets the value of nx2d
  CALL LIST_REAL(sRes, nx2d, A2D)
  
  CALL SCANINICHAR(bakfile,inifile,'Discrete','3DXPhi','0.0',sRes)
  Phix3D(:)=C_0_R; nx3d = MAX_FRC_FREC
  CALL LIST_REAL(sRes, nx3d, Phix3d)
  CALL SCANINICHAR(bakfile,inifile,'Discrete','3DZPhi','0.0',sRes)
  Phiz3D(:)=C_0_R; nz3d = MAX_FRC_FREC
  CALL LIST_REAL(sRes, nz3d, Phiz3D)
  CALL SCANINICHAR(bakfile,inifile,'Discrete','3DAmpl','0.0',sRes)
  A3D(:)=C_0_R; nx3d = MAX_FRC_FREC ! The amplitude sets the value of nx3d
  CALL LIST_REAL(sRes, nx3d, A3D)

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'Type', 'Varicose', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'varicose' ) THEN; ifrcdsc_mode = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'sinuous'  ) THEN; ifrcdsc_mode = 2
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'gaussian' ) THEN; ifrcdsc_mode = 3
  ELSE
     CALL IO_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Error in Discrete.Type.')
     CALL DNS_STOP(DNS_ERROR_INFDISCR)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', frc_delta)

  RETURN
END SUBROUTINE FLOW_READ_LOCAL
