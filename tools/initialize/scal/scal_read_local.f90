#include "types.h"
#include "dns_error.h"

SUBROUTINE SCAL_READ_LOCAL(inifile)

  USE DNS_CONSTANTS, ONLY : efile, lfile, wfile
  USE DNS_GLOBAL, ONLY : inb_scal, ycoor_i
  USE SCAL_LOCAL

  IMPLICIT NONE

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
! Scalar fluctuation field
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile,'#')
  CALL IO_WRITE_ASCII(bakfile,'#[IniFields]')
  CALL IO_WRITE_ASCII(bakfile,'#Scalar=<option>')
  CALL IO_WRITE_ASCII(bakfile,'#ThickIniS=<value>')
  CALL IO_WRITE_ASCII(bakfile,'#NormalizeS=<value>')
  CALL IO_WRITE_ASCII(bakfile,'#YCoorIni=<Relative Y reference point>')
  CALL IO_WRITE_ASCII(bakfile,'#Mixture=<string>')

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Scalar', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'           ) THEN; flag_s = 0
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'layerdiscrete'  ) THEN; flag_s = 1
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'layerbroadband' ) THEN; flag_s = 2
!  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'both'           ) THEN; flag_s = 3
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planebroadband' ) THEN; flag_s = 4
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'planediscrete'  ) THEN; flag_s = 5
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'deltabroadband' ) THEN; flag_s = 6
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'deltadiscrete'  ) THEN; flag_s = 7
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'fluxbroadband'  ) THEN; flag_s = 8
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'fluxdiscrete'   ) THEN; flag_s = 9; ENDIF

! Geometry and scaling of perturbation
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIniS', '1.0', sRes)
  thick_ini(:) = C_0_R; idummy = MAX_NSP
  CALL LIST_REAL(sRes, idummy, thick_ini)
  IF ( idummy .NE. inb_scal ) THEN ! Consistency check 
     IF ( idummy .EQ. 1 ) THEN
        thick_ini(2:) = thick_ini(1)     
        CALL IO_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using ThickIniS(1) for all scalars.')
     ELSE
        CALL IO_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. ThickIniS size does not match number of scalars.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'NormalizeS', '1.0', sRes)
  norm_ini_s(:) = C_0_R; idummy = MAX_NSP
  CALL LIST_REAL(sRes, idummy, norm_ini_s)
  IF ( idummy .NE. inb_scal ) THEN ! Consistency check 
     IF ( idummy .EQ. 1 ) THEN
        norm_ini_s(2:) = norm_ini_s(1)     
        CALL IO_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using NormalizeS(1) for all scalars.')
     ELSE
        CALL IO_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. NormalizeS size does not match number of scalars.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIni', 'dummy', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'dummy' ) THEN; ycoor_ini = ycoor_i;
  ELSE
     ycoor_ini(:) = C_0_R; idummy = MAX_NSP
     CALL LIST_REAL(sRes, idummy, ycoor_ini)
     IF ( idummy .NE. inb_scal ) THEN ! Consistency check
        IF ( idummy .EQ. 1 ) THEN
           ycoor_ini(2:) = ycoor_ini(1)     
           CALL IO_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Using YCoorIni(1) for all scalars.')
        ELSE
           CALL IO_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. YCoorIni size does not match number of scalars.')
           CALL DNS_STOP(DNS_ERROR_OPTION)
        ENDIF
     ENDIF

  ENDIF
  
  CALL SCANINIREAL(bakfile,inifile,'IniFields', 'NormalizeR', '0.0', norm_ini_radiation) ! Radiation field

! For backwards compatibility; to be removed
  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'Broadening', 'dummy', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .NE. 'dummy' )  THEN
     CALL IO_WRITE_ASCII(wfile, 'SCAL_READ_LOCAL. Broadening obsolete, use ThickIni instead.')
     CALL SCANINIREAL(bakfile, inifile, 'IniFields', 'Broadening', '1.0', thick_ini(1))
     thick_ini(2:) = thick_ini(1)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIni', 'dummy', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .NE. 'dummy' )  THEN
     CALL IO_WRITE_ASCII(wfile, 'FLOW_READ_LOCAL: ThickIni obsolete, use ThickIniS instead.')
     CALL SCANINIREAL(bakfile, inifile, 'IniFields', 'ThickIni',  '1.0', thick_ini(1) )
     thick_ini(2:) = thick_ini(1)
  ENDIF

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

! An initial effect of radiation as an accumulation during a certain interval of time
! is considered by means of rad_ini
  CALL SCANINIREAL(bakfile,inifile,'Discrete','RadStart','0.0',rad_ini)

  CALL SCANINICHAR(bakfile, inifile, 'Discrete', 'Type', 'Varicose', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'varicose' ) THEN; imode_discrete = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'sinuous'  ) THEN; imode_discrete = 2
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'gaussian' ) THEN; imode_discrete = 3
  ELSE
     CALL IO_WRITE_ASCII(efile, 'SCAL_READ_LOCAL. Error in Discrete.Type.')
     CALL DNS_STOP(DNS_ERROR_INFDISCR)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', delta_discrete)

  RETURN
END SUBROUTINE SCAL_READ_LOCAL
