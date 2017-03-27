#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2003/01/01 - J.P. Mellado
!#              Modified
!# 2007/05/07 - J.P. Mellado
!#              Formulation on temperature is added.
!#              New variables are added and input field names are
!#              modified for consistency => version 4.5
!# 2011/04/18 - A. Lozar
!#              Adding inertial effects
!#
!########################################################################
!# DESCRIPTION
!#
!# Reading general data from file dns.ini, setting up general parameters
!# and doing cross-check of these general data.
!#
!########################################################################
SUBROUTINE DNS_READ_GLOBAL(inifile)

  USE DNS_CONSTANTS, ONLY : lfile, efile, wfile
  USE DNS_GLOBAL
  USE THERMO_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) inifile

! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*64 lstr
  CHARACTER*32 bakfile
  TINTEGER iMajorVersion, iMinorVersion
  TINTEGER is, inb_scal_local1, inb_scal_local2, idummy
!  TINTEGER nspa_storage, nlin_storage
  TREAL dummy

! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  iMajorVersion = 6; iMinorVersion = 2

  CALL IO_WRITE_ASCII(lfile, 'Reading global input data.')

! ###################################################################
! Version Checking
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#[Version]')
  CALL IO_WRITE_ASCII(bakfile, '#Major=<mayor version number>')
  CALL IO_WRITE_ASCII(bakfile, '#Minor=<minor version number>')

  CALL SCANINIINT(bakfile, inifile, 'Version', 'Major', '0', idummy)
  IF ( iMajorVersion .NE. idummy ) THEN
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Major version error.')
     CALL DNS_STOP(DNS_ERROR_VERSION)
  ENDIF
  CALL SCANINIINT(bakfile, inifile, 'Version', 'Minor', '0', idummy)
  IF ( iMinorVersion .NE. idummy ) THEN
     WRITE(sRes,'(I5)') iMinorVersion
     CALL IO_WRITE_ASCII(wfile, 'DNS_REAL_GLOBAL. Minor version warning. Expected : '//sRes)
  ENDIF

! ###################################################################
! Global information
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Main]')
  CALL IO_WRITE_ASCII(bakfile, '#FileFormat=<RawArray/RawSplit/NetCDF>')
  CALL IO_WRITE_ASCII(bakfile, '#VerbosityLevel=<0/1/2>')
  CALL IO_WRITE_ASCII(bakfile, '#Type=<temporal/spatial>')
  CALL IO_WRITE_ASCII(bakfile, '#Flow=<shear/jet/isotropic>')
  CALL IO_WRITE_ASCII(bakfile, '#CalculateFlow=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#CalculateScalar=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#Equations=<total/internal/incompressible/anelastic>')
  CALL IO_WRITE_ASCII(bakfile, '#TermAdvection=<divergence/skewsymmetric>')
  CALL IO_WRITE_ASCII(bakfile, '#TermViscous=<divergence/explicit>')
  CALL IO_WRITE_ASCII(bakfile, '#TermDiffusion=<divergence/explicit>')
  CALL IO_WRITE_ASCII(bakfile, '#TermBodyForce=<none/Explicit/Linear/Bilinear/Quadratic>')
  CALL IO_WRITE_ASCII(bakfile, '#TermCoriolis=<none/explicit/normalized>')
  CALL IO_WRITE_ASCII(bakfile, '#TermRadiation=<none/Bulk1dGlobal/Bulk1dLocal>')
  CALL IO_WRITE_ASCII(bakfile, '#TermTransport=<constant/powerlaw/sutherland/Airwater/AirwaterSimplified>')
  CALL IO_WRITE_ASCII(bakfile, '#TermChemistry=<none/quadratic/layeredrelaxation/ozone>')
  CALL IO_WRITE_ASCII(bakfile, '#SpaceOrder=<CompactJacobian4/CompactJacobian6/CompactJacobian8/CompactDirect6>')

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'FileFormat', 'RawSplit', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'rawarray' ) THEN; imode_files = DNS_FILE_RAWARRAY
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'rawsplit' ) THEN; imode_files = DNS_FILE_RAWSPLIT
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'netcdf'   ) THEN; imode_files = DNS_FILE_NETCDF
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'none'   )   THEN; imode_files = DNS_NOFILE
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Wrong file format.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  CALL SCANINIINT(bakfile, inifile, 'Main', 'VerbosityLevel', '1', imode_verbosity)

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Type', 'temporal', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'temporal' ) THEN; imode_sim = DNS_MODE_TEMPORAL
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'spatial'  ) THEN; imode_sim = DNS_MODE_SPATIAL
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Simulation must be temporal/spatial')
     CALL DNS_STOP(DNS_ERROR_SIMTYPE)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Flow', 'shear', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'shear'       ) THEN; imode_flow = DNS_FLOW_SHEAR
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'jet'         ) THEN; imode_flow = DNS_FLOW_JET
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'isotropic'   ) THEN; imode_flow = DNS_FLOW_ISOTROPIC
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Flow type must be shear/jet/isotropic')
     CALL DNS_STOP(DNS_ERROR_SIMFLOW)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'CalculateFlow', 'yes', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_flow = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_flow = 0
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. CalculateFlow must be yes or no')
     CALL DNS_STOP(DNS_ERROR_CALCFLOW)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'CalculateScalar', 'yes', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_scal = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_scal = 0
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Scalar must be yes or no')
     CALL DNS_STOP(DNS_ERROR_CALCSCALAR)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'CalculateParticle', 'no', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_particle = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_particle = 0
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. CalculateParticle must be yes or no')
     CALL DNS_STOP(DNS_ERROR_CALCPARTICLE)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Equations', 'internal', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'total'          ) THEN; imode_eqns = DNS_EQNS_TOTAL
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'internal'       ) THEN; imode_eqns = DNS_EQNS_INTERNAL
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'incompressible' ) THEN; imode_eqns = DNS_EQNS_INCOMPRESSIBLE
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'anelastic'      ) THEN; imode_eqns = DNS_EQNS_ANELASTIC
  ELSE
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Wrong Equations option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN; ifourier = 1
  ELSE;                                         ifourier = 0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Mixture', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'          )     THEN; imixture = MIXT_TYPE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'airvapor'      ) THEN; imixture = MIXT_TYPE_AIRVAPOR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'airwater'      ) THEN; imixture = MIXT_TYPE_AIRWATER
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'airwaterlinear') THEN; imixture = MIXT_TYPE_AIRWATER_LINEAR
  ELSE 
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong multispecies model.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermAdvection', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'          ) THEN; iadvection = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'divergence'    ) THEN; iadvection = EQNS_DIVERGENCE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'skewsymmetric' ) THEN; iadvection = EQNS_SKEWSYMMETRIC
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'convective'    ) THEN; iadvection = EQNS_CONVECTIVE
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong TermAdvection option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermViscous', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'       ) THEN; iviscous = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'divergence' ) THEN; iviscous = EQNS_DIVERGENCE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit'   ) THEN; iviscous = EQNS_EXPLICIT
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong TermViscous option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermDiffusion', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'       ) THEN; idiffusion = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'divergence' ) THEN; idiffusion = EQNS_DIVERGENCE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit'   ) THEN; idiffusion = EQNS_EXPLICIT
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong TermDiffusion option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermBodyForce', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'        ) THEN; buoyancy%type = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'explicit'    ) THEN; buoyancy%type = EQNS_EXPLICIT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'homogeneous' ) THEN; buoyancy%type = EQNS_BOD_HOMOGENEOUS
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linear'      ) THEN; buoyancy%type = EQNS_BOD_LINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'bilinear'    ) THEN; buoyancy%type = EQNS_BOD_BILINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'quadratic'   ) THEN; buoyancy%type = EQNS_BOD_QUADRATIC
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong TermBodyForce option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermCoriolis', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'       ) THEN; coriolis%type = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit'   ) THEN; coriolis%type = EQNS_EXPLICIT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'normalized' ) THEN; coriolis%type = EQNS_COR_NORMALIZED
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong TermCoriolis option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermRadiation', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'          ) THEN; radiation%type = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'bulk1dglobal'  ) THEN; radiation%type = EQNS_RAD_BULK1D_GLOBAL
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'bulk1dlocal'   ) THEN; radiation%type = EQNS_RAD_BULK1D_LOCAL
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong TermRadiation option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermTransport', 'constant', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'sutherland'         ) THEN; transport%type = EQNS_TRANS_SUTHERLAND;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'powerlaw'           ) THEN; transport%type = EQNS_TRANS_POWERLAW; 
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'airwater'           ) THEN; transport%type = EQNS_TRANS_AIRWATER;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'airwatersimplified' ) THEN; transport%type = EQNS_TRANS_AIRWATERSIMPLIFIED;
  ELSE;                                                          transport%type = EQNS_NONE; ENDIF

  itransport = transport%type

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermChemistry', 'none', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'quadratic'        ) THEN; chemistry%type = EQNS_CHEM_QUADRATIC;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'layeredrelaxation') THEN; chemistry%type = EQNS_CHEM_LAYEREDRELAXATION; 
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'ozone'            ) THEN; chemistry%type = EQNS_CHEM_OZONE; 
  ELSE;                                                        chemistry%type = EQNS_NONE; ENDIF

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'SpaceOrder', 'void', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'compactjacobian4' ) THEN; imode_fdm = FDM_COM4_JACOBIAN; 
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'compactjacobian6' ) THEN; imode_fdm = FDM_COM6_JACOBIAN; 
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'compactjacobian8' ) THEN; imode_fdm = FDM_COM8_JACOBIAN; 
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'compactdirect6'   ) THEN; imode_fdm = FDM_COM6_DIRECT;
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong SpaceOrder option.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

  g(1:3)%mode_fdm = imode_fdm
  
! ###################################################################
! Iteration Section
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Iteration]')
  CALL SCANINIINT(bakfile, inifile, 'Statistics', 'StatSave', '10', nspa_rest)
  CALL SCANINIINT(bakfile, inifile, 'Statistics', 'StatStep', '10', nspa_step)

! ###################################################################
! Dimensionles parameters
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile,  '#')
  CALL IO_WRITE_ASCII(bakfile,  '#[Parameters]')
  CALL IO_WRITE_ASCII(bakfile,  '#Reynolds=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Prandtl=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Froude=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Rossby=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Mach=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Gama=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Schmidt=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Damkohler=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Stokes=<value>')
  CALL IO_WRITE_ASCII(bakfile,  '#Settling=<value>')

  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Reynolds', '100', reynolds  )
  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Gama',     '1.4', gama0     )
  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Prandtl',  '1.0', prandtl   )
  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Mach',     '1.0', mach      )
  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Froude',   '1.0', froude    )
  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Rossby',   '1.0', rossby    )
  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Stokes',   '0.0', stokes    )
  CALL SCANINIREAL(bakfile, inifile, 'Parameters', 'Settling', '0.0', settling  )

  CALL SCANINICHAR(bakfile, inifile, 'Parameters', 'Schmidt',  '1.0', sRes)
  schmidt(:) = C_0_R; inb_scal_local1 = MAX_NSP
  CALL LIST_REAL(sRes, inb_scal_local1, schmidt  )

  lstr='0.0'; DO is = 2,inb_scal_local1; lstr = TRIM(ADJUSTL(lstr))//',0.0'; ENDDO
  CALL SCANINICHAR(bakfile, inifile, 'Parameters', 'Damkohler', lstr, sRes)
  damkohler(:) = C_0_R; inb_scal_local2 = MAX_NSP
  CALL LIST_REAL(sRes, inb_scal_local2, damkohler)
  IF ( inb_scal_local1 .NE. inb_scal_local2 ) THEN ! Consistency check
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Schmidt and Damkholer sizes do not match.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF

! ###################################################################
! Buoyancy
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[BodyForce]')
  CALL IO_WRITE_ASCII(bakfile, '#Vector=<Gx,Gy,Gz>')
  CALL IO_WRITE_ASCII(bakfile, '#Parameters=<value>')

  buoyancy%vector = C_0_R; buoyancy%active = .FALSE.
  IF ( buoyancy%type .NE. EQNS_NONE ) THEN
     CALL SCANINICHAR(bakfile, inifile, 'BodyForce', 'Vector', '0.0,-1.0,0.0', sRes)
     idummy = 3
     CALL LIST_REAL(sRes, idummy, buoyancy%vector)

     IF ( ABS(buoyancy%vector(1)) .GT. C_0_R ) THEN; buoyancy%active(1) = .TRUE.; CALL IO_WRITE_ASCII(lfile, 'Body force along Ox.'); ENDIF
     IF ( ABS(buoyancy%vector(2)) .GT. C_0_R ) THEN; buoyancy%active(2) = .TRUE.; CALL IO_WRITE_ASCII(lfile, 'Body force along Oy.'); ENDIF
     IF ( ABS(buoyancy%vector(3)) .GT. C_0_R ) THEN; buoyancy%active(3) = .TRUE.; CALL IO_WRITE_ASCII(lfile, 'Body force along Oz.'); ENDIF

     IF ( froude .GT. C_0_R ) THEN
           buoyancy%vector(:) = buoyancy%vector(:) /froude ! adding the froude number into de vector g
     ELSE
        CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Froude number must be nonzero if buoyancy is retained.')
        CALL DNS_STOP(DNS_ERROR_OPTION)        
     ENDIF

     buoyancy%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'BodyForce', 'Parameters', '0.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, buoyancy%parameters)

  ENDIF

! ###################################################################
! Rotation
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Rotation]')
  CALL IO_WRITE_ASCII(bakfile, '#Vector=<Fx,Fy,Fz>')
  CALL IO_WRITE_ASCII(bakfile, '#Parameters=<value>')

  coriolis%vector(:) = C_0_R; coriolis%active = .FALSE.
  IF ( coriolis%type .NE. EQNS_NONE ) THEN
     CALL SCANINICHAR(bakfile, inifile, 'Rotation', 'Vector', '0.0,1.0,0.0', sRes)
     idummy = 3
     CALL LIST_REAL(sRes, idummy, coriolis%vector)  
     
     IF ( ABS(coriolis%vector(1)) .GT. C_0_R ) THEN; coriolis%active(2) = .TRUE.; coriolis%active(3) = .TRUE.; CALL IO_WRITE_ASCII(lfile, 'Angular velocity along Ox.'); ENDIF
     IF ( ABS(coriolis%vector(2)) .GT. C_0_R ) THEN; coriolis%active(3) = .TRUE.; coriolis%active(1) = .TRUE.; CALL IO_WRITE_ASCII(lfile, 'Angular velocity along Oy.'); ENDIF
     IF ( ABS(coriolis%vector(3)) .GT. C_0_R ) THEN; coriolis%active(1) = .TRUE.; coriolis%active(2) = .TRUE.; CALL IO_WRITE_ASCII(lfile, 'Angular velocity along Oz.'); ENDIF
              
     IF ( rossby .GT. C_0_R ) THEN
        coriolis%vector(:) = coriolis%vector(:) /rossby ! adding the rossby number into the vector
     ELSE
        CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Rossby number must be nonzero if coriolis is retained.')
        CALL DNS_STOP(DNS_ERROR_OPTION)        
     ENDIF
     
     coriolis%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Rotation', 'Parameters', '0.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, coriolis%parameters)

  ENDIF

! Consistency check
  IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
     IF ( coriolis%active(2) ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. TermCoriolis option only allows for angular velocity along Oy.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

! ###################################################################
! Radiation
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Radiation]')
  CALL IO_WRITE_ASCII(bakfile, '#Scalar=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Parameters=<value>')

  radiation%active = .FALSE.
  IF ( radiation%type .NE. EQNS_NONE ) THEN
     CALL SCANINIINT(bakfile, inifile, 'Radiation', 'Scalar', '1', idummy)   
     radiation%active(idummy) = .TRUE.
     
     radiation%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Radiation', 'Parameters', '1.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, radiation%parameters)

  ENDIF

! ###################################################################
! Transport
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Transport]')
  CALL IO_WRITE_ASCII(bakfile, '#Parameters=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Exponent=<value>')

  transport%active = .FALSE.
  IF ( transport%type .NE. EQNS_NONE ) THEN
     transport%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Transport', 'Parameters', '1.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, transport%parameters)

     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        transport%active = .TRUE. ! All scalars are affected

        CALL SCANINIREAL(bakfile, inifile, 'Transport', 'Exponent', '0.0', transport%auxiliar(1))
     ENDIF
     
  ENDIF
  
! ###################################################################
! Chemistry
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Chemistry]')
  CALL IO_WRITE_ASCII(bakfile, '#Parameters=<value>')

  IF      ( chemistry%type .NE. EQNS_NONE ) THEN
     chemistry%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Chemistry', 'Parameters', '1.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, chemistry%parameters)

  ENDIF

! Activating terms
  chemistry%active = .FALSE.
  DO is = 1,inb_scal_local1
     IF ( ABS(damkohler(is)) .GT. C_0_R ) chemistry%active(is) = .TRUE.
  ENDDO
  
! ###################################################################
! Thermodynamics
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Thermodynamics]')
  CALL IO_WRITE_ASCII(bakfile, '#Parameters=<value>')

  IF ( imixture .NE. EQNS_NONE ) THEN
     thermo_param(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Thermodynamics', 'Parameters', '1.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, thermo_param)

  ENDIF

! ###################################################################
! Grid Parameters
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Grid]')
  CALL IO_WRITE_ASCII(bakfile, '#Imax=<imax>')
  CALL IO_WRITE_ASCII(bakfile, '#Imax(*)=<imax_proc>')
  CALL IO_WRITE_ASCII(bakfile, '#Jmax=<jmax>')
  CALL IO_WRITE_ASCII(bakfile, '#Kmax=<kmax>')
  CALL IO_WRITE_ASCII(bakfile, '#Kmax(*)=<kmax_proc>')
  CALL IO_WRITE_ASCII(bakfile, '#XUniform=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#YUniform=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#ZUniform=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#XPeriodic=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#YPeriodic=<yes/no>')
  CALL IO_WRITE_ASCII(bakfile, '#ZPeriodic=<yes/no>')

  CALL SCANINIINT(bakfile, inifile, 'Grid', 'Imax', '0', imax_total)
  CALL SCANINIINT(bakfile, inifile, 'Grid', 'Jmax', '0', jmax_total)
  CALL SCANINIINT(bakfile, inifile, 'Grid', 'Kmax', '0', kmax_total)
  g(1)%size = imax_total
  g(2)%size = jmax_total
  g(3)%size = kmax_total
  
! default
  imax = imax_total
  jmax = jmax_total
  kmax = kmax_total

! -------------------------------------------------------------------
! Domain decomposition in parallel mode
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro .GT. 1 ) THEN
     CALL SCANINIINT(bakfile, inifile, 'Grid', 'Kmax(*)', '-1', kmax)
     IF ( kmax .GT. 0 .AND. MOD(kmax_total,kmax) .EQ. 0 ) THEN
        ims_npro_k = kmax_total/kmax
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Input kmax incorrect')
        CALL DNS_STOP(DNS_ERROR_KMAXTOTAL)
     ENDIF
     
     CALL SCANINIINT(bakfile, inifile, 'Grid', 'Imax(*)', '-1', imax)
     IF ( imax .GT. 0 .AND. MOD(imax_total,imax) .EQ. 0 ) THEN
        ims_npro_i = imax_total/imax
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Input imax incorrect')
        CALL DNS_STOP(DNS_ERROR_KMAXTOTAL)
     ENDIF
 
     IF ( ims_npro_i*ims_npro_k .EQ. ims_npro ) THEN ! check
        WRITE(lstr,*) ims_npro_i; WRITE(sRes,*) ims_npro_k
        lstr = TRIM(ADJUSTL(lstr))//'x'//TRIM(ADJUSTL(sRes))
        CALL IO_WRITE_ASCII(lfile, 'Initialize domain partition '//TRIM(ADJUSTL(lstr)))
     ELSE
        CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Inconsistency in total number of PEs')
        CALL DNS_STOP(DNS_ERROR_KMAXTOTAL)
     ENDIF

  ENDIF

#endif

! -------------------------------------------------------------------
! Uniform
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'XUniform', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; iunifx = 0; g(1)%uniform = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; iunifx = 1; g(1)%uniform = .FALSE.
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Error in Uniform X grid')
     CALL DNS_STOP(DNS_ERROR_UNIFORMX)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'YUniform', 'void', sRes)      
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; iunify = 0; g(2)%uniform = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; iunify = 1; g(2)%uniform = .FALSE.
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Error in Uniform Y grid')
     CALL DNS_STOP(DNS_ERROR_UNIFORMY)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'ZUniform', 'void', sRes)      
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; iunifz = 0; g(3)%uniform = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; iunifz = 1; g(3)%uniform = .FALSE.
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Error in Uniform Z grid')
     CALL DNS_STOP(DNS_ERROR_UNIFORMZ)
  ENDIF

! -------------------------------------------------------------------
! Periodic
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'XPeriodic', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; i1bc = 0; g(1)%periodic = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; i1bc = 1; g(1)%periodic = .FALSE.
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Error in Periodic X grid')
     CALL DNS_STOP(DNS_ERROR_IBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'YPeriodic', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; j1bc = 0; g(2)%periodic = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; j1bc = 1; g(2)%periodic = .FALSE.
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Error in Periodic Y grid')
     CALL DNS_STOP(DNS_ERROR_JBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'ZPeriodic', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; k1bc = 0; g(3)%periodic = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; k1bc = 1; g(3)%periodic = .FALSE.
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Error in Periodic Z grid')
     CALL DNS_STOP(DNS_ERROR_KBC)
  ENDIF


! ###################################################################
! Statistics Control   
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Statistics]')
  CALL IO_WRITE_ASCII(bakfile, '#IAvera=<plane1,plane2,...>')
  CALL IO_WRITE_ASCII(bakfile, '#ILines=<plane1,plane2,...>')
  CALL IO_WRITE_ASCII(bakfile, '#JLines=<plane1,plane2,...>')
  CALL IO_WRITE_ASCII(bakfile, '#IPlane=<plane1,plane2,...>')

  nstatavg = MAX_STATS_SPATIAL
  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'IAvera', '1', sRes)
  CALL LIST_INTEGER(sRes, nstatavg, statavg)

  nstatlin = MAX_STATS_SPATIAL
  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'ILines', '1', sRes)
  CALL LIST_INTEGER(sRes, nstatlin, statlin_i)
  idummy = MAX_STATS_SPATIAL
  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'JLines', '1', sRes)
  CALL LIST_INTEGER(sRes, idummy, statlin_j)
  IF ( idummy .NE. nstatlin ) THEN
     CALL IO_WRITE_ASCII(efile,'Mismatch in the number of i-j-lines')
     CALL DNS_STOP(DNS_ERROR_IJMISMATCH)
  ENDIF

  nstatpln = MAX_STATS_SPATIAL
  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'IPlane', '1', sRes)
  CALL LIST_INTEGER(sRes, nstatpln, statpln)

! ###################################################################
! Flow physical properties of the system
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Flow]')
  CALL IO_WRITE_ASCII(bakfile, '#VelocityX=<mean velocity X>')
  CALL IO_WRITE_ASCII(bakfile, '#VelocityY=<mean velocity Y>')
  CALL IO_WRITE_ASCII(bakfile, '#VelocityZ=<mean velocity Z>')
  CALL IO_WRITE_ASCII(bakfile, '#Pressure=<mean pressure>')
  CALL IO_WRITE_ASCII(bakfile, '#Density=<mean density>')
  CALL IO_WRITE_ASCII(bakfile, '#Temperature=<mean temperature>')
  CALL IO_WRITE_ASCII(bakfile, '#ProfileVelocity=<None/Linear/Tanh/Erf/Ekman/EkmanP/Parabolic>')
  CALL IO_WRITE_ASCII(bakfile, '#YCoorVelocity=<Relative Y reference point>')
  CALL IO_WRITE_ASCII(bakfile, '#DiamVelocity=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#ThickVelocity=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#DeltaVelocity=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#ProfileDensity=<None/Linear/Tanh/Erf>')
  CALL IO_WRITE_ASCII(bakfile, '#YCoorDensity=<Relative Y reference point>')
  CALL IO_WRITE_ASCII(bakfile, '#DiamDensity=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#ThickDensity=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#DeltaDensity=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#ProfileTemperature=<None/Linear/Tanh/Erf>')
  CALL IO_WRITE_ASCII(bakfile, '#YCoorTemperature=<Relative Y reference point>')
  CALL IO_WRITE_ASCII(bakfile, '#DiamTemperature=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#ThickTemperature=<value>') 
  CALL IO_WRITE_ASCII(bakfile, '#DeltaTemperature=<value>')

! -------------------------------------------------------------------
! Mean flow values
! -------------------------------------------------------------------
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'VelocityX',  '0.0', mean_u)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'VelocityY',  '0.0', mean_v)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'VelocityZ',  '0.0', mean_w)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'Pressure',   '0.0', pbg%mean)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'Density',    '0.0', rbg%mean)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'Temperature','0.0', tbg%mean)

! -------------------------------------------------------------------
! Shear flow values
! -------------------------------------------------------------------
! streamwise velocity
  CALL SCANINICHAR(bakfile, inifile, 'Flow', 'ProfileVelocity', 'Tanh', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'      ) THEN; iprof_u = PROFILE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linear'    ) THEN; iprof_u = PROFILE_LINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'tanh'      ) THEN; iprof_u = PROFILE_TANH
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'erf'       ) THEN; iprof_u = PROFILE_ERF
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'bickley'   ) THEN; iprof_u = PROFILE_BICKLEY
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'gaussian'  ) THEN; iprof_u = PROFILE_GAUSSIAN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'ekman'     ) THEN; iprof_u = PROFILE_EKMAN_U
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'ekmanp'    ) THEN; iprof_u = PROFILE_EKMAN_U_P
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'parabolic' ) THEN; iprof_u = PROFILE_PARABOLIC
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linearcrop') THEN; iprof_u = PROFILE_LINEAR_CROP
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'mixedlayer') THEN; iprof_u = PROFILE_MIXEDLAYER
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong velocity profile.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'YCoorVelocity', '0.5', ycoor_u)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'DiamVelocity',  '1.0', diam_u )
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickVelocity', '0.0', thick_u)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'DeltaVelocity', '0.0', delta_u)

! density
  CALL SCANINICHAR(bakfile, inifile, 'Flow', 'ProfileDensity', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'      ) THEN; rbg%type = PROFILE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linear'    ) THEN; rbg%type = PROFILE_LINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'tanh'      ) THEN; rbg%type = PROFILE_TANH
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'erf'       ) THEN; rbg%type = PROFILE_ERF
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'parabolic' ) THEN; rbg%type = PROFILE_PARABOLIC
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linearcrop') THEN; rbg%type = PROFILE_LINEAR_CROP
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'mixedlayer') THEN; rbg%type = PROFILE_MIXEDLAYER
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong density profile.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'YCoorDensity', '0.5', rbg%ymean)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'DiamDensity',  '1.0', rbg%diam )
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickDensity', '0.0', rbg%thick)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'DeltaDensity', '0.0', rbg%delta)

! temperature/enthalpy
  CALL SCANINICHAR(bakfile, inifile, 'Flow', 'ProfileTemperature', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'              ) THEN; tbg%type = PROFILE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linear'            ) THEN; tbg%type = PROFILE_LINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'tanh'              ) THEN; tbg%type = PROFILE_TANH
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'erf'               ) THEN; tbg%type = PROFILE_ERF
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linearerf'         ) THEN; tbg%type = PROFILE_LINEAR_ERF
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'enthalpyerf'       ) THEN; tbg%type =-PROFILE_ERF
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'enthalpylinearerf' ) THEN; tbg%type =-PROFILE_LINEAR_ERF
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'parabolic'         ) THEN; tbg%type = PROFILE_PARABOLIC
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linearcrop'        ) THEN; tbg%type = PROFILE_LINEAR_CROP
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'mixedlayer'        ) THEN; tbg%type = PROFILE_MIXEDLAYER
  ELSE
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong temperature profile.')
     CALL DNS_STOP(DNS_ERROR_OPTION)
  ENDIF
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'YCoorTemperature', '0.5', tbg%ymean)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'DiamTemperature',  '1.0', tbg%diam )
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickTemperature', '0.0', tbg%thick)
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'DeltaTemperature', '0.0', tbg%delta)

! pressure
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'YCoorPressure','0.5', pbg%ymean        )
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'RefPressure',  '0.0', pbg%reference    )
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ScaleHeight',  '0.0', pbg%parameters(1))
  
! additional specific data
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'BottomSlope', '0.0',  tbg%parameters(1))
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'UpperSlope',  '0.0',  tbg%parameters(2))

! consistency check
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     IF ( rbg%type .EQ. PROFILE_NONE .AND. tbg%type .EQ. PROFILE_NONE ) THEN
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Specify density or temperature.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
     
     IF ( rbg%type .NE. PROFILE_NONE .AND. tbg%type .NE. PROFILE_NONE ) THEN
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Specify only density or only temperature.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

! -------------------------------------------------------------------
! Spatial case
! Thickness evolutions delta_i/diam_i=a*(x/diam_i+b)
! -------------------------------------------------------------------
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     CALL IO_WRITE_ASCII(bakfile, '#ThickAVelocity=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#ThickBVelocity=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#FluxVelocity=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#ThickADensity=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#ThickBDensity=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#FluxDensity=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#ThickATemperature=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#ThickBTemperature=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#FluxTemperature=<value>')

! Bradbury profile is the default (x0=a*b)
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickAVelocity', '0.1235', jet_u(1))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBVelocity', '-0.873', jet_u(2))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'FluxVelocity',   '0.96',   jet_u(3))

! Ramaprian is the default (x0=a*b)
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickADensity', '0.14', rbg%parameters(1))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBDensity', '2.0',  rbg%parameters(2))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'FluxDensity',   '0.94', rbg%parameters(3))

     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickATemperature', '0.14', tbg%parameters(1))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBTemperature', '2.0',  tbg%parameters(2))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'FluxTemperature',   '0.94', tbg%parameters(3))

  ENDIF

! ###################################################################
! Scalars physical properties of the system
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')
  CALL IO_WRITE_ASCII(bakfile, '#[Scalar]')
  CALL IO_WRITE_ASCII(bakfile, '#Profile=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#YCoor=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Diam=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Thick=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Mean=<value>')
  CALL IO_WRITE_ASCII(bakfile, '#Delta=<value>')

! read scalars profiles
  DO is = 1,MAX_NSP
     WRITE(lstr,*) is; lstr='ProfileScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINICHAR(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), 'None', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'      ) THEN; iprof_i(is) = PROFILE_NONE
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linear'    ) THEN; iprof_i(is) = PROFILE_LINEAR
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'tanh'      ) THEN; iprof_i(is) = PROFILE_TANH
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'erf'       ) THEN; iprof_i(is) = PROFILE_ERF
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linearerf' ) THEN; iprof_i(is) = PROFILE_LINEAR_ERF
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'erfantisym') THEN; iprof_i(is) = PROFILE_ERF_ANTISYM
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'parabolic' ) THEN; iprof_i(is) = PROFILE_PARABOLIC
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linearcrop') THEN; iprof_i(is) = PROFILE_LINEAR_CROP
     ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'mixedlayer') THEN; iprof_i(is) = PROFILE_MIXEDLAYER
     ELSE
        CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Wrong species profile.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
     WRITE(lstr,*) is; lstr='MeanScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.0', mean_i(is) )
     WRITE(lstr,*) is; lstr='YCoorScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.5', ycoor_i(is))
     WRITE(lstr,*) is; lstr='DiamScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '1.0', diam_i(is) )
     WRITE(lstr,*) is; lstr='ThickScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.0', thick_i(is))
     WRITE(lstr,*) is; lstr='DeltaScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.0', delta_i(is))

! additional specific data     
     WRITE(lstr,*) is; lstr='BottomSlopeScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.0',  prof_i(1,is))
     WRITE(lstr,*) is; lstr='UpperSlopeScalar'//TRIM(ADJUSTL(lstr))
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.0',  prof_i(2,is))

     IF ( iprof_i(is) .EQ. PROFILE_ERF_ANTISYM ) THEN
        WRITE(lstr,*) is; lstr='YCoorSymmetry'//TRIM(ADJUSTL(lstr))
        CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.0',  prof_i(1,is))
     ENDIF
     
  ENDDO

! use chemkin for the thermodynamic data
  CALL SCANINICHAR(bakfile, inifile, 'Scalar', 'ChemkinFile', 'none', chemkin_file)
  IF ( TRIM(ADJUSTL(chemkin_file)) .EQ. 'none' ) THEN; iuse_chemkin = 0
  ELSE;                                                iuse_chemkin = 1; ENDIF

! Special data
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     CALL SCANINIREAL(bakfile, inifile, 'Scalar', 'SmoothFactor', '0.1', dsmooth)
  ENDIF

! -------------------------------------------------------------------
! Spatial case
! Thickness evolutions delta_i/diam_i=a*(x/diam_i+b)
! Ramaprian is the default (x0=a*b)
! -------------------------------------------------------------------
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     CALL IO_WRITE_ASCII(bakfile, '#ThickA=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#ThickB=<value>')
     CALL IO_WRITE_ASCII(bakfile, '#Flux=<value>')

     DO is = 1, MAX_NSP
        WRITE(lstr,*) is; lstr='ThickA'//TRIM(ADJUSTL(lstr))
        CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.14', jet_i(1,is))
        WRITE(lstr,*) is; lstr='ThickB'//TRIM(ADJUSTL(lstr))
        CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '2.0',  jet_i(2,is))
        WRITE(lstr,*) is; lstr='Flux'//TRIM(ADJUSTL(lstr))
        CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.94', jet_i(3,is))
     ENDDO

  ENDIF

! ###################################################################
! Final initialization
! ###################################################################
  CALL IO_WRITE_ASCII(bakfile, '#')

  IF ( iviscous .EQ. EQNS_NONE ) THEN; visc = C_0_R
  ELSE;                                visc = C_1_R/reynolds; ENDIF

! -------------------------------------------------------------------
! Initializing thermodynamic data of the mixture
! -------------------------------------------------------------------
  inb_scal       = inb_scal_local1 ! Default is general N scalars; gama0 has been already read above.
  inb_scal_array = inb_scal
  NSP            = inb_scal

  IF ( imixture .NE. MIXT_TYPE_NONE ) THEN ! particular mixture (requires implementation)
     CALL THERMO_INITIALIZE                ! gama0 is defined here
     IF ( inb_scal_local1 .NE. inb_scal )  THEN 
        CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Incorrect number of Schmidt numbers.')
        CALL DNS_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF
! Value of R_0/(C_{p,0}W_0) is called GRATIO
  IF ( gama0 .GT. C_0_R ) GRATIO = (gama0-C_1_R)/gama0

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     mach    = C_0_R
     MRATIO  = C_1_R
     prandtl = schmidt(1)
  ELSE
     MRATIO = gama0*mach*mach
  ENDIF

  IF ( buoyancy%type .EQ. EQNS_BOD_LINEAR   .OR. &
       buoyancy%type .EQ. EQNS_BOD_BILINEAR .OR. &
       buoyancy%type .EQ. EQNS_BOD_QUADRATIC ) THEN
     IF ( inb_scal .EQ. 0 ) THEN
        CALL IO_WRITE_ASCII(wfile,'DNS_READ_GLOBAL. Zero scalars; setting TermBodyForce equal to none.')
        buoyancy%type = EQNS_NONE
     ENDIF
  ENDIF
  
! mean_rho and delta_rho need to be defined, because of old version.
! Note that rho1 and rho2 are the values defined by equation of state,
! being then mean_rho=(rho1+rho2)/2.
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     IF ( rbg%type .EQ. PROFILE_NONE ) THEN
        dummy     = tbg%delta /tbg%mean
        rbg%mean  = MRATIO *pbg%mean / tbg%mean /(C_1_R-C_025_R*dummy*dummy)
        rbg%delta =-rbg%mean *dummy
        rbg%thick = tbg%thick
        rbg%diam  = tbg%diam
     ELSE
        dummy     = rbg%delta /rbg%mean
        tbg%mean  = MRATIO *pbg%mean /rbg%mean /(C_1_R-C_025_R*dummy*dummy)
        tbg%delta =-tbg%mean *dummy
        tbg%thick = rbg%thick
        tbg%diam  = rbg%diam
     ENDIF
  ENDIF

! Consistency check
  IF ( imixture .GT. 0 ) THEN

     IF      ( imixture .EQ. MIXT_TYPE_BS          &
          .OR. imixture .EQ. MIXT_TYPE_BSZELDOVICH ) THEN
! These cases force Sc_i=Sc_Z=Pr (Lewis unity)
        schmidt(inb_scal) = prandtl

     ELSE IF ( imixture .EQ. MIXT_TYPE_QUASIBS     ) THEN
! These cases force Sc_i=Sc_Z, already read

     ELSE IF ( imixture .EQ. MIXT_TYPE_AIRWATER    ) THEN
        IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
           schmidt(2:3) = schmidt(1) ! used in diffusion eqns, though should be fixed
        ENDIF

        IF ( damkohler(1) .EQ. C_0_R .AND. damkohler(2) .EQ. C_0_R ) THEN
           damkohler(1:2) = damkohler(3)
        ELSE
           CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. AirWater requires at least first 2 Damkholer numbers zero.')
           CALL DNS_STOP(DNS_ERROR_OPTION)
        ENDIF
        
     ENDIF

  ENDIF

! -------------------------------------------------------------------
! Arrays size
! -------------------------------------------------------------------
  isize_field = imax*jmax*kmax

  IF ( imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE .OR. imode_eqns .EQ. DNS_EQNS_ANELASTIC ) THEN
     inb_flow       = 3
     inb_flow_array = inb_flow
  ELSE
     inb_flow       = 5
     inb_flow_array = inb_flow + 2                                ! space for p, T
     IF ( transport%type .EQ. EQNS_TRANS_SUTHERLAND .OR. transport%type .EQ. EQNS_TRANS_POWERLAW ) inb_flow_array = inb_flow_array + 1 ! space for vis
  ENDIF
  inb_vars = inb_flow + inb_scal

  inb_wrk1d = 16
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN; inb_wrk2d = 11
  ELSE;                                        inb_wrk2d =  2; ENDIF

  isize_wrk1d = MAX(imax_total,MAX(jmax_total,kmax_total))
  isize_wrk2d = MAX(imax*jmax, MAX(imax*kmax,jmax*kmax)  )

! grid array
  DO is = 1,3
     g(is)%inb_grid = 1                  ! Nodes
     g(is)%inb_grid = g(is)%inb_grid &
                    + 2                  ! Jacobians of first- and second-order derivatives

     IF ( g(is)%periodic ) THEN
        g(is)%inb_grid = g(is)%inb_grid  &
                       + 5               & ! LU decomposition 1. order
                       + 5               & ! LU decomposition 2. order
                       + 5 *(1+inb_scal) & ! LU decomposition 2. order with diffusivities
                       + 2                 ! modified wavenumbers
     ELSE
        g(is)%inb_grid = g(is)%inb_grid  &
                       + 3 *4            & ! LU decomposition 1. order, 4 bcs
                       + 3 *4            & ! LU decomposition 2. order, 4 bcs
                       + 3 *(1+inb_scal)   ! LU decomposition 2. order w/ diffusivities, 1 bcs
! In Direct mode, we only need 10 instead of 3*4 because only 1 bcs is considered
     ENDIF
  END DO
     
! auxiliar array txc
  isize_txc_field = imax*jmax*kmax
  IF ( ifourier .EQ. 1 ) THEN
     isize_txc_dimz  = (imax+2)*(jmax+2)
     isize_txc_dimx  =  kmax   *(jmax+2)
     isize_txc_field = isize_txc_dimz*kmax ! space for FFTW lib
#ifdef USE_MPI
     IF ( ims_npro_k .GT. 1 ) THEN
        IF ( MOD(isize_txc_dimz,(2*ims_npro_k)) .NE. 0 ) THEN ! add space for MPI transposition
           isize_txc_dimz =  isize_txc_dimz     /(2*ims_npro_k)
           isize_txc_dimz = (isize_txc_dimz + 1)*(2*ims_npro_k)
        ENDIF
        isize_txc_field = MAX(isize_txc_field,isize_txc_dimz*kmax)
     ENDIF
     IF ( ims_npro_i .GT. 1 ) THEN
        IF ( MOD(isize_txc_dimx,(2*ims_npro_i)) .NE. 0 ) THEN ! add space for MPI transposition
           isize_txc_dimx =  isize_txc_dimx     /(2*ims_npro_i)
           isize_txc_dimx = (isize_txc_dimx + 1)*(2*ims_npro_i)
        ENDIF
        isize_txc_field = MAX(isize_txc_field,isize_txc_dimx*(imax+2))
     ENDIF
#endif
     IF ( MOD(imax,2) .NE. 0 ) THEN
        CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Imax must be a multiple of 2 for the FFT operations.')
        CALL DNS_STOP(DNS_ERROR_DIMGRID)
     ENDIF
  ENDIF

! loop counters over the whole domain are integer*4
  IF ( isize_field .GT. HUGE(imax) ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Integer model of 4 bytes not big enough.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)     
  ENDIF

! -------------------------------------------------------------------
! Test periodicity constrains
! -------------------------------------------------------------------
  IF ( i1bc .EQ. 0 .AND. iunifx .EQ. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Grid must be uniform in periodic direction X')
     CALL DNS_STOP(DNS_ERROR_CHECKUNIFX)
  ENDIF

  IF ( j1bc .EQ. 0 .AND. iunify .EQ. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Grid must be uniform in periodic direction Y')
     CALL DNS_STOP(DNS_ERROR_CHECKUNIFY)
  ENDIF

  IF ( k1bc .EQ. 0 .AND. iunifz .EQ. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Grid must be uniform in periodic direction Z')
     CALL DNS_STOP(DNS_ERROR_CHECKUNIFZ)
  ENDIF

! -------------------------------------------------------------------
! Other parameters
! -------------------------------------------------------------------
! By default, transport and radiation are caused by last scalar
  transport%scalar = inb_scal_array
  radiation%scalar = inb_scal_array
  
  IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
     IF ( transport%type .NE. EQNS_NONE ) THEN 
        transport%scalar = inb_scal_array                ! Transport is caused by liquid
        transport%parameters(inb_scal_array    ) = C_1_R ! liquid
        transport%parameters(inb_scal_array + 1) = C_1_R ! buoyancy
! Adding the settling number in the parameter definitions
        transport%parameters = transport%parameters *settling
     ENDIF

     IF ( radiation%type .NE. EQNS_NONE ) THEN 
        radiation%scalar = inb_scal_array             ! Radiation is caused by liquid
        radiation%active(inb_scal_array    ) = .TRUE. ! liquid
        radiation%active(inb_scal_array + 1) = .TRUE. ! buoyancy
     ENDIF
     
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. i1bc .NE. 0 ) THEN
     CALL IO_WRITE_ASCII(efile,'DNS_READ_GLOBAL. Grid must be uniform and periodic in direction X for temporal simulation')
     CALL DNS_STOP(DNS_ERROR_CHECKUNIFX)
  ENDIF

  IF ( iMajorVersion*10+iMinorVersion .LE. 42 ) THEN; nstatplnextra = 0
  ELSE;                                               nstatplnextra = inb_scal + inb_flow; ENDIF
! Add TOTAL number of variables + extra + temperature
  nstatplnvars = inb_vars + nstatplnextra + 1 

  IF ( inb_vars .GT. MAX_VARS ) THEN
     CALL IO_WRITE_ASCII(efile, 'DNS_READ_GLOBAL. Error MAX_VARS < inb_vars')
     CALL DNS_STOP(DNS_ERROR_TOTALVARS)
  ENDIF

! Minimum number of plane saves
  nspa_rest = MAX(i1, nspa_rest)

!  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
!     IF ( frunplane .EQ. 1 ) THEN
!        nspa_storage = nstatplnvars*nstatpln*nspa_rest/imax
!        WRITE(lstr,*) nspa_storage
!        CALL IO_WRITE_ASCII(bakfile, 'Total Plane Storage :'//lstr)
!     ENDIF
!     IF ( frunline .EQ. 1 ) THEN
!        nlin_storage = inb_vars*nstatlin*nspa_rest/(imax*jmax)
!        WRITE(lstr,*) nlin_storage
!        CALL IO_WRITE_ASCII(bakfile, 'Total Line Storage :'//lstr)
!     ENDIF
!  ENDIF

! Verification of averages planes
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     CALL SORT_INTEGER(nstatavg,statavg)
! include of last plane
     IF ( statavg(nstatpln) .LT. imax ) THEN
        nstatavg = nstatavg + 1; statavg(nstatavg) = imax
     ENDIF
! include first plane
     IF ( statavg(1) .NE. 1 ) THEN
        nstatavg = nstatavg + 1; statavg(nstatavg) = 1
     ENDIF
     CALL SORT_INTEGER(nstatavg,statavg)
  ENDIF

  RETURN
END SUBROUTINE DNS_READ_GLOBAL
