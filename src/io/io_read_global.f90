#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "IO_READ_LOCAL"

!########################################################################
!# DESCRIPTION
!#
!# Reading general data from file dns.ini, setting up general parameters
!# and doing cross-check of these general data.
!#
!########################################################################
SUBROUTINE IO_READ_GLOBAL(inifile)

  USE TLAB_CONSTANTS, ONLY : lfile, efile, wfile, MajorVersion, MinorVersion
  USE TLAB_VARS
  USE TLAB_PROCS
  USE THERMO_VARS
#ifdef USE_MPI
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"

  CHARACTER*(*) inifile

! -------------------------------------------------------------------
  CHARACTER*512 sRes
  CHARACTER*64 lstr, default
  CHARACTER*32 bakfile
  TINTEGER is, ig, inb_scal_local1, inb_scal_local2, idummy
  TREAL dummy

! ###################################################################
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL TLAB_WRITE_ASCII(lfile, 'Reading global input data.')

! ###################################################################
! Version Checking
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#[Version]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Major=<mayor version number>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Minor=<minor version number>')

  CALL SCANINIINT(bakfile, inifile, 'Version', 'Major', '0', idummy)
  IF ( MajorVersion .NE. idummy ) THEN
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Major version error.')
     CALL TLAB_STOP(DNS_ERROR_VERSION)
  ENDIF
  CALL SCANINIINT(bakfile, inifile, 'Version', 'Minor', '0', idummy)
  IF ( MinorVersion .NE. idummy ) THEN
     WRITE(sRes,'(I5)') MinorVersion
     CALL TLAB_WRITE_ASCII(wfile, 'DNS_REAL_GLOBAL. Minor version warning. Expected : '//sRes)
  ENDIF

! ###################################################################
! Global information
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Main]')
  CALL TLAB_WRITE_ASCII(bakfile, '#FileFormat=<mpiio/RawSplit/NetCDF>')
  CALL TLAB_WRITE_ASCII(bakfile, '#VerbosityLevel=<0/1/2>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Type=<temporal/spatial>')
  CALL TLAB_WRITE_ASCII(bakfile, '#CalculateFlow=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#CalculateScalar=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Equations=<total/internal/incompressible/anelastic>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermAdvection=<divergence/skewsymmetric>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermViscous=<divergence/explicit>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermDiffusion=<divergence/explicit>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermBodyForce=<none/Explicit/Homogeneous/Linear/Bilinear/Quadratic>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermCoriolis=<none/explicit/normalized>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermRadiation=<none/Bulk1dGlobal/Bulk1dLocal>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermSubsidence=<none/ConstantDivergenceLocal/ConstantDivergenceGlobal>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermTransport=<constant/powerlaw/sutherland/Airwater/AirwaterSimplified>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermChemistry=<none/quadratic/layeredrelaxation/ozone>')
  CALL TLAB_WRITE_ASCII(bakfile, '#TermRandom=<value>') 
  CALL TLAB_WRITE_ASCII(bakfile, '#SpaceOrder=<CompactJacobian4/CompactJacobian6/CompactJacpenta6/CompactJacobian8/CompactDirect6>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ComModeITranspose=<none,asynchronous,sendrecv>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ComModeKTranspose=<none,asynchronous,sendrecv>')

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'FileFormat', 'MpiIO', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'mpiio'    ) THEN; imode_files = IO_MPIIO
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'netcdf'   ) THEN; imode_files = IO_NETCDF
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'none'     ) THEN; imode_files = IO_NOFILE
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Wrong Main.FileFormat.')
     CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'FileType', 'Double', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'double'   ) THEN; imode_precision_files = IO_TYPE_DOUBLE
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'single'   ) THEN; imode_precision_files = IO_TYPE_SINGLE
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Wrong Main.FileType.')
     CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

  CALL SCANINIINT(bakfile, inifile, 'Main', 'VerbosityLevel', '1', imode_verbosity)

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Type', 'temporal', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'temporal' ) THEN; imode_sim = DNS_MODE_TEMPORAL
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'spatial'  ) THEN; imode_sim = DNS_MODE_SPATIAL
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Entry Main.Type must be temporal or spatial')
     CALL TLAB_STOP(DNS_ERROR_SIMTYPE)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'CalculateFlow', 'yes', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_flow = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_flow = 0
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Entry Main.CalculateFlow must be yes or no')
     CALL TLAB_STOP(DNS_ERROR_CALCFLOW)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'CalculateScalar', 'yes', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_scal = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_scal = 0
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Entry Main.CalculateScalar must be yes or no')
     CALL TLAB_STOP(DNS_ERROR_CALCSCALAR)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'CalculateParticle', 'no', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; icalc_part = 1
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; icalc_part = 0
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Entry Main.CalculateParticle must be yes or no')
     CALL TLAB_STOP(DNS_ERROR_CALCPARTICLE)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Equations', 'internal', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'total'          ) THEN; imode_eqns = DNS_EQNS_TOTAL
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'internal'       ) THEN; imode_eqns = DNS_EQNS_INTERNAL
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'incompressible' ) THEN; imode_eqns = DNS_EQNS_INCOMPRESSIBLE
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'anelastic'      ) THEN; imode_eqns = DNS_EQNS_ANELASTIC
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Wrong entry Main.Equations option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL ) THEN; ifourier = 1
  ELSE;                                         ifourier = 0; ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'Mixture', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'          ) THEN; imixture = MIXT_TYPE_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'air'           ) THEN; imixture = MIXT_TYPE_AIR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'airvapor'      ) THEN; imixture = MIXT_TYPE_AIRVAPOR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'airwater'      ) THEN; imixture = MIXT_TYPE_AIRWATER
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'airwaterlinear') THEN; imixture = MIXT_TYPE_AIRWATER_LINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'chemkin'       ) THEN; imixture = MIXT_TYPE_CHEMKIN ! use chemkin file to read the data
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong entry Main.Mixture model.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermAdvection', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'          ) THEN; iadvection = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'divergence'    ) THEN; iadvection = EQNS_DIVERGENCE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'skewsymmetric' ) THEN; iadvection = EQNS_SKEWSYMMETRIC
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'convective'    ) THEN; iadvection = EQNS_CONVECTIVE
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermAdvection option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermViscous', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'       ) THEN; iviscous = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'divergence' ) THEN; iviscous = EQNS_DIVERGENCE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit'   ) THEN; iviscous = EQNS_EXPLICIT
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermViscous option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermDiffusion', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'       ) THEN; idiffusion = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'divergence' ) THEN; idiffusion = EQNS_DIVERGENCE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit'   ) THEN; idiffusion = EQNS_EXPLICIT
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermDiffusion option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermBodyForce', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .EQ. 'none'        ) THEN; buoyancy%type = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'explicit'    ) THEN; buoyancy%type = EQNS_EXPLICIT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'homogeneous' ) THEN; buoyancy%type = EQNS_BOD_HOMOGENEOUS
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'linear'      ) THEN; buoyancy%type = EQNS_BOD_LINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'bilinear'    ) THEN; buoyancy%type = EQNS_BOD_BILINEAR
  ELSE IF ( TRIM(ADJUSTL(sRes)) .EQ. 'quadratic'   ) THEN; buoyancy%type = EQNS_BOD_QUADRATIC
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermBodyForce option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermCoriolis', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'       ) THEN; coriolis%type = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit'   ) THEN; coriolis%type = EQNS_EXPLICIT
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'normalized' ) THEN; coriolis%type = EQNS_COR_NORMALIZED
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermCoriolis option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermRadiation', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'          ) THEN; radiation%type = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'bulk1dglobal'  ) THEN; radiation%type = EQNS_RAD_BULK1D_GLOBAL
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'bulk1dlocal'   ) THEN; radiation%type = EQNS_RAD_BULK1D_LOCAL
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermRadiation option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Main', 'TermSubsidence', 'None', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'                    ) THEN; subsidence%type = EQNS_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'constantdivergencelocal' ) THEN; subsidence%type = EQNS_SUB_CONSTANT_LOCAL
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'constantdivergenceglobal') THEN; subsidence%type = EQNS_SUB_CONSTANT_GLOBAL
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermSubsidence option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
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
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'quadratic3'       ) THEN; chemistry%type = EQNS_CHEM_QUADRATIC3;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'layeredrelaxation') THEN; chemistry%type = EQNS_CHEM_LAYEREDRELAXATION;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'ozone'            ) THEN; chemistry%type = EQNS_CHEM_OZONE;
  ELSE;                                                        chemistry%type = EQNS_NONE; ENDIF

! -------------------------------------------------------------------
  CALL SCANINIREAL(bakfile, inifile, 'Main', 'TermRandom', '0.0', dummy)
  IF ( ABS(dummy) .GT. 0.0 ) THEN
     random%type = EQNS_RAND_MULTIPLY
     random%parameters=dummy 
     random%active(1:3)=.TRUE.
  ELSE
     random%type = EQNS_NONE
     random%active(1:3)=.FALSE. 
  ENDIF
     
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Main', 'SpaceOrder', 'void', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .EQ. 'compactjacobian4' ) THEN; imode_fdm = FDM_COM4_JACOBIAN;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'compactjacobian6' ) THEN; imode_fdm = FDM_COM6_JACOBIAN;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'compactjacpenta6' ) THEN; imode_fdm = FDM_COM6_JACPENTA;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'compactjacobian8' ) THEN; imode_fdm = FDM_COM8_JACOBIAN;
  ELSEIF ( TRIM(ADJUSTL(sRes)) .EQ. 'compactdirect6'   ) THEN; imode_fdm = FDM_COM6_DIRECT;
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong SpaceOrder option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  g(1:3)%mode_fdm = imode_fdm

! -------------------------------------------------------------------
#ifdef USE_MPI
  CALL SCANINICHAR(bakfile,inifile, 'Main', 'ComModeITranspose', 'asynchronous',sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'none')         THEN; ims_trp_mode_i = TLAB_MPI_TRP_NONE
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'asynchronous') THEN; ims_trp_mode_i = TLAB_MPI_TRP_ASYNCHRONOUS
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'sendrecv'    ) THEN; ims_trp_mode_i = TLAB_MPI_TRP_SENDRECV
  ELSE
     CALL TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong ComModeITranspose option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile,inifile, 'Main', 'ComModeKTranspose', 'asynchronous',sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'none')         THEN; ims_trp_mode_k = TLAB_MPI_TRP_NONE
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'asynchronous') THEN; ims_trp_mode_k = TLAB_MPI_TRP_ASYNCHRONOUS
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'sendrecv'    ) THEN; ims_trp_mode_k = TLAB_MPI_TRP_SENDRECV
  ELSE
     CALL TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong ComModeKTranspose option.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF
#endif

! ###################################################################
! Pressure staggering 
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Staggering]')
  CALL TLAB_WRITE_ASCII(bakfile, '#StaggerHorizontalPressure=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#FilterVerticalPressure=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#FilterParameter=<value>')

  CALL SCANINICHAR(bakfile, inifile, 'Staggering', 'StaggerHorizontalPressure', 'no', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; istagger = 1; CALL TLAB_WRITE_ASCII(lfile, 'Horizontal staggering of the pressure along Ox and Oz.')
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; istagger = 0
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Entry Main. StaggerHorizontalPressure must be yes or no')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Staggering', 'FilterVerticalPressure', 'no', sRes)
  IF     ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; ivfilter = 1;  CALL TLAB_WRITE_ASCII(lfile, 'Vertical filtering of the pressure and dpdy along Oy.')
  ELSEIF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; ivfilter = 0
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Entry Main. FilterVerticalPressure must be yes or no')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINIREAL(bakfile, inifile, 'Staggering', 'Filterparameter', '4.0', vfilter_param)

! Consistency check
  IF ( istagger .EQ. 1 ) THEN
     IF ( .NOT. ( (imode_eqns .EQ. DNS_EQNS_INCOMPRESSIBLE) .OR. (imode_eqns .EQ. DNS_EQNS_ANELASTIC) ) ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Horizontal pressure staggering only implemented for anelastic or incompressible mode.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
     IF ( .NOT. ( (iadvection .EQ. EQNS_CONVECTIVE) .OR. (iadvection .EQ. EQNS_SKEWSYMMETRIC) ) ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Horizontal pressure staggering not implemented for current advection scheme.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
  ENDIF
  IF ( ivfilter .EQ. 1 ) THEN
     IF ( .NOT. (istagger .EQ. 1) ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Vertical pressure filtering only in combination with horizontal pressure staggering.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

! ###################################################################
! Dimensionles parameters
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile,  '#')
  CALL TLAB_WRITE_ASCII(bakfile,  '#[Parameters]')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Reynolds=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Prandtl=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Froude=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Rossby=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Mach=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Gama=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Schmidt=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Damkohler=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Stokes=<value>')
  CALL TLAB_WRITE_ASCII(bakfile,  '#Settling=<value>')

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
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Schmidt and Damkholer sizes do not match.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

! ###################################################################
! Buoyancy
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[BodyForce]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Vector=<Gx,Gy,Gz>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

  buoyancy%vector = C_0_R; buoyancy%active = .FALSE.
  IF ( buoyancy%type .NE. EQNS_NONE ) THEN
     CALL SCANINICHAR(bakfile, inifile, 'BodyForce', 'Vector', '0.0,-1.0,0.0', sRes)
     idummy = 3
     CALL LIST_REAL(sRes, idummy, buoyancy%vector)

     IF ( ABS(buoyancy%vector(1)) .GT. C_0_R ) THEN; buoyancy%active(1) = .TRUE.; CALL TLAB_WRITE_ASCII(lfile, 'Body force along Ox.'); ENDIF
     IF ( ABS(buoyancy%vector(2)) .GT. C_0_R ) THEN; buoyancy%active(2) = .TRUE.; CALL TLAB_WRITE_ASCII(lfile, 'Body force along Oy.'); ENDIF
     IF ( ABS(buoyancy%vector(3)) .GT. C_0_R ) THEN; buoyancy%active(3) = .TRUE.; CALL TLAB_WRITE_ASCII(lfile, 'Body force along Oz.'); ENDIF

     IF ( froude .GT. C_0_R ) THEN
           buoyancy%vector(:) = buoyancy%vector(:) /froude ! adding the froude number into de vector g
     ELSE
        CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Froude number must be nonzero if buoyancy is retained.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
     ENDIF

     buoyancy%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'BodyForce', 'Parameters', '0.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, buoyancy%parameters)

  ENDIF

! ###################################################################
! Rotation
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Rotation]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Vector=<Fx,Fy,Fz>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

  coriolis%vector(:) = C_0_R; coriolis%active = .FALSE.
  IF ( coriolis%type .NE. EQNS_NONE ) THEN
     CALL SCANINICHAR(bakfile, inifile, 'Rotation', 'Vector', '0.0,1.0,0.0', sRes)
     idummy = 3
     CALL LIST_REAL(sRes, idummy, coriolis%vector)

     IF ( ABS(coriolis%vector(1)) .GT. C_0_R ) THEN; coriolis%active(2) = .TRUE.; coriolis%active(3) = .TRUE.; CALL TLAB_WRITE_ASCII(lfile, 'Angular velocity along Ox.'); ENDIF
     IF ( ABS(coriolis%vector(2)) .GT. C_0_R ) THEN; coriolis%active(3) = .TRUE.; coriolis%active(1) = .TRUE.; CALL TLAB_WRITE_ASCII(lfile, 'Angular velocity along Oy.'); ENDIF
     IF ( ABS(coriolis%vector(3)) .GT. C_0_R ) THEN; coriolis%active(1) = .TRUE.; coriolis%active(2) = .TRUE.; CALL TLAB_WRITE_ASCII(lfile, 'Angular velocity along Oz.'); ENDIF

     IF ( rossby .GT. C_0_R ) THEN
        coriolis%vector(:) = coriolis%vector(:) /rossby ! adding the rossby number into the vector
     ELSE
        CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Rossby number must be nonzero if coriolis is retained.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
     ENDIF

     coriolis%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Rotation', 'Parameters', '0.0,1.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, coriolis%parameters)

     IF ( coriolis%parameters(2) .EQ. C_0_R ) THEN
        CALL TLAB_WRITE_ASCII(lfile,C_FILE_LOC//'. Default normalized geostrophic velocity set to one.')
        coriolis%parameters(2) = C_1_R
     ENDIF

  ENDIF

! Consistency check
  IF ( coriolis%type .EQ. EQNS_COR_NORMALIZED ) THEN
     IF ( coriolis%active(2) ) THEN
        CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. TermCoriolis option only allows for angular velocity along Oy.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

! ###################################################################
! Radiation
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Radiation]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

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
! Subsidence
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Subsidence]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

  subsidence%active = .FALSE.
  IF ( subsidence%type .NE. EQNS_NONE ) THEN
     subsidence%active = .TRUE.

     subsidence%parameters(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Subsidence', 'Parameters', '0.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, subsidence%parameters)

  ENDIF

! This subsidence type is implemented in opr_burgers_y only
! to speed up calculation
  IF ( subsidence%type .EQ. EQNS_SUB_CONSTANT_LOCAL ) subsidence%active = .FALSE.

! ###################################################################
! Transport
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Transport]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Exponent=<value>')

  transport%active = .FALSE.
  IF ( transport%type .NE. EQNS_NONE ) THEN
     transport%parameters(:) = C_1_R ! default values
     CALL SCANINICHAR(bakfile, inifile, 'Transport', 'Parameters', 'void', sRes)
     IF ( TRIM(ADJUSTL(sRes)) .NE. 'void' ) THEN
        idummy = MAX_PROF
        CALL LIST_REAL(sRes, idummy, transport%parameters)
     ENDIF

     IF ( settling .GT. C_0_R ) THEN
        transport%parameters = transport%parameters *settling ! adding the settling number in the parameter definitions
     ELSE
        CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Settling number must be nonzero if transport is retained.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
     ENDIF

     IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
        transport%active = .TRUE. ! All scalars are affected

        CALL SCANINIREAL(bakfile, inifile, 'Transport', 'Exponent', '0.0', transport%auxiliar(1))
     ENDIF

  ENDIF

! ###################################################################
! Chemistry
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Chemistry]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

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
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Thermodynamics]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Nondimensional=<yes,no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

  CALL SCANINICHAR(bakfile, inifile, 'Thermodynamics', 'Nondimensional', 'yes', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; nondimensional = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; nondimensional = .FALSE.
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Thermodynamics.Nondimensional')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  IF ( imixture .NE. EQNS_NONE ) THEN
     thermo_param(:) = C_0_R
     CALL SCANINICHAR(bakfile, inifile, 'Thermodynamics', 'Parameters', '1.0', sRes)
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, thermo_param)

  ENDIF

  IF ( imixture .NE. MIXT_TYPE_CHEMKIN ) THEN
    CALL SCANINICHAR(bakfile, inifile, 'Thermodynamics', 'ChemkinFile', 'none', chemkin_file)
  ENDIF

  IF ( imixture .EQ. MIXT_TYPE_AIRWATER ) THEN
     CALL SCANINIREAL(bakfile, inifile, 'Thermodynamics', 'SmoothFactor', '0.1', dsmooth)
  ENDIF


! ###################################################################
! Grid Parameters
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Grid]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Imax=<imax>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Imax(*)=<imax_proc>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Jmax=<jmax>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Kmax=<kmax>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Kmax(*)=<kmax_proc>')
  CALL TLAB_WRITE_ASCII(bakfile, '#XUniform=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#YUniform=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ZUniform=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#XPeriodic=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#YPeriodic=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ZPeriodic=<yes/no>')

  CALL SCANINIINT(bakfile, inifile, 'Grid', 'Imax', '0', g(1)%size)
  CALL SCANINIINT(bakfile, inifile, 'Grid', 'Jmax', '0', g(2)%size)
  CALL SCANINIINT(bakfile, inifile, 'Grid', 'Kmax', '0', g(3)%size)

! default
  imax = g(1)%size
  jmax = g(2)%size
  kmax = g(3)%size

  g(1)%name = 'x'
  g(2)%name = 'y'
  g(3)%name = 'z'

! -------------------------------------------------------------------
! Domain decomposition in parallel mode
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro .GT. 1 ) THEN
     CALL SCANINIINT(bakfile, inifile, 'Grid', 'Kmax(*)', '-1', kmax)
     IF ( kmax .GT. 0 .AND. MOD(g(3)%size,kmax) .EQ. 0 ) THEN
        ims_npro_k = g(3)%size/kmax
     ELSE
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Input kmax incorrect')
        CALL TLAB_STOP(DNS_ERROR_KMAXTOTAL)
     ENDIF

     CALL SCANINIINT(bakfile, inifile, 'Grid', 'Imax(*)', '-1', imax)
     IF ( imax .GT. 0 .AND. MOD(g(1)%size,imax) .EQ. 0 ) THEN
        ims_npro_i = g(1)%size/imax
     ELSE
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Input imax incorrect')
        CALL TLAB_STOP(DNS_ERROR_KMAXTOTAL)
     ENDIF

     IF ( ims_npro_i*ims_npro_k .EQ. ims_npro ) THEN ! check
        WRITE(lstr,*) ims_npro_i; WRITE(sRes,*) ims_npro_k
        lstr = TRIM(ADJUSTL(lstr))//'x'//TRIM(ADJUSTL(sRes))
        CALL TLAB_WRITE_ASCII(lfile, 'Initializing domain partition '//TRIM(ADJUSTL(lstr)))
     ELSE
        CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Inconsistency in total number of PEs')
        CALL TLAB_STOP(DNS_ERROR_KMAXTOTAL)
     ENDIF

  ENDIF

#endif

! -------------------------------------------------------------------
! Uniform
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'XUniform', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; g(1)%uniform = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; g(1)%uniform = .FALSE.
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Uniform X grid')
     CALL TLAB_STOP(DNS_ERROR_UNIFORMX)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'YUniform', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; g(2)%uniform = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; g(2)%uniform = .FALSE.
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Uniform Y grid')
     CALL TLAB_STOP(DNS_ERROR_UNIFORMY)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'ZUniform', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; g(3)%uniform = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; g(3)%uniform = .FALSE.
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Uniform Z grid')
     CALL TLAB_STOP(DNS_ERROR_UNIFORMZ)
  ENDIF

! -------------------------------------------------------------------
! Periodic
! -------------------------------------------------------------------
  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'XPeriodic', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; g(1)%periodic = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; g(1)%periodic = .FALSE.
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Periodic X grid')
     CALL TLAB_STOP(DNS_ERROR_IBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'YPeriodic', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; g(2)%periodic = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; g(2)%periodic = .FALSE.
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Periodic Y grid')
     CALL TLAB_STOP(DNS_ERROR_JBC)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Grid', 'ZPeriodic', 'void', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'yes' ) THEN; g(3)%periodic = .TRUE.
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'no'  ) THEN; g(3)%periodic = .FALSE.
  ELSE
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Periodic Z grid')
     CALL TLAB_STOP(DNS_ERROR_KBC)
  ENDIF

! ###################################################################
! Dealising (a filter type)
! ###################################################################
    CALL TLAB_WRITE_ASCII(bakfile, '#')
    CALL TLAB_WRITE_ASCII(bakfile, '#[Dealising]')
    CALL TLAB_WRITE_ASCII(bakfile, '#Type=<none/compact>')
    CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<values>')
    CALL TLAB_WRITE_ASCII(bakfile, '#ActiveX=<yes/no>')
    CALL TLAB_WRITE_ASCII(bakfile, '#ActiveY=<yes/no>')
    CALL TLAB_WRITE_ASCII(bakfile, '#ActiveZ=<yes/no>')

    Dealiasing(:)%size       = g(:)%size
    Dealiasing(:)%periodic   = g(:)%periodic
    Dealiasing(:)%uniform    = g(:)%uniform
    CALL SCANINICHAR(bakfile, inifile, 'Dealiasing', 'Type', 'none', sRes)
    IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; Dealiasing(:)%type = DNS_FILTER_NONE
    ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'compact'   ) THEN; Dealiasing(:)%type = DNS_FILTER_COMPACT
        Dealiasing(:)%parameters(1) = 0.49 ! default alpha value
        Dealiasing(:)%inb_filter    = 10
        Dealiasing(:)%BcsMin = DNS_FILTER_BCS_BIASED
        Dealiasing(:)%BcsMax = DNS_FILTER_BCS_BIASED
    ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'compactcutoff'   ) THEN; Dealiasing(:)%type = DNS_FILTER_COMPACT_CUTOFF
        Dealiasing(:)%inb_filter    = 7
        Dealiasing(2)%type = DNS_FILTER_COMPACT ! nonuniform version not yet implemented; fall back to compact
        Dealiasing(2)%parameters(1) = 0.49
        Dealiasing(2)%inb_filter    = 10
        Dealiasing(2)%BcsMin = DNS_FILTER_BCS_BIASED
        Dealiasing(2)%BcsMax = DNS_FILTER_BCS_BIASED
    ENDIF

    ! Boundary conditions correction
    DO ig = 1,3
        IF ( FilterDomain(ig)%periodic ) THEN
            Dealiasing(ig)%BcsMin = DNS_FILTER_BCS_PERIODIC
            Dealiasing(ig)%BcsMax = DNS_FILTER_BCS_PERIODIC
        ENDIF
    ENDDO

    CALL SCANINICHAR(bakfile, inifile, 'Dealiasing', 'Parameters', 'void', sRes)
    IF ( TRIM(ADJUSTL(sRes)) .NE. 'void' ) THEN
        idummy = MAX_PROF
        CALL LIST_REAL(sRes, idummy, Dealiasing(1)%parameters(:) )
        IF ( idummy .LT. 3 ) & ! Fill 3 directions; if global, filled information is unused
        Dealiasing(1)%parameters(idummy+1:3) = FilterDomain(1)%parameters(idummy)
        DO ig = 2,3
            Dealiasing(ig)%parameters(1) = FilterDomain(1)%parameters(ig)
        ENDDO
    ENDIF

    CALL SCANINICHAR(bakfile, inifile, 'Dealiasing', 'ActiveX', 'yes', sRes)
    IF ( TRIM(ADJUSTL(sRes)) .EQ. 'no' ) Dealiasing(1)%type = DNS_FILTER_NONE
    CALL SCANINICHAR(bakfile, inifile, 'Dealiasing', 'ActiveY', 'yes', sRes)
    IF ( TRIM(ADJUSTL(sRes)) .EQ. 'no' ) Dealiasing(2)%type = DNS_FILTER_NONE
    CALL SCANINICHAR(bakfile, inifile, 'Dealiasing', 'ActiveZ', 'yes', sRes)
    IF ( TRIM(ADJUSTL(sRes)) .EQ. 'no' ) Dealiasing(3)%type = DNS_FILTER_NONE

    ! Further control
    DO ig = 1,3
        IF ( Dealiasing(ig)%size .EQ. 1 ) Dealiasing(ig)%type = DNS_FILTER_NONE
    ENDDO
   
#ifdef USE_MPI
    Dealiasing(1)%mpitype = TLAB_MPI_I_PARTIAL
    Dealiasing(3)%mpitype = TLAB_MPI_K_PARTIAL
#endif

! ###################################################################
! Filter
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Filter]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Type=<none/compact/explicit6/explicit4/adm/helmholtz/SpectralBand/SpectralErf/tophat>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<values>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ActiveX=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ActiveY=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ActiveZ=<yes/no>')
  CALL TLAB_WRITE_ASCII(bakfile, '#BcsJmin=<free,solid,zero>')
  CALL TLAB_WRITE_ASCII(bakfile, '#BcsJmax=<free,solid,zero>')

  FilterDomain(:)%size       = g(:)%size
  FilterDomain(:)%periodic   = g(:)%periodic
  FilterDomain(:)%uniform    = g(:)%uniform
  FilterDomain(:)%inb_filter = 0        ! default array size
  default                    = 'biased' ! default boundary condition

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'Type', 'none', sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'none'      ) THEN; FilterDomain(:)%type = DNS_FILTER_NONE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'compact'   ) THEN; FilterDomain(:)%type = DNS_FILTER_COMPACT
     FilterDomain(:)%parameters(1) = 0.49 ! default alpha value
     FilterDomain(:)%inb_filter    = 10
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit6' ) THEN; FilterDomain(:)%type = DNS_FILTER_6E
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'explicit4' ) THEN; FilterDomain(:)%type = DNS_FILTER_4E
     FilterDomain(:)%inb_filter = 5
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'adm'       ) THEN; FilterDomain(:)%type = DNS_FILTER_ADM
     FilterDomain(:)%inb_filter = 5
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'tophat'    ) THEN; FilterDomain(:)%type = DNS_FILTER_TOPHAT
     FilterDomain(:)%parameters(1) = 2    ! default filter size (in grid-step units)
     default                    = 'free'
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'spectralcutoff' ) THEN; FilterDomain(:)%type = DNS_FILTER_BAND
     ! The frequency interval is (Parameter1, Parameter2)
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'spectralerf'    ) THEN; FilterDomain(:)%type = DNS_FILTER_ERF
     ! Parameter1 is the transition wavenumber in physical units:
     ! >0: high-pass filter
     ! <0; low-pass filter
     ! Parameter2 is the characteristic width--in log units (relative to domain size)'
     FilterDomain(:)%parameters(3) = C_1_R    ! used to normalise wavenumbers in z-direction
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'helmholtz'      ) THEN; FilterDomain(:)%type = DNS_FILTER_HELMHOLTZ
     FilterDomain(:)%parameters(1) = C_1_R    ! default filter size
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Wrong Filter.Type.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  ! Boundary conditions
  DO ig = 1,3
     IF ( FilterDomain(ig)%periodic ) THEN
        FilterDomain(ig)%BcsMin = DNS_FILTER_BCS_PERIODIC
        FilterDomain(ig)%BcsMax = DNS_FILTER_BCS_PERIODIC
     ENDIF
  ENDDO

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'BcsJmin', TRIM(ADJUSTL(default)), sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'periodic'  ) THEN; FilterDomain(2)%BcsMin = DNS_FILTER_BCS_PERIODIC
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'biased'    ) THEN; FilterDomain(2)%BcsMin = DNS_FILTER_BCS_BIASED
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'free'      ) THEN; FilterDomain(2)%BcsMin = DNS_FILTER_BCS_FREE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'solid'     ) THEN; FilterDomain(2)%BcsMin = DNS_FILTER_BCS_SOLID
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; FilterDomain(2)%BcsMin = DNS_FILTER_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; FilterDomain(2)%BcsMin = DNS_FILTER_BCS_NEUMANN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'zero'      ) THEN; FilterDomain(2)%BcsMin = DNS_FILTER_BCS_ZERO
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Wrong Filter.BcsJmin.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'BcsJmax', TRIM(ADJUSTL(default)), sRes)
  IF      ( TRIM(ADJUSTL(sRes)) .eq. 'periodic'  ) THEN; FilterDomain(2)%BcsMax = DNS_FILTER_BCS_PERIODIC
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'biased'    ) THEN; FilterDomain(2)%BcsMax = DNS_FILTER_BCS_BIASED
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'free'      ) THEN; FilterDomain(2)%BcsMax = DNS_FILTER_BCS_FREE
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'solid'     ) THEN; FilterDomain(2)%BcsMax = DNS_FILTER_BCS_SOLID
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; FilterDomain(2)%BcsMax = DNS_FILTER_BCS_DIRICHLET
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; FilterDomain(2)%BcsMax = DNS_FILTER_BCS_NEUMANN
  ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'zero'      ) THEN; FilterDomain(2)%BcsMax = DNS_FILTER_BCS_ZERO
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Wrong Filter.BcsJmax.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'Parameters', 'void', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .NE. 'void' ) THEN
     idummy = MAX_PROF
     CALL LIST_REAL(sRes, idummy, FilterDomain(1)%parameters(:) )
     IF ( idummy .LT. 3 ) & ! Fill 3 directions; if global, filled information is unused
          FilterDomain(1)%parameters(idummy+1:3) = FilterDomain(1)%parameters(idummy)
     DO ig = 2,3
        FilterDomain(ig)%parameters(1) = FilterDomain(1)%parameters(ig)
     ENDDO
  ENDIF

  CALL SCANINIINT(bakfile, inifile, 'Filter', 'Repeat', '1', idummy)
  IF ( idummy .GT. 0 ) THEN
     FilterDomain(:)%repeat = idummy
  ELSE
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Entry Filter.Repeat must be positive.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  ENDIF

  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'ActiveX', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'no' ) FilterDomain(1)%type = DNS_FILTER_NONE
  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'ActiveY', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'no' ) FilterDomain(2)%type = DNS_FILTER_NONE
  CALL SCANINICHAR(bakfile, inifile, 'Filter', 'ActiveZ', 'yes', sRes)
  IF ( TRIM(ADJUSTL(sRes)) .EQ. 'no' ) FilterDomain(3)%type = DNS_FILTER_NONE

! To eventually allow for control field by field
  FilterDomainActive(:) = .TRUE.

! Further control
  DO ig = 1,3
     IF ( FilterDomain(ig)%size .EQ. 1 ) FilterDomain(ig)%type = DNS_FILTER_NONE
     IF ( FilterDomain(ig)%type .EQ. DNS_FILTER_TOPHAT .AND. &
          FilterDomain(ig)%parameters(1) .EQ. 0              ) FilterDomain(ig)%type = DNS_FILTER_NONE

     IF ( FilterDomain(ig)%type .EQ. DNS_FILTER_TOPHAT ) THEN
        IF ( MOD(INT(FilterDomain(is)%parameters(1)),2) .NE. 0 ) THEN
           CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Tophat filter size must be even.')
           CALL TLAB_STOP(DNS_ERROR_PARAMETER)
        ENDIF
        FilterDomain(ig)%inb_filter = INT(FilterDomain(ig)%parameters(1)) +1
     ENDIF

  ENDDO

#ifdef USE_MPI
  FilterDomain(1)%mpitype = TLAB_MPI_I_PARTIAL
  FilterDomain(3)%mpitype = TLAB_MPI_K_PARTIAL
#endif

! ###################################################################
! Statistics Control
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Statistics]')
  CALL TLAB_WRITE_ASCII(bakfile, '#IAvera=<plane1,plane2,...>')

  nstatavg = MAX_STATS_SPATIAL
  CALL SCANINICHAR(bakfile, inifile, 'Statistics', 'IAvera', '1', sRes)
  CALL LIST_INTEGER(sRes, nstatavg, statavg)

! ###################################################################
! Flow physical properties of the system
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Flow]')
  CALL TLAB_WRITE_ASCII(bakfile, '#VelocityX=<mean velocity X>')
  CALL TLAB_WRITE_ASCII(bakfile, '#VelocityY=<mean velocity Y>')
  CALL TLAB_WRITE_ASCII(bakfile, '#VelocityZ=<mean velocity Z>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Pressure=<mean pressure>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Density=<mean density>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Temperature=<mean temperature>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ProfileVelocity=<None/Linear/Tanh/Erf/Ekman/EkmanP/Parabolic>')
  CALL TLAB_WRITE_ASCII(bakfile, '#YCoorVelocity=<Relative Y reference point>')
  CALL TLAB_WRITE_ASCII(bakfile, '#DiamVelocity=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ThickVelocity=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#DeltaVelocity=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ProfileDensity=<None/Linear/Tanh/Erf>')
  CALL TLAB_WRITE_ASCII(bakfile, '#YCoorDensity=<Relative Y reference point>')
  CALL TLAB_WRITE_ASCII(bakfile, '#DiamDensity=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ThickDensity=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#DeltaDensity=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ProfileTemperature=<None/Linear/Tanh/Erf>')
  CALL TLAB_WRITE_ASCII(bakfile, '#YCoorTemperature=<Relative Y reference point>')
  CALL TLAB_WRITE_ASCII(bakfile, '#DiamTemperature=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#ThickTemperature=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#DeltaTemperature=<value>')


  call PROFILE_READBLOCK(bakfile, inifile, 'Flow', 'Velocity', qbg(1))  
! originally, all velocity data was contained in block 'Velocity' except for the mean value
  CALL SCANINIREAL(bakfile, inifile, 'Flow', 'VelocityX',     '0.0', qbg(1)%mean )
  call PROFILE_READBLOCK(bakfile, inifile, 'Flow', 'VelocityY', qbg(2))
  call PROFILE_READBLOCK(bakfile, inifile, 'Flow', 'VelocityZ', qbg(3))

  ! Consistency check
 IF ( qbg(1)%type .EQ. PROFILE_EKMAN_U .OR. qbg(1)%type .EQ. PROFILE_EKMAN_U_P ) THEN
     qbg(3)%type  = PROFILE_EKMAN_V
     qbg(3)%ymean = qbg(1)%ymean
     qbg(3)%thick = qbg(1)%thick
     qbg(3)%delta = qbg(1)%delta
  ENDIF

  call PROFILE_READBLOCK(bakfile, inifile, 'Flow', 'Density', rbg)
  call PROFILE_READBLOCK(bakfile, inifile, 'Flow', 'Temperature', tbg)  ! temperature/enthalpy
  call PROFILE_READBLOCK(bakfile, inifile, 'Flow', 'Pressure', pbg)

! consistency check
  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     IF ( rbg%type .EQ. PROFILE_NONE .AND. tbg%type .EQ. PROFILE_NONE ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Specify density or temperature.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
     ENDIF

     IF ( rbg%type .NE. PROFILE_NONE .AND. tbg%type .NE. PROFILE_NONE ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Specify only density or only temperature.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
     ENDIF
  ENDIF

! -------------------------------------------------------------------
! Spatial case
! Thickness evolutions delta_i/diam_i=a*(x/diam_i+b)
! -------------------------------------------------------------------
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickAVelocity=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickBVelocity=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#FluxVelocity=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickADensity=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickBDensity=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#FluxDensity=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickATemperature=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickBTemperature=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#FluxTemperature=<value>')

! Bradbury profile is the default (x0=a*b)
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickAVelocity', '0.1235', qbg(1)%parameters(2))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBVelocity', '-0.873', qbg(1)%parameters(3))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'FluxVelocity',   '0.96',   qbg(1)%parameters(4))

! Ramaprian is the default (x0=a*b)
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickADensity', '0.14', rbg%parameters(2))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBDensity', '2.0',  rbg%parameters(3))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'FluxDensity',   '0.94', rbg%parameters(4))

     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickATemperature', '0.14', tbg%parameters(2))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBTemperature', '2.0',  tbg%parameters(3))
     CALL SCANINIREAL(bakfile, inifile, 'Flow', 'FluxTemperature',   '0.94', tbg%parameters(4))

  ENDIF

! ###################################################################
! Scalars physical properties of the system
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')
  CALL TLAB_WRITE_ASCII(bakfile, '#[Scalar]')
  CALL TLAB_WRITE_ASCII(bakfile, '#Profile=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#YCoor=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Diam=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Thick=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Mean=<value>')
  CALL TLAB_WRITE_ASCII(bakfile, '#Delta=<value>')

  DO is = 1,MAX_NSP
     WRITE(lstr,*) is
     call PROFILE_READBLOCK(bakfile, inifile, 'Scalar', 'Scalar'//TRIM(ADJUSTL(lstr)), sbg(is))
  ENDDO

! -------------------------------------------------------------------
! Spatial case
! Thickness evolutions delta_i/diam_i=a*(x/diam_i+b)
! Ramaprian is the default (x0=a*b)
! -------------------------------------------------------------------
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickA=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#ThickB=<value>')
     CALL TLAB_WRITE_ASCII(bakfile, '#Flux=<value>')

     DO is = 1, MAX_NSP
        WRITE(lstr,*) is; lstr='ThickA'//TRIM(ADJUSTL(lstr))
        CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.14', sbg(is)%parameters(2))
        WRITE(lstr,*) is; lstr='ThickB'//TRIM(ADJUSTL(lstr))
        CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '2.0',  sbg(is)%parameters(3))
        WRITE(lstr,*) is; lstr='Flux'//TRIM(ADJUSTL(lstr))
        CALL SCANINIREAL(bakfile, inifile, 'Scalar', TRIM(ADJUSTL(lstr)), '0.94', sbg(is)%parameters(4))
     ENDDO

  ENDIF

! ###################################################################
! Final consistency check and initialization
! ###################################################################
  CALL TLAB_WRITE_ASCII(bakfile, '#')

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
        CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Incorrect number of Schmidt numbers.')
        CALL TLAB_STOP(DNS_ERROR_OPTION)
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

  if ( .not. nondimensional ) then
    MRATIO = C_1_R
    GRATIO = C_1_R
  end if
  
  IF ( imode_eqns .EQ. DNS_EQNS_ANELASTIC .AND. &
       imixture .NE. MIXT_TYPE_AIR .AND. imixture .NE. MIXT_TYPE_AIRVAPOR .AND. imixture .NE. MIXT_TYPE_AIRWATER ) THEN
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Incorrect mixture type.')
     CALL TLAB_STOP(DNS_ERROR_OPTION)
  END IF

  IF ( buoyancy%type .EQ. EQNS_BOD_LINEAR   .OR. &
       buoyancy%type .EQ. EQNS_BOD_BILINEAR .OR. &
       buoyancy%type .EQ. EQNS_BOD_QUADRATIC ) THEN
     IF ( inb_scal .EQ. 0 ) THEN
        CALL TLAB_WRITE_ASCII(wfile,C_FILE_LOC//'. Zero scalars; setting TermBodyForce equal to none.')
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
           CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. AirWater requires at least first 2 Damkholer numbers zero.')
           CALL TLAB_STOP(DNS_ERROR_OPTION)
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

  inb_wrk1d = 18
  IF ( imode_sim .EQ. DNS_MODE_SPATIAL ) THEN; inb_wrk2d = 11
  ELSE;                                        inb_wrk2d =  2; ENDIF

  isize_wrk1d = MAX(g(1)%size,MAX(g(2)%size,g(3)%size))
  isize_wrk2d = MAX(imax*jmax, MAX(imax*kmax,jmax*kmax)  )

! grid array
  DO is = 1,3
     g(is)%inb_grid = 1                  ! Nodes
     g(is)%inb_grid = g(is)%inb_grid &
                    + 2              &   ! Jacobians of first- and second-order derivatives
                    + 2                  ! 1/dx and 1/dx**2 used in time-step stability constraint

     IF ( g(is)%periodic ) THEN
        g(is)%inb_grid = g(is)%inb_grid  &
                       + 7               & ! LU decomposition 1. order
                       + 5               & ! LU decomposition 2. order
                       + 5 *(1+inb_scal) & ! LU decomposition 2. order with diffusivities
                       + 2                 ! modified wavenumbers
     ELSE
        g(is)%inb_grid = g(is)%inb_grid  &
                       + 5 *4            & ! LU decomposition 1. order, 4 bcs
                       + 3 *4            & ! LU decomposition 2. order, 4 bcs
                       + 3 *(1+inb_scal)   ! LU decomposition 2. order w/ diffusivities, 1 bcs
! In Direct mode, we only need 10 instead of 3*4 because only 1 bcs is considered
     ENDIF
     g(is)%inb_grid = g(is)%inb_grid  &
                    + 1                    ! Density correction in anelastic mode
     IF ( (istagger .EQ. 1) .AND. g(is)%periodic ) THEN 
        g(is)%inb_grid = g(is)%inb_grid  &
                       + 5               & ! LU decomposition interpolation
                       + 5                 ! LU decomposition 1. order interpolatory
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
        CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Imax must be a multiple of 2 for the FFT operations.')
        CALL TLAB_STOP(DNS_ERROR_DIMGRID)
     ENDIF
  ENDIF

! loop counters over the whole domain are integer*4
  IF ( isize_field .GT. HUGE(imax) ) THEN
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Integer model of 4 bytes not big enough.')
     CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF

! -------------------------------------------------------------------
! Test periodicity constrains
! -------------------------------------------------------------------
  IF ( g(1)%periodic .AND. ( .NOT. g(1)%uniform ) ) THEN
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Grid must be uniform in periodic direction X')
     CALL TLAB_STOP(DNS_ERROR_CHECKUNIFX)
  ENDIF

  IF ( g(2)%periodic .AND. ( .NOT. g(2)%uniform ) ) THEN
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Grid must be uniform in periodic direction Y')
     CALL TLAB_STOP(DNS_ERROR_CHECKUNIFY)
  ENDIF

  IF ( g(3)%periodic .AND. ( .NOT. g(3)%uniform ) ) THEN
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Grid must be uniform in periodic direction Z')
     CALL TLAB_STOP(DNS_ERROR_CHECKUNIFZ)
  ENDIF

! -------------------------------------------------------------------
! Helmholtz filter that maintains prognostic bcs; I need inb_{flow,scal}
! -------------------------------------------------------------------
  FilterDomainBcsFlow(:) = FilterDomain(2)%BcsMin
  FilterDomainBcsScal(:) = FilterDomain(2)%BcsMin

  IF ( FilterDomain(1)%type   .EQ. DNS_FILTER_HELMHOLTZ     .AND. &
       FilterDomain(2)%BcsMin .NE. DNS_FILTER_BCS_DIRICHLET .AND. &
       FilterDomain(2)%BcsMin .NE. DNS_FILTER_BCS_SOLID     .AND. &
       FilterDomain(2)%BcsMin .NE. DNS_FILTER_BCS_NEUMANN ) THEN
     CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmin', 'void', sRes)
     IF      ( TRIM(ADJUSTL(sRes)) .eq. 'noslip'   ) THEN; FilterDomainBcsFlow(1:3) = DNS_FILTER_BCS_DIRICHLET
     ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'freeslip' ) THEN; FilterDomainBcsFlow(1:3) = DNS_FILTER_BCS_NEUMANN
     ELSE
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. BoundaryConditions.VelocityJmin.')
        CALL TLAB_STOP(DNS_ERROR_IBC)
     ENDIF
     FilterDomainBcsFlow(2) = DNS_FILTER_BCS_DIRICHLET ! Normal velocity is always Dirichlet
     DO is = 1,inb_scal
        WRITE(lstr,*) is; lstr='Scalar'//TRIM(ADJUSTL(lstr))//'Jmin'
        CALL SCANINICHAR(bakfile, inifile, 'BoundaryConditions', TRIM(ADJUSTL(lstr)), 'void', sRes)
        IF      ( TRIM(ADJUSTL(sRes)) .eq. 'dirichlet' ) THEN; FilterDomainBcsScal(is) = DNS_FILTER_BCS_DIRICHLET
        ELSE IF ( TRIM(ADJUSTL(sRes)) .eq. 'neumann'   ) THEN; FilterDomainBcsScal(is) = DNS_FILTER_BCS_NEUMANN
        ELSE
           CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. BoundaryConditions.'//TRIM(ADJUSTL(lstr)))
           CALL TLAB_STOP(DNS_ERROR_IBC)
        ENDIF
     ENDDO
  ENDIF

! -------------------------------------------------------------------
! Other parameters
! -------------------------------------------------------------------
! By default, transport and radiation are caused by last scalar
! The variable inb_scal_array is only available at the end of this routine
  transport%scalar = inb_scal_array
  radiation%scalar = inb_scal_array

  IF ( imixture .EQ. MIXT_TYPE_AIRWATER .OR. imixture .EQ. MIXT_TYPE_AIRWATER_LINEAR ) THEN
     IF ( radiation%type .NE. EQNS_NONE ) THEN
        radiation%active(inb_scal_array    ) = .TRUE. ! liquid
        radiation%active(inb_scal_array + 1) = .TRUE. ! buoyancy
     ENDIF

  ENDIF

  IF ( imode_sim .EQ. DNS_MODE_TEMPORAL .AND. ( .NOT. g(1)%periodic ) ) THEN
     CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Grid must be uniform and periodic in direction X for temporal simulation')
     CALL TLAB_STOP(DNS_ERROR_CHECKUNIFX)
  ENDIF

  IF ( inb_flow + inb_scal .GT. MAX_VARS ) THEN
     CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error MAX_VARS < inb_flow + inb_scal')
     CALL TLAB_STOP(DNS_ERROR_TOTALVARS)
  ENDIF

  SELECT CASE ( imode_eqns )
  CASE( DNS_EQNS_INCOMPRESSIBLE,DNS_EQNS_ANELASTIC )
    IF ( iviscous /= EQNS_EXPLICIT ) THEN
      CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Main.TermViscous undeveloped.')
      CALL TLAB_STOP(DNS_ERROR_OPTION)
    END IF
    IF ( idiffusion /= EQNS_EXPLICIT ) THEN
      CALL TLAB_WRITE_ASCII(efile,C_FILE_LOC//'. Main.TermDiffusion undeveloped.')
      CALL TLAB_STOP(DNS_ERROR_OPTION)
    END IF
  CASE( DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL )
  END SELECT

  RETURN
END SUBROUTINE IO_READ_GLOBAL