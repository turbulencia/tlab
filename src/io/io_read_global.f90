#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "IO_READ_GLOBAL"

!########################################################################
!# Reading general data from file dns.ini, setting up general parameters
!# and doing cross-check of these general data.
!########################################################################
subroutine IO_READ_GLOBAL(inifile)

    use TLAB_CONSTANTS, only: wp, wi, lfile, efile, wfile, MajorVersion, MinorVersion
    use TLAB_VARS
    use TLAB_PROCS
    use THERMO_VARS
    use PROFILES
#ifdef USE_MPI
    use TLAB_MPI_VARS
#endif

    implicit none

    character(len=*), intent(in) :: inifile

! -------------------------------------------------------------------
    character*512 sRes
    character*64 lstr
    character*32 bakfile
    integer(wi) is, inb_scal_local1, inb_scal_local2, idummy
    real(wp) dummy

! ###################################################################
    bakfile = trim(adjustl(inifile))//'.bak'

    call TLAB_WRITE_ASCII(lfile, 'Reading global input data.')

! ###################################################################
! Version Checking
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#[Version]')
    call TLAB_WRITE_ASCII(bakfile, '#Major=<mayor version number>')
    call TLAB_WRITE_ASCII(bakfile, '#Minor=<minor version number>')

    call SCANINIINT(bakfile, inifile, 'Version', 'Major', '0', idummy)
    if (MajorVersion /= idummy) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Major version error.')
        call TLAB_STOP(DNS_ERROR_VERSION)
    end if
    call SCANINIINT(bakfile, inifile, 'Version', 'Minor', '0', idummy)
    if (MinorVersion /= idummy) then
        write (sRes, '(I5)') MinorVersion
        call TLAB_WRITE_ASCII(wfile, 'DNS_REAL_GLOBAL. Minor version warning. Expected : '//sRes)
    end if

! ###################################################################
! Global information
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Main]')
    call TLAB_WRITE_ASCII(bakfile, '#FileFormat=<mpiio/RawSplit/NetCDF>')
    call TLAB_WRITE_ASCII(bakfile, '#VerbosityLevel=<0/1/2>')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<temporal/spatial>')
    call TLAB_WRITE_ASCII(bakfile, '#CalculateFlow=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#CalculateScalar=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#Equations=<total/internal/incompressible/anelastic>')
    call TLAB_WRITE_ASCII(bakfile, '#TermAdvection=<divergence/skewsymmetric>')
    call TLAB_WRITE_ASCII(bakfile, '#TermViscous=<divergence/explicit>')
    call TLAB_WRITE_ASCII(bakfile, '#TermDiffusion=<divergence/explicit>')
    call TLAB_WRITE_ASCII(bakfile, '#TermBodyForce=<none/Explicit/Homogeneous/Linear/Bilinear/Quadratic>')
    call TLAB_WRITE_ASCII(bakfile, '#TermCoriolis=<none/explicit/normalized>')
    call TLAB_WRITE_ASCII(bakfile, '#TermRadiation=<none/Bulk1dGlobal/Bulk1dLocal>')
    call TLAB_WRITE_ASCII(bakfile, '#TermSubsidence=<none/ConstantDivergenceLocal/ConstantDivergenceGlobal>')
    call TLAB_WRITE_ASCII(bakfile, '#TermTransport=<constant/powerlaw/sutherland/Airwater/AirwaterSimplified>')
    call TLAB_WRITE_ASCII(bakfile, '#TermChemistry=<none/quadratic/layeredrelaxation/ozone>')
    call TLAB_WRITE_ASCII(bakfile, '#TermRandom=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#SpaceOrder=<CompactJacobian4/CompactJacobian6/CompactJacpenta6/CompactDirect6>')
    call TLAB_WRITE_ASCII(bakfile, '#ComModeITranspose=<none,asynchronous,sendrecv>')
    call TLAB_WRITE_ASCII(bakfile, '#ComModeKTranspose=<none,asynchronous,sendrecv>')

    call SCANINICHAR(bakfile, inifile, 'Main', 'FileFormat', 'MpiIO', sRes)
    if (trim(adjustl(sRes)) == 'mpiio') then; imode_files = IO_MPIIO
    elseif (trim(adjustl(sRes)) == 'netcdf') then; imode_files = IO_NETCDF
    elseif (trim(adjustl(sRes)) == 'none') then; imode_files = IO_NOFILE
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong Main.FileFormat.')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'FileType', 'Double', sRes)
    if (trim(adjustl(sRes)) == 'double') then; imode_precision_files = IO_TYPE_DOUBLE
    elseif (trim(adjustl(sRes)) == 'single') then; imode_precision_files = IO_TYPE_SINGLE
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong Main.FileType.')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    end if

    call SCANINIINT(bakfile, inifile, 'Main', 'VerbosityLevel', '1', imode_verbosity)

    call SCANINICHAR(bakfile, inifile, 'Main', 'Type', 'temporal', sRes)
    if (trim(adjustl(sRes)) == 'temporal') then; imode_sim = DNS_MODE_TEMPORAL
    elseif (trim(adjustl(sRes)) == 'spatial') then; imode_sim = DNS_MODE_SPATIAL
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Entry Main.Type must be temporal or spatial')
        call TLAB_STOP(DNS_ERROR_SIMTYPE)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'CalculateFlow', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; flow_on = .true.
    elseif (trim(adjustl(sRes)) == 'no') then; flow_on = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Entry Main.CalculateFlow must be yes or no')
        call TLAB_STOP(DNS_ERROR_CALCFLOW)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'CalculateScalar', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; scal_on = .true.
    elseif (trim(adjustl(sRes)) == 'no') then; scal_on = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Entry Main.CalculateScalar must be yes or no')
        call TLAB_STOP(DNS_ERROR_CALCSCALAR)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'Equations', 'internal', sRes)
    if (trim(adjustl(sRes)) == 'total') then; imode_eqns = DNS_EQNS_TOTAL
    elseif (trim(adjustl(sRes)) == 'internal') then; imode_eqns = DNS_EQNS_INTERNAL
    elseif (trim(adjustl(sRes)) == 'incompressible') then; imode_eqns = DNS_EQNS_INCOMPRESSIBLE
    elseif (trim(adjustl(sRes)) == 'anelastic') then; imode_eqns = DNS_EQNS_ANELASTIC
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong entry Main.Equations option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    if (imode_sim == DNS_MODE_TEMPORAL) fourier_on = .true.

    call SCANINICHAR(bakfile, inifile, 'Main', 'Mixture', 'None', sRes)
    if (trim(adjustl(sRes)) == 'none') then; imixture = MIXT_TYPE_NONE
    else if (trim(adjustl(sRes)) == 'air') then; imixture = MIXT_TYPE_AIR
    else if (trim(adjustl(sRes)) == 'airvapor') then; imixture = MIXT_TYPE_AIRVAPOR
    else if (trim(adjustl(sRes)) == 'airwater') then; imixture = MIXT_TYPE_AIRWATER
    else if (trim(adjustl(sRes)) == 'airwaterlinear') then; imixture = MIXT_TYPE_AIRWATER_LINEAR
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong entry Main.Mixture model.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'Main', 'TermAdvection', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; iadvection = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'divergence') then; iadvection = EQNS_DIVERGENCE
    else if (trim(adjustl(sRes)) == 'skewsymmetric') then; iadvection = EQNS_SKEWSYMMETRIC
    else if (trim(adjustl(sRes)) == 'convective') then; iadvection = EQNS_CONVECTIVE
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermAdvection option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'TermViscous', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; iviscous = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'divergence') then; iviscous = EQNS_DIVERGENCE
    else if (trim(adjustl(sRes)) == 'explicit') then; iviscous = EQNS_EXPLICIT
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermViscous option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'TermDiffusion', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; idiffusion = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'divergence') then; idiffusion = EQNS_DIVERGENCE
    else if (trim(adjustl(sRes)) == 'explicit') then; idiffusion = EQNS_EXPLICIT
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermDiffusion option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'TermBodyForce', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; buoyancy%type = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'explicit') then; buoyancy%type = EQNS_EXPLICIT
    else if (trim(adjustl(sRes)) == 'homogeneous') then; buoyancy%type = EQNS_BOD_HOMOGENEOUS
    else if (trim(adjustl(sRes)) == 'linear') then; buoyancy%type = EQNS_BOD_LINEAR
    else if (trim(adjustl(sRes)) == 'bilinear') then; buoyancy%type = EQNS_BOD_BILINEAR
    else if (trim(adjustl(sRes)) == 'quadratic') then; buoyancy%type = EQNS_BOD_QUADRATIC
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermBodyForce option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'TermCoriolis', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; coriolis%type = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'explicit') then; coriolis%type = EQNS_EXPLICIT
    else if (trim(adjustl(sRes)) == 'normalized') then; coriolis%type = EQNS_COR_NORMALIZED
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermCoriolis option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'TermRadiation', 'None', sRes)
    if (trim(adjustl(sRes)) == 'none') then; radiation%type = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'bulk1dglobal') then; radiation%type = EQNS_RAD_BULK1D_GLOBAL
    else if (trim(adjustl(sRes)) == 'bulk1dlocal') then; radiation%type = EQNS_RAD_BULK1D_LOCAL
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermRadiation option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'TermSubsidence', 'None', sRes)
    if (trim(adjustl(sRes)) == 'none') then; subsidence%type = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'constantdivergencelocal') then; subsidence%type = EQNS_SUB_CONSTANT_LOCAL
    else if (trim(adjustl(sRes)) == 'constantdivergenceglobal') then; subsidence%type = EQNS_SUB_CONSTANT_GLOBAL
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong TermSubsidence option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'Main', 'TermTransport', 'constant', sRes)
    if (trim(adjustl(sRes)) == 'sutherland') then; transport%type = EQNS_TRANS_SUTHERLAND; 
    elseif (trim(adjustl(sRes)) == 'powerlaw') then; transport%type = EQNS_TRANS_POWERLAW; 
    elseif (trim(adjustl(sRes)) == 'airwater') then; transport%type = EQNS_TRANS_AIRWATER; 
    elseif (trim(adjustl(sRes)) == 'airwatersimplified') then; transport%type = EQNS_TRANS_AIRWATERSIMPLIFIED; 
    else; transport%type = EQNS_NONE; end if

    itransport = transport%type

! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'Main', 'TermChemistry', 'none', sRes)
    if (trim(adjustl(sRes)) == 'quadratic') then; chemistry%type = EQNS_CHEM_QUADRATIC; 
    elseif (trim(adjustl(sRes)) == 'quadratic3') then; chemistry%type = EQNS_CHEM_QUADRATIC3; 
    elseif (trim(adjustl(sRes)) == 'layeredrelaxation') then; chemistry%type = EQNS_CHEM_LAYEREDRELAXATION; 
    elseif (trim(adjustl(sRes)) == 'ozone') then; chemistry%type = EQNS_CHEM_OZONE; 
    else; chemistry%type = EQNS_NONE; end if

! -------------------------------------------------------------------
    call SCANINIREAL(bakfile, inifile, 'Main', 'TermRandom', '0.0', dummy)
    if (abs(dummy) > 0.0) then
        random%type = EQNS_RAND_MULTIPLY
        random%parameters = dummy
        random%active(1:3) = .true.
    else
        random%type = EQNS_NONE
        random%active(1:3) = .false.
    end if

! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'Main', 'SpaceOrder', 'void', sRes)
    if (trim(adjustl(sRes)) == 'compactjacobian4') then; imode_fdm = FDM_COM4_JACOBIAN; 
    elseif (trim(adjustl(sRes)) == 'compactjacobian6') then; imode_fdm = FDM_COM6_JACOBIAN; 
    elseif (trim(adjustl(sRes)) == 'compactjacpenta6') then; imode_fdm = FDM_COM6_JACPENTA; 
    elseif (trim(adjustl(sRes)) == 'compactjacobian8') then; imode_fdm = FDM_COM8_JACOBIAN; 
    elseif (trim(adjustl(sRes)) == 'compactdirect6') then; imode_fdm = FDM_COM6_DIRECT; 
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Wrong SpaceOrder option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    g(1:3)%mode_fdm = imode_fdm

! -------------------------------------------------------------------
#ifdef USE_MPI
    call SCANINICHAR(bakfile, inifile, 'Main', 'ComModeITranspose', 'asynchronous', sRes)
    if (trim(adjustl(sRes)) == 'none') then; ims_trp_mode_i = TLAB_MPI_TRP_NONE
    elseif (trim(adjustl(sRes)) == 'asynchronous') then; ims_trp_mode_i = TLAB_MPI_TRP_ASYNCHRONOUS
    elseif (trim(adjustl(sRes)) == 'sendrecv') then; ims_trp_mode_i = TLAB_MPI_TRP_SENDRECV
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong ComModeITranspose option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'Main', 'ComModeKTranspose', 'asynchronous', sRes)
    if (trim(adjustl(sRes)) == 'none') then; ims_trp_mode_k = TLAB_MPI_TRP_NONE
    elseif (trim(adjustl(sRes)) == 'asynchronous') then; ims_trp_mode_k = TLAB_MPI_TRP_ASYNCHRONOUS
    elseif (trim(adjustl(sRes)) == 'sendrecv') then; ims_trp_mode_k = TLAB_MPI_TRP_SENDRECV
    else
        call TLAB_WRITE_ASCII(efile, 'DNS_READ_LOCAL. Wrong ComModeKTranspose option.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if
#endif

! -------------------------------------------------------------------
    select case (imode_eqns)
    case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
        if (itransport == EQNS_TRANS_POWERLAW) then
            call TLAB_WRITE_ASCII(efile, 'RHS_SCAL_GLOBAL_2. Only constant viscosity.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        if (imode_eqns == DNS_EQNS_TOTAL) then
            call TLAB_WRITE_ASCII(efile, 'RHS_SCAL_GLOBAL_2. No total energy formulation.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

    end select

! ###################################################################
! Pressure staggering
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Staggering]')
    call TLAB_WRITE_ASCII(bakfile, '#StaggerHorizontalPressure=<yes/no>')

    call SCANINICHAR(bakfile, inifile, 'Staggering', 'StaggerHorizontalPressure', 'no', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; stagger_on = .true.; call TLAB_WRITE_ASCII(lfile, 'Horizontal staggering of the pressure along Ox and Oz.')
    elseif (trim(adjustl(sRes)) == 'no') then; stagger_on = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Entry Main. StaggerHorizontalPressure must be yes or no')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! Consistency check
    if (stagger_on) then
        if (.not. ((imode_eqns == DNS_EQNS_INCOMPRESSIBLE) .or. (imode_eqns == DNS_EQNS_ANELASTIC))) then
 call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Horizontal pressure staggering only implemented for anelastic or incompressible mode.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. ((iadvection == EQNS_CONVECTIVE) .or. (iadvection == EQNS_SKEWSYMMETRIC))) then
          call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Horizontal pressure staggering not implemented for current advection scheme.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. (imode_fdm == FDM_COM6_JACOBIAN)) then
call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Horizontal pressure staggering only implemented for compact jacobian 6th-order scheme.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
    end if

! ###################################################################
! Dimensionles parameters
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Parameters]')
    call TLAB_WRITE_ASCII(bakfile, '#Reynolds=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Prandtl=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Froude=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Rossby=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Mach=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Gama=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Schmidt=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Damkohler=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Stokes=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Settling=<value>')

    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Reynolds', '100', reynolds)
    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Gama', '1.4', gama0)
    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Prandtl', '1.0', prandtl)
    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Mach', '1.0', mach)
    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Froude', '1.0', froude)
    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Rossby', '1.0', rossby)
    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Stokes', '0.0', stokes)
    call SCANINIREAL(bakfile, inifile, 'Parameters', 'Settling', '0.0', settling)

    call SCANINICHAR(bakfile, inifile, 'Parameters', 'Schmidt', '1.0', sRes)
    schmidt(:) = 0.0_wp; inb_scal_local1 = MAX_NSP
    call LIST_REAL(sRes, inb_scal_local1, schmidt)

    lstr = '0.0'; do is = 2, inb_scal_local1; lstr = trim(adjustl(lstr))//',0.0'; end do
    call SCANINICHAR(bakfile, inifile, 'Parameters', 'Damkohler', lstr, sRes)
    damkohler(:) = 0.0_wp; inb_scal_local2 = MAX_NSP
    call LIST_REAL(sRes, inb_scal_local2, damkohler)
    if (inb_scal_local1 /= inb_scal_local2) then ! Consistency check
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Schmidt and Damkholer sizes do not match.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! ###################################################################
! Buoyancy
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[BodyForce]')
    call TLAB_WRITE_ASCII(bakfile, '#Vector=<Gx,Gy,Gz>')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

    buoyancy%vector = 0.0_wp; buoyancy%active = .false.
    if (buoyancy%type /= EQNS_NONE) then
        call SCANINICHAR(bakfile, inifile, 'BodyForce', 'Vector', '0.0,-1.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, buoyancy%vector)

        if (abs(buoyancy%vector(1)) > 0.0_wp) then; buoyancy%active(1) = .true.; call TLAB_WRITE_ASCII(lfile, 'Body force along Ox.'); end if
        if (abs(buoyancy%vector(2)) > 0.0_wp) then; buoyancy%active(2) = .true.; call TLAB_WRITE_ASCII(lfile, 'Body force along Oy.'); end if
        if (abs(buoyancy%vector(3)) > 0.0_wp) then; buoyancy%active(3) = .true.; call TLAB_WRITE_ASCII(lfile, 'Body force along Oz.'); end if

        if (froude > 0.0_wp) then
            buoyancy%vector(:) = buoyancy%vector(:)/froude ! adding the froude number into de vector g
        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Froude number must be nonzero if buoyancy is retained.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        buoyancy%parameters(:) = 0.0_wp
        call SCANINICHAR(bakfile, inifile, 'BodyForce', 'Parameters', '0.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, buoyancy%parameters)

    end if

! ###################################################################
! Rotation
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Rotation]')
    call TLAB_WRITE_ASCII(bakfile, '#Vector=<Fx,Fy,Fz>')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

    coriolis%vector(:) = 0.0_wp; coriolis%active = .false.
    if (coriolis%type /= EQNS_NONE) then
        call SCANINICHAR(bakfile, inifile, 'Rotation', 'Vector', '0.0,1.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, coriolis%vector)

        if (abs(coriolis%vector(1)) > 0.0_wp) then; coriolis%active(2) = .true.; coriolis%active(3) = .true.; call TLAB_WRITE_ASCII(lfile, 'Angular velocity along Ox.'); end if
        if (abs(coriolis%vector(2)) > 0.0_wp) then; coriolis%active(3) = .true.; coriolis%active(1) = .true.; call TLAB_WRITE_ASCII(lfile, 'Angular velocity along Oy.'); end if
        if (abs(coriolis%vector(3)) > 0.0_wp) then; coriolis%active(1) = .true.; coriolis%active(2) = .true.; call TLAB_WRITE_ASCII(lfile, 'Angular velocity along Oz.'); end if

        if (rossby > 0.0_wp) then
            coriolis%vector(:) = coriolis%vector(:)/rossby ! adding the rossby number into the vector
        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Rossby number must be nonzero if coriolis is retained.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        coriolis%parameters(:) = 0.0_wp
        call SCANINICHAR(bakfile, inifile, 'Rotation', 'Parameters', '0.0,1.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, coriolis%parameters)

        if (coriolis%parameters(2) == 0.0_wp) then
            call TLAB_WRITE_ASCII(lfile, C_FILE_LOC//'. Default normalized geostrophic velocity set to one.')
            coriolis%parameters(2) = 1.0_wp
        end if

    end if

! Consistency check
    if (coriolis%type == EQNS_COR_NORMALIZED) then
        if (coriolis%active(2)) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. TermCoriolis option only allows for angular velocity along Oy.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if
    end if

! ###################################################################
! Radiation
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Radiation]')
    call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

    radiation%active = .false.
    if (radiation%type /= EQNS_NONE) then
        call SCANINIINT(bakfile, inifile, 'Radiation', 'Scalar', '1', idummy)
        radiation%active(idummy) = .true.

        radiation%parameters(:) = 0.0_wp
        call SCANINICHAR(bakfile, inifile, 'Radiation', 'Parameters', '1.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, radiation%parameters)

    end if

! ###################################################################
! Subsidence
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Subsidence]')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

    subsidence%active = .false.
    if (subsidence%type /= EQNS_NONE) then
        subsidence%active = .true.

        subsidence%parameters(:) = 0.0_wp
        call SCANINICHAR(bakfile, inifile, 'Subsidence', 'Parameters', '0.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, subsidence%parameters)

    end if

! This subsidence type is implemented in opr_burgers_y only
! to speed up calculation
    if (subsidence%type == EQNS_SUB_CONSTANT_LOCAL) subsidence%active = .false.

! ###################################################################
! Transport
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Transport]')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Exponent=<value>')

    transport%active = .false.
    if (transport%type /= EQNS_NONE) then
        transport%parameters(:) = 1.0_wp ! default values
        call SCANINICHAR(bakfile, inifile, 'Transport', 'Parameters', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, transport%parameters)
        end if

        if (settling > 0.0_wp) then
            transport%parameters = transport%parameters*settling ! adding the settling number in the parameter definitions
        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Settling number must be nonzero if transport is retained.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
            transport%active = .true. ! All scalars are affected

            call SCANINIREAL(bakfile, inifile, 'Transport', 'Exponent', '0.0', transport%auxiliar(1))
        end if

    end if

! ###################################################################
! Chemistry
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Chemistry]')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

    if (chemistry%type /= EQNS_NONE) then
        chemistry%parameters(:) = 0.0_wp
        call SCANINICHAR(bakfile, inifile, 'Chemistry', 'Parameters', '1.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, chemistry%parameters)

    end if

! Activating terms
    chemistry%active = .false.
    do is = 1, inb_scal_local1
        if (abs(damkohler(is)) > 0.0_wp) chemistry%active(is) = .true.
    end do

! ###################################################################
! Thermodynamics
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Thermodynamics]')
    call TLAB_WRITE_ASCII(bakfile, '#Nondimensional=<yes,no>')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

    call SCANINICHAR(bakfile, inifile, 'Thermodynamics', 'Nondimensional', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; nondimensional = .true.
    else if (trim(adjustl(sRes)) == 'no') then; nondimensional = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Thermodynamics.Nondimensional')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    if (imixture /= EQNS_NONE) then
        thermo_param(:) = 0.0_wp
        call SCANINICHAR(bakfile, inifile, 'Thermodynamics', 'Parameters', '1.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, thermo_param)

    end if

    call SCANINIREAL(bakfile, inifile, 'Thermodynamics', 'ScaleHeight', '0.0', scaleheight)

    if (imixture == MIXT_TYPE_AIRWATER) then
        call SCANINIREAL(bakfile, inifile, 'Thermodynamics', 'SmoothFactor', '0.1', dsmooth)
    end if

! ###################################################################
! Grid Parameters
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Grid]')
    call TLAB_WRITE_ASCII(bakfile, '#Imax=<imax>')
    call TLAB_WRITE_ASCII(bakfile, '#Imax(*)=<imax_proc>')
    call TLAB_WRITE_ASCII(bakfile, '#Jmax=<jmax>')
    call TLAB_WRITE_ASCII(bakfile, '#Kmax=<kmax>')
    call TLAB_WRITE_ASCII(bakfile, '#Kmax(*)=<kmax_proc>')
    call TLAB_WRITE_ASCII(bakfile, '#XUniform=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#YUniform=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#ZUniform=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#XPeriodic=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#YPeriodic=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#ZPeriodic=<yes/no>')

    call SCANINIINT(bakfile, inifile, 'Grid', 'Imax', '0', g(1)%size)
    call SCANINIINT(bakfile, inifile, 'Grid', 'Jmax', '0', g(2)%size)
    call SCANINIINT(bakfile, inifile, 'Grid', 'Kmax', '0', g(3)%size)

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
    if (ims_npro > 1) then
        call SCANINIINT(bakfile, inifile, 'Grid', 'Kmax(*)', '-1', kmax)
        if (kmax > 0 .and. mod(g(3)%size, kmax) == 0) then
            ims_npro_k = g(3)%size/kmax
        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Input kmax incorrect')
            call TLAB_STOP(DNS_ERROR_KMAXTOTAL)
        end if

        call SCANINIINT(bakfile, inifile, 'Grid', 'Imax(*)', '-1', imax)
        if (imax > 0 .and. mod(g(1)%size, imax) == 0) then
            ims_npro_i = g(1)%size/imax
        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Input imax incorrect')
            call TLAB_STOP(DNS_ERROR_KMAXTOTAL)
        end if

        if (ims_npro_i*ims_npro_k == ims_npro) then ! check
            write (lstr, *) ims_npro_i; write (sRes, *) ims_npro_k
            lstr = trim(adjustl(lstr))//'x'//trim(adjustl(sRes))
            call TLAB_WRITE_ASCII(lfile, 'Initializing domain partition '//trim(adjustl(lstr)))
        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Inconsistency in total number of PEs')
            call TLAB_STOP(DNS_ERROR_KMAXTOTAL)
        end if

    end if

#endif

! -------------------------------------------------------------------
! Uniform
! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'Grid', 'XUniform', 'void', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; g(1)%uniform = .true.
    else if (trim(adjustl(sRes)) == 'no') then; g(1)%uniform = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Uniform X grid')
        call TLAB_STOP(DNS_ERROR_UNIFORMX)
    end if

    call SCANINICHAR(bakfile, inifile, 'Grid', 'YUniform', 'void', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; g(2)%uniform = .true.
    else if (trim(adjustl(sRes)) == 'no') then; g(2)%uniform = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Uniform Y grid')
        call TLAB_STOP(DNS_ERROR_UNIFORMY)
    end if

    call SCANINICHAR(bakfile, inifile, 'Grid', 'ZUniform', 'void', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; g(3)%uniform = .true.
    else if (trim(adjustl(sRes)) == 'no') then; g(3)%uniform = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Uniform Z grid')
        call TLAB_STOP(DNS_ERROR_UNIFORMZ)
    end if

! -------------------------------------------------------------------
! Periodic
! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, inifile, 'Grid', 'XPeriodic', 'void', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; g(1)%periodic = .true.
    else if (trim(adjustl(sRes)) == 'no') then; g(1)%periodic = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Periodic X grid')
        call TLAB_STOP(DNS_ERROR_IBC)
    end if

    call SCANINICHAR(bakfile, inifile, 'Grid', 'YPeriodic', 'void', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; g(2)%periodic = .true.
    else if (trim(adjustl(sRes)) == 'no') then; g(2)%periodic = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Periodic Y grid')
        call TLAB_STOP(DNS_ERROR_JBC)
    end if

    call SCANINICHAR(bakfile, inifile, 'Grid', 'ZPeriodic', 'void', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; g(3)%periodic = .true.
    else if (trim(adjustl(sRes)) == 'no') then; g(3)%periodic = .false.
    else
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error in Periodic Z grid')
        call TLAB_STOP(DNS_ERROR_KBC)
    end if

! ###################################################################
! Dealising (a filter type)
! ###################################################################
    call FILTER_READBLOCK(bakfile, inifile, 'Dealiasing', Dealiasing)

! ###################################################################
! Pressure Filter (a filter type)
! ###################################################################
    call FILTER_READBLOCK(bakfile, inifile, 'PressureFilter', PressureFilter)

    ! Consistency check
    if (PressureFilter(1)%type /= DNS_FILTER_NONE) call TLAB_WRITE_ASCII(lfile, 'Pressure and dpdy filter along Ox.')
    if (PressureFilter(2)%type /= DNS_FILTER_NONE) call TLAB_WRITE_ASCII(lfile, 'Pressure and dpdy filter along Oy.')
    if (PressureFilter(3)%type /= DNS_FILTER_NONE) call TLAB_WRITE_ASCII(lfile, 'Pressure and dpdy filter along Oz.')

    if (any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        if (.not. ((imode_eqns == DNS_EQNS_INCOMPRESSIBLE) .or. (imode_eqns == DNS_EQNS_ANELASTIC))) then
       call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Pressure and dpdy filter only implemented for anelastic or incompressible mode.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. (iadvection == EQNS_CONVECTIVE)) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Pressure and dpdy filter not implemented for current advection scheme.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
    end if

! ###################################################################
! Domain Filter (a filter type)
! ###################################################################
    call FILTER_READBLOCK(bakfile, inifile, 'Filter', FilterDomain)

! To eventually allow for control field by field
    FilterDomainActive(:) = .true.

! ###################################################################
! Statistics Control
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Statistics]')
    call TLAB_WRITE_ASCII(bakfile, '#IAvera=<plane1,plane2,...>')

    nstatavg = MAX_STATS_SPATIAL
    call SCANINICHAR(bakfile, inifile, 'Statistics', 'IAvera', '1', sRes)
    call LIST_INTEGER(sRes, nstatavg, statavg)

! ###################################################################
! Flow physical properties of the system
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Flow]')

    call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'VelocityX', qbg(1))
    call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'VelocityY', qbg(2))
    call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'VelocityZ', qbg(3))

    ! backwards compatilibity; originally, all velocity data was contained in block 'Velocity' except for the mean value
    call SCANINICHAR(bakfile, inifile, 'Flow', 'ProfileVelocity', 'void', sRes)
    if (trim(adjustl(sRes)) /= 'void') then
        call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'Velocity', qbg(1))
        call TLAB_WRITE_ASCII(wfile, 'Update tag Flow.Velocity to Flow.VelocityX.')
    end if

    ! Consistency check
    if (qbg(1)%type == PROFILE_EKMAN_U .or. qbg(1)%type == PROFILE_EKMAN_U_P) then
        qbg(3)%type = PROFILE_EKMAN_V
        qbg(3)%ymean = qbg(1)%ymean
        qbg(3)%ymean_rel = qbg(1)%ymean_rel
        qbg(3)%thick = qbg(1)%thick
        qbg(3)%delta = qbg(1)%delta
    end if

    call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'Pressure', pbg)
    call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'Density', rbg)
    call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'Temperature', tbg)
    call PROFILES_READBLOCK(bakfile, inifile, 'Flow', 'Enthalpy', hbg)

! ! consistency check; two and only two are givem TO BE CHECKED BECAUSE PROFILE_NONE is used as constant profile
    ! if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
    !     idummy=0
    !     if (pbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (rbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (tbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (hbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (idummy /= 2) then
    !         call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Specify only 2 thermodynamic profiles.')
    !         call TLAB_STOP(DNS_ERROR_OPTION)
    !     end if
    ! end if

! Scalars
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Scalar]')

    do is = 1, MAX_NSP
        write (lstr, *) is
        call PROFILES_READBLOCK(bakfile, inifile, 'Scalar', 'Scalar'//trim(adjustl(lstr)), sbg(is))
    end do

! -------------------------------------------------------------------
! Spatial case
! Thickness evolutions delta_i/diam_i=a*(x/diam_i+b)
! -------------------------------------------------------------------
    if (imode_sim == DNS_MODE_SPATIAL) then
        call TLAB_WRITE_ASCII(bakfile, '#ThickAVelocity=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#ThickBVelocity=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#FluxVelocity=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#ThickADensity=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#ThickBDensity=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#FluxDensity=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#ThickATemperature=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#ThickBTemperature=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#FluxTemperature=<value>')

! Bradbury profile is the default (x0=a*b)
        call SCANINIREAL(bakfile, inifile, 'Flow', 'ThickAVelocity', '0.1235', qbg(1)%parameters(2))
        call SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBVelocity', '-0.873', qbg(1)%parameters(3))
        call SCANINIREAL(bakfile, inifile, 'Flow', 'FluxVelocity', '0.96', qbg(1)%parameters(4))

! Ramaprian is the default (x0=a*b)
        call SCANINIREAL(bakfile, inifile, 'Flow', 'ThickADensity', '0.14', rbg%parameters(2))
        call SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBDensity', '2.0', rbg%parameters(3))
        call SCANINIREAL(bakfile, inifile, 'Flow', 'FluxDensity', '0.94', rbg%parameters(4))

        call SCANINIREAL(bakfile, inifile, 'Flow', 'ThickATemperature', '0.14', tbg%parameters(2))
        call SCANINIREAL(bakfile, inifile, 'Flow', 'ThickBTemperature', '2.0', tbg%parameters(3))
        call SCANINIREAL(bakfile, inifile, 'Flow', 'FluxTemperature', '0.94', tbg%parameters(4))

        ! Scalars
        call TLAB_WRITE_ASCII(bakfile, '#ThickA=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#ThickB=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Flux=<value>')

        do is = 1, MAX_NSP
            write (lstr, *) is; lstr = 'ThickA'//trim(adjustl(lstr))
            call SCANINIREAL(bakfile, inifile, 'Scalar', trim(adjustl(lstr)), '0.14', sbg(is)%parameters(2))
            write (lstr, *) is; lstr = 'ThickB'//trim(adjustl(lstr))
            call SCANINIREAL(bakfile, inifile, 'Scalar', trim(adjustl(lstr)), '2.0', sbg(is)%parameters(3))
            write (lstr, *) is; lstr = 'Flux'//trim(adjustl(lstr))
            call SCANINIREAL(bakfile, inifile, 'Scalar', trim(adjustl(lstr)), '0.94', sbg(is)%parameters(4))
        end do

    end if

! ###################################################################
! Final consistency check and initialization
! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')

    if (iviscous == EQNS_NONE) then; visc = 0.0_wp
    else; visc = 1.0_wp/reynolds; end if

! -------------------------------------------------------------------
! Initializing thermodynamic data of the mixture
! -------------------------------------------------------------------
    inb_scal = inb_scal_local1 ! Default is general N scalars; gama0 has been already read above.
    inb_scal_array = inb_scal
    NSP = inb_scal

    if (imixture /= MIXT_TYPE_NONE) then ! particular mixture (requires implementation)
        call THERMO_INITIALIZE                ! gama0 is defined here
        if (inb_scal_local1 /= inb_scal) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Incorrect number of Schmidt numbers.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if
    end if
! Value of R_0/(C_{p,0}W_0) is called GRATIO
    if (gama0 > 0.0_wp) GRATIO = (gama0 - 1.0_wp)/gama0

    if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
        mach = 0.0_wp
        MRATIO = 1.0_wp
        prandtl = schmidt(1)
    else
        MRATIO = gama0*mach*mach
    end if

    if (.not. nondimensional) then
        MRATIO = 1.0_wp
        GRATIO = 1.0_wp
    end if

    if (imode_eqns == DNS_EQNS_ANELASTIC .and. &
        imixture /= MIXT_TYPE_AIR .and. imixture /= MIXT_TYPE_AIRVAPOR .and. imixture /= MIXT_TYPE_AIRWATER) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Incorrect mixture type.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    if (buoyancy%type == EQNS_BOD_LINEAR .or. &
        buoyancy%type == EQNS_BOD_BILINEAR .or. &
        buoyancy%type == EQNS_BOD_QUADRATIC) then
        if (inb_scal == 0) then
            call TLAB_WRITE_ASCII(wfile, C_FILE_LOC//'. Zero scalars; setting TermBodyForce equal to none.')
            buoyancy%type = EQNS_NONE
        end if
    end if

! mean_rho and delta_rho need to be defined, because of old version.
! Note that rho1 and rho2 are the values defined by equation of state,
! being then mean_rho=(rho1+rho2)/2.
    if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
        if (rbg%type == PROFILE_NONE) then
            dummy = tbg%delta/tbg%mean
            rbg%mean = MRATIO*pbg%mean/tbg%mean/(1.0_wp - 0.25_wp*dummy*dummy)
            rbg%delta = -rbg%mean*dummy
            rbg%thick = tbg%thick
            rbg%diam = tbg%diam
        else
            dummy = rbg%delta/rbg%mean
            tbg%mean = MRATIO*pbg%mean/rbg%mean/(1.0_wp - 0.25_wp*dummy*dummy)
            tbg%delta = -tbg%mean*dummy
            tbg%thick = rbg%thick
            tbg%diam = rbg%diam
        end if
    end if

! Consistency check
    if (imixture > 0) then

        if (imixture == MIXT_TYPE_BS &
            .or. imixture == MIXT_TYPE_BSZELDOVICH) then
! These cases force Sc_i=Sc_Z=Pr (Lewis unity)
            schmidt(inb_scal) = prandtl

        else if (imixture == MIXT_TYPE_QUASIBS) then
! These cases force Sc_i=Sc_Z, already read

        else if (imixture == MIXT_TYPE_AIRWATER) then
            if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
                schmidt(2:3) = schmidt(1) ! used in diffusion eqns, though should be fixed
            end if

            if (damkohler(1) == 0.0_wp .and. damkohler(2) == 0.0_wp) then
                damkohler(1:2) = damkohler(3)
            else
                call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. AirWater requires at least first 2 Damkholer numbers zero.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if

        end if

    end if

! -------------------------------------------------------------------
! Arrays size
! -------------------------------------------------------------------
    isize_field = imax*jmax*kmax

    if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
        inb_flow = 3
        inb_flow_array = inb_flow
    else
        inb_flow = 5
        inb_flow_array = inb_flow + 2                                ! space for p, T
        if (transport%type == EQNS_TRANS_SUTHERLAND .or. transport%type == EQNS_TRANS_POWERLAW) inb_flow_array = inb_flow_array + 1 ! space for vis
    end if

    inb_wrk1d = 18
    if (imode_sim == DNS_MODE_SPATIAL) then; inb_wrk2d = 11
    else; inb_wrk2d = 2; end if

    isize_wrk1d = max(g(1)%size, max(g(2)%size, g(3)%size))
    isize_wrk2d = max(imax*jmax, max(imax*kmax, jmax*kmax))

! grid array
    do is = 1, 3
        g(is)%inb_grid = 1                  ! Nodes
        g(is)%inb_grid = g(is)%inb_grid &
                         + 2 &   ! Jacobians of first- and second-order derivatives
                         + 2                  ! 1/dx and 1/dx**2 used in time-step stability constraint

        if (g(is)%periodic) then
            g(is)%inb_grid = g(is)%inb_grid &
                             + 7 & ! LU decomposition 1. order
                             + 5 & ! LU decomposition 2. order
                             + 5*(1 + inb_scal) & ! LU decomposition 2. order with diffusivities
                             + 2                 ! modified wavenumbers
        else
            g(is)%inb_grid = g(is)%inb_grid &
                             + 5*4 & ! LU decomposition 1. order, 4 bcs
                             + 3*4 & ! LU decomposition 2. order, 4 bcs
                             + 3*(1 + inb_scal)   ! LU decomposition 2. order w/ diffusivities, 1 bcs
! In Direct mode, we only need 10 instead of 3*4 because only 1 bcs is considered
        end if
        g(is)%inb_grid = g(is)%inb_grid &
                         + 1                    ! Density correction in anelastic mode
        if ((stagger_on) .and. g(is)%periodic) then
            g(is)%inb_grid = g(is)%inb_grid &
                             + 5 & ! LU decomposition interpolation
                             + 5                 ! LU decomposition 1. order interpolatory
        end if
    end do

! auxiliar array txc
    isize_txc_field = imax*jmax*kmax
    if (fourier_on) then
        isize_txc_dimz = (imax + 2)*(jmax + 2)
        isize_txc_dimx = kmax*(jmax + 2)
        isize_txc_field = isize_txc_dimz*kmax ! space for FFTW lib
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            if (mod(isize_txc_dimz, (2*ims_npro_k)) /= 0) then ! add space for MPI transposition
                isize_txc_dimz = isize_txc_dimz/(2*ims_npro_k)
                isize_txc_dimz = (isize_txc_dimz + 1)*(2*ims_npro_k)
            end if
            isize_txc_field = max(isize_txc_field, isize_txc_dimz*kmax)
        end if
        if (ims_npro_i > 1) then
            if (mod(isize_txc_dimx, (2*ims_npro_i)) /= 0) then ! add space for MPI transposition
                isize_txc_dimx = isize_txc_dimx/(2*ims_npro_i)
                isize_txc_dimx = (isize_txc_dimx + 1)*(2*ims_npro_i)
            end if
            isize_txc_field = max(isize_txc_field, isize_txc_dimx*(imax + 2))
        end if
#endif
        if (mod(imax, 2) /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Imax must be a multiple of 2 for the FFT operations.')
            call TLAB_STOP(DNS_ERROR_DIMGRID)
        end if
    end if

! loop counters over the whole domain are integer*4
    if (isize_field > huge(imax)) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Integer model of 4 bytes not big enough.')
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)
    end if

! -------------------------------------------------------------------
! Test periodicity constrains
! -------------------------------------------------------------------
    if (g(1)%periodic .and. (.not. g(1)%uniform)) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Grid must be uniform in periodic direction X')
        call TLAB_STOP(DNS_ERROR_CHECKUNIFX)
    end if

    if (g(2)%periodic .and. (.not. g(2)%uniform)) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Grid must be uniform in periodic direction Y')
        call TLAB_STOP(DNS_ERROR_CHECKUNIFY)
    end if

    if (g(3)%periodic .and. (.not. g(3)%uniform)) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Grid must be uniform in periodic direction Z')
        call TLAB_STOP(DNS_ERROR_CHECKUNIFZ)
    end if

! -------------------------------------------------------------------
! Helmholtz filter that maintains prognostic bcs; I need inb_{flow,scal}
! -------------------------------------------------------------------
    FilterDomainBcsFlow(:) = FilterDomain(2)%BcsMin
    FilterDomainBcsScal(:) = FilterDomain(2)%BcsMin

    if (FilterDomain(1)%type == DNS_FILTER_HELMHOLTZ .and. &
        FilterDomain(2)%BcsMin /= DNS_FILTER_BCS_DIRICHLET .and. &
        FilterDomain(2)%BcsMin /= DNS_FILTER_BCS_SOLID .and. &
        FilterDomain(2)%BcsMin /= DNS_FILTER_BCS_NEUMANN) then
        call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmin', 'void', sRes)
        if (trim(adjustl(sRes)) == 'noslip') then; FilterDomainBcsFlow(1:3) = DNS_FILTER_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; FilterDomainBcsFlow(1:3) = DNS_FILTER_BCS_NEUMANN
        else
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. BoundaryConditions.VelocityJmin.')
            call TLAB_STOP(DNS_ERROR_IBC)
        end if
        FilterDomainBcsFlow(2) = DNS_FILTER_BCS_DIRICHLET ! Normal velocity is always Dirichlet
        do is = 1, inb_scal
            write (lstr, *) is; lstr = 'Scalar'//trim(adjustl(lstr))//'Jmin'
            call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'dirichlet') then; FilterDomainBcsScal(is) = DNS_FILTER_BCS_DIRICHLET
            else if (trim(adjustl(sRes)) == 'neumann') then; FilterDomainBcsScal(is) = DNS_FILTER_BCS_NEUMANN
            else
                call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. BoundaryConditions.'//trim(adjustl(lstr)))
                call TLAB_STOP(DNS_ERROR_IBC)
            end if
        end do
    end if

! -------------------------------------------------------------------
! Other parameters
! -------------------------------------------------------------------
! By default, transport and radiation are caused by last scalar
! The variable inb_scal_array is only available at the end of this routine
    transport%scalar = inb_scal_array
    radiation%scalar = inb_scal_array

    if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
        if (radiation%type /= EQNS_NONE) then
            radiation%active(inb_scal_array) = .true. ! liquid
            radiation%active(inb_scal_array + 1) = .true. ! buoyancy
        end if

    end if

    if (imode_sim == DNS_MODE_TEMPORAL .and. (.not. g(1)%periodic)) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Grid must be uniform and periodic in direction X for temporal simulation')
        call TLAB_STOP(DNS_ERROR_CHECKUNIFX)
    end if

    if (inb_flow + inb_scal > MAX_VARS) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error MAX_VARS < inb_flow + inb_scal')
        call TLAB_STOP(DNS_ERROR_TOTALVARS)
    end if

    select case (imode_eqns)
    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        if (iviscous /= EQNS_EXPLICIT) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Main.TermViscous undeveloped.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if
        if (idiffusion /= EQNS_EXPLICIT) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Main.TermDiffusion undeveloped.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if
    case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
    end select

    if (imode_fdm == FDM_COM6_JACPENTA) then ! CFL_max depends on max[g%mwn(:,1)]
        call TLAB_WRITE_ASCII(wfile, C_FILE_LOC//'. Main.SpaceOrder.CompactJacpenta6 requires adjusted CFL-number depending on C1N6M_ALPHA, C1N6M_BETA values.')
    end if

    return
end subroutine IO_READ_GLOBAL

!########################################################################
! Should be in module opr_filters but cannot because of the dependencies...
!########################################################################
subroutine FILTER_READBLOCK(bakfile, inifile, tag, variable)
    use TLAB_CONSTANTS, only: efile, wp, MAX_PARS
    use TLAB_TYPES, only: filter_dt
    use TLAB_VARS, only: g
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    implicit none

    character(len=*) bakfile, inifile, tag
    type(filter_dt) variable(3)

! -------------------------------------------------------------------
    character(len=512) sRes
    character(len=10) default
    integer idummy, ig

!########################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(tag))//']')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<none/compact/helmholtz/SpectralBand/SpectralErf/tophat>')
    call TLAB_WRITE_ASCII(bakfile, '#Parameters=<values>')
    call TLAB_WRITE_ASCII(bakfile, '#ActiveX=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#ActiveY=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#ActiveZ=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#BcsJmin=<free,solid,zero>')
    call TLAB_WRITE_ASCII(bakfile, '#BcsJmax=<free,solid,zero>')

    variable(:)%size = g(:)%size
    variable(:)%periodic = g(:)%periodic
    variable(:)%uniform = g(:)%uniform
    variable(:)%inb_filter = 0          ! default array size
    default = 'biased'                  ! default boundary condition

    call SCANINICHAR(bakfile, inifile, trim(adjustl(tag)), 'Type', 'none', sRes)
    if (trim(adjustl(sRes)) == 'none') then; variable(:)%type = DNS_FILTER_NONE
    else if (trim(adjustl(sRes)) == 'compact') then; variable(:)%type = DNS_FILTER_COMPACT
        variable(:)%parameters(1) = 0.49 ! default alpha value
        variable(:)%inb_filter = 10
    else if (trim(adjustl(sRes)) == 'explicit6') then; variable(:)%type = DNS_FILTER_6E
    else if (trim(adjustl(sRes)) == 'explicit4') then; variable(:)%type = DNS_FILTER_4E
        variable(:)%inb_filter = 5
    else if (trim(adjustl(sRes)) == 'adm') then; variable(:)%type = DNS_FILTER_ADM
        variable(:)%inb_filter = 5
    else if (trim(adjustl(sRes)) == 'tophat') then; variable(:)%type = DNS_FILTER_TOPHAT
        variable(:)%parameters(1) = 2    ! default filter size (in grid-step units)
        default = 'free'
    else if (trim(adjustl(sRes)) == 'compactcutoff') then; variable(:)%type = DNS_FILTER_COMPACT_CUTOFF
        variable(:)%inb_filter = 7
        variable(2)%type = DNS_FILTER_COMPACT ! nonuniform version not yet implemented; fall back to compact
        variable(2)%parameters(1) = 0.49
        variable(2)%inb_filter = 10
    else if (trim(adjustl(sRes)) == 'spectralcutoff') then; variable(:)%type = DNS_FILTER_BAND
        ! The frequency interval is (Parameter1, Parameter2)
    else if (trim(adjustl(sRes)) == 'spectralerf') then; variable(:)%type = DNS_FILTER_ERF
        ! Parameter1 is the transition wavenumber in physical units:
        ! >0: high-pass filter
        ! <0; low-pass filter
        ! Parameter2 is the characteristic width--in log units (relative to domain size)'
        variable(:)%parameters(3) = 1.0_wp    ! used to normalise wavenumbers in z-direction
    else if (trim(adjustl(sRes)) == 'helmholtz') then; variable(:)%type = DNS_FILTER_HELMHOLTZ
        variable(:)%parameters(1) = 1.0_wp    ! default filter size
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong '//trim(adjustl(tag))//'Type.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    ! Boundary conditions correction
    do ig = 1, 3
        if (variable(ig)%periodic) then
            variable(ig)%BcsMin = DNS_FILTER_BCS_PERIODIC
            variable(ig)%BcsMax = DNS_FILTER_BCS_PERIODIC
        end if
    end do

    call SCANINICHAR(bakfile, inifile, trim(adjustl(tag)), 'BcsJmin', trim(adjustl(default)), sRes)
    if (trim(adjustl(sRes)) == 'periodic') then; variable(2)%BcsMin = DNS_FILTER_BCS_PERIODIC
    else if (trim(adjustl(sRes)) == 'biased') then; variable(2)%BcsMin = DNS_FILTER_BCS_BIASED
    else if (trim(adjustl(sRes)) == 'free') then; variable(2)%BcsMin = DNS_FILTER_BCS_FREE
    else if (trim(adjustl(sRes)) == 'solid') then; variable(2)%BcsMin = DNS_FILTER_BCS_SOLID
    else if (trim(adjustl(sRes)) == 'dirichlet') then; variable(2)%BcsMin = DNS_FILTER_BCS_DIRICHLET
    else if (trim(adjustl(sRes)) == 'neumann') then; variable(2)%BcsMin = DNS_FILTER_BCS_NEUMANN
    else if (trim(adjustl(sRes)) == 'zero') then; variable(2)%BcsMin = DNS_FILTER_BCS_ZERO
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Filter.BcsJmin.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, trim(adjustl(tag)), 'BcsJmax', trim(adjustl(default)), sRes)
    if (trim(adjustl(sRes)) == 'periodic') then; variable(2)%BcsMax = DNS_FILTER_BCS_PERIODIC
    else if (trim(adjustl(sRes)) == 'biased') then; variable(2)%BcsMax = DNS_FILTER_BCS_BIASED
    else if (trim(adjustl(sRes)) == 'free') then; variable(2)%BcsMax = DNS_FILTER_BCS_FREE
    else if (trim(adjustl(sRes)) == 'solid') then; variable(2)%BcsMax = DNS_FILTER_BCS_SOLID
    else if (trim(adjustl(sRes)) == 'dirichlet') then; variable(2)%BcsMax = DNS_FILTER_BCS_DIRICHLET
    else if (trim(adjustl(sRes)) == 'neumann') then; variable(2)%BcsMax = DNS_FILTER_BCS_NEUMANN
    else if (trim(adjustl(sRes)) == 'zero') then; variable(2)%BcsMax = DNS_FILTER_BCS_ZERO
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Filter.BcsJmax.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, trim(adjustl(tag)), 'Parameters', 'void', sRes)
    if (trim(adjustl(sRes)) /= 'void') then
        idummy = MAX_PARS
        call LIST_REAL(sRes, idummy, variable(1)%parameters(:))
        if (idummy < 3) & ! Fill 3 directions; if global, filled information is unused
            variable(1)%parameters(idummy + 1:3) = variable(1)%parameters(idummy)
        do ig = 2, 3
            variable(ig)%parameters(1) = variable(1)%parameters(ig)
        end do
    end if

    call SCANINIINT(bakfile, inifile, trim(adjustl(tag)), 'Repeat', '1', idummy)
    if (idummy > 0) then
        variable(:)%repeat = idummy
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Entry Filter.Repeat must be positive.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, trim(adjustl(tag)), 'ActiveX', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'no') variable(1)%type = DNS_FILTER_NONE
    call SCANINICHAR(bakfile, inifile, trim(adjustl(tag)), 'ActiveY', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'no') variable(2)%type = DNS_FILTER_NONE
    call SCANINICHAR(bakfile, inifile, trim(adjustl(tag)), 'ActiveZ', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'no') variable(3)%type = DNS_FILTER_NONE

    ! Further control
    do ig = 1, 3
        if (variable(ig)%size == 1) variable(ig)%type = DNS_FILTER_NONE

        if (variable(ig)%type == DNS_FILTER_TOPHAT .and. &
            variable(ig)%parameters(1) == 0) variable(ig)%type = DNS_FILTER_NONE

        if (variable(ig)%type == DNS_FILTER_TOPHAT) then
            if (mod(int(variable(ig)%parameters(1)), 2) /= 0) then
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Tophat filter size must be even.')
                call TLAB_STOP(DNS_ERROR_PARAMETER)
            end if
            variable(ig)%inb_filter = int(variable(ig)%parameters(1)) + 1
        end if

    end do

#ifdef USE_MPI
    variable(1)%mpitype = TLAB_MPI_I_PARTIAL
    variable(3)%mpitype = TLAB_MPI_K_PARTIAL
#endif

    return
end subroutine FILTER_READBLOCK
