#include "dns_error.h"
#include "dns_const.h"

subroutine NavierStokes_Initialize_Parameters(inifile)
    use TLab_Constants, only: wp, wi, lfile, efile, wfile, MAX_PROF, MAX_VARS
    use TLAB_VARS, only: imode_eqns, iadvection, iviscous, idiffusion, itransport
    use TLAB_VARS, only: inb_flow, inb_flow_array, inb_scal, inb_scal_array
    use TLAB_VARS, only: inb_wrk1d, inb_wrk2d
    use TLAB_VARS, only: qbg, sbg, pbg, rbg, tbg, hbg
    use TLAB_VARS, only: buoyancy, coriolis, subsidence
    use TLAB_VARS, only: visc, prandtl, schmidt, mach, damkohler, froude, rossby, stokes, settling
    use TLAB_VARS, only: imode_sim, stagger_on
    use TLAB_VARS, only: g
    use TLAB_VARS, only: FilterDomain, FilterDomainBcsFlow, FilterDomainBcsScal
    use Thermodynamics, only: gama0
    use TLab_Spatial
    use TLab_WorkFlow
    use Profiles, only: Profiles_ReadBlock, PROFILE_EKMAN_U, PROFILE_EKMAN_U_P, PROFILE_EKMAN_V
    implicit none

    character(len=*), intent(in) :: inifile

! -------------------------------------------------------------------
    character(len=32) bakfile, block
    character(len=512) sRes
    character*64 lstr
    integer(wi) is, idummy
    real(wp) dummy, reynolds

! ###################################################################
    bakfile = trim(adjustl(inifile))//'.bak'

! ###################################################################
    block = 'Main'

    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    call TLab_Write_ASCII(bakfile, '#Equations=<total/internal/incompressible/anelastic>')
    call TLab_Write_ASCII(bakfile, '#TermAdvection=<divergence/skewsymmetric>')
    call TLab_Write_ASCII(bakfile, '#TermViscous=<divergence/explicit>')
    call TLab_Write_ASCII(bakfile, '#TermDiffusion=<divergence/explicit>')
    call TLab_Write_ASCII(bakfile, '#TermTransport=<constant/powerlaw/sutherland>')
    !
    call TLab_Write_ASCII(bakfile, '#TermBodyForce=<none/Explicit/Homogeneous/Linear/Bilinear/Quadratic>')
    call TLab_Write_ASCII(bakfile, '#TermCoriolis=<none/explicit/normalized>')
    call TLab_Write_ASCII(bakfile, '#TermSubsidence=<none/ConstantDivergenceLocal/ConstantDivergenceGlobal>')

! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, inifile, 'Main', 'Equations', 'internal', sRes)
    if (trim(adjustl(sRes)) == 'total') then; imode_eqns = DNS_EQNS_TOTAL
    elseif (trim(adjustl(sRes)) == 'internal') then; imode_eqns = DNS_EQNS_INTERNAL
    elseif (trim(adjustl(sRes)) == 'incompressible') then; imode_eqns = DNS_EQNS_INCOMPRESSIBLE
    elseif (trim(adjustl(sRes)) == 'anelastic') then; imode_eqns = DNS_EQNS_ANELASTIC
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong entry Main.Equations option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'TermAdvection', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; iadvection = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'divergence') then; iadvection = EQNS_DIVERGENCE
    else if (trim(adjustl(sRes)) == 'skewsymmetric') then; iadvection = EQNS_SKEWSYMMETRIC
    else if (trim(adjustl(sRes)) == 'convective') then; iadvection = EQNS_CONVECTIVE
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermAdvection option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'TermViscous', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; iviscous = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'divergence') then; iviscous = EQNS_DIVERGENCE
    else if (trim(adjustl(sRes)) == 'explicit') then; iviscous = EQNS_EXPLICIT
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermViscous option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'TermDiffusion', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; idiffusion = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'divergence') then; idiffusion = EQNS_DIVERGENCE
    else if (trim(adjustl(sRes)) == 'explicit') then; idiffusion = EQNS_EXPLICIT
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermDiffusion option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'TermTransport', 'constant', sRes)
    if (trim(adjustl(sRes)) == 'sutherland') then; itransport = EQNS_TRANS_SUTHERLAND; 
    elseif (trim(adjustl(sRes)) == 'powerlaw') then; itransport = EQNS_TRANS_POWERLAW; 
    else; itransport = EQNS_NONE; end if

    ! consistency check
    select case (imode_eqns)
    case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
        if (itransport == EQNS_TRANS_POWERLAW) then
            call TLab_Write_ASCII(efile, 'RHS_SCAL_GLOBAL_2. Only constant viscosity.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        if (imode_eqns == DNS_EQNS_TOTAL) then
            call TLab_Write_ASCII(efile, 'RHS_SCAL_GLOBAL_2. No total energy formulation.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        if (iviscous /= EQNS_EXPLICIT) then
            call TLab_Write_ASCII(efile, __FILE__//'. Main.TermViscous undeveloped.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if
        if (idiffusion /= EQNS_EXPLICIT) then
            call TLab_Write_ASCII(efile, __FILE__//'. Main.TermDiffusion undeveloped.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

    end select

! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, inifile, 'Main', 'TermCoriolis', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; coriolis%type = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'explicit') then; coriolis%type = EQNS_EXPLICIT
    else if (trim(adjustl(sRes)) == 'normalized') then; coriolis%type = EQNS_COR_NORMALIZED
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermCoriolis option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'TermSubsidence', 'None', sRes)
    if (trim(adjustl(sRes)) == 'none') then; subsidence%type = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'constantdivergencelocal') then; subsidence%type = EQNS_SUB_CONSTANT_LOCAL
    else if (trim(adjustl(sRes)) == 'constantdivergenceglobal') then; subsidence%type = EQNS_SUB_CONSTANT_GLOBAL
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermSubsidence option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

! ###################################################################
    block = 'Parameters'

    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    call TLab_Write_ASCII(bakfile, '#Reynolds=<value>')
    call TLab_Write_ASCII(bakfile, '#Schmidt=<value>')
    call TLab_Write_ASCII(bakfile, '#Froude=<value>')
    call TLab_Write_ASCII(bakfile, '#Rossby=<value>')
    call TLab_Write_ASCII(bakfile, '#Damkohler=<value>')
    call TLab_Write_ASCII(bakfile, '#Stokes=<value>')
    call TLab_Write_ASCII(bakfile, '#Settling=<value>')
    call TLab_Write_ASCII(bakfile, '#Mach=<value>')
    call TLab_Write_ASCII(bakfile, '#Gama=<value>')
    call TLab_Write_ASCII(bakfile, '#Prandtl=<value>')

    ! Molecular transport
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Reynolds', '-1.0', reynolds)
    if (reynolds <= 0.0) then
        call ScanFile_Real(bakfile, inifile, 'Parameters', 'Viscosity', '-1.0', dummy)
        if (dummy <= 0.0) then
            call TLab_Write_ASCII(efile, __FILE__//'. Molecular transport coefficients need to be positive.')
            call TLab_Stop(DNS_ERROR_OPTION)
        else
            reynolds = 1.0_wp/dummy
        end if
    end if

    call ScanFile_Char(bakfile, inifile, 'Parameters', 'Schmidt', '1.0', sRes)
    schmidt(:) = 0.0_wp; inb_scal = MAX_VARS
    call LIST_REAL(sRes, inb_scal, schmidt)

    ! Gravity
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Froude', '-1.0', froude)
    if (froude <= 0.0) then
        call ScanFile_Real(bakfile, inifile, 'Parameters', 'Gravity', '1.0', dummy)   ! default value
        froude = 1.0_wp/dummy
    end if

    ! Coriolis
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Rossby', '-1.0', rossby)
    if (rossby <= 0.0) then
        call ScanFile_Real(bakfile, inifile, 'Parameters', 'Coriolis', '1.0', dummy)   ! default value
        rossby = 1.0_wp/dummy
    end if

    ! Chemistry
    lstr = '0.0'
    do is = 2, inb_scal
        lstr = trim(adjustl(lstr))//',0.0'
    end do
    call ScanFile_Char(bakfile, inifile, 'Parameters', 'Damkohler', lstr, sRes)
    damkohler(:) = 0.0_wp; idummy = MAX_VARS
    call LIST_REAL(sRes, idummy, damkohler)
    if (inb_scal /= idummy) then ! Consistency check
        call TLab_Write_ASCII(efile, __FILE__//'. Schmidt and Damkholer sizes do not match.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    ! Compressible flows
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Gama', '1.4', gama0)
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Prandtl', '1.0', prandtl)
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Mach', '1.0', mach)

    ! Particle-laden flows
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Stokes', '0.0', stokes)
    call ScanFile_Real(bakfile, inifile, 'Parameters', 'Settling', '0.0', settling)

    ! consistency check
    if (iviscous == EQNS_NONE) then
        visc = 0.0_wp
    else
        visc = 1.0_wp/reynolds
    end if

    select case (imode_eqns)
    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        prandtl = schmidt(1)

    end select

! ###################################################################
! Buoyancy
! ###################################################################
    ! I wonder if this should be part of [BodyForce] instead of [Main]
    call ScanFile_Char(bakfile, inifile, 'Main', 'TermBodyForce', 'void', sRes)
    if (trim(adjustl(sRes)) == 'none') then; buoyancy%type = EQNS_NONE
    else if (trim(adjustl(sRes)) == 'explicit') then; buoyancy%type = EQNS_EXPLICIT
    else if (trim(adjustl(sRes)) == 'homogeneous') then; buoyancy%type = EQNS_BOD_HOMOGENEOUS
    else if (trim(adjustl(sRes)) == 'linear') then; buoyancy%type = EQNS_BOD_LINEAR
    else if (trim(adjustl(sRes)) == 'bilinear') then; buoyancy%type = EQNS_BOD_BILINEAR
    else if (trim(adjustl(sRes)) == 'quadratic') then; buoyancy%type = EQNS_BOD_QUADRATIC
    else if (trim(adjustl(sRes)) == 'normalizedmean') then; buoyancy%type = EQNS_BOD_NORMALIZEDMEAN
    else if (trim(adjustl(sRes)) == 'subtractmean') then; buoyancy%type = EQNS_BOD_SUBTRACTMEAN
    else
        call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermBodyForce option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    if (any([EQNS_BOD_LINEAR, EQNS_BOD_BILINEAR, EQNS_BOD_QUADRATIC] == buoyancy%type) .and. inb_scal == 0) then
        call TLab_Write_ASCII(wfile, __FILE__//'. Zero scalars; setting TermBodyForce equal to none.')
        buoyancy%type = EQNS_NONE
    end if

    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[BodyForce]')
    call TLab_Write_ASCII(bakfile, '#Vector=<Gx,Gy,Gz>')
    call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

    buoyancy%vector = 0.0_wp; buoyancy%active = .false.
    if (buoyancy%type /= EQNS_NONE) then
        call ScanFile_Char(bakfile, inifile, 'BodyForce', 'Vector', '0.0,-1.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, buoyancy%vector)

        if (abs(buoyancy%vector(1)) > 0.0_wp) then; buoyancy%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Body force along Ox.'); end if
        if (abs(buoyancy%vector(2)) > 0.0_wp) then; buoyancy%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Body force along Oy.'); end if
        if (abs(buoyancy%vector(3)) > 0.0_wp) then; buoyancy%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Body force along Oz.'); end if

        if (froude > 0.0_wp) then
            buoyancy%vector(:) = buoyancy%vector(:)/froude ! adding the froude number into the vector g
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Froude number must be nonzero if buoyancy is retained.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        buoyancy%parameters(:) = 0.0_wp
        call ScanFile_Char(bakfile, inifile, 'BodyForce', 'Parameters', '0.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, buoyancy%parameters)
        buoyancy%scalar(1) = idummy

    end if

! ###################################################################
    block = 'Rotation'

    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    call TLab_Write_ASCII(bakfile, '#Vector=<Fx,Fy,Fz>')
    call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

    coriolis%vector(:) = 0.0_wp; coriolis%active = .false.
    if (coriolis%type /= EQNS_NONE) then
        call ScanFile_Char(bakfile, inifile, 'Rotation', 'Vector', '0.0,1.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, coriolis%vector)

        if (abs(coriolis%vector(1)) > 0.0_wp) then; coriolis%active(2) = .true.; coriolis%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Ox.'); end if
        if (abs(coriolis%vector(2)) > 0.0_wp) then; coriolis%active(3) = .true.; coriolis%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Oy.'); end if
        if (abs(coriolis%vector(3)) > 0.0_wp) then; coriolis%active(1) = .true.; coriolis%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Oz.'); end if

        if (rossby > 0.0_wp) then
            coriolis%vector(:) = coriolis%vector(:)/rossby ! adding the rossby number into the vector
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Rossby number must be nonzero if coriolis is retained.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        coriolis%parameters(:) = 0.0_wp
        call ScanFile_Char(bakfile, inifile, 'Rotation', 'Parameters', '0.0,1.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, coriolis%parameters)

        if (coriolis%parameters(2) == 0.0_wp) then
            call TLab_Write_ASCII(lfile, __FILE__//'. Default normalized geostrophic velocity set to one.')
            coriolis%parameters(2) = 1.0_wp
        end if

    end if

! Consistency check
    if (coriolis%type == EQNS_COR_NORMALIZED) then
        if (coriolis%active(2)) then
            call TLab_Write_ASCII(efile, __FILE__//'. TermCoriolis option only allows for angular velocity along Oy.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if
    end if

! ###################################################################
    block = 'Subsidence'

    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

    subsidence%active = .false.
    if (subsidence%type /= EQNS_NONE) then
        subsidence%active = .true.

        subsidence%parameters(:) = 0.0_wp
        call ScanFile_Char(bakfile, inifile, 'Subsidence', 'Parameters', '0.0', sRes)
        idummy = MAX_PROF
        call LIST_REAL(sRes, idummy, subsidence%parameters)

    end if

! This subsidence type is implemented in opr_burgers_y only
! to speed up calculation
    if (subsidence%type == EQNS_SUB_CONSTANT_LOCAL) subsidence%active = .false.

! ###################################################################
    block = 'Scalar'

    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    do is = 1, MAX_VARS
        write (lstr, *) is
        call Profiles_ReadBlock(bakfile, inifile, 'Scalar', 'Scalar'//trim(adjustl(lstr)), sbg(is))
    end do

! ###################################################################
    block = 'Flow'

    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityX', qbg(1))
    call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityY', qbg(2))
    call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityZ', qbg(3))

    ! backwards compatilibity; originally, all velocity data was contained in block 'Velocity' except for the mean value
    call ScanFile_Char(bakfile, inifile, 'Flow', 'ProfileVelocity', 'void', sRes)
    if (trim(adjustl(sRes)) /= 'void') then
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Velocity', qbg(1))
        call TLab_Write_ASCII(wfile, 'Update tag Flow.Velocity to Flow.VelocityX.')
    end if

    ! Consistency check
    if (any([PROFILE_EKMAN_U, PROFILE_EKMAN_U_P] == qbg(1)%type)) then
        qbg(3) = qbg(1)
        qbg(3)%type = PROFILE_EKMAN_V
    end if

    call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Pressure', pbg)
    call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Density', rbg)
    call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Temperature', tbg)
    call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Enthalpy', hbg)

! ! consistency check; two and only two are givem TO BE CHECKED BECAUSE PROFILE_NONE is used as constant profile
    ! if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
    !     idummy=0
    !     if (pbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (rbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (tbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (hbg%type == PROFILE_NONE) idummy=idummy+1
    !     if (idummy /= 2) then
    !         call TLab_Write_ASCII(efile, __FILE__//'. Specify only 2 thermodynamic profiles.')
    !         call TLab_Stop(DNS_ERROR_OPTION)
    !     end if
    ! end if

! -------------------------------------------------------------------
! Spatial case
! Thickness evolutions delta_i/diam_i=a*(x/diam_i+b)
! -------------------------------------------------------------------
    if (imode_sim == DNS_MODE_SPATIAL) then
        call TLab_Write_ASCII(bakfile, '#ThickAVelocity=<value>')
        call TLab_Write_ASCII(bakfile, '#ThickBVelocity=<value>')
        call TLab_Write_ASCII(bakfile, '#FluxVelocity=<value>')
        call TLab_Write_ASCII(bakfile, '#ThickADensity=<value>')
        call TLab_Write_ASCII(bakfile, '#ThickBDensity=<value>')
        call TLab_Write_ASCII(bakfile, '#FluxDensity=<value>')
        call TLab_Write_ASCII(bakfile, '#ThickATemperature=<value>')
        call TLab_Write_ASCII(bakfile, '#ThickBTemperature=<value>')
        call TLab_Write_ASCII(bakfile, '#FluxTemperature=<value>')

! Bradbury profile is the default (x0=a*b)
        call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickAVelocity', '0.1235', qbg(1)%parameters(2))
        call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickBVelocity', '-0.873', qbg(1)%parameters(3))
        call ScanFile_Real(bakfile, inifile, 'Flow', 'FluxVelocity', '0.96', qbg(1)%parameters(4))

! Ramaprian is the default (x0=a*b)
        call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickADensity', '0.14', rbg%parameters(2))
        call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickBDensity', '2.0', rbg%parameters(3))
        call ScanFile_Real(bakfile, inifile, 'Flow', 'FluxDensity', '0.94', rbg%parameters(4))

        call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickATemperature', '0.14', tbg%parameters(2))
        call ScanFile_Real(bakfile, inifile, 'Flow', 'ThickBTemperature', '2.0', tbg%parameters(3))
        call ScanFile_Real(bakfile, inifile, 'Flow', 'FluxTemperature', '0.94', tbg%parameters(4))

        ! Scalars
        call TLab_Write_ASCII(bakfile, '#ThickA=<value>')
        call TLab_Write_ASCII(bakfile, '#ThickB=<value>')
        call TLab_Write_ASCII(bakfile, '#Flux=<value>')

        do is = 1, MAX_VARS
            write (lstr, *) is; lstr = 'ThickA'//trim(adjustl(lstr))
            call ScanFile_Real(bakfile, inifile, 'Scalar', trim(adjustl(lstr)), '0.14', sbg(is)%parameters(2))
            write (lstr, *) is; lstr = 'ThickB'//trim(adjustl(lstr))
            call ScanFile_Real(bakfile, inifile, 'Scalar', trim(adjustl(lstr)), '2.0', sbg(is)%parameters(3))
            write (lstr, *) is; lstr = 'Flux'//trim(adjustl(lstr))
            call ScanFile_Real(bakfile, inifile, 'Scalar', trim(adjustl(lstr)), '0.94', sbg(is)%parameters(4))
        end do

    end if

! ###################################################################
    if (FilterDomain(1)%type == DNS_FILTER_HELMHOLTZ .and. &
        all([DNS_FILTER_BCS_DIRICHLET, DNS_FILTER_BCS_SOLID, DNS_FILTER_BCS_NEUMANN] /= FilterDomain(2)%BcsMin)) then
        call ScanFile_Char(bakfile, inifile, 'BoundaryConditions', 'VelocityJmin', 'void', sRes)
        if (trim(adjustl(sRes)) == 'noslip') then; FilterDomainBcsFlow(1:3) = DNS_FILTER_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; FilterDomainBcsFlow(1:3) = DNS_FILTER_BCS_NEUMANN
        else
            call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.VelocityJmin.')
            call TLab_Stop(DNS_ERROR_IBC)
        end if
        FilterDomainBcsFlow(2) = DNS_FILTER_BCS_DIRICHLET ! Normal velocity is always Dirichlet

        do is = 1, inb_scal
            write (lstr, *) is; lstr = 'Scalar'//trim(adjustl(lstr))//'Jmin'
            call ScanFile_Char(bakfile, inifile, 'BoundaryConditions', trim(adjustl(lstr)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'dirichlet') then; FilterDomainBcsScal(is) = DNS_FILTER_BCS_DIRICHLET
            else if (trim(adjustl(sRes)) == 'neumann') then; FilterDomainBcsScal(is) = DNS_FILTER_BCS_NEUMANN
            else
                call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr)))
                call TLab_Stop(DNS_ERROR_IBC)
            end if
        end do

    end if

! ###################################################################
! Initialization of array sizes
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')

! prognostic and diagnostic variables
    select case (imode_eqns)
    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        inb_flow = 3                            ! space for u, v, w
        inb_flow_array = inb_flow

    case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
        inb_flow = 5                            ! space for u, v, w, e, rho
        inb_flow_array = inb_flow + 2           ! space for p, T
        if (any([EQNS_TRANS_SUTHERLAND, EQNS_TRANS_POWERLAW] == itransport)) inb_flow_array = inb_flow_array + 1    ! space for viscosity

    end select

    inb_scal_array = inb_scal ! Default is that array contains only the prognostic variables;
    !                           can be changed in Thermodynamics_Initialize_Parameters(ifile)

! scratch arrays
    inb_wrk1d = 18

    inb_wrk2d = 2
    if (imode_sim == DNS_MODE_SPATIAL) inb_wrk2d = max(11, inb_wrk2d)

! grid array
    do is = 1, 3
        g(is)%inb_grid = 1                          ! Nodes
        g(is)%inb_grid = g(is)%inb_grid &
                         + 2 &                      ! Jacobians of first- and second-order derivatives
                         + 2                        ! 1/dx and 1/dx**2 used in time-step stability constraint

        g(is)%inb_grid = g(is)%inb_grid &
                         + 5 &                      ! max # of diagonals in LHS for 1. order derivative
                         + 7 &                      ! max # of diagonals in RHS for 1. order derivative
                         + 5 &                      ! max # of diagonals in LHS for 2. order derivative
                         + 7 + 5                    ! max # of diagonals in RHS for 2. order + diagonals for Jacobian case
        g(is)%inb_grid = g(is)%inb_grid &
                         + 5*2 &                    ! max # of diagonals in LHS for 1. integral, 2 bcs
                         + 7*2                      ! max # of diagonals in RHS for 1. integral, 2 bcs
        if (g(is)%periodic) then
            g(is)%inb_grid = g(is)%inb_grid &
                             + 5 + 2 &                      ! LU decomposition 1. order
                             + 5 + 2 &                      ! LU decomposition 2. order
                             + (5 + 2)*(1 + inb_scal) &     ! LU decomposition 2. order with diffusivities
                             + 2                            ! modified wavenumbers
        else
            g(is)%inb_grid = g(is)%inb_grid &
                             + 5*4 &                ! LU decomposition 1. order, 4 bcs
                             + 5 &                  ! LU decomposition 2. order, 1bcs
                             + 5*(1 + inb_scal)     ! LU decomposition 2. order w/ diffusivities, 1 bcs
        end if
        g(is)%inb_grid = g(is)%inb_grid &
                         + 1                        ! Density correction in anelastic mode
        if ((stagger_on) .and. g(is)%periodic) then
            g(is)%inb_grid = g(is)%inb_grid &
                             + 5 &                  ! LU decomposition interpolation
                             + 5                    ! LU decomposition 1. order interpolatory
        end if

    end do

    return
end subroutine NavierStokes_Initialize_Parameters
