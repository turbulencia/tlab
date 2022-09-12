#include "dns_error.h"
#include "dns_const.h"

subroutine FLOW_READ_LOCAL(inifile)
    use TLAB_CONSTANTS, only: efile, lfile, wfile
    use TLAB_VARS, only: qbg
    use TLAB_PROCS
    use FLOW_LOCAL

    implicit none

    integer(ci) :: bcs_flow_jmin, bcs_flow_jmax ! Boundary conditions

    character*(*) inifile

    ! -------------------------------------------------------------------
    real(cp) dummy(1)
    integer(ci) idummy
    character*512 sRes
    character*32 bakfile

    ! ###################################################################
    bakfile = TRIM(ADJUSTL(inifile))//'.bak'

    call TLAB_WRITE_ASCII(lfile, 'Reading local input data')

    ! ###################################################################
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[IniFields]')
    call TLAB_WRITE_ASCII(bakfile, '#Velocity=<option>')
    call TLAB_WRITE_ASCII(bakfile, '#Temperature=<option>')
    call TLAB_WRITE_ASCII(bakfile, '#ForceDilatation=<yes/no>')
    call TLAB_WRITE_ASCII(bakfile, '#ProfileIniK=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#ThickIniK=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#YCoorIniK=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#NormalizeK=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#NormalizeP=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Mixture=<string>')

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'Velocity', 'None', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; flag_u = 0
    else if (TRIM(ADJUSTL(sRes)) == 'velocitydiscrete') then; flag_u = 1
    else if (TRIM(ADJUSTL(sRes)) == 'velocitybroadband') then; flag_u = 2
    else if (TRIM(ADJUSTL(sRes)) == 'vorticitybroadband') then; flag_u = 3
    else if (TRIM(ADJUSTL(sRes)) == 'potentialbroadband') then; flag_u = 4
    else
        call TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Velocity forcing type unknown')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'ForceDilatation', 'yes', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'no') then; flag_dilatation = 0
    else; flag_dilatation = 1; end if

    Kini = qbg(1) ! default geometry and scaling of perturbation

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'ProfileIniK', 'GaussianSurface', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; Kini%type = PROFILE_NONE
    else if (TRIM(ADJUSTL(sRes)) == 'gaussian') then; Kini%type = PROFILE_GAUSSIAN
    else if (TRIM(ADJUSTL(sRes)) == 'gaussianvaricose') then; Kini%type = PROFILE_GAUSSIAN_ANTISYM
    else if (TRIM(ADJUSTL(sRes)) == 'gaussiansinuous') then; Kini%type = PROFILE_GAUSSIAN_SYM
    else if (TRIM(ADJUSTL(sRes)) == 'gaussiansurface') then; Kini%type = PROFILE_GAUSSIAN_SURFACE
    else if (TRIM(ADJUSTL(sRes)) == 'parabolicsurface') then; Kini%type = PROFILE_PARABOLIC_SURFACE
    else
        call TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Wrong ProfileIni parameter.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIniK', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') & ! backwards compatilibity
        call SCANINICHAR(bakfile, inifile, 'IniFields', 'ThickIni', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) /= 'void') then
        dummy(1) = 1.0_cp; idummy = 1
        call LIST_REAL(sRes, idummy, dummy)
        Kini%thick = dummy(1)
    end if

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIniK', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') & ! backwards compatilibity
        call SCANINICHAR(bakfile, inifile, 'IniFields', 'YCoorIni', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) /= 'void') then
        dummy(1) = 1.0_cp; idummy = 1
        call LIST_REAL(sRes, idummy, dummy)
        Kini%ymean_rel = dummy(1)
    end if

    call SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeK', '-1.0', norm_ini_u)
    call SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeP', '-1.0', norm_ini_p)

    call SCANINICHAR(bakfile, inifile, 'IniFields', 'Temperature', 'None', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; flag_t = 0
    else if (TRIM(ADJUSTL(sRes)) == 'planebroadband') then; flag_t = 4
    else if (TRIM(ADJUSTL(sRes)) == 'planediscrete') then; flag_t = 5
    else
        call TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. Temperature forcing type unknown')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    ! Additional parameters
    call SCANINICHAR(bakfile, inifile, 'IniFields', 'Mixture', 'None', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; flag_mixture = 0
    else if (TRIM(ADJUSTL(sRes)) == 'equilibrium') then; flag_mixture = 1
    else if (TRIM(ADJUSTL(sRes)) == 'loadfields') then; flag_mixture = 2
    end if

    ! Boundary conditions
    flag_wall = 0
    call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmin', 'freeslip', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; bcs_flow_jmin = DNS_BCS_NONE
    else if (TRIM(ADJUSTL(sRes)) == 'noslip') then; bcs_flow_jmin = DNS_BCS_DIRICHLET; flag_wall = flag_wall + 1
    else if (TRIM(ADJUSTL(sRes)) == 'freeslip') then; bcs_flow_jmin = DNS_BCS_NEUMANN
    else
        call TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. BoundaryConditions.VelocityJmin.')
        call TLAB_STOP(DNS_ERROR_IBC)
    end if
    call SCANINICHAR(bakfile, inifile, 'BoundaryConditions', 'VelocityJmax', 'freeslip', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; bcs_flow_jmax = DNS_BCS_NONE
    else if (TRIM(ADJUSTL(sRes)) == 'noslip') then; bcs_flow_jmax = DNS_BCS_DIRICHLET; flag_wall = flag_wall + 2
    else if (TRIM(ADJUSTL(sRes)) == 'freeslip') then; bcs_flow_jmax = DNS_BCS_NEUMANN
    else
        call TLAB_WRITE_ASCII(efile, 'FLOW_READ_LOCAL. BoundaryConditions.VelocityJmax.')
        call TLAB_STOP(DNS_ERROR_IBC)
    end if

    ! ###################################################################
    ! Discrete Forcing
    ! ###################################################################
    call DISCRETE_READBLOCK(bakfile, inifile, 'Discrete', fp)
    ! Modulation type in fp%type

!   specific for this tool
    call SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', fp%parameters(1))
    call SCANINIREAL(bakfile, inifile, 'Discrete', 'ThickStep', '-1.0', fp%parameters(2))

    return
end subroutine FLOW_READ_LOCAL
