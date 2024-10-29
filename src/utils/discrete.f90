#include "dns_const.h"
#include "dns_error.h"

subroutine DISCRETE_READBLOCK(bakfile, inifile, block, var)
    use TLab_Types, only: discrete_dt, MAX_MODES, MAX_PARS, wp, wi
    use TLab_Constants, only: efile
    use TLab_WorkFlow
    use Profiles
    implicit none

    character(LEN=*), intent(in) :: bakfile, inifile, block
    type(discrete_dt), intent(out) :: var

! -------------------------------------------------------------------
    integer(wi) idummy
    character(LEN=512) sRes

! -------------------------------------------------------------------
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Discrete]')
    call TLab_Write_ASCII(bakfile, '#Type=<Varicose/Sinuous/Gaussian/Step>')
    call TLab_Write_ASCII(bakfile, '#Aplitude=<value>')
    call TLab_Write_ASCII(bakfile, '#ModeX=<value>')
    call TLab_Write_ASCII(bakfile, '#ModeZ=<value>')
    call TLab_Write_ASCII(bakfile, '#PhaseX=<value>')
    call TLab_Write_ASCII(bakfile, '#PhaseZ=<value>')
    call TLab_Write_ASCII(bakfile, '#Broadening=<value>')
!    CALL TLab_Write_ASCII(bakfile, '#Parameters=<values>')

    call ScanFile_Char(bakfile, inifile, block, 'Amplitude', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') &        ! backwards compatilibity
        call ScanFile_Char(bakfile, inifile, block, '2DAmpl', '0.0', sRes)
    var%amplitude(:) = 0.0_wp; var%size = MAX_MODES
    call LIST_REAL(sRes, var%size, var%amplitude)

    call ScanFile_Char(bakfile, inifile, block, 'ModeX', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        do idummy = 1, var%size; var%modex(idummy) = idummy; end do
    else
        idummy = MAX_MODES
        call LIST_INTEGER(sRes, idummy, var%modex)
        if (idummy /= var%size) then
            call TLab_Write_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeX.')
            call TLab_Stop(DNS_ERROR_INFDISCR)
        end if
    end if

    call ScanFile_Char(bakfile, inifile, block, 'ModeZ', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        var%modez = 0
    else
        idummy = MAX_MODES
        call LIST_INTEGER(sRes, idummy, var%modez)
        if (idummy /= var%size) then
            call TLab_Write_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeZ.')
            call TLab_Stop(DNS_ERROR_INFDISCR)
        end if
    end if

    call ScanFile_Char(bakfile, inifile, block, 'PhaseX', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') &        ! backwards compatilibity
        call ScanFile_Char(bakfile, inifile, block, '2DPhi', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        var%phasex = 0.0_wp
    else
        idummy = MAX_MODES
        call LIST_REAL(sRes, idummy, var%phasex)
        if (idummy /= var%size) then
            call TLab_Write_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseX.')
            call TLab_Stop(DNS_ERROR_INFDISCR)
        end if
    end if

    call ScanFile_Char(bakfile, inifile, block, 'PhaseZ', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        var%phasez = 0.0_wp
    else
        idummy = MAX_MODES
        call LIST_REAL(sRes, idummy, var%phasez)
        if (idummy /= var%size) then
            call TLab_Write_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseZ.')
            call TLab_Stop(DNS_ERROR_INFDISCR)
        end if
    end if

    call ScanFile_Char(bakfile, inifile, block, 'Type', 'none', sRes) ! Modulation type
    if (TRIM(ADJUSTL(sRes)) == 'none') then; var%type = PROFILE_NONE
    elseif (TRIM(ADJUSTL(sRes)) == 'varicose') then; var%type = PROFILE_GAUSSIAN_ANTISYM
    elseif (TRIM(ADJUSTL(sRes)) == 'sinuous') then; var%type = PROFILE_GAUSSIAN_SYM
    elseif (TRIM(ADJUSTL(sRes)) == 'gaussian') then; var%type = PROFILE_GAUSSIAN
    elseif (TRIM(ADJUSTL(sRes)) == 'step') then; var%type = PROFILE_TANH_COS
    else
        call TLab_Write_ASCII(efile, 'FLOW_READ_GLOBAL. Error in Discrete.Type.')
        call TLab_Stop(DNS_ERROR_INFDISCR)
    end if

    var%parameters(:) = 0.0_wp
    call ScanFile_Char(bakfile, inifile, 'Discrete', 'Parameters', '-1.0,-1.0', sRes)
    idummy = MAX_PARS
    call LIST_REAL(sRes, idummy, var%parameters)

    return
end subroutine DISCRETE_READBLOCK
