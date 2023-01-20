#include "dns_const.h"
#include "dns_error.h"

subroutine DISCRETE_READBLOCK(bakfile, inifile, block, var)
    use TLAB_TYPES, only: discrete_dt, MAX_MODES, MAX_PARS, wp, wi
    use TLAB_CONSTANTS, only: efile
    use TLAB_PROCS
    implicit none

    character(LEN=*), intent(in) :: bakfile, inifile, block
    type(discrete_dt), intent(out) :: var

! -------------------------------------------------------------------
    integer(wi) idummy
    character(LEN=512) sRes

! -------------------------------------------------------------------
    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#[Discrete]')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<Varicose/Sinuous/Gaussian/Step>')
    call TLAB_WRITE_ASCII(bakfile, '#Aplitude=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#ModeX=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#ModeZ=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#PhaseX=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#PhaseZ=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Broadening=<value>')
!    CALL TLAB_WRITE_ASCII(bakfile, '#Parameters=<values>')

    call SCANINICHAR(bakfile, inifile, block, 'Amplitude', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') &        ! backwards compatilibity
        call SCANINICHAR(bakfile, inifile, block, '2DAmpl', '0.0', sRes)
    var%amplitude(:) = 0.0_wp; var%size = MAX_MODES
    call LIST_REAL(sRes, var%size, var%amplitude)

    call SCANINICHAR(bakfile, inifile, block, 'ModeX', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        do idummy = 1, var%size; var%modex(idummy) = idummy; end do
    else
        idummy = MAX_MODES
        call LIST_INTEGER(sRes, idummy, var%modex)
        if (idummy /= var%size) then
            call TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeX.')
            call TLAB_STOP(DNS_ERROR_INFDISCR)
        end if
    end if

    call SCANINICHAR(bakfile, inifile, block, 'ModeZ', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        var%modez = 0
    else
        idummy = MAX_MODES
        call LIST_INTEGER(sRes, idummy, var%modez)
        if (idummy /= var%size) then
            call TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.ModeZ.')
            call TLAB_STOP(DNS_ERROR_INFDISCR)
        end if
    end if

    call SCANINICHAR(bakfile, inifile, block, 'PhaseX', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') &        ! backwards compatilibity
        call SCANINICHAR(bakfile, inifile, block, '2DPhi', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        var%phasex = 0
    else
        idummy = MAX_MODES
        call LIST_REAL(sRes, idummy, var%phasex)
        if (idummy /= var%size) then
            call TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseX.')
            call TLAB_STOP(DNS_ERROR_INFDISCR)
        end if
    end if

    call SCANINICHAR(bakfile, inifile, block, 'PhaseZ', 'void', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'void') then     ! Default
        var%phasez = 0
    else
        idummy = MAX_MODES
        call LIST_REAL(sRes, idummy, var%phasez)
        if (idummy /= var%size) then
            call TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Inconsistent Discrete.PhaseZ.')
            call TLAB_STOP(DNS_ERROR_INFDISCR)
        end if
    end if

    call SCANINICHAR(bakfile, inifile, block, 'Type', 'none', sRes) ! Modulation type
    if (TRIM(ADJUSTL(sRes)) == 'none') then; var%type = PROFILE_NONE
    elseif (TRIM(ADJUSTL(sRes)) == 'varicose') then; var%type = PROFILE_GAUSSIAN_ANTISYM
    elseif (TRIM(ADJUSTL(sRes)) == 'sinuous') then; var%type = PROFILE_GAUSSIAN_SYM
    elseif (TRIM(ADJUSTL(sRes)) == 'gaussian') then; var%type = PROFILE_GAUSSIAN
    elseif (TRIM(ADJUSTL(sRes)) == 'step') then; var%type = PROFILE_TANH_COS
    else
        call TLAB_WRITE_ASCII(efile, 'FLOW_READ_GLOBAL. Error in Discrete.Type.')
        call TLAB_STOP(DNS_ERROR_INFDISCR)
    end if

    var%parameters(:) = 0.0_wp
    call SCANINICHAR(bakfile, inifile, 'Discrete', 'Parameters', '-1.0,-1.0', sRes)
    idummy = MAX_PARS
    call LIST_REAL(sRes, idummy, var%parameters)

    return
end subroutine DISCRETE_READBLOCK
