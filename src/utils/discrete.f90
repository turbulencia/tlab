#include "dns_const.h"
#include "dns_error.h"

module Discrete
    use TLab_Constants, only: wp, wi, efile, MAX_PARS, MAX_MODES
    implicit none
    private

    type, public :: discrete_dt
        sequence
        integer type, size
        integer, dimension(MAX_MODES) :: modex, modez
        real(wp), dimension(MAX_MODES) :: amplitude, phasex, phasez
        real(wp), dimension(MAX_PARS) :: parameters
    end type discrete_dt

    public :: Discrete_ReadBlock

contains
    subroutine Discrete_ReadBlock(bakfile, inifile, block, var)
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        use Profiles, only: PROFILE_NONE, PROFILE_TANH_COS, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_ANTISYM

        character(LEN=*), intent(in) :: bakfile, inifile, block
        type(discrete_dt), intent(out) :: var

! -------------------------------------------------------------------
        integer(wi) idummy
        character(LEN=512) sRes

! -------------------------------------------------------------------
        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<Varicose/Sinuous/Gaussian/Step>')
        call TLab_Write_ASCII(bakfile, '#Aplitude=<value>')
        call TLab_Write_ASCII(bakfile, '#ModeX=<value>')
        call TLab_Write_ASCII(bakfile, '#ModeZ=<value>')
        call TLab_Write_ASCII(bakfile, '#PhaseX=<value>')
        call TLab_Write_ASCII(bakfile, '#PhaseZ=<value>')
        call TLab_Write_ASCII(bakfile, '#Broadening=<value>')
!    CALL TLab_Write_ASCII(bakfile, '#Parameters=<values>')

        call ScanFile_Char(bakfile, inifile, block, 'Amplitude', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &        ! backwards compatilibity
            call ScanFile_Char(bakfile, inifile, block, '2DAmpl', '0.0', sRes)
        var%amplitude(:) = 0.0_wp; var%size = MAX_MODES
        call LIST_REAL(sRes, var%size, var%amplitude)

        call ScanFile_Char(bakfile, inifile, block, 'ModeX', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') then     ! Default
            do idummy = 1, var%size; var%modex(idummy) = idummy; end do
        else
            idummy = MAX_MODES
            call LIST_INTEGER(sRes, idummy, var%modex)
            if (idummy /= var%size) then
                call TLab_Write_ASCII(efile, __FILE__//'. Inconsistent Discrete.ModeX.')
                call TLab_Stop(DNS_ERROR_INFDISCR)
            end if
        end if

        call ScanFile_Char(bakfile, inifile, block, 'ModeZ', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') then     ! Default
            var%modez = 0
        else
            idummy = MAX_MODES
            call LIST_INTEGER(sRes, idummy, var%modez)
            if (idummy /= var%size) then
                call TLab_Write_ASCII(efile, __FILE__//'. Inconsistent Discrete.ModeZ.')
                call TLab_Stop(DNS_ERROR_INFDISCR)
            end if
        end if

        call ScanFile_Char(bakfile, inifile, block, 'PhaseX', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &        ! backwards compatilibity
            call ScanFile_Char(bakfile, inifile, block, '2DPhi', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') then     ! Default
            var%phasex = 0.0_wp
        else
            idummy = MAX_MODES
            call LIST_REAL(sRes, idummy, var%phasex)
            if (idummy /= var%size) then
                call TLab_Write_ASCII(efile, __FILE__//'. Inconsistent Discrete.PhaseX.')
                call TLab_Stop(DNS_ERROR_INFDISCR)
            end if
        end if

        call ScanFile_Char(bakfile, inifile, block, 'PhaseZ', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') then     ! Default
            var%phasez = 0.0_wp
        else
            idummy = MAX_MODES
            call LIST_REAL(sRes, idummy, var%phasez)
            if (idummy /= var%size) then
                call TLab_Write_ASCII(efile, __FILE__//'. Inconsistent Discrete.PhaseZ.')
                call TLab_Stop(DNS_ERROR_INFDISCR)
            end if
        end if

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'none', sRes) ! Modulation type
        if (trim(adjustl(sRes)) == 'none') then; var%type = PROFILE_NONE
        elseif (trim(adjustl(sRes)) == 'varicose') then; var%type = PROFILE_GAUSSIAN_ANTISYM
        elseif (trim(adjustl(sRes)) == 'sinuous') then; var%type = PROFILE_GAUSSIAN_SYM
        elseif (trim(adjustl(sRes)) == 'gaussian') then; var%type = PROFILE_GAUSSIAN
        elseif (trim(adjustl(sRes)) == 'step') then; var%type = PROFILE_TANH_COS
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Error in Discrete.Type.')
            call TLab_Stop(DNS_ERROR_INFDISCR)
        end if

        var%parameters(:) = 0.0_wp
        call ScanFile_Char(bakfile, inifile, 'Discrete', 'Parameters', '-1.0,-1.0', sRes)
        idummy = MAX_PARS
        call LIST_REAL(sRes, idummy, var%parameters)

        return
    end subroutine Discrete_ReadBlock

end module Discrete
