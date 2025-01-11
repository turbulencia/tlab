#include "dns_const.h"
#include "dns_error.h"

module Microphysics
    use TLab_Constants, only: wp, wi, pi_wp, efile, MAX_PARS
    use TLab_Types, only: term_dt
    use FDM, only: grid_dt
    use TLAB_VARS, only: imode_eqns, inb_scal_array
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thermodynamics, only: imixture
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
    use OPR_ODES
    implicit none
    private
   
    type(term_dt), public, protected :: sedimentationProps          ! Microphysics parameters
    ! type(term_dt), public, protected :: evaporationProps            ! Microphysics parameters

    public :: Microphysics_Initialize
    public :: Microphysics_Sedimentation

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_SED_AIRWATER = 1
    integer, parameter :: TYPE_SED_AIRWATERSIMPLIFIED = 2

contains
!########################################################################
!########################################################################
    subroutine Microphysics_Initialize(inifile)
        use TLAB_VARS, only: settling
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile
        character(len=512) sRes
        integer(wi) idummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#[Sedimentation]')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')
        call TLab_Write_ASCII(bakfile, '#Exponent=<value>')

        call ScanFile_Char(bakfile, inifile, 'Transport', 'Parameters', 'void', sRes)                 ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) /= 'void') then
            call TLab_Write_ASCII(efile, __FILE__//'. Deprecated block [Transport]. Update to [Sedimentation].')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, 'Sedimentation', 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'TermTransport', 'None', sRes)               ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) == 'none') then; sedimentationProps%type = TYPE_NONE
        elseif (trim(adjustl(sRes)) == 'airwater') then; sedimentationProps%type = TYPE_SED_AIRWATER
        elseif (trim(adjustl(sRes)) == 'airwatersimplified') then; sedimentationProps%type = TYPE_SED_AIRWATERSIMPLIFIED
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Error in Sedimentation.Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        sedimentationProps%active = .false.
        if (sedimentationProps%type /= EQNS_NONE) then
            if (any([MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR] == imixture)) then
                sedimentationProps%active = .true.           ! All scalars are affected
            end if

            sedimentationProps%parameters(:) = 1.0_wp        ! default values
            call ScanFile_Char(bakfile, inifile, 'Sedimentation', 'Parameters', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                idummy = MAX_PARS
                call LIST_REAL(sRes, idummy, sedimentationProps%parameters)
            end if

            if (any([MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR] == imixture)) then
                call ScanFile_Real(bakfile, inifile, 'Sedimentation', 'Exponent', '0.0', sedimentationProps%auxiliar(1))
            end if

        end if

        ! -------------------------------------------------------------------
        ! By default, sedimentation is caused by last scalar
        sedimentationProps%scalar = inb_scal_array

        if (sedimentationProps%type /= EQNS_NONE) then
            if (settling > 0.0_wp) then
                sedimentationProps%parameters = sedimentationProps%parameters*settling ! adding the settling number in the parameter definitions
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Settling number must be nonzero if sedimentation is retained.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        return
    end subroutine Microphysics_Initialize

!########################################################################
!########################################################################
    subroutine Microphysics_Sedimentation(locProps, nx, ny, nz, is, g, s, source, tmp1, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz, is
        type(grid_dt), intent(in) :: g
        real(wp), intent(in) :: s(nx*ny*nz, inb_scal_array)
        real(wp), intent(out) :: source(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)
        real(wp), intent(out), optional :: flux(nx*ny*nz)

        target s, source

        ! -----------------------------------------------------------------------
        real(wp) dummy, exponent
        integer(wi) bcs(2, 2)

        real(wp), pointer :: s_active(:) => null()

        !########################################################################
        bcs = 0

        exponent = locProps%auxiliar(1)
        dummy = 1.0_wp + exponent

        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx, ny, nz, rbackground, s(:, locProps%scalar(is)), source)
            s_active => source
        else
            s_active => s(:, locProps%scalar(is))
        end if

        select case (locProps%type)
        case (TYPE_SED_AIRWATER)
            select case (is)
            case (2, 3)         ! q_t, q_l
                if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                    tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*(s_active**dummy)
                else
                    tmp1 = locProps%parameters(is)*(1.0_wp - s(:, is))*s_active
                end if

            case default        ! energy variables
                call THERMO_ANELASTIC_STATIC_L(nx, ny, nz, s, tmp1)
                if (exponent > 0.0_wp) then
                    tmp1 = locProps%parameters(is)*tmp1*(s_active**dummy)
                else
                    tmp1 = locProps%parameters(is)*tmp1*s_active
                end if

            end select

            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g, tmp1, source)
            if (present(flux)) flux = -tmp1

        case (TYPE_SED_AIRWATERSIMPLIFIED)
            ! if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
            !     tmp1 = locProps%parameters(is)*(s_active**dummy)
            ! else
            !     tmp1 = locProps%parameters(is)*s_active
            ! end if

            ! call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g, tmp1, source)
            ! if (present(flux)) flux = -tmp1

            ! the previous formulation yields oscillations at sharp gradients
            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g, s_active, tmp1)
            if (exponent > 0.0_wp) tmp1 = tmp1*(s_active**exponent)
            source = locProps%parameters(is)*dummy*tmp1

            if (present(flux)) flux = -locProps%parameters(is)*(s_active**dummy)

        end select

! ###################################################################
        nullify (s_active)

    end subroutine Microphysics_Sedimentation

! !########################################################################
! !########################################################################
!     subroutine Microphysics_Evaporation()

!         return
!     end subroutine Microphysics_Evaporation

end module
