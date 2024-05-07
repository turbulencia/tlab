#include "dns_const.h"
#include "dns_error.h"

module Microphysics
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, efile, MAX_PROF
    use TLAB_TYPES, only: term_dt, grid_dt
    use TLAB_VARS, only: imode_eqns, inb_scal_array
    use TLAB_PROCS, only: TLAB_WRITE_ASCII, TLAB_STOP
    use THERMO_VARS, only: imixture
    use OPR_PARTIAL, only: OPR_PARTIAL_Y
    use OPR_ODES
    implicit none

    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_SED_AIRWATER = 1
    integer, parameter :: TYPE_SED_AIRWATERSIMPLIFIED = 2

    public :: Microphysics_Initialize
    public :: Microphysics_Sedimentation

contains
!########################################################################
!########################################################################
    subroutine Microphysics_Initialize(inifile)
        use TLAB_VARS, only: sedimentation, settling
        character(len=*), intent(in) :: inifile

        character(len=32) bakfile
        character(len=512) sRes
        integer(wi) idummy

        ! -------------------------------------------------------------------
        bakfile = trim(adjustl(inifile))//'.bak'

        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#[Sedimentation]')
        call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<value>')
        call TLAB_WRITE_ASCII(bakfile, '#Parameters=<value>')

        call SCANINICHAR(bakfile, inifile, 'Sedimentation', 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') &
            call SCANINICHAR(bakfile, inifile, 'Main', 'TermTransport', 'None', sRes)               ! backwards compatibility, to be removed
        if (trim(adjustl(sRes)) == 'none') then; sedimentation%type = TYPE_NONE
        elseif (trim(adjustl(sRes)) == 'airwater') then; sedimentation%type = TYPE_SED_AIRWATER
        elseif (trim(adjustl(sRes)) == 'airwatersimplified') then; sedimentation%type = TYPE_SED_AIRWATERSIMPLIFIED
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Sedimentation option.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        sedimentation%active = .false.
        if (sedimentation%type /= EQNS_NONE) then
            if (any([MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR] == imixture)) then
                sedimentation%active = .true. ! All scalars are affected
            end if

            sedimentation%parameters(:) = 1.0_wp ! default values
            call SCANINICHAR(bakfile, inifile, 'Sedimentation', 'Parameters', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                idummy = MAX_PROF
                call LIST_REAL(sRes, idummy, sedimentation%parameters)
            end if

            if (any([MIXT_TYPE_AIRWATER, MIXT_TYPE_AIRWATER_LINEAR] == imixture)) then
                call SCANINIREAL(bakfile, inifile, 'Sedimentation', 'Exponent', '0.0', sedimentation%auxiliar(1))
            end if

        end if

        ! -------------------------------------------------------------------
        ! By default, transport and radiation are caused by last scalar
        sedimentation%scalar = inb_scal_array

        if (sedimentation%type /= EQNS_NONE) then
            if (settling > 0.0_wp) then
                sedimentation%parameters = sedimentation%parameters*settling ! adding the settling number in the parameter definitions
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Settling number must be nonzero if sedimentation is retained.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end if

        return
    end subroutine Microphysics_Initialize

!########################################################################
!########################################################################
    subroutine Microphysics_Sedimentation(sedimentation, nx, ny, nz, is, g, s, source, tmp1, flux)
        use THERMO_ANELASTIC
        type(term_dt), intent(in) :: sedimentation
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

        exponent = sedimentation%auxiliar(1)
        dummy = 1.0_wp + exponent

        if (imode_eqns == DNS_EQNS_ANELASTIC) then
            call THERMO_ANELASTIC_WEIGHT_OUTPLACE(nx, ny, nz, rbackground, s(:, sedimentation%scalar(is)), source)
            s_active => source
        else
            s_active => s(:, sedimentation%scalar(is))
        end if

        select case (sedimentation%type)
        case (TYPE_SED_AIRWATER)
            select case (is)
            case (2, 3)         ! q_t, q_l
                if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                    tmp1 = sedimentation%parameters(is)*(1.0_wp - s(:, is))*(s_active**dummy)
                else
                    tmp1 = sedimentation%parameters(is)*(1.0_wp - s(:, is))*s_active
                end if

            case default        ! energy variables
                call THERMO_ANELASTIC_STATIC_L(nx, ny, nz, s, tmp1)
                if (exponent > 0.0_wp) then
                    tmp1 = sedimentation%parameters(is)*tmp1*(s_active**dummy)
                else
                    tmp1 = sedimentation%parameters(is)*tmp1*s_active
                end if

            end select

        case (TYPE_SED_AIRWATERSIMPLIFIED)
            if (exponent > 0.0_wp) then ! to avoid the calculation of a power, if not necessary
                tmp1 = sedimentation%parameters(is)*(s_active**dummy)
            else
                tmp1 = sedimentation%parameters(is)*s_active
            end if

        end select

        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g, tmp1, source)
        if (present(flux)) flux = -tmp1

    end subroutine Microphysics_Sedimentation

! !########################################################################
! !########################################################################
!     subroutine Microphysics_Evaporation()

!         return
!     end subroutine Microphysics_Evaporation

end module
