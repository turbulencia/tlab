#include "dns_const.h"
#include "dns_error.h"

! Everything has been read from input file.
! Check for cross dependencies and undeveloped options.
subroutine TLab_Consistency_Check()
    use TLab_Constants, only: efile, lfile, MAX_VARS
    use TLAB_VARS, only: imode_eqns, iadvection, subsidence
    use TLAB_VARS, only: inb_flow, inb_scal
    use TLAB_VARS, only: stagger_on, PressureFilter
    use TLAB_VARS, only: g
    use IBM_VARS, only: imode_ibm
    use Radiation
    use Microphysics
    use Chemistry
    use SpecialForcing
    use TLab_WorkFlow
    implicit none

! ###################################################################
    if (stagger_on) then
        if (.not. ((imode_eqns == DNS_EQNS_INCOMPRESSIBLE) .or. (imode_eqns == DNS_EQNS_ANELASTIC))) then
            call TLab_Write_ASCII(efile, __FILE__//'. Horizontal pressure staggering only implemented for anelastic or incompressible mode.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. ((iadvection == EQNS_CONVECTIVE) .or. (iadvection == EQNS_SKEWSYMMETRIC))) then
            call TLab_Write_ASCII(efile, __FILE__//'. Horizontal pressure staggering not implemented for current advection scheme.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (any([g(1)%mode_fdm1, g(2)%mode_fdm1, g(3)%mode_fdm1] /= FDM_COM6_JACOBIAN)) then
            call TLab_Write_ASCII(efile, __FILE__//'. Horizontal pressure staggering only implemented for compact jacobian 6th-order scheme.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
    end if

    ! ###################################################################
    if (imode_ibm == 1) then
        if (.not. (imode_eqns == DNS_EQNS_INCOMPRESSIBLE)) then
            call TLab_Write_ASCII(efile, __FILE__//'. IBM only implemented for incompressible mode.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. ((iadvection == EQNS_CONVECTIVE) .or. (iadvection == EQNS_SKEWSYMMETRIC))) then
            call TLab_Write_ASCII(efile, __FILE__//'. IBM only implemented for convective advection scheme.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        if ((infraredProps%type /= EQNS_NONE) .or. &
            (sedimentationProps%type /= EQNS_NONE) .or. &
            (infraredProps%type /= EQNS_NONE) .or. &
            (chemistryProps%type /= EQNS_NONE) .or. &
            (subsidence%type /= EQNS_NONE)) then
            call TLab_Write_ASCII(efile, __FILE__//'. IBM not implemented for infraredProps, sedimentationProps, chemistry, subsidence.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
    end if

    ! ###################################################################
    if (max(inb_flow, inb_scal) > MAX_VARS) then
        call TLab_Write_ASCII(efile, __FILE__//'. Error MAX_VARS should be larger than or equal to inb_flow and inb_scal')
        call TLab_Stop(DNS_ERROR_TOTALVARS)
    end if

    ! ###################################################################
    ! if (PressureFilter(1)%type /= DNS_FILTER_NONE) call TLab_Write_ASCII(lfile, 'Pressure and dpdy filter along Ox.')
    ! if (PressureFilter(2)%type /= DNS_FILTER_NONE) call TLab_Write_ASCII(lfile, 'Pressure and dpdy filter along Oy.')
    ! if (PressureFilter(3)%type /= DNS_FILTER_NONE) call TLab_Write_ASCII(lfile, 'Pressure and dpdy filter along Oz.')

    if (any(PressureFilter(:)%type /= DNS_FILTER_NONE)) then
        if (.not. ((imode_eqns == DNS_EQNS_INCOMPRESSIBLE) .or. (imode_eqns == DNS_EQNS_ANELASTIC))) then
            call TLab_Write_ASCII(efile, __FILE__//'. Pressure and dpdy filter only implemented for anelastic or incompressible mode.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. (iadvection == EQNS_CONVECTIVE)) then
            call TLab_Write_ASCII(efile, __FILE__//'. Pressure and dpdy filter not implemented for current advection scheme.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
    end if

    return
end subroutine TLab_Consistency_Check
