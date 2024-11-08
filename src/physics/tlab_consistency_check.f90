#include "dns_const.h"
#include "dns_error.h"

! Everything has been read from input file.
! Check for cross dependencies and undeveloped options.
subroutine TLab_Consistency_Check()
    use TLab_Constants, only: efile, lfile, MAX_VARS
    use TLAB_VARS, only: imode_eqns, iadvection, subsidence
    use TLAB_VARS, only: inb_flow, inb_flow_array, inb_scal
    use TLAB_VARS, only: stagger_on, PressureFilter
    use TLAB_VARS, only: schmidt
    use FDM, only: g
    use IBM_VARS, only: imode_ibm
    use Thermodynamics, only: imixture, itransport
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

    ! ###################################################################
    if (any([EQNS_TRANS_SUTHERLAND, EQNS_TRANS_POWERLAW] == itransport)) inb_flow_array = inb_flow_array + 1    ! space for viscosity

    ! ###################################################################
    if (imode_eqns == DNS_EQNS_ANELASTIC .and. all([MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER] /= imixture)) then
        call TLab_Write_ASCII(efile, __FILE__//'. Incorrect mixture type.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    select case (imixture)
        ! case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH)
        !     schmidt(inb_scal) = prandtl ! These cases force Sc_i=Sc_Z=Pr (Lewis unity)

    case (MIXT_TYPE_AIRWATER)
        if (any([DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL] == imode_eqns)) schmidt(2:3) = schmidt(1) ! used in diffusion eqns, though should be fixed

        ! if (all([damkohler(1:2)] == 0.0_wp)) then
        !     damkohler(1:2) = damkohler(3)
        ! else
        !     call TLab_Write_ASCII(efile, __FILE__//'. AirWater requires at least first 2 Damkholer numbers zero.')
        !     call TLab_Stop(DNS_ERROR_OPTION)
        ! end if

    end select

    return
end subroutine TLab_Consistency_Check
