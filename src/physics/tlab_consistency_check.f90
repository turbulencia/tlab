#include "dns_const.h"
#include "dns_error.h"

! Everything has been read from input file.
! Check for cross dependencies and undeveloped options.
subroutine TLab_Consistency_Check()
    use TLab_Constants, only: efile, lfile, MAX_VARS
    use TLab_Memory, only: inb_flow, inb_flow_array, inb_scal
    use TLab_WorkFlow, only: stagger_on
    use TLab_Time, only: rtime
    use OPR_Filters, only: PressureFilter
    use FDM, only: g
    use IBM_VARS, only: imode_ibm
    use IO_FIELDS
    use NavierStokes
    use Thermodynamics
    use Radiation
    use Microphysics
    use Chemistry
    use SpecialForcing
    use LargeScaleForcing, only: subsidenceProps
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none

    ! -------------------------------------------------------------------
    integer is

    ! ###################################################################
    if (stagger_on) then
        if (.not. ((nse_eqns == DNS_EQNS_INCOMPRESSIBLE) .or. (nse_eqns == DNS_EQNS_ANELASTIC))) then
            call TLab_Write_ASCII(efile, __FILE__//'. Horizontal pressure staggering only implemented for anelastic or incompressible mode.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. ((nse_advection == EQNS_CONVECTIVE) .or. (nse_advection == EQNS_SKEWSYMMETRIC))) then
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
        if (.not. (nse_eqns == DNS_EQNS_INCOMPRESSIBLE)) then
            call TLab_Write_ASCII(efile, __FILE__//'. IBM only implemented for incompressible mode.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. ((nse_advection == EQNS_CONVECTIVE) .or. (nse_advection == EQNS_SKEWSYMMETRIC))) then
            call TLab_Write_ASCII(efile, __FILE__//'. IBM only implemented for convective advection scheme.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        if ((infraredProps%type /= EQNS_NONE) .or. &
            (sedimentationProps%type /= EQNS_NONE) .or. &
            (infraredProps%type /= EQNS_NONE) .or. &
            (chemistryProps%type /= EQNS_NONE) .or. &
            (subsidenceProps%type /= EQNS_NONE)) then
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
        if (.not. ((nse_eqns == DNS_EQNS_INCOMPRESSIBLE) .or. (nse_eqns == DNS_EQNS_ANELASTIC))) then
            call TLab_Write_ASCII(efile, __FILE__//'. Pressure and dpdy filter only implemented for anelastic or incompressible mode.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (.not. (nse_advection == EQNS_CONVECTIVE)) then
            call TLab_Write_ASCII(efile, __FILE__//'. Pressure and dpdy filter not implemented for current advection scheme.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
    end if

    ! ###################################################################
    if (any([EQNS_TRANS_SUTHERLAND, EQNS_TRANS_POWERLAW] == itransport)) inb_flow_array = inb_flow_array + 1    ! space for viscosity

    ! ###################################################################
    if (any([DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL] == nse_eqns) .and. imode_thermo /= THERMO_TYPE_COMPRESSIBLE) then
        call TLab_Write_ASCII(efile, __FILE__//'. Incorrect combination of thermodynamics and type of evolution equations.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    if (nse_eqns == DNS_EQNS_ANELASTIC .and. imode_thermo /= THERMO_TYPE_ANELASTIC) then
        call TLab_Write_ASCII(efile, __FILE__//'. Incorrect combination of thermodynamics and type of evolution equations.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    select case (imixture)
        ! case (MIXT_TYPE_BS, MIXT_TYPE_BSZELDOVICH)
        !     schmidt(inb_scal) = prandtl ! These cases force Sc_i=Sc_Z=Pr (Lewis unity)

    case (MIXT_TYPE_AIRWATER)
        if (any([DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL] == nse_eqns)) schmidt(2:3) = schmidt(1) ! used in diffusion eqns, though should be fixed

        ! if (all([damkohler(1:2)] == 0.0_wp)) then
        !     damkohler(1:2) = damkohler(3)
        ! else
        !     call TLab_Write_ASCII(efile, __FILE__//'. AirWater requires at least first 2 Damkholer numbers zero.')
        !     call TLab_Stop(DNS_ERROR_OPTION)
        ! end if

    end select

    ! ###################################################################
    ! preparing headers of restart files
    io_header_q(1)%size = 0
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = rtime
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = visc
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = froude
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = rossby
    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
        io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = gama0
        io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = prandtl
        io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = mach
    end if

    do is = 1, inb_scal
        io_header_s(is)%size = 0
        io_header_s(is)%size = io_header_s(is)%size + 1; io_header_s(is)%params(io_header_s(is)%size) = rtime
        io_header_s(is)%size = io_header_s(is)%size + 1; io_header_s(is)%params(io_header_s(is)%size) = visc
        io_header_s(is)%size = io_header_s(is)%size + 1; io_header_s(is)%params(io_header_s(is)%size) = schmidt(is)
        io_header_s(is)%size = io_header_s(is)%size + 1; io_header_s(is)%params(io_header_s(is)%size) = damkohler(is)
    end do

    return
end subroutine TLab_Consistency_Check
