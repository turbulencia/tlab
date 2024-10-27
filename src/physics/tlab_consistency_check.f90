#include "dns_const.h"
#include "dns_error.h"

! Everything has been read from input file.
! Check for cross dependencies and undeveloped options.
subroutine TLab_Consistency_Check()
    use TLab_Constants, only: efile
    use TLAB_VARS, only: imode_eqns, iadvection, subsidence
    use IBM_VARS, only: imode_ibm
    use Radiation
    use Microphysics
    use Chemistry
    use SpecialForcing
    use TLab_WorkFlow
    implicit none

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

    return
end subroutine TLab_Consistency_Check
