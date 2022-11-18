#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY
!#
!# 2015/06/25 - A de Lozar
!#              Created
!# 2016/06/25 - JP Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the transport terms due to settling in an airwater mixture.
!#
!########################################################################
subroutine FI_TRANSPORT(transport, flag_grad, nx, ny, nz, is, s, trans, tmp, wrk2d, wrk3d)

    use TLAB_TYPES, only: term_dt
    use TLAB_VARS, only: g, epbackground, inb_scal_array

    implicit none

    type(term_dt), intent(IN) :: transport
    TINTEGER, intent(IN) :: nx, ny, nz, flag_grad
    TINTEGER, intent(IN) :: is
    TREAL, dimension(nx*ny*nz, *), intent(IN) :: s
    TREAL, dimension(nx*ny*nz, 1), intent(OUT) :: trans ! Transport component. It could have eventually three directions
    TREAL, dimension(nx*ny*nz, 1), intent(INOUT) :: tmp   ! To avoid re-calculations when repetedly calling this routine
    TREAL, dimension(*), intent(INOUT) :: wrk2d, wrk3d

! -----------------------------------------------------------------------
    TREAL dummy, exponent
    TINTEGER is_ref, bcs(2, 2)

!########################################################################
    bcs = 0

    exponent = transport%auxiliar(1)
    is_ref = transport%scalar(1)

    if (transport%type == EQNS_TRANS_AIRWATERSIMPLIFIED) then
        if (flag_grad == 1) then
            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), s(1, is_ref), tmp, wrk3d, wrk2d, wrk3d)
            if (exponent > C_0_R) tmp(:, 1) = tmp(:, 1)*(s(:, is_ref)**exponent)
        end if

        dummy = transport%parameters(is)*(C_1_R + exponent)
        trans(:, 1) = dummy*tmp(:, 1)

    elseif (transport%type == EQNS_TRANS_AIRWATER) then
        dummy = C_1_R + exponent

        select case (is)
        case (2, 3)         ! q_t, q_l
            if (exponent > C_0_R) then
                tmp(:, 1) = (transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*(s(:, is_ref)**dummy)
            else
                tmp(:, 1) = (transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*s(:, is_ref)
            end if

        case default        ! energy variables
            call THERMO_ANELASTIC_STATIC_L(nx, ny, nz, s, epbackground, tmp(:, 1))
            if (exponent > C_0_R) then
                tmp(:, 1) = transport%parameters(is)*tmp(:, 1)*(s(:, is_ref)**dummy)
            else
                tmp(:, 1) = transport%parameters(is)*tmp(:, 1)*s(:, is_ref)
            end if

        end select

        call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp(1, 1), trans(1, 1), wrk3d, wrk2d, wrk3d)

    end if

    return
end subroutine FI_TRANSPORT

!########################################################################
!########################################################################
subroutine FI_TRANSPORT_FLUX(transport, nx, ny, nz, is, s, trans)

    use TLAB_TYPES, only: term_dt
    use TLAB_VARS, only: epbackground, inb_scal_array

    implicit none

    type(term_dt), intent(IN) :: transport
    TINTEGER, intent(IN) :: nx, ny, nz
    TINTEGER, intent(IN) :: is
    TREAL, dimension(nx*ny*nz, *), intent(IN) :: s
    TREAL, dimension(nx*ny*nz, 1), intent(OUT) :: trans ! Transport component. It could have eventually three directions

! -----------------------------------------------------------------------
    TREAL dummy, exponent
    TINTEGER is_ref

!########################################################################
    exponent = transport%auxiliar(1)
    is_ref = transport%scalar(1)

    if (transport%type == EQNS_TRANS_AIRWATERSIMPLIFIED) then
        dummy = C_1_R + exponent

        trans(:, 1) = -transport%parameters(is)*(s(:, is_ref)**dummy)

    elseif (transport%type == EQNS_TRANS_AIRWATER) then
        dummy = C_1_R + exponent

        select case (is)
        case (2, 3)         ! q_t, q_l
            if (exponent > C_0_R) then
                trans(:, 1) = -(transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*(s(:, is_ref)**dummy)
            else
                trans(:, 1) = -(transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*s(:, is_ref)
            end if

        case default        ! energy variables
            call THERMO_ANELASTIC_STATIC_L(nx, ny, nz, s, epbackground, trans(:, 1))
            if (exponent > C_0_R) then
                trans(:, 1) = -transport%parameters(is)*trans(:, 1)*(s(:, is_ref)**dummy)
            else
                trans(:, 1) = -transport%parameters(is)*trans(:, 1)*s(:, is_ref)
            end if

        end select

    end if

    return
end subroutine FI_TRANSPORT_FLUX
