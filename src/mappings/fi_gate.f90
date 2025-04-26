#include "dns_const.h"

subroutine FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
                   nx, ny, nz, igate_size, gate_threshold, q, s, txc, gate)
    use TLab_Constants, only: wp, wi, small_wp
    use FDM, only: g
    use FI_VECTORCALCULUS, only: FI_CURL
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN
    use OPR_Partial
    implicit none

    integer(wi), intent(IN) :: opt_cond, opt_cond_relative, opt_cond_scal
    integer(wi), intent(IN) :: nx, ny, nz, igate_size
    real(wp), intent(INOUT) :: gate_threshold(igate_size)
    real(wp), dimension(nx*ny*nz, *), intent(IN) :: q, s
    real(wp), dimension(nx*ny*nz, 5), intent(INOUT) :: txc
    integer(1), dimension(nx*ny*nz), intent(OUT) :: gate

! -----------------------------------------------------------------------
    integer(wi) ij, n
    real(wp) umin, umax, gate_threshold_loc(igate_size)
    integer(wi) bcs(2, 2)

! #######################################################################
! Preparing indicator field

    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    if (opt_cond == 2) then ! Based on scalar
        txc(:, 1) = s(:, opt_cond_scal)

    else if (opt_cond == 3) then ! Based on vorticity
        call FI_VORTICITY(nx, ny, nz, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))

    else if (opt_cond == 4) then ! Based on scalar gradient
        call FI_GRADIENT(nx, ny, nz, s, txc(1, 1), txc(1, 2))

    else if (opt_cond == 5) then ! Based on vertical velocity
        txc(:, 1) = q(:, 2)

    else if (opt_cond == 6 .or. opt_cond == 7) then ! Based on scalar fluctuation
        txc(:, 1) = s(:, opt_cond_scal)
        call FI_FLUCTUATION_INPLACE(nx, ny, nz, txc(1, 1))

    else if (opt_cond == 8) then ! Based on potential vorticity
        call FI_CURL(nx, ny, nz, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))

        txc(:, 4) = s(:, opt_cond_scal)
        call OPR_Partial_X(OPR_P1, nx, ny, nz, bcs, g(1), txc(1, 4), txc(1, 5))
        txc(:, 1) = txc(:, 1)*txc(:, 5)
        call OPR_Partial_Y(OPR_P1, nx, ny, nz, bcs, g(2), txc(1, 4), txc(1, 5))
        txc(:, 1) = txc(:, 1) + txc(:, 2)*txc(:, 5)
        call OPR_Partial_Z(OPR_P1, nx, ny, nz, bcs, g(3), txc(1, 4), txc(1, 5))
        txc(:, 1) = txc(:, 1) + txc(:, 3)*txc(:, 5)

        txc(:, 1) = txc(:, 1)*txc(:, 1)
        txc(:, 1) = log(txc(:, 1) + small_wp)

    end if

! #######################################################################
! define gate field

    if (opt_cond == 7) then ! double conditioning; flux
        do ij = 1, nx*ny*nz
            if (txc(ij, 1) > 0.0_wp .and. q(ij, 2) >= 0.0_wp) then; gate(ij) = 1; 
            else if (txc(ij, 1) <= 0.0_wp .and. q(ij, 2) > 0.0_wp) then; gate(ij) = 2; 
            else if (txc(ij, 1) < 0.0_wp .and. q(ij, 2) <= 0.0_wp) then; gate(ij) = 3; 
            else; gate(ij) = 4; end if
        end do

    else                             ! Local file
        if (opt_cond_relative == 1) then ! case of threshold relative to maximum
            call MINMAX(nx, ny, nz, txc(1, 1), umin, umax)
            gate_threshold_loc = gate_threshold*(umax - umin)

        else
            gate_threshold_loc = gate_threshold

        end if

! Loop over all thresholds (# thresholds is 1 less than # gate levels)
! Thresholds are in ascending order
        do ij = 1, nx*ny*nz
            do n = 1, igate_size - 1
                if (txc(ij, 1) < gate_threshold_loc(n)) exit
            end do
            gate(ij) = int(n, KIND=1) ! note that gate can get -- correctly -- the value igate_size
        end do

    end if

    return
end subroutine FI_GATE
