#include "dns_const.h"

!#######################################################################
!#######################################################################
subroutine RHS_PART_1()

    use TLAB_TYPES, only: pointers_dt, pointers3d_dt
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: g
    use TLAB_VARS, only: visc, radiation, stokes
    use TLAB_ARRAYS
    use DNS_ARRAYS
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_INTERPOLATE
    use THERMO_VARS, only: thermo_param
    use OPR_PARTIAL

    implicit none

! -------------------------------------------------------------------
    real(wp) dummy, dummy2
    integer(wi) bcs(2, 2), nvar
    integer(wi) i
    real(wp) delta_inv0, delta_inv2, delta_inv4

    type(pointers3d_dt), dimension(7) :: data
    type(pointers_dt), dimension(7) :: data_out

! #####################################################################
    bcs = 0

! #####################################################################
! Setting pointers to Eulerian fields that need to be interpolated into particle positions
! #####################################################################
    nvar = 0

    select case (part%type)

    case (PART_TYPE_TRACER)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 1); data_out(nvar)%field => l_hq(:, 1)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 2); data_out(nvar)%field => l_hq(:, 2)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 3); data_out(nvar)%field => l_hq(:, 3)
    
    case (PART_TYPE_INERTIA)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 1); data_out(nvar)%field => l_txc(:, 1)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 2); data_out(nvar)%field => l_txc(:, 2)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 3); data_out(nvar)%field => l_txc(:, 3)
        l_txc(:, 1:3) = 0.0_wp  ! The interpolation routine adds on top of previous data

    case (PART_TYPE_BIL_CLOUD_3, PART_TYPE_BIL_CLOUD_4)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 1); data_out(nvar)%field => l_hq(:, 1)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 2); data_out(nvar)%field => l_hq(:, 2)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, 3); data_out(nvar)%field => l_hq(:, 3)

        dummy2 = -thermo_param(2)
        dummy = -thermo_param(1)

!LAPLACE Xi
        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), s(1, 1), txc(1, 6), txc(1, 3), wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), s(1, 1), txc(1, 5), txc(1, 3), wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), s(1, 1), txc(1, 4), txc(1, 3), wrk2d, wrk3d)

        txc(:, 1) = visc*dummy*(txc(:, 4) + txc(:, 5) + txc(:, 6))

!LAPLACE delta
        call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), s(1, 2), txc(1, 6), txc(1, 3), wrk2d, wrk3d)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), s(1, 2), txc(1, 5), txc(1, 3), wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), s(1, 2), txc(1, 4), txc(1, 3), wrk2d, wrk3d)

        txc(:, 1) = txc(:, 1) + (visc*dummy2*(txc(:, 4) + txc(:, 5) + txc(:, 6))) !first eq. without ds/dxi
        txc(:, 2) = 1.0_wp - dummy*s(:, 1) - dummy2*s(:, 2) !xi field in txc(1,2)

        call FI_GRADIENT(imax, jmax, kmax, txc(1, 2), txc(1, 3), txc(1, 4), wrk2d, wrk3d) ! square of chi gradient in txc(1,3)
        txc(:, 3) = visc*txc(:, 3)

        call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), s(1, radiation%scalar(1)), txc(1, 4), wrk1d, wrk3d)
! Radiation *** ATTENTION RADIATION IS MINUS
        txc(:, 1) = txc(:, 1) + dummy2*txc(:, 4)

! Radiation SECOND FORMULATION *** ATTENTION RADIATION IS MINUS
        txc(:, 4) = dummy2*txc(:, 4)

! Setting pointers
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 1); data_out(nvar)%field => l_txc(:, 1)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 2); data_out(nvar)%field => l_txc(:, 2)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 3); data_out(nvar)%field => l_txc(:, 3)
        nvar = nvar + 1; data(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 4); data_out(nvar)%field => l_txc(:, 4)
        l_txc(:, 1:4) = 0.0_wp

    end select

! -------------------------------------------------------------------
! Interpolating field data into particles
! The interpolated data is added to the existing data, which
! consitutes already the evolution equation for particle position
! -------------------------------------------------------------------
    call FIELD_TO_PARTICLE(data(1:nvar), data_out(1:nvar), l_g, l_q)

! #####################################################################
! Completing evolution equations
! #####################################################################
    select case (part%type)

    case (PART_TYPE_TRACER)
        ! equation dx_p/dt = v_p = v(x_p) already calculated via interpolation step

    case (PART_TYPE_INERTIA)
        dummy = 1.0_wp /stokes
        ! equation dx_p/dt = v_p 
        l_hq(1:l_g%np,1:3) = l_hq(1:l_g%np,1:3) + l_q(1:l_g%np,1:3)
        ! equation dv_p/dt = (v(x_p)-v_p)/stokes
        l_hq(1:l_g%np,4:6) = l_hq(1:l_g%np,4:6) + dummy*(l_txc(1:l_g%np,1:3)-l_q(1:l_g%np,4:6))

    case (PART_TYPE_BIL_CLOUD_3, PART_TYPE_BIL_CLOUD_4)
! l_txc(1) = equation without ds/dxi
! l_txc(2) = xi
! l_txc(3) = evaporation/condensation term without d2s/dxi2
! l_txc(4) = radiation term without ds/dxi

        delta_inv0 = 1.0_wp/thermo_param(1)/thermo_param(3)
        delta_inv2 = -0.5_wp/thermo_param(1)/thermo_param(3)
        delta_inv4 = -0.25_wp/thermo_param(1)/thermo_param(3)

        do i = 1, l_g%np
            l_hq(i, 4) = l_hq(i, 4) - l_txc(i, 1)/(1.0_wp + exp(l_txc(i, 2)*delta_inv0))

            l_hq(i, 5) = l_hq(i, 5) - l_txc(i, 4)/(1.0_wp + exp(l_txc(i, 2)*delta_inv0)) &
                         - l_txc(i, 3)*delta_inv4/(cosh(l_txc(i, 2)*delta_inv2)**2)
        end do

    case (PART_TYPE_SIMPLE_SETT)
        l_hq(1:l_g%np, 2) = l_hq(1:l_g%np, 2) - part%parameters(1)

    end select

! -------------------------------------------------------------------
    do i = 1, nvar
        nullify (data(i)%field, data_out(i)%field)
    end do

    return
end subroutine RHS_PART_1
