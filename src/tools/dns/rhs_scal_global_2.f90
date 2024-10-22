#include "dns_const.h"

!########################################################################
!#
!# Modified from RHS_SCAL_EULER_SKEWSYMMETRIC to include diffusion terms
!# from RHS_SCAL_DIFFUSION_EXPLICIT and avoid duplication of derivatives
!# in routines OPR_PARTIAL_XX, OPR_PARTIAL_YY, OPR_PARTIAL_ZZ.
!# Internal energy formulation only.
!# Additional convective part due to skewsymmetric formulation Y_i d(\rho u_k)/dx_k
!# done in RHS_FLOW_GLOBAL_2
!#
!########################################################################
subroutine RHS_SCAL_GLOBAL_2(is)

    use TLab_Constants, only: efile, wp, wi
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
    use TLab_WorkFlow, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_VARS, only: g
    use TLAB_VARS, only: idiffusion, visc, prandtl, schmidt
    use TLAB_ARRAYS, only: s
    use TLAB_POINTERS, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, T, rho
    use DNS_ARRAYS, only: hs, hq
    use Thermodynamics, only: imixture, THERMO_AI, THERMO_TLIM, NSP, NCP
    use BOUNDARY_BCS
    use OPR_PARTIAL

#ifdef USE_OPENMP
    use OMP_LIB
#endif

    implicit none

    integer(wi), intent(in) :: is

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i, im, icp
    real(wp) diff, cond, dummy

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_SCAL_GLOBAL_2')
#endif

    if (idiffusion == EQNS_NONE) then; diff = 0.0_wp; cond = 0.0_wp
    else; diff = visc/schmidt(is); cond = visc/prandtl; end if

! ###################################################################
! divergence terms
! ###################################################################
!$omp parallel default( shared ) private( i, dummy )
!$omp do
    do i = 1, imax*jmax*kmax
        dummy = 0.5_wp*rho(i)*s(i, is)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
    end do
!$omp end do
!$omp end parallel
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hs(i, is) = hs(i, is) - (tmp2(i) + tmp3(i) + tmp4(i))
    end do
!$omp end do
!$omp end parallel

! ###################################################################
! convective part + diffusion
! ###################################################################
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), s(:, is), tmp6, tmp3)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), s(:, is), tmp5, tmp2)
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), s(:, is), tmp4, tmp1)

!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hs(i, is) = hs(i, is) - 0.5_wp*rho(i)*(u(i)*tmp1(i) + v(i)*tmp2(i) + w(i)*tmp3(i)) &
                    + diff*(tmp4(i) + tmp5(i) + tmp6(i))
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! enthalpy transport by diffusion velocities
! -------------------------------------------------------------------
! enthalpy difference (h_is-h_NSP) is calculated in array tmp4
    if (imixture > 0 .and. is < NSP .and. schmidt(is) /= prandtl) then
        tmp4(:) = 0.0_wp
        do i = 1, imax*jmax*kmax
            if (T(i) < THERMO_TLIM(3, is)) then; im = 2
            else; im = 1; end if
            do icp = NCP, 1, -1
                tmp4(i) = tmp4(i)*T(i) &
                          + (THERMO_AI(icp, im, is) - THERMO_AI(icp, im, NSP))/real(icp, wp)
            end do
! factor (diff-cond) added now
            tmp4(i) = (diff - cond)*(tmp4(i)*T(i) + THERMO_AI(6, im, is) - THERMO_AI(6, im, NSP))
        end do
!$omp parallel default( shared ) private( i )
!$omp do
        do i = 1, imax*jmax*kmax
            hq(i, 4) = hq(i, 4) + tmp4(i)*(tmp1(i) + tmp2(i) + tmp3(i))
        end do
!$omp end do
!$omp end parallel

! cross-gradients
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp1)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp2)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(:, is), tmp4)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(:, is), tmp5)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(:, is), tmp6)
!$omp parallel default( shared ) private( i )
!$omp do
        do i = 1, imax*jmax*kmax
            hq(i, 4) = hq(i, 4) + (tmp1(i)*tmp4(i) + tmp2(i)*tmp5(i) + tmp3(i)*tmp6(i))
        end do
!$omp end do
!$omp end parallel

    end if

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_SCAL_GLOBAL_2')
#endif

    return
end subroutine RHS_SCAL_GLOBAL_2
