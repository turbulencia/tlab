#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!#
!# Implementation of skewsymmetric formulation and viscous explicit
!# in one routine, to avoid duplication of derivatives. Internal energy
!# formulation only.
!# 30 first derivative operations and 9 second derivative operations.
!# Viscosity needs to be homogeneous, because of the explicit treatment
!# of diffusion terms.
!#
!########################################################################
subroutine RHS_FLOW_GLOBAL_2()
    use TLAB_CONSTANTS, only: efile, wp, wi
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
    use TLAB_PROCS, only: TLAB_WRITE_ASCII
#endif
    use TLAB_VARS, only: imax, jmax, kmax, inb_scal
    use TLAB_VARS, only: g, buoyancy
    use TLAB_VARS, only: idiffusion, visc, prandtl, mach
    use TLAB_ARRAYS, only: wrk2d, wrk3d
    use TLAB_POINTERS
    use DNS_ARRAYS
    use THERMO_VARS, only: gama0
    use BOUNDARY_BCS
    use OPR_PARTIAL

#ifdef USE_OPENMP
    use OMP_LIB
#endif

    implicit none

! -------------------------------------------------------------------
    integer(wi) bcs(2, 1), i, is
    real(wp) g1, g2, g3, prefactor, cond, dummy, c13, dum1, dum2, dum3

! ###################################################################
#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'ENTERING RHS_FLOW_GLOBAL_2')
#endif

    bcs = 0

    g1 = buoyancy%vector(1)
    g2 = buoyancy%vector(2)
    g3 = buoyancy%vector(3)
    prefactor = (gama0 - 1.0_wp)*mach*mach
    c13 = 1.0_wp/3.0_wp

! ###################################################################
! Terms \rho u in mass, u-momentum equations
! ###################################################################
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*u(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i) + p(i)
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp4, tmp6, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 5) = hq(i, 5) - 2.0_wp*tmp6(i)
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! u-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) - (tmp2(i) + tmp3(i) + tmp4(i)) + g1*rho(i)
    end do
!$omp end do
!$omp end parallel

! ###################################################################
! Terms \rho v in mass, v-momentum equations
! ###################################################################
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*v(i)
        tmp3(i) = tmp4(i)*w(i)
        tmp2(i) = tmp4(i)*v(i) + p(i)
        tmp1(i) = tmp4(i)*u(i)
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp4, tmp5, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 5) = hq(i, 5) - 2.0_wp*tmp5(i)
        tmp6(i) = tmp6(i) + tmp5(i)
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! v-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 2) = hq(i, 2) - (tmp2(i) + tmp3(i) + tmp4(i)) + g2*rho(i)
    end do
!$omp end do
!$omp end parallel

! ###################################################################
! Terms \rho w in mass, w-momentum equations
! ###################################################################
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        tmp4(i) = 0.5_wp*rho(i)*w(i)
        tmp3(i) = tmp4(i)*w(i) + p(i)
        tmp2(i) = tmp4(i)*v(i)
        tmp1(i) = tmp4(i)*u(i)
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! mass
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp5, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 5) = hq(i, 5) - 2.0_wp*tmp5(i)
        tmp6(i) = tmp6(i) + tmp5(i)
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! w-momentum
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 3) = hq(i, 3) - (tmp2(i) + tmp3(i) + tmp4(i)) + g3*rho(i)
    end do
!$omp end do
!$omp end parallel

! ###################################################################
! Term \rho e in energy equation
! ###################################################################
!$omp parallel default( shared ) private( i, dummy )
!$omp do
    do i = 1, imax*jmax*kmax
        dummy = 0.5_wp*rho(i)*e(i)
        tmp3(i) = dummy*w(i)
        tmp2(i) = dummy*v(i)
        tmp1(i) = dummy*u(i)
    end do
!$omp end do
!$omp end parallel
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp3, tmp4, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp3, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 4) = hq(i, 4) - (tmp2(i) + tmp3(i) + tmp4(i))
    end do
!$omp end do
!$omp end parallel

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms u_i d(\rho u_k)/dx_k
! ###################################################################
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) - u(i)*tmp6(i)
        hq(i, 2) = hq(i, 2) - v(i)*tmp6(i)
        hq(i, 3) = hq(i, 3) - w(i)*tmp6(i)
        hq(i, 4) = hq(i, 4) - e(i)*tmp6(i)
    end do
!$omp end do
!$omp end parallel

    do is = 1, inb_scal
!$omp parallel default( shared ) private( i )
!$omp do
        do i = 1, imax*jmax*kmax
            hs(i, is) = hs(i, is) - s(i, is)*tmp6(i)
        end do
!$omp end do
!$omp end parallel
    end do

! ###################################################################
! Additional convective part due to skewsymmetric formulation
! Terms \rho u_k du_i/dx_k
! First derivatives in internal energy equation from dissipation
! Second derivative terms in the momentum equation
! Note that second derivative routines give also first derivatives !
! ###################################################################
! -------------------------------------------------------------------
! Cross derivatives
! -------------------------------------------------------------------
! energy equation
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), e, tmp4, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), e, tmp3, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), e, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 4) = hq(i, 4) - 0.5_wp*rho(i)*(u(i)*tmp2(i) + v(i)*tmp3(i) + w(i)*tmp4(i))
    end do
!$omp end do
!$omp end parallel

! momentum equations
    dummy = prefactor*visc

    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), u, tmp5, tmp1, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), u, tmp6, tmp2, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) + visc*(tmp5(i) + tmp6(i)) - 0.5_wp*rho(i)*(v(i)*tmp1(i) + w(i)*tmp2(i))
    end do
!$omp end do
!$omp end parallel
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), v, tmp5, tmp3, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 3), g(3), v, tmp6, tmp4, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i, dum1 )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 2) = hq(i, 2) + visc*(tmp5(i) + tmp6(i)) - 0.5_wp*rho(i)*(u(i)*tmp3(i) + w(i)*tmp4(i))
        dum1 = tmp1(i) + tmp3(i)
        hq(i, 4) = hq(i, 4) + dummy*dum1*dum1
    end do
!$omp end do
!$omp end parallel

! arrays tmp1 and tmp3 can be reused
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 1), g(1), w, tmp5, tmp1, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs_out(:, :, 2), g(2), w, tmp6, tmp3, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i, dum2, dum3 )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 3) = hq(i, 3) + visc*(tmp5(i) + tmp6(i)) - 0.5_wp*rho(i)*(u(i)*tmp1(i) + v(i)*tmp3(i))
        dum2 = tmp4(i) + tmp3(i)
        dum3 = tmp1(i) + tmp2(i)
        hq(i, 4) = hq(i, 4) + dummy*(dum2*dum2 + dum3*dum3)
    end do
!$omp end do
!$omp end parallel

! -------------------------------------------------------------------
! Dilatation part
! -------------------------------------------------------------------
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), w, tmp6, tmp3, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), v, tmp5, tmp2, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs_inf(:, :, 1), g(1), u, tmp4, tmp1, wrk2d, wrk3d)
    dummy = 2.0_wp*visc
!$omp parallel default( shared ) private( i, dum1 )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) - 0.5_wp*rho(i)*u(i)*tmp1(i)
        hq(i, 2) = hq(i, 2) - 0.5_wp*rho(i)*v(i)*tmp2(i)
        hq(i, 3) = hq(i, 3) - 0.5_wp*rho(i)*w(i)*tmp3(i)

        dum1 = tmp1(i) + tmp2(i) + tmp3(i)
        hq(i, 4) = hq(i, 4) + prefactor*( &
                   dummy*(tmp3(i)*tmp3(i) + tmp2(i)*tmp2(i) + tmp1(i)*tmp1(i) - c13*dum1*dum1) - p(i)*dum1)
! array tmp1 no longer needed
        tmp1(i) = c13*dum1
    end do
!$omp end do
!$omp end parallel

! Second derivative terms in the momentum equation
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 1), g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 1) = hq(i, 1) + visc*(tmp4(i) + tmp2(i))
    end do
!$omp end do
!$omp end parallel

    call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 2), g(2), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 2) = hq(i, 2) + visc*(tmp5(i) + tmp2(i))
    end do
!$omp end do
!$omp end parallel

    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs_inf(:, :, 3), g(3), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 3) = hq(i, 3) + visc*(tmp6(i) + tmp2(i))
    end do
!$omp end do
!$omp end parallel

! ###################################################################
! Enthalpy diffusion in energy equation
! ###################################################################
    if (idiffusion == EQNS_NONE) then; cond = 0.0_wp
    else; cond = visc/prandtl; end if

! calculate the enthalpy
    call THERMO_CALORIC_ENTHALPY(imax, jmax, kmax, s, T, tmp4)

! total flux
    call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 3), g(3), tmp4, tmp3, tmp5, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 2), g(2), tmp4, tmp2, tmp5, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs_out(:, :, 1), g(1), tmp4, tmp1, tmp5, wrk2d, wrk3d)
!$omp parallel default( shared ) private( i )
!$omp do
    do i = 1, imax*jmax*kmax
        hq(i, 4) = hq(i, 4) + cond*(tmp1(i) + tmp2(i) + tmp3(i))
    end do
!$omp end do
!$omp end parallel

#ifdef TRACE_ON
    call TLAB_WRITE_ASCII(tfile, 'LEAVING RHS_FLOW_GLOBAL_2')
#endif

    return
end subroutine RHS_FLOW_GLOBAL_2
