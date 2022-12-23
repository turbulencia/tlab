#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!# DESCRIPTION
!#
!# Momentum equations, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!# BCs need 2 1st order derivatives in Oy
!# Scalar needed for the buoyancy term
!#
!########################################################################
subroutine RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1 &
    (u, v, w, h1, h2, h3, q, hq, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, &
     wrk1d, wrk2d, wrk3d)

#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_wrk1d, inb_flow
    use TLAB_VARS, only: g
    use TLAB_VARS, only: visc
    use TIME, only: dte
    use DNS_LOCAL, only: remove_divergence
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use OPR_PARTIAL
    use OPR_BURGERS
    
    implicit none

#include "integers.h"

    TREAL, dimension(isize_field) :: u, v, w, h1, h2, h3
    TREAL, dimension(isize_field, *) :: q, hq
    TREAL, dimension(isize_field) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    TREAL, dimension(isize_wrk1d, *) :: wrk1d
    TREAL, dimension(*) :: wrk2d, wrk3d

    target tmp2, h2

! -----------------------------------------------------------------------
    TINTEGER iq, ij, k, nxy, ip_b, ip_t
    TINTEGER ibc, bcs(2, 2)
    TREAL dummy

    TINTEGER siz, srt, end    !  Variables for OpenMP Partitioning

    TREAL, dimension(:), pointer :: p_bcs

#ifdef USE_ESSL
    integer ilen
#endif

! #######################################################################
    nxy = imax*jmax

    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

#ifdef USE_ESSL
    ilen = isize_field
#endif

! #######################################################################
! Diffusion and convection terms in Ox momentum eqn
! #######################################################################
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), u, tmp6, tmp3, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), u, tmp5, tmp2, wrk2d, wrk3d)
    call OPR_BURGERS_X(i0, i0, imax, jmax, kmax, bcs, g(1), u, u, u, tmp4, tmp1)

!$omp parallel default( shared ) private( ij, srt,end,siz )
    call DNS_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        h1(ij) = h1(ij) + tmp4(ij) + visc*(tmp6(ij) + tmp5(ij)) &
                 - (w(ij)*tmp3(ij) + v(ij)*tmp2(ij))
    end do
!$omp end parallel

! #######################################################################
! Diffusion and convection terms in Oz momentum eqn
! #######################################################################
    if (g(3)%size > 1) then
        call OPR_BURGERS_Z(i0, i0, imax, jmax, kmax, bcs, g(3), w, w, w, tmp6, tmp3)
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), w, tmp5, tmp2, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), w, tmp4, tmp1, wrk2d, wrk3d)

!$omp parallel default( shared ) private( ij, srt,end,siz )
        call DNS_OMP_PARTITION(isize_field, srt, end, siz)
        do ij = srt, end
            h3(ij) = h3(ij) + tmp6(ij) + visc*(tmp5(ij) + tmp4(ij)) &
                     - (v(ij)*tmp2(ij) + u(ij)*tmp1(ij))
        end do
!$omp end parallel

    end if

! #######################################################################
! Diffusion and convection terms in Oy momentum eqn
! #######################################################################
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), v, tmp6, tmp3, wrk2d, wrk3d)
    call OPR_BURGERS_Y(i0, i0, imax, jmax, kmax, bcs, g(2), v, v, v, tmp5, tmp2)
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), v, tmp4, tmp1, wrk2d, wrk3d)

#ifdef USE_ESSL
!$omp parallel default( shared ) &
!$omp private( ilen, srt,end,siz)
#else
!$omp parallel default( shared ) &
!$omp private( ij,   srt,end,siz)
#endif
    call DNS_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        h2(ij) = h2(ij) + tmp5(ij) + visc*(tmp6(ij) + tmp4(ij)) &
                 - (w(ij)*tmp3(ij) + u(ij)*tmp1(ij))
    end do
!$omp end parallel

! #######################################################################
! Impose buffer zone as relaxation terms
! #######################################################################
    if (BuffType == DNS_BUFFER_RELAX .or. BuffType == DNS_BUFFER_BOTH) then
        call BOUNDARY_BUFFER_RELAX_FLOW(q, hq)
    end if

! #######################################################################
! Pressure term
! #######################################################################
    if (remove_divergence) then ! remove residual divergence

#ifdef USE_ESSL
!$omp parallel default( shared )&
!$omp private( ilen, dummy, srt,end,siz )
#else
!$omp parallel default( shared )&
!$omp private( ij,   dummy, srt,end,siz )
#endif

        call DNS_OMP_PARTITION(isize_field, srt, end, siz)
        dummy = C_1_R/dte

#ifdef USE_ESSL
        ilen = siz
        call DZAXPY(ilen, dummy, v(srt), 1, h2(srt), 1, tmp2(srt), 1)
        call DZAXPY(ilen, dummy, u(srt), 1, h1(srt), 1, tmp3(srt), 1)
        call DZAXPY(ilen, dummy, w(srt), 1, h3(srt), 1, tmp4(srt), 1)

#else
        do ij = srt, end
            tmp2(ij) = h2(ij) + v(ij)*dummy
            tmp3(ij) = h1(ij) + u(ij)*dummy
            tmp4(ij) = h3(ij) + w(ij)*dummy
        end do

#endif
!$omp end parallel

        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp3, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp4, tmp3, wrk3d, wrk2d, wrk3d)

    else
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), h2, tmp1, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), h1, tmp2, wrk3d, wrk2d, wrk3d)
        call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), h3, tmp3, wrk3d, wrk2d, wrk3d)

    end if

! -----------------------------------------------------------------------
!$omp parallel default( shared ) private( ij,srt,end,siz )
    call DNS_OMP_PARTITION(isize_field, srt, end, siz)
    do ij = srt, end
        tmp1(ij) = tmp1(ij) + tmp2(ij) + tmp3(ij) ! forcing term in tmp1
    end do
!$omp end parallel

! -----------------------------------------------------------------------
! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    ip_b = 1
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        p_bcs => h2(ip_b:); BcsFlowJmin%ref(1:imax, k, 2) = p_bcs(1:imax); ip_b = ip_b + nxy ! bottom
        p_bcs => h2(ip_t:); BcsFlowJmax%ref(1:imax, k, 2) = p_bcs(1:imax); ip_t = ip_t + nxy ! top
    end do

! pressure in tmp1, Oy derivative in tmp3
    ibc = 3
    call OPR_POISSON_FXZ(.true., imax, jmax, kmax, g, ibc, &
                         tmp1, tmp3, tmp2, tmp4, BcsFlowJmin%ref(1, 1, 2), BcsFlowJmax%ref(1, 1, 2), wrk1d, wrk1d(1, 5), wrk3d)

! horizontal derivatives
    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2, wrk3d, wrk2d, wrk3d)
    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), tmp1, tmp4, wrk3d, wrk2d, wrk3d)

! -----------------------------------------------------------------------
! Add pressure gradient
! -----------------------------------------------------------------------
#ifdef USE_ESSL
!$omp parallel default( shared ) &
!$omp private( ilen, srt,end,siz,dummy )
#else
!$omp parallel default( shared ) &
!$omp private( ij,   srt,end,siz,dummy )
#endif
    call DNS_OMP_PARTITION(isize_field, srt, end, siz)

#ifdef USE_ESSL
    ilen = siz
    dummy = -C_1_R
    call DAXPY(ilen, dummy, tmp2(srt), 1, h1(srt), 1)
    call DAXPY(ilen, dummy, tmp3(srt), 1, h2(srt), 1)
    call DAXPY(ilen, dummy, tmp4(srt), 1, h3(srt), 1)
#else
    do ij = srt, end
        h1(ij) = h1(ij) - tmp2(ij)
        h2(ij) = h2(ij) - tmp3(ij)
        h3(ij) = h3(ij) - tmp4(ij)
    end do
#endif

!$omp end parallel

! #######################################################################
! Boundary conditions
! #######################################################################
! -----------------------------------------------------------------------
! Preliminaries
! -----------------------------------------------------------------------
    BcsFlowJmin%ref = C_0_R ! default is no-slip (dirichlet)
    BcsFlowJmax%ref = C_0_R

    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowJmin%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 1
        if (BcsFlowJmax%type(iq) == DNS_BCS_NEUMANN) ibc = ibc + 2
        if (ibc > 0) then
            call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hq(1, iq), &
                                        BcsFlowJmin%ref(1, 1, iq), BcsFlowJmax%ref(1, 1, iq), wrk1d, tmp1, wrk3d)
        end if
    end do

! -----------------------------------------------------------------------
! Impose bottom BCs at Jmin
! -----------------------------------------------------------------------
    ip_b = 1
    do k = 1, kmax
        h1(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 1)
        h2(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 2)
        h3(ip_b:ip_b + imax - 1) = BcsFlowJmin%ref(1:imax, k, 3); ip_b = ip_b + nxy
    end do

! -----------------------------------------------------------------------
! Impose top BCs at Jmax
! -----------------------------------------------------------------------
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        h1(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 1)
        h2(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 2)
        h3(ip_t:ip_t + imax - 1) = BcsFlowJmax%ref(1:imax, k, 3); ip_t = ip_t + nxy
    end do

    return
end subroutine RHS_FLOW_GLOBAL_INCOMPRESSIBLE_1
