#include "dns_const.h"

!########################################################################
!#
!# Scalar equation, nonlinear term in convective form and the
!# diffusion term explicit: 3 2nd order + 3 1st order derivatives.
!# BCs need 1 1st order derivatives in Oy
!#
!########################################################################
subroutine RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1(is)
    use TLAB_CONSTANTS, only: wp, wi
#ifdef USE_OPENMP
    use OMP_LIB
#endif
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: g
    use TLAB_VARS, only: idiffusion, visc, schmidt
    use BOUNDARY_BCS, only: BcsScalJmin, BcsScalJmax
    use TLAB_ARRAYS, only: s, wrk1d, wrk2d, wrk3d
    use TLAB_POINTERS, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    use DNS_ARRAYS, only: hs
    use OPR_PARTIAL

    implicit none

    integer(wi), intent(in) :: is

    real(wp), dimension(:), pointer :: p_bcs

! -----------------------------------------------------------------------
    integer(wi) ij, k, nxy, ip, ip_b, ip_t, ibc, bcs(2, 2)
    real(wp) diff

! #######################################################################
    nxy = imax*jmax

    bcs = 0

    if (idiffusion == EQNS_NONE) then; diff = 0.0_wp
    else; diff = visc/schmidt(is); end if

! #######################################################################
! Preliminaries for Scalar BC
! (flow BCs initialized below as they are used for pressure in between)
! #######################################################################
! Default is zero (Dirichlet)
    BcsScalJmin%ref(:, :, is) = 0.0_wp
    BcsScalJmax%ref(:, :, is) = 0.0_wp

! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
    ip_b = 1
    ip_t = imax*(jmax - 1) + 1
    do k = 1, kmax
        if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) then
            p_bcs => hs(ip_b:,is); BcsScalJmin%ref(1:imax, k, is) = p_bcs(1:imax); end if
        if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) then
            p_bcs => hs(ip_t:,is); BcsScalJmax%ref(1:imax, k, is) = p_bcs(1:imax); end if
        ip_b = ip_b + nxy ! bottom BC address
        ip_t = ip_t + nxy ! top BC address
    end do

! #######################################################################
! Diffusion and convection terms in scalar equations
! #######################################################################
    call OPR_PARTIAL_Z(OPR_P2_P1, imax, jmax, kmax, bcs, g(3), s(:,is), tmp6, tmp3, wrk2d, wrk3d)
    call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), s(:,is), tmp5, tmp2, wrk2d, wrk3d)
    call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), s(:,is), tmp4, tmp1, wrk2d, wrk3d)

!$omp parallel default( shared ) private( ij )
!$omp do
    do ij = 1, isize_field
        hs(ij,is) = hs(ij,is) + diff*(tmp6(ij) + tmp5(ij) + tmp4(ij)) &
                    - (w(ij)*tmp3(ij) + v(ij)*tmp2(ij) + u(ij)*tmp1(ij))
    end do
!$omp end do
!$omp end parallel

! #######################################################################
! Boundary conditions
! #######################################################################
    ibc = 0
    if (BcsScalJmin%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 1
    if (BcsScalJmax%type(is) == DNS_BCS_NEUMANN) ibc = ibc + 2
    if (ibc > 0) then
        call BOUNDARY_BCS_NEUMANN_Y(ibc, imax, jmax, kmax, g(2), hs(:,is), &
                                    BcsScalJmin%ref(1, 1, is), BcsScalJmax%ref(1, 1, is), wrk1d, tmp1, wrk3d)
    end if

    if (BcsScalJmin%type(is) /= DNS_SFC_STATIC .or. &
        BcsScalJmax%type(is) /= DNS_SFC_STATIC) then
        call BOUNDARY_SURFACE_J(is, bcs, s(:,is), hs(:,is), tmp1, tmp2, tmp3, wrk1d, wrk2d, wrk3d)
    end if

! -----------------------------------------------------------------------
! Impose BC at Jmin/max
! -----------------------------------------------------------------------
    ip = 1
    ip_t = 1 + imax*(jmax - 1)
    do k = 1, kmax
        hs(ip_b:ip_b + imax - 1,is) = BcsScalJmin%ref(1:imax, k, is); ip_b = ip + nxy
        hs(ip_t:ip_t + imax - 1,is) = BcsScalJmax%ref(1:imax, k, is); ip_t = ip + nxy
    end do

    return
end subroutine RHS_SCAL_GLOBAL_INCOMPRESSIBLE_1
