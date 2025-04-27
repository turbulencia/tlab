#include "dns_const.h"

!########################################################################
!########################################################################
subroutine DNS_FILTER()

    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_flow, inb_scal
    use NavierStokes, only: nse_eqns, DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL
    use TLab_WorkFlow, only: imode_sim
    use NavierStokes, only: visc
    use TLab_Time, only: itime, rtime
    use FDM, only: g
    use TLab_Arrays
    use OPR_FILTERS
    use DNS_LOCAL, only: DNS_BOUNDS_LIMIT
    use DNS_LOCAL, only: nitera_stats_spa, nitera_first, nitera_stats
    use DNS_STATISTICS, only: mean_flow, mean
    use Averages, only: AVG_IK_V

    implicit none

    ! -----------------------------------------------------------------------
    integer iq, is
    integer, parameter :: i1 = 1
    character*64 fname, varnames(1), groupnames(1)

    ! #######################################################################
    ! Statistics
#define Tke0(j)   mean(j,1)
#define Eps0(j)   mean(j,2)
#define Tke1(j)   mean(j,3)
#define Eps1(j)   mean(j,4)

    varnames(1) = 'TkeBefore EpsBefore TkeAfter EpsAfter'
    groupnames(1) = ''

    if (mod(itime - nitera_first, nitera_stats_spa) == 0) then   ! Accumulate statistics in spatially evolving cases
        call AVG_TKE_ZT_REDUCE(q(1, 5), q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), &
                               txc(1, 5), txc(1, 6), txc(1, 7), mean_flow)
    end if

    if (imode_sim == DNS_MODE_TEMPORAL .and. mod(itime - nitera_first, nitera_stats) == 0) then
        call FI_RTKE(imax, jmax, kmax, q, wrk3d)
        call AVG_IK_V(imax, jmax, kmax, wrk3d, Tke0(1), wrk1d)
        call FI_DISSIPATION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
        txc(1:isize_field, 1) = txc(1:isize_field, 1)*visc
        call AVG_IK_V(imax, jmax, kmax, txc(:, 1), Eps0(1), wrk1d)
    end if

    ! -------------------------------------------------------------------
    ! filtering

    ! Might be better to filter the pressure instead of the energy in compressible flows
    ! IF ( nse_eqns .EQ. DNS_EQNS_TOTAL .OR. nse_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
    !    iq_loc = (/ 5,1,2,3,6 /) ! Filtered variables: rho, u,v,w, p
    ! ELSE
    !    iq_loc = (/ 1,2,3 /)
    ! ENDIF

    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then ! contruct fields per unit volume
        do iq = 1, inb_flow - 1
            q(:, iq) = q(:, iq)*q(:, inb_flow)
        end do
        do is = 1, inb_scal
            s(:, is) = s(:, is)*q(:, inb_flow)
        end do
    end if

    do iq = 1, inb_flow
        call OPR_FILTER(imax, jmax, kmax, FilterDomain, q(1, iq), txc)
    end do
    do is = 1, inb_scal
        call OPR_FILTER(imax, jmax, kmax, FilterDomain, s(1, is), txc)
    end do

    if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then ! re-contruct fields per unit mass
        do iq = 1, inb_flow - 1
            q(:, iq) = q(:, iq)/q(:, inb_flow)
        end do
        do is = 1, inb_scal
            s(:, is) = s(:, is)/q(:, inb_flow)
        end do
    end if

    ! -------------------------------------------------------------------
    call DNS_BOUNDS_LIMIT()

    ! -------------------------------------------------------------------
    ! statistics
    if (imode_sim == DNS_MODE_TEMPORAL .and. mod(itime - nitera_first, nitera_stats) == 0) then
        call FI_RTKE(imax, jmax, kmax, q, wrk3d)
        call AVG_IK_V(imax, jmax, kmax, wrk3d, Tke1(1), wrk1d)
        call FI_DISSIPATION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
        txc(1:isize_field, 1) = txc(1:isize_field, 1)*visc
        call AVG_IK_V(imax, jmax, kmax, txc, Eps1(1), wrk1d)

        write (fname, *) itime; fname = 'kin'//trim(adjustl(fname))
        call IO_WRITE_AVERAGES(fname, itime, rtime, jmax, 4, 1, g(2)%nodes, varnames, groupnames, mean)

    end if

    ! -------------------------------------------------------------------
    ! recalculation of diagnostic variables
    call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)

    return
end subroutine DNS_FILTER

! #######################################################################
! Calculate kinetic energy of fluctuating field per unit volume
! #######################################################################
#define rR(j)     wrk1d(j,1)
#define fU(j)     wrk1d(j,2)
#define fV(j)     wrk1d(j,3)
#define fW(j)     wrk1d(j,4)
#define aux(j)    wrk1d(j,5)

subroutine FI_RTKE(nx, ny, nz, q, ke)
    use TLab_Constants, only: wp, wi
    use NavierStokes, only: nse_eqns, DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL, DNS_EQNS_ANELASTIC, DNS_EQNS_INCOMPRESSIBLE
    use TLab_Memory, only: inb_flow
    use TLab_Arrays, only: wrk1d
    use THERMO_ANELASTIC, only : rbackground
    use Averages, only: AVG_IK_V

    implicit none

    integer(wi) nx, ny, nz
    real(wp), intent(in) :: q(nx, ny, nz, inb_flow)
    real(wp), intent(out) :: ke(nx, ny, nz)

    ! -----------------------------------------------------------------------
    integer(wi) j

    ! #######################################################################
    select case (nse_eqns)
    case (DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL)
        call AVG_IK_V(nx, ny, nz, q(1, 1, 1, 5), rR(1), aux(1))

        ke = q(:, :, :, 5)*q(:, :, :, 1)
        call AVG_IK_V(nx, ny, nz, ke, fU(1), aux(1))
        fU(:) = fU(:)/rR(:)

        ke = q(:, :, :, 5)*q(:, :, :, 2)
        call AVG_IK_V(nx, ny, nz, ke, fV(1), aux(1))
        fV(:) = fV(:)/rR(:)

        ke = q(:, :, :, 5)*q(:, :, :, 3)
        call AVG_IK_V(nx, ny, nz, ke, fW(1), aux(1))
        fW(:) = fW(:)/rR(:)

    case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
        rR(:) = rbackground(:)
        call AVG_IK_V(nx, ny, nz, q(1, 1, 1, 1), fU(1), aux(1))
        call AVG_IK_V(nx, ny, nz, q(1, 1, 1, 2), fV(1), aux(1))
        call AVG_IK_V(nx, ny, nz, q(1, 1, 1, 3), fW(1), aux(1))

    end select

    do j = 1, ny
        ke(:, j, :) = 0.5_wp*rR(j)*((q(:, j, :, 1) - fU(j))**2 + (q(:, j, :, 2) - fV(j))**2 + (q(:, j, :, 3) - fW(j))**2)
    end do

    return
end subroutine FI_RTKE

