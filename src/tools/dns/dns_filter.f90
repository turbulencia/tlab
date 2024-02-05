#include "dns_const.h"

!########################################################################
!########################################################################
subroutine DNS_FILTER()

    use TLAB_VARS, only: imax, jmax, kmax, inb_flow, inb_scal
    use TLAB_VARS, only: imode_eqns, imode_sim
    use TLAB_VARS, only: itime, rtime
    use TLAB_VARS, only: FilterDomain
    use TLAB_VARS, only: g, area
    use TLAB_ARRAYS
    use OPR_FILTERS
    use DNS_LOCAL, only: DNS_BOUNDS_LIMIT
    use DNS_LOCAL, only: nitera_stats_spa, nitera_first, nitera_stats
    use STATISTICS
    use AVGS, only: AVG_IK_V

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
        call AVG_IK_V(imax, jmax, kmax, jmax, wrk3d, g(1)%jac, g(3)%jac, Tke0(1), wrk1d, area)
        call FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5))
        call AVG_IK_V(imax, jmax, kmax, jmax, txc, g(1)%jac, g(3)%jac, Eps0(1), wrk1d, area)
    end if

    ! -------------------------------------------------------------------
    ! filtering

    ! Might be better to filter the pressure instead of the energy in compressible flows
    ! IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
    !    iq_loc = (/ 5,1,2,3,6 /) ! Filtered variables: rho, u,v,w, p
    ! ELSE
    !    iq_loc = (/ 1,2,3 /)
    ! ENDIF

    if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then ! contruct fields per unit volume
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

    if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then ! re-contruct fields per unit mass
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
        call AVG_IK_V(imax, jmax, kmax, jmax, wrk3d, g(1)%jac, g(3)%jac, Tke1(1), wrk1d, area)
        call FI_DISSIPATION(i1, imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1), txc(1,2),txc(1,3),txc(1,4),txc(1,5))
        call AVG_IK_V(imax, jmax, kmax, jmax, txc, g(1)%jac, g(3)%jac, Eps1(1), wrk1d, area)

        write (fname, *) itime; fname = 'kin'//trim(adjustl(fname))
        call IO_WRITE_AVERAGES(fname, itime, rtime, jmax, 4, 1, g(2)%nodes, varnames, groupnames, mean)

    end if

    ! -------------------------------------------------------------------
    ! recalculation of diagnostic variables
    call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)

    return
end subroutine DNS_FILTER
