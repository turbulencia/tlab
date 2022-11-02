#include "types.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!#
!# Initialize data of reference thermodynamic profiles
!#
!########################################################################
subroutine FI_BACKGROUND_INITIALIZE(wrk1d)

    use TLAB_CONSTANTS, only: lfile
    use TLAB_VARS, only: inb_scal, inb_scal_array, imax, jmax, kmax, imode_eqns
    use TLAB_VARS, only: g
    use TLAB_VARS, only: qbg, pbg, rbg, tbg, hbg, sbg
    use TLAB_VARS, only: damkohler, froude, schmidt
    use TLAB_VARS, only: rbackground, ribackground, bbackground, pbackground, tbackground, epbackground
    use TLAB_VARS, only: buoyancy
    use TLAB_PROCS
    use THERMO_VARS, only: imixture, GRATIO
    use PROFILES
#ifdef USE_MPI
    use TLAB_MPI_VARS
#endif

    implicit none

#include "integers.h"

    TREAL, dimension(g(2)%size, *), intent(INOUT) :: wrk1d

! -----------------------------------------------------------------------
    TINTEGER is, j, ip, nlines, offset

! #######################################################################
    do is = 1,size(qbg)
        if (qbg(is)%relative) qbg(is)%ymean = g(2)%nodes(1) + g(2)%scale*qbg(is)%ymean_rel
    enddo
    if (pbg%relative) pbg%ymean = g(2)%nodes(1) + g(2)%scale*pbg%ymean_rel
    if (rbg%relative) rbg%ymean = g(2)%nodes(1) + g(2)%scale*rbg%ymean_rel
    if (tbg%relative) tbg%ymean = g(2)%nodes(1) + g(2)%scale*tbg%ymean_rel
    if (hbg%relative) hbg%ymean = g(2)%nodes(1) + g(2)%scale*hbg%ymean_rel
    do is = 1,size(sbg)
        if (sbg(is)%relative) sbg(is)%ymean = g(2)%nodes(1) + g(2)%scale*sbg(is)%ymean_rel
    enddo

! #######################################################################
    allocate (pbackground(g(2)%size))
    allocate (rbackground(g(2)%size))
    allocate (ribackground(g(2)%size))
    allocate (bbackground(g(2)%size))
    allocate (tbackground(g(2)%size))
    allocate (epbackground(g(2)%size))

! #######################################################################
! Thermodynamic background profiles
! #######################################################################
    rbackground = C_1_R ! defaults
    ribackground = C_1_R
    pbackground = C_1_R
    tbackground = C_1_R
    epbackground = C_0_R

! Construct given thermodynamic profiles
    do is = 1, inb_scal
        do j = 1, g(2)%size
            wrk1d(j, is) = PROFILES_CALCULATE(sbg(is), g(2)%nodes(j))
        end do
!     wrk1d(:,is) = sbg(is)%reference
    end do

    if (pbg%parameters(5) > C_0_R) then
! Calculate derived thermodynamic profiles
        epbackground = (g(2)%nodes - pbg%ymean)*GRATIO/pbg%parameters(5)

        if (buoyancy%active(2)) then
!        CALL FI_HYDROSTATIC_H_OLD(g(2)%size, g(2)%nodes, wrk1d, epbackground, tbackground, pbackground, wrk1d(1,4))
            call FI_HYDROSTATIC_H(g(2), wrk1d, epbackground, tbackground, pbackground, wrk1d(1, inb_scal_array + 1))
        end if

    end if

    if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= C_0_R) then ! Calculate q_l
        call THERMO_AIRWATER_PH(i1, g(2)%size, i1, wrk1d(1, 2), wrk1d(1, 1), epbackground, pbackground)
    else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
        call THERMO_AIRWATER_LINEAR(i1, g(2)%size, i1, wrk1d, wrk1d(1, inb_scal_array))
    end if

    if (pbg%parameters(5) > C_0_R) then
        call THERMO_ANELASTIC_DENSITY(i1, g(2)%size, i1, wrk1d, epbackground, pbackground, rbackground)
        ribackground = C_1_R/rbackground
    end if

! Calculate buoyancy profile
    if (buoyancy%type == EQNS_EXPLICIT) then
        call THERMO_ANELASTIC_BUOYANCY(i1, g(2)%size, i1, wrk1d, epbackground, pbackground, rbackground, bbackground)
    else
        wrk1d(:, inb_scal_array + 1) = C_0_R
        call FI_BUOYANCY(buoyancy, i1, g(2)%size, i1, wrk1d, bbackground, wrk1d(1, inb_scal_array + 1))
    end if

! #######################################################################
! Add diagnostic fields to reference profile data, if any
! #######################################################################
    do is = inb_scal + 1, inb_scal_array ! Add diagnostic fields, if any
        sbg(is) = sbg(1)
        schmidt(is) = schmidt(1)
    end do
! Buoyancy as next scalar, current value of counter is=inb_scal_array+1
    sbg(is) = sbg(1)
    sbg(is)%mean = (bbackground(1) + bbackground(g(2)%size))/froude
    sbg(is)%delta = ABS(bbackground(1) - bbackground(g(2)%size))/froude
    schmidt(is) = schmidt(1)

    if (imixture == MIXT_TYPE_AIRWATER) then
        is = is + 1
        call THERMO_ANELASTIC_THETA_L(i1, g(2)%size, i1, wrk1d, epbackground, pbackground, wrk1d(1, inb_scal_array + 1))
        sbg(is) = sbg(1)
        sbg(is)%mean = (wrk1d(1, inb_scal_array + 1) + wrk1d(g(2)%size, inb_scal_array + 1))*C_05_R
        sbg(is)%delta = ABS(wrk1d(1, inb_scal_array + 1) - wrk1d(g(2)%size, inb_scal_array + 1))
        schmidt(is) = schmidt(1)
    end if

! #######################################################################
! Anelastic density correction term in burgers operator
! #######################################################################
    if (imode_eqns == DNS_EQNS_ANELASTIC) then
        call TLAB_WRITE_ASCII(lfile, 'Initialize anelastic density correction in burgers operator.')

! Density correction term in the burgers operator along X
        g(1)%anelastic = .true.
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            nlines = ims_size_i(TLAB_MPI_I_PARTIAL)
            offset = nlines*ims_pro_i
        else
#endif
            nlines = jmax*kmax
            offset = 0
#ifdef USE_MPI
        end if
#endif
        allocate (g(1)%rhoinv(nlines))
        do j = 1, nlines
            ip = MOD(offset + j - 1, g(2)%size) + 1
            g(1)%rhoinv(j) = ribackground(ip)
        end do

! Density correction term in the burgers operator along Y; see fdm_initialize
! implemented directly in the tridiagonal system
        ip = 0
        do is = 0, inb_scal ! case 0 for the velocity
            g(2)%lu2d(:, ip + 2) = g(2)%lu2d(:, ip + 2)*ribackground(:)  ! matrix U; 1/diagonal
            g(2)%lu2d(:g(2)%size - 1, ip + 3) = g(2)%lu2d(:, ip + 3)*rbackground(2:) ! matrix U; 1. superdiagonal
            ip = ip + 3
        end do

! Density correction term in the burgers operator along Z
        g(3)%anelastic = .true.
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            nlines = ims_size_k(TLAB_MPI_K_PARTIAL)
            offset = nlines*ims_pro_k
        else
#endif
            nlines = imax*jmax
            offset = 0
#ifdef USE_MPI
        end if
#endif
        allocate (g(3)%rhoinv(nlines))
        do j = 1, nlines
            ip = (offset + j - 1)/imax + 1
            g(3)%rhoinv(j) = ribackground(ip)
        end do

    end if

    return
end subroutine FI_BACKGROUND_INITIALIZE
