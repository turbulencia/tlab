#include "dns_const.h"

!########################################################################
!#
!# Sources from the evolution equations.
!#
!########################################################################
subroutine FI_SOURCES_FLOW(q, s, hq, tmp1)
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: buoyancy, coriolis, subsidence, random
    use TLAB_VARS, only: bbackground, pbackground, rbackground, epbackground
    use TLAB_ARRAYS, only: wrk1d, wrk2d, wrk3d

    implicit none

    real(wp), dimension(isize_field, *), intent(IN) :: q, s
    real(wp), dimension(isize_field, *), intent(OUT) :: hq
    real(wp), dimension(isize_field), intent(INOUT) :: tmp1

    ! -----------------------------------------------------------------------
    integer(wi) ij, iq
    real(wp) dummy, u_geo, w_geo, dtr1, dtr3

    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning

#ifdef USE_BLAS
    integer ILEN
#endif

    ! #######################################################################
#ifdef USE_BLAS
    ILEN = isize_field
#endif

    ! -----------------------------------------------------------------------
    ! Coriolis. Remember that coriolis%vector already contains the Rossby #.
    ! -----------------------------------------------------------------------
    select case (coriolis%type)
    case (EQNS_EXPLICIT)
        do ij = 1, isize_field
            hq(ij, 1) = hq(ij, 1) + coriolis%vector(3)*q(ij, 2) - coriolis%vector(2)*q(ij, 3)
            hq(ij, 2) = hq(ij, 2) + coriolis%vector(1)*q(ij, 3) - coriolis%vector(3)*q(ij, 1)
            hq(ij, 3) = hq(ij, 3) + coriolis%vector(2)*q(ij, 1) - coriolis%vector(1)*q(ij, 2)
        end do

    case (EQNS_COR_NORMALIZED) ! So far, rotation only in the Oy direction.
        u_geo = cos(coriolis%parameters(1))*coriolis%parameters(2)
        w_geo = -sin(coriolis%parameters(1))*coriolis%parameters(2)

!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
        call DNS_OMP_PARTITION(isize_field, srt, end, siz)

        dummy = coriolis%vector(2)
        dtr3 = 0.0_wp; dtr1 = 0.0_wp
        do ij = srt, end
            hq(ij, 1) = hq(ij, 1) + dummy*(w_geo - q(ij, 3))
            hq(ij, 3) = hq(ij, 3) + dummy*(q(ij, 1) - u_geo)
        end do
!$omp end parallel

    end select

    ! -----------------------------------------------------------------------
    do iq = 1, 3

        ! -----------------------------------------------------------------------
        ! Buoyancy. Remember that buoyancy%vector already contains the Froude #.
        ! -----------------------------------------------------------------------
        if (buoyancy%active(iq)) then

            if (buoyancy%type == EQNS_EXPLICIT) then
                call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, epbackground, pbackground, rbackground, wrk3d)

            else
                if (iq == 2) then
                    call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, bbackground)
                else
                    wrk1d(:, 1) = 0.0_wp
                    call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, wrk1d)
                end if
            end if

#ifdef USE_BLAS
!$omp parallel default( shared ) &
!$omp private( ilen, dummy, srt,end,siz)
#else
!$omp parallel default( shared ) &
!$omp private( ij,   dummy, srt,end,siz )
#endif
            call DNS_OMP_PARTITION(isize_field, srt, end, siz)

            dummy = buoyancy%vector(iq)
#ifdef USE_BLAS
            ILEN = siz
            call DAXPY(ILEN, dummy, wrk3d(srt), 1, hq(srt, iq), 1)
#else
            do ij = srt, end
                hq(ij, iq) = hq(ij, iq) + dummy*wrk3d(ij)
            end do
#endif
!$omp end parallel

        end if

        ! -----------------------------------------------------------------------
        ! Subsidence
        ! -----------------------------------------------------------------------
        if (subsidence%active(iq)) then
            call FI_SUBSIDENCE(subsidence, imax, jmax, kmax, q(1, iq), tmp1, wrk1d, wrk2d, wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
            call DNS_OMP_PARTITION(isize_field, srt, end, siz)

            do ij = srt, end
                hq(ij, iq) = hq(ij, iq) + tmp1(ij)
            end do
!$omp end parallel

        end if

        ! -----------------------------------------------------------------------
        ! Random forcing
        ! -----------------------------------------------------------------------
        if (random%active(iq)) then
            call FI_RANDOM(random, imax, jmax, kmax, hq(1, iq), tmp1)
        end if
    end do

    return
end subroutine FI_SOURCES_FLOW

! #######################################################################
! #######################################################################
subroutine FI_SOURCES_SCAL(s, hs, tmp1, tmp2)
    use TLAB_CONSTANTS, only: wp, wi
    use TLAB_VARS, only: imax, jmax, kmax, inb_scal, isize_field
    use TLAB_VARS, only: imode_eqns
    use TLAB_VARS, only: g
    use TLAB_VARS, only: radiation, transport, chemistry, subsidence
    use TLAB_VARS, only: rbackground, ribackground
    use TLAB_ARRAYS, only: wrk1d, wrk2d, wrk3d

    implicit none

    real(wp), dimension(isize_field, *), intent(IN) :: s
    real(wp), dimension(isize_field, *), intent(OUT) :: hs
    real(wp), dimension(isize_field), intent(INOUT) :: tmp1, tmp2

    ! -----------------------------------------------------------------------
    integer(wi) ij, is, flag_grad

    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning

#ifdef USE_BLAS
    integer ILEN
#endif

    ! #######################################################################
#ifdef USE_BLAS
    ILEN = isize_field
#endif

    do is = 1, inb_scal ! Start loop over the N scalars

        ! -----------------------------------------------------------------------
        ! Radiation
        ! -----------------------------------------------------------------------
        if (radiation%active(is)) then
            if (imode_eqns == DNS_EQNS_ANELASTIC) then
                call THERMO_ANELASTIC_WEIGHT_OUTPLACE(imax, jmax, kmax, rbackground, s(1, radiation%scalar(is)), tmp2)
                call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), tmp2, tmp1, wrk1d, wrk3d)
                call THERMO_ANELASTIC_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(1, is))

            else
                call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), s(1, radiation%scalar(is)), tmp1, wrk1d, wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
                call DNS_OMP_PARTITION(isize_field, srt, end, siz)

                do ij = srt, end
                    hs(ij, is) = hs(ij, is) + tmp1(ij)
                end do
!$omp end parallel

            end if

        end if

        ! -----------------------------------------------------------------------
        ! Transport, such as settling
        ! array tmp2 should not be used inside the loop on is
        ! -----------------------------------------------------------------------
        if (transport%active(is)) then
            if (is == 1) then; flag_grad = 1; 
            else; flag_grad = 0
            end if
            call FI_TRANSPORT(transport, flag_grad, imax, jmax, kmax, is, s, tmp1, tmp2, wrk2d, wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
            call DNS_OMP_PARTITION(isize_field, srt, end, siz)

            do ij = srt, end
                hs(ij, is) = hs(ij, is) + tmp1(ij)
            end do
!$omp end parallel

        end if

        ! -----------------------------------------------------------------------
        ! Chemistry
        ! -----------------------------------------------------------------------
        if (chemistry%active(is)) then
            call FI_CHEM(chemistry, imax, jmax, kmax, is, s, tmp1)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
            call DNS_OMP_PARTITION(isize_field, srt, end, siz)

            do ij = srt, end
                hs(ij, is) = hs(ij, is) + tmp1(ij)
            end do
!$omp end parallel

        end if

        ! -----------------------------------------------------------------------
        ! Subsidence
        ! -----------------------------------------------------------------------
        if (subsidence%active(is)) then
            call FI_SUBSIDENCE(subsidence, imax, jmax, kmax, s(1, is), tmp1, wrk1d, wrk2d, wrk3d)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
            call DNS_OMP_PARTITION(isize_field, srt, end, siz)

            do ij = srt, end
                hs(ij, is) = hs(ij, is) + tmp1(ij)
            end do
!$omp end parallel

        end if

    end do

    return
end subroutine FI_SOURCES_SCAL
