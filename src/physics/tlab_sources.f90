#include "dns_const.h"

!########################################################################
!# Sources (processes) in the evolution equations.
!########################################################################
module TLab_Sources
    use TLab_Constants, only: wp, wi, small_wp
    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_scal, inb_scal_array
    use NavierStokes, only: nse_eqns
    use FDM, only: g
    use FDM_Integral, only: fdm_Int0
    use TLab_OpenMP
    use Thermo_Anelastic
    use Gravity, only: buoyancy, bbackground, Gravity_Buoyancy
    use Rotation, only: coriolis, Rotation_Coriolis
    use Radiation
    use Microphysics
    use Chemistry
    use SpecialForcing
    use LargeScaleForcing
    implicit none
    private

    integer(wi) ij, iq, is
    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning
#ifdef USE_BLAS
    integer ILEN
#endif

    public :: TLab_Sources_Flow
    public :: TLab_Sources_Scal

contains
! #######################################################################
! #######################################################################
    subroutine TLab_Sources_Flow(q, s, hq, tmp1)
        use TLab_Time, only: rtime
        real(wp), intent(in) :: q(isize_field, *), s(isize_field, *)
        real(wp), intent(out) :: hq(isize_field, *)
        real(wp), intent(inout) :: tmp1(isize_field)

        ! -----------------------------------------------------------------------
        real(wp) dummy

        ! #######################################################################
#ifdef USE_BLAS
        ILEN = isize_field
#endif

        ! -----------------------------------------------------------------------
        ! Coriolis. Remember that coriolis%vector already contains the Rossby #.
        ! -----------------------------------------------------------------------

        call Rotation_Coriolis(coriolis, imax, jmax, kmax, q, hq)

        ! -----------------------------------------------------------------------
        do iq = 1, 3

            ! -----------------------------------------------------------------------
            ! Buoyancy. Remember that buoyancy%vector already contains the Froude #.
            ! -----------------------------------------------------------------------
            if (buoyancy%active(iq)) then

                if (buoyancy%type == EQNS_EXPLICIT) then
                    call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, tmp1)

                else
                    call Gravity_Buoyancy(buoyancy, imax, jmax, kmax, s, tmp1, bbackground)

                end if

#ifdef USE_BLAS
!$omp parallel default( shared ) &
!$omp private( ilen, dummy, srt,end,siz)
#else
!$omp parallel default( shared ) &
!$omp private( ij,   dummy, srt,end,siz )
#endif
                call TLab_OMP_PARTITION(isize_field, srt, end, siz)

                dummy = buoyancy%vector(iq)
#ifdef USE_BLAS
                ILEN = siz
                call DAXPY(ILEN, dummy, tmp1(srt), 1, hq(srt, iq), 1)
#else
                do ij = srt, end
                    hq(ij, iq) = hq(ij, iq) + dummy*tmp1(ij)
                end do
#endif
!$omp end parallel

            end if

            ! -----------------------------------------------------------------------
            ! Subsidence
            ! -----------------------------------------------------------------------
            if (subsidenceProps%active(iq)) then
                call LargeScaleForcing_Subsidence(subsidenceProps, imax, jmax, kmax, q(:, iq), tmp1)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
                call TLab_OMP_PARTITION(isize_field, srt, end, siz)

                do ij = srt, end
                    hq(ij, iq) = hq(ij, iq) + tmp1(ij)
                end do
!$omp end parallel

            end if

            ! -----------------------------------------------------------------------
            ! special forcing
            ! -----------------------------------------------------------------------
            if (forcingProps%active(iq)) then
                call SpecialForcing_Source(forcingProps, imax, jmax, kmax, iq, rtime, q(:,iq), hq(:, iq), tmp1)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
                call TLab_OMP_PARTITION(isize_field, srt, end, siz)

                do ij = srt, end
                    ! hq(ij, iq) = hq(ij, iq) + tmp1(ij)*forcingProps%vector(iq)
                    hq(ij, iq) = hq(ij, iq) + tmp1(ij)
                end do
!$omp end parallel

            end if
        end do

        return
    end subroutine TLab_Sources_Flow

    ! #######################################################################
    ! #######################################################################
    subroutine TLab_Sources_Scal(s, hs, tmp1, tmp2, tmp3, tmp4)
        real(wp), intent(in) :: s(isize_field, *)
        real(wp), intent(out) :: hs(isize_field, *)
        real(wp), intent(inout) :: tmp1(isize_field), tmp2(isize_field), tmp3(isize_field), tmp4(isize_field)

        ! -----------------------------------------------------------------------

        ! #######################################################################
#ifdef USE_BLAS
        ILEN = isize_field
#endif

        do is = 1, inb_scal ! Start loop over the N scalars

            ! -----------------------------------------------------------------------
            ! Radiation
            ! -----------------------------------------------------------------------
            if (infraredProps%active(is)) then
                call Radiation_Infrared_Y(infraredProps, imax, jmax, kmax, fdm_Int0, s, tmp1, tmp2, tmp3, tmp4)

                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(:, is))
                else
!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
                    call TLab_OMP_PARTITION(isize_field, srt, end, siz)

                    do ij = srt, end
                        hs(ij, is) = hs(ij, is) + tmp1(ij)
                    end do
!$omp end parallel
                end if

            end if

            ! -----------------------------------------------------------------------
            ! Microphysics
            ! -----------------------------------------------------------------------
            if (sedimentationProps%active(is)) then
                call Microphysics_Sedimentation(sedimentationProps, imax, jmax, kmax, is, g(2), s, tmp1, tmp2)

                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(:, is))
                else
!$omp parallel default( shared ) &forcingProps%vector
!$omp private( ij, srt,end,siz )
                    call TLab_OMP_PARTITION(isize_field, srt, end, siz)

                    do ij = srt, end
                        hs(ij, is) = hs(ij, is) + tmp1(ij)
                    end do
!$omp end parallel
                end if

            end if

            ! -----------------------------------------------------------------------
            ! Chemistry
            ! -----------------------------------------------------------------------
            if (chemistryProps%active(is)) then
                call Chemistry_Source(chemistryProps, imax, jmax, kmax, is, s, tmp1)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
                call TLab_OMP_PARTITION(isize_field, srt, end, siz)

                do ij = srt, end
                    hs(ij, is) = hs(ij, is) + tmp1(ij)
                end do
!$omp end parallel

            end if

            ! -----------------------------------------------------------------------
            ! Subsidence
            ! -----------------------------------------------------------------------
            if (subsidenceProps%active(is)) then
                call LargeScaleForcing_Subsidence(subsidenceProps, imax, jmax, kmax, s(:, is), tmp1)

!$omp parallel default( shared ) &
!$omp private( ij, srt,end,siz )
                call TLab_OMP_PARTITION(isize_field, srt, end, siz)

                do ij = srt, end
                    hs(ij, is) = hs(ij, is) + tmp1(ij)
                end do
!$omp end parallel

            end if

        end do

        return
    end subroutine TLab_Sources_Scal

end module TLab_Sources
