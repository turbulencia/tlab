#include "dns_const.h"

!########################################################################
!# Sources (processes) in the evolution equations.
!########################################################################
module TLab_Sources
    use TLab_Constants, only: wp, wi, small_wp
    use TLab_Types, only: term_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal, inb_scal_array
    use NavierStokes, only: nse_eqns
    use FDM, only: g
    use TLAB_VARS, only: coriolis
    use TLab_OpenMP
    use THERMO_ANELASTIC
    use Gravity, only: buoyancy, bbackground, Gravity_Buoyancy
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
    public :: FI_CORIOLIS

contains
! #######################################################################
! #######################################################################
    subroutine TLab_Sources_Flow(q, s, hq, tmp1)
        use TLAB_VARS, only: rtime
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

        call FI_CORIOLIS(coriolis, imax, jmax, kmax, q, hq)

        ! -----------------------------------------------------------------------
        do iq = 1, 3

            ! -----------------------------------------------------------------------
            ! Buoyancy. Remember that buoyancy%vector already contains the Froude #.
            ! -----------------------------------------------------------------------
            if (buoyancy%active(iq)) then

                if (buoyancy%type == EQNS_EXPLICIT) then
                    call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, tmp1)

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
                call Radiation_Infrared_Y(infraredProps, imax, jmax, kmax, g(2), s, tmp1, tmp2, tmp3, tmp4)

                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call THERMO_ANELASTIC_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(:, is))
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
                    call THERMO_ANELASTIC_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(:, is))
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

!########################################################################
!# Calculating the coriolis term
!########################################################################
    subroutine FI_CORIOLIS(coriolis, nx, ny, nz, u, r)
        type(term_dt), intent(in) :: coriolis
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz, *)
        real(wp), intent(out) :: r(nx*ny*nz, *)

        ! -----------------------------------------------------------------------
        integer(wi) ii, field_sz
        real(wp) dummy, dtr3, dtr1, geo_u, geo_w
        field_sz = nx*ny*nz
        ! -----------------------------------------------------------------------
        ! Coriolis. Remember that coriolis%vector already contains the Rossby #.
        ! -----------------------------------------------------------------------
        select case (coriolis%type)
        case (EQNS_EXPLICIT)
            do ii = 1, field_sz
                r(ii, 1) = r(ii, 1) + coriolis%vector(3)*u(ii, 2) - coriolis%vector(2)*u(ii, 3)
                r(ii, 2) = r(ii, 2) + coriolis%vector(1)*u(ii, 3) - coriolis%vector(3)*u(ii, 1)
                r(ii, 3) = r(ii, 3) + coriolis%vector(2)*u(ii, 1) - coriolis%vector(1)*u(ii, 2)
            end do

        case (EQNS_COR_NORMALIZED)
            geo_u = cos(coriolis%parameters(1))*coriolis%parameters(2)
            geo_w = -sin(coriolis%parameters(1))*coriolis%parameters(2)

!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
            call TLab_OMP_PARTITION(field_sz, srt, end, siz)

            dummy = coriolis%vector(2)
            dtr3 = 0.0_wp; dtr1 = 0.0_wp
            do ii = srt, end
                r(ii, 1) = r(ii, 1) + dummy*(geo_w - u(ii, 3))
                r(ii, 3) = r(ii, 3) + dummy*(u(ii, 1) - geo_u)
            end do
!$omp end parallel
        end select

    end subroutine FI_CORIOLIS

end module TLab_Sources
