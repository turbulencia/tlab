#include "dns_const.h"

!########################################################################
!# Sources (processes) in the evolution equations.
!########################################################################
module FI_SOURCES
    use TLab_Constants, only: wp, wi, small_wp
    use TLab_Types, only: term_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal, inb_scal_array
    use TLAB_VARS, only: imode_eqns
    use FDM, only: g
    use TLAB_VARS, only: buoyancy, coriolis, subsidence
    use TLab_OpenMP
    use THERMO_ANELASTIC
    use Radiation
    use Microphysics
    use Chemistry
    use SpecialForcing
    implicit none
    private

    integer(wi) ij, iq, is
    integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning
#ifdef USE_BLAS
    integer ILEN
#endif

    public :: FI_SOURCES_FLOW
    public :: FI_SOURCES_SCAL
    public :: FI_BUOYANCY, FI_BUOYANCY_SOURCE
    public :: FI_CORIOLIS

    real(wp), allocatable, public :: bbackground(:)     ! Buoyancy

contains
! #######################################################################
! #######################################################################
    subroutine FI_SOURCES_FLOW(q, s, hq, tmp1)
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
                    call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, tmp1, bbackground)

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
            if (subsidence%active(iq)) then
                call FI_SUBSIDENCE(subsidence, imax, jmax, kmax, q(:, iq), tmp1)

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
    end subroutine FI_SOURCES_FLOW

    ! #######################################################################
    ! #######################################################################
    subroutine FI_SOURCES_SCAL(s, hs, tmp1, tmp2, tmp3, tmp4)
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

                if (imode_eqns == DNS_EQNS_ANELASTIC) then
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

                if (imode_eqns == DNS_EQNS_ANELASTIC) then
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
            if (subsidence%active(is)) then
                call FI_SUBSIDENCE(subsidence, imax, jmax, kmax, s(:, is), tmp1)

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
    end subroutine FI_SOURCES_SCAL

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

!########################################################################
!# Determine the buoyancy term (density difference rho -rho_0
!# when it is a function of a scalar
!########################################################################
    subroutine FI_BUOYANCY(buoyancy, nx, ny, nz, s, b, ref)
        type(term_dt), intent(in) :: buoyancy
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx, ny, nz, inb_scal_array)
        real(wp), intent(out) :: b(nx, ny, nz)
        real(wp), intent(inout) :: ref(ny)         ! reference profile

        ! -----------------------------------------------------------------------
        integer(wi) j, k
        real(wp) c0_loc, c1_loc, c2_loc, c3_loc, dummy

        ! #######################################################################
        select case (buoyancy%type)

        case (EQNS_BOD_HOMOGENEOUS)
            b = buoyancy%parameters(1)

        case (EQNS_BOD_LINEAR)
            c1_loc = buoyancy%parameters(1); 
            c2_loc = buoyancy%parameters(2); 
            c3_loc = buoyancy%parameters(3)                     ! proportionality factors
            c0_loc = buoyancy%parameters(inb_scal_array + 1)    ! independent term

            if (buoyancy%scalar(1) == 1) then
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) - dummy
                    end do
                end do

            else if (buoyancy%scalar(1) == 2) then
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) + c2_loc*s(1:nx, j, k, 2) - dummy
                    end do
                end do

            else if (buoyancy%scalar(1) == 3) then
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) + c2_loc*s(1:nx, j, k, 2) + c3_loc*s(1:nx, j, k, 3) - dummy
                    end do
                end do

            else ! general
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j)
                        b(1:nx, j, k) = c0_loc - dummy
                    end do
                end do

                do is = 1, buoyancy%scalar(1)
                    if (abs(buoyancy%parameters(is)) > small_wp) b(:, :, :) = b(:, :, :) + buoyancy%parameters(is)*s(:, :, :, is)
                end do

            end if

        case (EQNS_BOD_BILINEAR)
            c0_loc = buoyancy%parameters(1); 
            c1_loc = buoyancy%parameters(2); 
            c2_loc = buoyancy%parameters(3)

            do k = 1, nz
                do j = 1, ny
                    dummy = ref(j)
                    b(1:nx, j, k) = c0_loc*s(1:nx, j, k, 1) + c1_loc*s(1:nx, j, k, 2) + c2_loc*s(1:nx, j, k, 1)*s(1:nx, j, k, 2) - dummy
                end do
            end do

        case (EQNS_BOD_QUADRATIC)
            c0_loc = -buoyancy%parameters(1)/(buoyancy%parameters(2)/2.0_wp)**2
            c1_loc = buoyancy%parameters(2)

            do k = 1, nz
                do j = 1, ny
                    dummy = ref(j)
                    b(1:nx, j, k) = c0_loc*s(1:nx, j, k, 1)*(s(1:nx, j, k, 1) - c1_loc) - dummy
                end do
            end do

        case (EQNS_BOD_NORMALIZEDMEAN)
            c1_loc = buoyancy%parameters(1)

            do k = 1, nz
                do j = 1, ny
                    dummy = 1.0_wp/bbackground(j)
                    b(1:nx, j, k) = c1_loc*(dummy*s(1:nx, j, k, 1) - 1.0_wp)
                end do
            end do

        case (EQNS_BOD_SUBTRACTMEAN)
            c1_loc = buoyancy%parameters(1)

            do k = 1, nz
                do j = 1, ny
                    dummy = bbackground(j)
                    b(1:nx, j, k) = c1_loc*(s(1:nx, j, k, 1) - dummy)
                end do
            end do

        case default
            b = 0.0_wp

        end select

        return
    end subroutine FI_BUOYANCY

    !########################################################################
    !########################################################################
    subroutine FI_BUOYANCY_SOURCE(buoyancy, nx, ny, nz, s, gradient, source)
        type(term_dt), intent(in) :: buoyancy
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx, ny, nz, inb_scal_array)
        real(wp), intent(inout) :: gradient(nx, ny, nz)         ! gradient magnitude
        real(wp), intent(out) :: source(nx, ny, nz)

        ! -----------------------------------------------------------------------
        real(wp) c0_loc

        ! #######################################################################
        select case (buoyancy%type)

        case (EQNS_BOD_HOMOGENEOUS)
            source = 0.0_wp

        case (EQNS_BOD_LINEAR)
            source = 0.0_wp

        case (EQNS_BOD_BILINEAR)
            source = 0.0_wp

        case (EQNS_BOD_QUADRATIC)
            c0_loc = -buoyancy%parameters(1)/(buoyancy%parameters(2)/2.0_wp)**2; c0_loc = c0_loc*2.0_wp

            source = c0_loc*gradient

        end select

        return
    end subroutine FI_BUOYANCY_SOURCE

! #######################################################################
! #######################################################################
    subroutine FI_SUBSIDENCE(subsidence, nx, ny, nz, a, source)
        use TLab_Arrays, only: wrk1d
        use OPR_PARTIAL, only: OPR_PARTIAL_Y
        use Averages, only: AVG1V2D_V
        type(term_dt), intent(in) :: subsidence
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: a(nx, ny, nz)
        real(wp), intent(out) :: source(nx, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) bcs(2, 2), j

        !########################################################################
        bcs = 0

        select case (subsidence%type)

        case (EQNS_SUB_CONSTANT_LOCAL)
            wrk1d(1:ny, 1) = g(2)%nodes(1:ny)*subsidence%parameters(1)      ! subsidence velocity

            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), a, source)

            do j = 1, ny
                source(:, j, :) = source(:, j, :)*wrk1d(j, 1)
            end do

        case (EQNS_SUB_CONSTANT_GLOBAL)
            wrk1d(1:ny, 1) = g(2)%nodes(1:ny)*subsidence%parameters(1)      ! subsidence velocity

            call AVG1V2D_V(nx, ny, nz, 1, a, wrk1d(:, 2), wrk1d(1, 3))     ! Calculate averaged scalar into wrk1d(:,2)
            call OPR_PARTIAL_Y(OPR_P1, 1, ny, 1, bcs, g(2), wrk1d(1, 2), wrk1d(1, 3))

            do j = 1, ny
                source(:, j, :) = wrk1d(j, 3)*wrk1d(j, 1)
            end do

        end select

        return
    end subroutine FI_SUBSIDENCE

end module FI_SOURCES
