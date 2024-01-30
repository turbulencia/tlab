#include "dns_const.h"

!########################################################################
!# Sources (processes) in the evolution equations.
!########################################################################
module FI_SOURCES
    use TLAB_CONSTANTS, only: wp, wi, small_wp
    use TLAB_TYPES, only: term_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal, inb_scal_array
    use TLAB_VARS, only: imode_eqns
    use TLAB_VARS, only: g
    use TLAB_VARS, only: buoyancy, coriolis, subsidence, random
    use TLAB_VARS, only: radiation, transport, chemistry, subsidence
    use TLAB_VARS, only: bbackground
    use THERMO_ANELASTIC
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
    public :: FI_TRANSPORT, FI_TRANSPORT_FLUX
    public :: FI_FORCING_0, FI_FORCING_1

contains
! #######################################################################
! #######################################################################
    subroutine FI_SOURCES_FLOW(q, s, hq, tmp1)
        use TLAB_ARRAYS, only: wrk1d
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
                    call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, epbackground, pbackground, rbackground, tmp1)

                else
                    if (iq == 2) then
                        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, tmp1, bbackground)
                    else
                        wrk1d(:, 1) = 0.0_wp
                        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, tmp1, wrk1d)
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
        real(wp), intent(in) :: s(isize_field, *)
        real(wp), intent(out) :: hs(isize_field, *)
        real(wp), intent(inout) :: tmp1(isize_field), tmp2(isize_field)

        ! -----------------------------------------------------------------------
        integer flag_grad

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
                    call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), tmp2, tmp1)
                    call THERMO_ANELASTIC_WEIGHT_ADD(imax, jmax, kmax, ribackground, tmp1, hs(1, is))

                else
                    call OPR_RADIATION(radiation, imax, jmax, kmax, g(2), s(1, radiation%scalar(is)), tmp1)

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
                call FI_TRANSPORT(transport, flag_grad, imax, jmax, kmax, is, s, tmp1, tmp2)

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
                call FI_SUBSIDENCE(subsidence, imax, jmax, kmax, s(:, is), tmp1)

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

!########################################################################
!# Calculating the coriolis term
!########################################################################
    subroutine FI_CORIOLIS(coriolis, nx, ny, nz, u, r)
        type(term_dt), intent(in) :: coriolis
        integer(wi), intent(in)   :: nx, ny, nz
        real(wp), intent(in)      :: u(nx*ny*nz,*)
        real(wp), intent(out)     :: r(nx*ny*nz,*)

        ! -----------------------------------------------------------------------
        integer(wi) ii, field_sz
        real(wp) dummy, dtr3, dtr1, geo_u, geo_w
        field_sz = nx*ny*nz
        ! -----------------------------------------------------------------------
        ! Coriolis. Remember that coriolis%vector already contains the Rossby #.
        ! -----------------------------------------------------------------------
        select case (coriolis%type)
        case (EQNS_EXPLICIT)
            do  ii = 1, field_sz
                r(ii, 1) = r(ii, 1) + coriolis%vector(3)*u(ii, 2) - coriolis%vector(2)*u(ii, 3)
                r(ii, 2) = r(ii, 2) + coriolis%vector(1)*u(ii, 3) - coriolis%vector(3)*u(ii, 1)
                r(ii, 3) = r(ii, 3) + coriolis%vector(2)*u(ii, 1) - coriolis%vector(1)*u(ii, 2)
            end do

        case (EQNS_COR_NORMALIZED)
            geo_u = cos(coriolis%parameters(1))*coriolis%parameters(2)
            geo_w = -sin(coriolis%parameters(1))*coriolis%parameters(2)

!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
            call DNS_OMP_PARTITION(field_sz, srt, end, siz)

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
!# Determine the buoyancy term (density difference \rho-\rho_0)
!# when it is a function of a scalar
!########################################################################
    subroutine FI_BUOYANCY(buoyancy, nx, ny, nz, s, b, ref)
        type(term_dt), intent(in) :: buoyancy
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx, ny, nz, inb_scal_array)
        real(wp), intent(out) :: b(nx, ny, nz)
        real(wp), intent(in) :: ref(ny)         ! reference profile

        ! -----------------------------------------------------------------------
        integer(wi) j, k
        real(wp) c0_loc, c1_loc, c2_loc, c3_loc, dummy

        ! #######################################################################
        select case (buoyancy%type)

        case (EQNS_BOD_HOMOGENEOUS)
            b = buoyancy%parameters(1)

        case (EQNS_BOD_LINEAR)
            c1_loc = buoyancy%parameters(1); c2_loc = buoyancy%parameters(2); c3_loc = buoyancy%parameters(3) ! proportionality factors
            c0_loc = buoyancy%parameters(inb_scal_array + 1)                                                    ! independent term

            if (inb_scal_array == 1 .or. &
                (inb_scal_array == 2 .and. abs(buoyancy%parameters(2)) < small_wp)) then ! avoid mem call for one common case is=2
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) - dummy
                    end do
                end do

            else if (inb_scal_array == 2) then
                do k = 1, nz
                    do j = 1, ny
                        dummy = ref(j) - c0_loc
                        b(1:nx, j, k) = c1_loc*s(1:nx, j, k, 1) + c2_loc*s(1:nx, j, k, 2) - dummy
                    end do
                end do

            else if (inb_scal_array == 3) then
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

                do is = 1, inb_scal_array
                    if (abs(buoyancy%parameters(is)) > small_wp) b(:, :, :) = b(:, :, :) + buoyancy%parameters(is)*s(:, :, :, is)
                end do

            end if

        case (EQNS_BOD_BILINEAR)
            c0_loc = buoyancy%parameters(1); c1_loc = buoyancy%parameters(2); c2_loc = buoyancy%parameters(3)

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
        use TLAB_ARRAYS, only: wrk1d
        use OPR_PARTIAL, only: OPR_PARTIAL_Y
        use AVGS, only: AVG1V2D_V
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
            wrk1d(1:ny, 1) = g(2)%nodes(1:ny)*subsidence%parameters(1)      ! sedimentation velocity

            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), a, source)

            do j = 1, ny
                source(:, j, :) = source(:, j, :)*wrk1d(j, 1)
            end do

        case (EQNS_SUB_CONSTANT_GLOBAL)
            wrk1d(1:ny, 1) = g(2)%nodes(1:ny)*subsidence%parameters(1)      ! sedimentation velocity

            call AVG1V2D_V(nx, ny, nz, 1, a, wrk1d(:, 2), wrk1d(1, 3))     ! Calculate averaged scalar into wrk1d(:,2)
            call OPR_PARTIAL_Y(OPR_P1, 1, ny, 1, bcs, g(2), wrk1d(1, 2), wrk1d(1, 3))

            do j = 1, ny
                source(:, j, :) = wrk1d(j, 3)*wrk1d(j, 1)
            end do

        end select

        return
    end subroutine FI_SUBSIDENCE

    !########################################################################
    !# Calculate the transport terms due to settling in an airwater mixture.
    !########################################################################
    subroutine FI_TRANSPORT(transport, flag_grad, nx, ny, nz, is, s, trans, tmp)
        use OPR_PARTIAL, only: OPR_PARTIAL_Y
        type(term_dt), intent(IN) :: transport
        integer(wi), intent(IN) :: nx, ny, nz, flag_grad
        integer(wi), intent(IN) :: is
        real(wp), intent(IN) :: s(nx*ny*nz, *)
        real(wp), intent(OUT) :: trans(nx*ny*nz, 1) ! Transport component. It could have eventually three directions
        real(wp), intent(INOUT) :: tmp(nx*ny*nz, 1)   ! To avoid re-calculations when repetedly calling this routine

! -----------------------------------------------------------------------
        real(wp) dummy, exponent
        integer(wi) is_ref, bcs(2, 2)

!########################################################################
        bcs = 0

        exponent = transport%auxiliar(1)
        is_ref = transport%scalar(1)

        if (transport%type == EQNS_TRANS_AIRWATERSIMPLIFIED) then
            if (flag_grad == 1) then
                call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), s(1, is_ref), tmp)
                if (exponent > 0.0_wp) tmp(:, 1) = tmp(:, 1)*(s(:, is_ref)**exponent)
            end if

            dummy = transport%parameters(is)*(1.0_wp + exponent)
            trans(:, 1) = dummy*tmp(:, 1)

        elseif (transport%type == EQNS_TRANS_AIRWATER) then
            dummy = 1.0_wp + exponent

            select case (is)
            case (2, 3)         ! q_t, q_l
                if (exponent > 0.0_wp) then
                    tmp(:, 1) = (transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*(s(:, is_ref)**dummy)
                else
                    tmp(:, 1) = (transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*s(:, is_ref)
                end if

            case default        ! energy variables
                call THERMO_ANELASTIC_STATIC_L(nx, ny, nz, s, epbackground, tmp(:, 1))
                if (exponent > 0.0_wp) then
                    tmp(:, 1) = transport%parameters(is)*tmp(:, 1)*(s(:, is_ref)**dummy)
                else
                    tmp(:, 1) = transport%parameters(is)*tmp(:, 1)*s(:, is_ref)
                end if

            end select

            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), tmp(1, 1), trans(1, 1))

        end if

        return
    end subroutine FI_TRANSPORT

!########################################################################
!########################################################################
    subroutine FI_TRANSPORT_FLUX(transport, nx, ny, nz, is, s, trans)
        type(term_dt), intent(in) :: transport
        integer(wi), intent(in) :: nx, ny, nz
        integer(wi), intent(in) :: is
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: trans(nx*ny*nz, 1) ! Transport component. It could have eventually three directions

! -----------------------------------------------------------------------
        real(wp) dummy, exponent
        integer(wi) is_ref

!########################################################################
        exponent = transport%auxiliar(1)
        is_ref = transport%scalar(1)

        if (transport%type == EQNS_TRANS_AIRWATERSIMPLIFIED) then
            dummy = 1.0_wp + exponent

            trans(:, 1) = -transport%parameters(is)*(s(:, is_ref)**dummy)

        elseif (transport%type == EQNS_TRANS_AIRWATER) then
            dummy = 1.0_wp + exponent

            select case (is)
            case (2, 3)         ! q_t, q_l
                if (exponent > 0.0_wp) then
                 trans(:, 1) = -(transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*(s(:, is_ref)**dummy)
                else
                    trans(:, 1) = -(transport%parameters(is) - transport%parameters(inb_scal_array + 3)*s(:, is))*s(:, is_ref)
                end if

            case default        ! energy variables
                call THERMO_ANELASTIC_STATIC_L(nx, ny, nz, s, epbackground, trans(:, 1))
                if (exponent > 0.0_wp) then
                    trans(:, 1) = -transport%parameters(is)*trans(:, 1)*(s(:, is_ref)**dummy)
                else
                    trans(:, 1) = -transport%parameters(is)*trans(:, 1)*s(:, is_ref)
                end if

            end select

        end if

        return
    end subroutine FI_TRANSPORT_FLUX

    ! #######################################################################
    ! #######################################################################
    subroutine FI_CHEM(chemistry, nx, ny, nz, is, s, source)
        use TLAB_TYPES, only: profiles_dt
        use TLAB_VARS, only: sbg, damkohler
        use PROFILES

        type(term_dt), intent(IN) :: chemistry
        integer(wi), intent(IN) :: nx, ny, nz, is
        real(wp), intent(IN) :: s(nx, ny, nz, inb_scal)
        real(wp), intent(OUT) :: source(nx, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy, dummy2
        type(profiles_dt) prof_loc

        !########################################################################
        select case (chemistry%type)

        case (EQNS_CHEM_LAYEREDRELAXATION)
            prof_loc%type = PROFILE_TANH
            prof_loc%ymean = sbg(is)%ymean
            prof_loc%thick = -chemistry%parameters(3)*0.5_wp
            prof_loc%mean = 0.5_wp
            prof_loc%delta = 1.0_wp
            prof_loc%lslope = 0.0_wp
            prof_loc%uslope = 0.0_wp

            do j = 1, ny
                source(:, j, :) = PROFILES_CALCULATE(prof_loc, g(2)%nodes(j) - chemistry%parameters(2)) ! strength constant
            end do

            dummy = -damkohler(is)/chemistry%parameters(1)
            source = dummy*source*s(:, :, :, is)

        case (EQNS_CHEM_QUADRATIC)
            dummy = damkohler(is)*chemistry%parameters(is)
            source = dummy*s(:, :, :, 2)*s(:, :, :, 3)

        case (EQNS_CHEM_QUADRATIC3)
            dummy = damkohler(is)*chemistry%parameters(is)

            if (is >= 1 .and. is <= 3) then
                source = dummy*s(:, :, :, 2)*s(:, :, :, 3)
            else if (is >= 4 .and. is <= 6) then
                source = dummy*s(:, :, :, 4)*s(:, :, :, 5)
            else if (is >= 7 .and. is <= 9) then
                source = dummy*s(:, :, :, 7)*s(:, :, :, 8)
            end if

        case (EQNS_CHEM_OZONE)
            dummy = damkohler(is)
            if (is == 4) dummy = -dummy

            source = -chemistry%parameters(1)/(1.0_wp + chemistry%parameters(2)*s(:, :, :, 1))
            source = exp(source)

            if (is == 4) then
                dummy2 = 1.0_wp + chemistry%parameters(3)
                source = dummy*(dummy2*s(:, :, :, 4) - source*s(:, :, :, 2)*s(:, :, :, 3))
            else
                source = dummy*(s(:, :, :, 4) - source*s(:, :, :, 2)*s(:, :, :, 3))
            end if

        end select

        return
    end subroutine FI_CHEM

    ! #######################################################################
    ! #######################################################################
    subroutine FI_RANDOM(random, nx, ny, nz, h, tmp)
        type(term_dt), intent(in) :: random
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(inout) :: h(nx, ny, nz)
        real(wp), intent(inout) :: tmp(nx, ny, nz)

        ! -----------------------------------------------------------------------
        real(wp) dummy

        dummy = random%parameters(1)

        call random_number(tmp)

        tmp = dummy*(tmp*2.0 - 1.0)
        h = h*(1 + tmp)

        return

    end subroutine FI_RANDOM

    !########################################################################
    ! Sinusoidal forcing
    !########################################################################
    subroutine FI_FORCING_0(imax, jmax, kmax, time, visc, u, v, h_u, h_v)
        integer(wi) imax, jmax, kmax
        real(wp) time, visc
        real(wp), dimension(imax, jmax, kmax) :: u, v, h_u, h_v

        ! -----------------------------------------------------------------------
        real(wp) omega, sigma, amplitude
        !  integer(wi) i,j

        !########################################################################
        omega = 2.0_wp*acos(-1.0_wp)
        sigma = 2.0_wp*omega*omega*visc

        !  amplitude =-( 1.0_wp + (sigma/omega)**2 )*omega
        !  amplitude = amplitude*sin(omega*time)
        !  DO j = 1,jmax; DO i = 1,imax
        !     u(i,j,1) = SIN(x(i)*omega)*COS(g(2)%nodes(j)*omega)
        !     v(i,j,1) =-COS(x(i)*omega)*SIN(g(2)%nodes(j)*omega)
        !  ENDDO; ENDDO

        amplitude = sin(omega*time)

        h_u = h_u + amplitude*u
        h_v = h_v + amplitude*v

        return
    end subroutine FI_FORCING_0

    !########################################################################
    ! Velocity field with no-slip
    !########################################################################
    subroutine FI_FORCING_1(imax, jmax, kmax, time, visc, h1, h2, tmp1, tmp2, tmp3, tmp4)
        use OPR_PARTIAL, only: OPR_PARTIAL_X, OPR_PARTIAL_Y
        integer(wi) imax, jmax, kmax
        real(wp) time, visc

        real(wp), dimension(imax*jmax*kmax) :: h1, h2
        real(wp), dimension(imax*jmax*kmax) :: tmp1, tmp2, tmp3, tmp4

        ! -----------------------------------------------------------------------
        integer(wi) ij, i, j, bcs(2, 2)
        real(wp) pi_loc

        bcs = 0

        pi_loc = acos(-1.0_wp)

        ! #######################################################################
        do j = 1, jmax; do i = 1, imax; ij = imax*(j - 1) + i
                !     tmp1(ij) = sin(2.0_wp*pi_loc*g(1)%nodes(i))*       sin(C_4_R*pi_loc*g(2)%nodes(j))
                !     tmp2(ij) =-cos(2.0_wp*pi_loc*g(1)%nodes(i))*(1.0_wp-cos(C_4_R*pi_loc*g(2)%nodes(j)))*0.5_wp
                tmp1(ij) = sin(pi_loc*g(1)%nodes(i))*sin(pi_loc*g(1)%nodes(i))*sin(2.0_wp*pi_loc*g(2)%nodes(j))
                tmp2(ij) = -sin(2.0_wp*pi_loc*g(1)%nodes(i))*sin(pi_loc*g(2)%nodes(j))*sin(pi_loc*g(2)%nodes(j))
            end do; end do

        ! Time terms
        do j = 1, jmax; do i = 1, imax; ij = imax*(j - 1) + i
                h1(ij) = h1(ij) - tmp1(ij)*2.0_wp*pi_loc*sin(2.0_wp*pi_loc*time)
                h2(ij) = h2(ij) - tmp2(ij)*2.0_wp*pi_loc*sin(2.0_wp*pi_loc*time)
            end do; end do

        ! velocities
        do j = 1, jmax; do i = 1, imax; ij = imax*(j - 1) + i
                tmp1(ij) = tmp1(ij)*cos(2.0_wp*pi_loc*time)
                tmp2(ij) = tmp2(ij)*cos(2.0_wp*pi_loc*time)
            end do; end do

        ! Diffusion and convection terms in Ox momentum eqn
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), tmp1, tmp4, tmp3)
        do ij = 1, imax*jmax*kmax
            h1(ij) = h1(ij) - visc*(tmp4(ij)) + (tmp3(ij)*tmp2(ij))
        end do

        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp4, tmp3)
        do ij = 1, imax*jmax*kmax
            h1(ij) = h1(ij) - visc*(tmp4(ij)) + (tmp3(ij)*tmp1(ij))
        end do

        ! Diffusion and convection terms in Oy momentum eqn
        call OPR_PARTIAL_Y(OPR_P2_P1, imax, jmax, kmax, bcs, g(2), tmp2, tmp4, tmp3)
        do ij = 1, imax*jmax*kmax
            h2(ij) = h2(ij) - visc*(tmp4(ij)) + (tmp3(ij)*tmp2(ij))
        end do

        call OPR_PARTIAL_X(OPR_P2_P1, imax, jmax, kmax, bcs, g(1), tmp2, tmp4, tmp3)
        do ij = 1, imax*jmax*kmax
            h2(ij) = h2(ij) - visc*(tmp4(ij)) + (tmp3(ij)*tmp1(ij))
        end do

        ! #######################################################################
        do j = 1, jmax; do i = 1, imax; ij = imax*(j - 1) + i
                !     tmp1(ij) = cos(C_4_R*pi_loc*g(1)%nodes(i))*(2.0_wp-cos(C_4_R*pi_loc*g(2)%nodes(j)))/C_8_R &
                !          - 0.5_wp*(sin(2.0_wp*pi_loc*g(2)%nodes(j)))**4
                tmp1(ij) = sin(2.0_wp*pi_loc*g(1)%nodes(i))*sin(2.0_wp*pi_loc*g(2)%nodes(j))
                tmp1(ij) = tmp1(ij)*(cos(2.0_wp*pi_loc*time))**2
            end do; end do

        ! Pressure gradient
        call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), tmp1, tmp2)
        call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), tmp1, tmp3)

        do ij = 1, imax*jmax*kmax
            h1(ij) = h1(ij) + tmp2(ij)
            h2(ij) = h2(ij) + tmp3(ij)
        end do

        return
    end subroutine FI_FORCING_1

end module FI_SOURCES
