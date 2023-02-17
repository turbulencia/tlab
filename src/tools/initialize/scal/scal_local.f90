#include "dns_const.h"

module SCAL_LOCAL
    use TLAB_CONSTANTS, only: wp, wi, pi_wp, big_wp
    use TLAB_TYPES, only: profiles_dt, discrete_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal, MAX_NSP
    use TLAB_VARS, only: g, sbg
    use TLAB_VARS, only: rtime ! rtime is overwritten in io_read_fields
    use TLAB_ARRAYS, only: wrk1d
    use TLAB_POINTERS_3D, only: p_wrk2d, p_wrk3d
    use IO_FIELDS
    use PROFILES
    use AVGS, only: AVG1V2D
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    implicit none
    save
    private

    public :: Sini, flag_s, flag_mixture, norm_ini_s, norm_ini_radiation, fp
    public :: SCAL_FLUCTUATION_PLANE, SCAL_FLUCTUATION_VOLUME

    ! -------------------------------------------------------------------
    integer(wi) :: flag_s, flag_mixture

    type(profiles_dt) :: Sini(MAX_NSP)                        ! Geometry of perturbation of initial boundary condition
    type(profiles_dt) :: prof_loc
    real(wp) :: norm_ini_s(MAX_NSP), norm_ini_radiation         ! Scaling of perturbation
    type(discrete_dt) :: fp                                     ! Discrete perturbation

    ! -------------------------------------------------------------------
    integer(wi) i, j, k

    integer(wi) im, idsp, kdsp
    real(wp) wx, wz, wx_1, wz_1

    real(wp), dimension(:), pointer :: xn, zn

contains

! ###################################################################
    subroutine SCAL_SHAPE(is, prof)
        integer(wi) is
        real(wp), dimension(jmax, 1), intent(out) :: prof

        ! -------------------------------------------------------------------
        real(wp) yr

        real(wp), dimension(:), pointer :: yn

        ! ###################################################################
        yn => g(2)%nodes

        prof_loc = Sini(is)
        prof_loc%delta = 1.0_wp
        prof_loc%mean = 0.0_wp
        do j = 1, jmax
            prof(j, 1) = PROFILES_CALCULATE(prof_loc, yn(j))
        end do

        select case (Sini(is)%type)
        case (PROFILE_GAUSSIAN_SURFACE) ! set perturbation and its normal derivative to zero at the boundaries
            do j = 1, jmax
                yr = 0.5_wp*(yn(j) - yn(1))/Sini(is)%thick
                prof(j, 1) = prof(j, 1)*TANH(yr)**2
                yr = -0.5_wp*(yn(j) - yn(jmax))/Sini(is)%thick
                prof(j, 1) = prof(j, 1)*TANH(yr)**2
            end do

        end select

        return
    end subroutine SCAL_SHAPE

! ###################################################################
    subroutine SCAL_FLUCTUATION_VOLUME(is, s, tmp)
        integer(wi) is
        real(wp), dimension(imax, jmax, kmax), intent(out) :: s
        real(wp), dimension(imax, jmax, kmax), intent(inout) :: tmp

        ! -------------------------------------------------------------------
        real(wp) dummy, amplify

! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; kdsp = ims_offset_k
#else
        idsp = 0; kdsp = 0
#endif

        xn => g(1)%nodes
        zn => g(3)%nodes

        call SCAL_SHAPE(is, wrk1d)

        select case (flag_s)
        case (1)   ! Broadband case
            dummy = rtime   ! rtime is overwritten in io_read_fields
            call IO_READ_FIELDS('scal.rand', IO_SCAL, imax, jmax, kmax, inb_scal, is, tmp)
            rtime = dummy

            amplify = 0.0_wp
            do j = 1, jmax
                dummy = AVG1V2D(imax, jmax, kmax, j, 1, tmp)       ! Calculate mean
                p_wrk3d(:, j, :) = (tmp(:, j, :) - dummy)*wrk1d(j, 1)  ! Remove mean and apply shape function
            end do

        case (2)   ! Discrete case
            wx_1 = 2.0_wp*pi_wp/g(1)%scale ! Fundamental wavelengths
            wz_1 = 2.0_wp*pi_wp/g(3)%scale

            p_wrk2d = 0.0_wp
            do im = 1, fp%size
                wx = real(fp%modex(im), wp)*wx_1
                wz = real(fp%modez(im), wp)*wz_1
                do k = 1, kmax
                    p_wrk2d(:, k, 1) = p_wrk2d(:, k, 1) &
                            + fp%amplitude(im)*COS(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*COS(wz*zn(kdsp + k) + fp%phasez(im))
                end do
            end do

            do k = 1, kmax
                do j = 1, jmax
                    p_wrk3d(:, j, k) = p_wrk2d(:, k, 1)*wrk1d(j, 1)
                end do
            end do

        end select

        if (norm_ini_s(is) > 0.0_wp) call SCAL_NORMALIZE(is, p_wrk3d)

        s = s + p_wrk3d

        return
    end subroutine SCAL_FLUCTUATION_VOLUME

! ###################################################################
    subroutine SCAL_FLUCTUATION_PLANE(is, s)
        integer(wi) is
        real(wp), dimension(imax, jmax, kmax) :: s

        ! -------------------------------------------------------------------
        real(wp) dummy
        real(wp) xcenter, zcenter, rcenter, amplify

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; kdsp = ims_offset_k
#else
        idsp = 0; kdsp = 0
#endif

        xn => g(1)%nodes
        zn => g(3)%nodes

        ! ###################################################################
#define disp(i,k) p_wrk2d(i,k,1)

        disp(:, :) = 0.0_wp
        select case (flag_s)
        case (4, 6, 8)   ! Broadband case
            dummy = rtime   ! rtime is overwritten in io_read_fields
            call IO_READ_FIELDS('scal.rand', IO_SCAL, imax, 1, kmax, 1, 1, disp(:, :))
            rtime = dummy
            dummy = AVG1V2D(imax, 1, kmax, 1, 1, disp(:, :))     ! remove mean
            disp(:, :) = disp(:, :) - dummy

        case (5, 7, 9)   ! Discrete case
            wx_1 = 2.0_wp*pi_wp/g(1)%scale ! Fundamental wavelengths
            wz_1 = 2.0_wp*pi_wp/g(3)%scale

            do im = 1, fp%size
                wx = real(fp%modex(im), wp)*wx_1
                wz = real(fp%modez(im), wp)*wz_1

                if (fp%type == PROFILE_TANH_COS) then ! Smoothed step funtion Tanh(a*Cos(\xi/b))
                    if (fp%parameters(2) <= 0.0_wp) then; dummy = big_wp; 
                    else; dummy = 0.5_wp/(wx*fp%parameters(2)); end if
                    do k = 1, kmax
                        disp(:, k) = disp(:, k) &
                + fp%amplitude(im)*TANH(dummy*COS(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*COS(wz*zn(kdsp + k) + fp%phasez(im)))
                    end do

                else
                    do k = 1, kmax
                        disp(:, k) = disp(:, k) &
                            + fp%amplitude(im)*COS(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*COS(wz*zn(kdsp + k) + fp%phasez(im))
                    end do

                end if

            end do

        end select

        ! Modulation
        if (fp%type == PROFILE_GAUSSIAN .and. fp%parameters(1) > 0.0_wp) then
            do k = 1, kmax
                do i = 1, imax
                    xcenter = g(1)%nodes(i + idsp) - g(1)%scale*fp%phasex(1) - g(1)%nodes(1)
                    if (g(3)%size > 1) then; zcenter = g(3)%nodes(k + kdsp) - g(3)%scale*fp%phasez(1) - g(3)%nodes(1)
                    else; zcenter = 0.0_wp; end if
                    rcenter = SQRT(xcenter**2 + zcenter**2)
                    amplify = EXP(-0.5_wp*(rcenter/fp%parameters(1))**2)
                    disp(i, k) = disp(i, k)*amplify
                end do
            end do
        end if

        ! ###################################################################
        select case (flag_s)
        case (4, 5)           ! Perturbation in the centerplane
            do k = 1, kmax
                do i = 1, imax
                    do j = 1, jmax
                        s(i, j, k) = PROFILES_CALCULATE(sbg(is), g(2)%nodes(j) - disp(i, k))
                    end do
                end do
            end do

        case (6, 7)           ! Perturbation in the thickness
            prof_loc = sbg(is)

            do k = 1, kmax
                do i = 1, imax
                    prof_loc%thick = sbg(is)%thick + disp(i, k)

                    do j = 1, jmax
                        s(i, j, k) = PROFILES_CALCULATE(prof_loc, g(2)%nodes(j))
                    end do

                end do
            end do

        case (8, 9)           ! Perturbation in the magnitude (constant derivative)
            prof_loc = sbg(is)

            do k = 1, kmax
                do i = 1, imax
                    prof_loc%delta = sbg(is)%delta + disp(i, k)
                    prof_loc%mean = (prof_loc%delta)*0.5_wp
                    if (sbg(is)%delta > 0) prof_loc%thick = prof_loc%delta/sbg(is)%delta*sbg(is)%thick

                    do j = 1, jmax
                        s(i, j, k) = PROFILES_CALCULATE(prof_loc, g(2)%nodes(j))
                    end do

                end do
            end do

        end select

#undef disp

        return
    end subroutine SCAL_FLUCTUATION_PLANE

! ###################################################################
    subroutine SCAL_NORMALIZE(is, s)
        integer(wi) is
        real(wp), dimension(imax, jmax, kmax), intent(inout) :: s

        ! -------------------------------------------------------------------
        real(wp) dummy, amplify

        ! ###################################################################
        amplify = 0.0_wp                                      ! Maximum across the layer
        do j = 1, jmax
            dummy = AVG1V2D(imax, jmax, kmax, j, 2, s)
            amplify = MAX(dummy, amplify)
        end do

        amplify = norm_ini_s(is)/SQRT(amplify)           ! Scaling factor to normalize to maximum rms

        s = s*amplify

        return
    end subroutine SCAL_NORMALIZE

end module SCAL_LOCAL
