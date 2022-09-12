#include "types.h"
#include "dns_const.h"

module SCAL_LOCAL

    use TLAB_TYPES, only: cp, ci, profiles_dt, discrete_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal, MAX_NSP
    use TLAB_VARS, only: g, sbg
    use TLAB_VARS, only: rtime ! rtime is overwritten in io_read_fields
    use IO_FIELDS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    implicit none
    save
    private
    public :: Sini, flag_s, flag_mixture, norm_ini_s, norm_ini_radiation, fp
    public :: SCAL_FLUCTUATION_PLANE, SCAL_FLUCTUATION_VOLUME

    ! -------------------------------------------------------------------
    integer(ci) :: flag_s, flag_mixture

    type(profiles_dt) :: Sini(MAX_NSP)                        ! Geometry of perturbation of initial boundary condition
    type(profiles_dt) :: prof_loc
    real(cp) :: norm_ini_s(MAX_NSP), norm_ini_radiation         ! Scaling of perturbation
    type(discrete_dt) :: fp                                     ! Discrete perturbation

    ! -------------------------------------------------------------------
    integer(ci) i, j, k

    integer(ci) im, idsp, kdsp
    real(cp) wx, wz, wx_1, wz_1

    real(cp), dimension(:), pointer :: xn, zn

contains

! ###################################################################
    subroutine SCAL_SHAPE(is, wrk1d)

        integer(ci) is
        real(cp), dimension(jmax, 1), intent(inout) :: wrk1d

        ! -------------------------------------------------------------------
        real(cp) PROFILES, ycenter, yr
        external PROFILES

        real(cp), dimension(:), pointer :: yn

        ! ###################################################################
        yn => g(2)%nodes

        prof_loc = Sini(is)
        prof_loc%delta = C_1_R
        prof_loc%mean = C_0_R
        ycenter = g(2)%nodes(1) + g(2)%scale*prof_loc%ymean_rel
        do j = 1, jmax
            wrk1d(j, 1) = PROFILES(prof_loc, ycenter, yn(j))
        end do

        select case (Sini(is)%type)
        case (PROFILE_GAUSSIAN_SURFACE) ! set perturbation and its normal derivative to zero at the boundaries
            do j = 1, jmax
                yr = C_05_R*(yn(j) - yn(1))/Sini(is)%thick
                wrk1d(j, 1) = wrk1d(j, 1)*TANH(yr)**2
                yr = -C_05_R*(yn(j) - yn(jmax))/Sini(is)%thick
                wrk1d(j, 1) = wrk1d(j, 1)*TANH(yr)**2
            end do

        end select

        return
    end subroutine SCAL_SHAPE

! ###################################################################
    subroutine SCAL_FLUCTUATION_VOLUME(is, s, tmp, wrk1d, wrk2d, wrk3d)

        integer(ci) is
        real(cp), dimension(imax, jmax, kmax), intent(out) :: s
        real(cp), dimension(jmax, 1), intent(inout) :: wrk1d
        real(cp), dimension(imax, kmax, 1), intent(inout) :: wrk2d
        real(cp), dimension(imax, jmax, kmax), intent(inout) :: tmp, wrk3d

        ! -------------------------------------------------------------------
        real(cp) AVG1V2D, dummy, amplify
        external AVG1V2D

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
            call IO_READ_FIELDS('scal.rand', IO_SCAL, imax, jmax, kmax, inb_scal, is, tmp, wrk3d)
            rtime = dummy

            amplify = C_0_R
            do j = 1, jmax
                dummy = AVG1V2D(imax, jmax, kmax, j, 1, tmp)       ! Calculate mean
                wrk3d(:, j, :) = (tmp(:, j, :) - dummy)*wrk1d(j, 1)  ! Remove mean and apply shape function
            end do

        case (2)   ! Discrete case
            wx_1 = C_2_R*C_PI_R/g(1)%scale ! Fundamental wavelengths
            wz_1 = C_2_R*C_PI_R/g(3)%scale

            wrk2d = C_0_R
            do im = 1, fp%size
                wx = M_REAL(fp%modex(im))*wx_1
                wz = M_REAL(fp%modez(im))*wz_1
                do k = 1, kmax
wrk2d(:,k,1) = wrk2d(:,k,1) + fp%amplitude(im) *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) )
                end do
            end do

            do k = 1, kmax
                do j = 1, jmax
                    wrk3d(:, j, k) = wrk2d(:, k, 1)*wrk1d(j, 1)
                end do
            end do

        end select

        if (norm_ini_s(is) > C_0_R) call SCAL_NORMALIZE(is, wrk3d)

        s = s + wrk3d

        return
    end subroutine SCAL_FLUCTUATION_VOLUME

! ###################################################################
    subroutine SCAL_FLUCTUATION_PLANE(is, s, disp)

        integer(ci) is
        real(cp), dimension(imax, jmax, kmax) :: s
        real(cp), dimension(imax, kmax) :: disp

        ! -------------------------------------------------------------------
        real(cp) dummy, ycenter
        real(cp) AVG1V2D, PROFILES
        real(cp) xcenter, zcenter, rcenter, amplify

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; kdsp = ims_offset_k
#else
        idsp = 0; kdsp = 0
#endif

        xn => g(1)%nodes
        zn => g(3)%nodes

        ! ###################################################################
        disp = C_0_R
        select case (flag_s)
        case (4, 6, 8)   ! Broadband case
            dummy = rtime   ! rtime is overwritten in io_read_fields
            call IO_READ_FIELDS('scal.rand', IO_SCAL, imax, 1, kmax, 1, 1, disp, s) ! using array s as aux array
            rtime = dummy
            dummy = AVG1V2D(imax, 1, kmax, 1, 1, disp)     ! remove mean
            disp = disp - dummy

        case (5, 7, 9)   ! Discrete case
            wx_1 = C_2_R*C_PI_R/g(1)%scale ! Fundamental wavelengths
            wz_1 = C_2_R*C_PI_R/g(3)%scale

            do im = 1, fp%size
                wx = M_REAL(fp%modex(im))*wx_1
                wz = M_REAL(fp%modez(im))*wz_1

                if (fp%type == PROFILE_TANH_COS) then ! Smoothed step funtion Tanh(a*Cos(\xi/b))
                    if (fp%parameters(2) <= C_0_R) then; dummy = C_BIG_R; 
                    else; dummy = C_05_R/(wx*fp%parameters(2)); end if
                    do k = 1, kmax
            disp(:,k) = disp(:,k) + fp%amplitude(im) *TANH( dummy *COS( wx *xn(idsp+1:idsp+imax) +fp%phasex(im) ) *COS( wz *zn(kdsp+k) +fp%phasez(im) ) )
                    end do

                else
                    do k = 1, kmax
    disp(:, k) = disp(:, k) + fp%amplitude(im)*COS(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*COS(wz*zn(kdsp + k) + fp%phasez(im))
                    end do

                end if

            end do

        end select

        ! Modulation
        if (fp%type == PROFILE_GAUSSIAN .and. fp%parameters(1) > C_0_R) then
            do k = 1, kmax
                do i = 1, imax
                    xcenter = g(1)%nodes(i + idsp) - g(1)%scale*fp%phasex(1) - g(1)%nodes(1)
                    if (g(3)%size > 1) then; zcenter = g(3)%nodes(k + kdsp) - g(3)%scale*fp%phasez(1) - g(3)%nodes(1)
                    else; zcenter = C_0_R; end if
                    rcenter = SQRT(xcenter**2 + zcenter**2)
                    amplify = EXP(-C_05_R*(rcenter/fp%parameters(1))**2)
                    disp(i, k) = disp(i, k)*amplify
                end do
            end do
        end if

        ! ###################################################################
        select case (flag_s)
        case (4, 5)           ! Perturbation in the centerplane
            do k = 1, kmax
                do i = 1, imax
                    ycenter = g(2)%nodes(1) + g(2)%scale*sbg(is)%ymean_rel + disp(i, k)
                    do j = 1, jmax
                        s(i, j, k) = PROFILES(sbg(is), ycenter, g(2)%nodes(j))
                    end do
                end do
            end do

        case (6, 7)           ! Perturbation in the thickness
            prof_loc = sbg(is)

            do k = 1, kmax
                do i = 1, imax
                    ycenter = g(2)%nodes(1) + g(2)%scale*sbg(is)%ymean_rel
                    prof_loc%thick = sbg(is)%thick + disp(i, k)

                    do j = 1, jmax
                        s(i, j, k) = PROFILES(prof_loc, ycenter, g(2)%nodes(j))
                    end do
                    
                end do
            end do

        case (8, 9)           ! Perturbation in the magnitude (constant derivative)
            prof_loc = sbg(is)

            do k = 1, kmax
                do i = 1, imax
                    ycenter = g(2)%nodes(1) + g(2)%scale*sbg(is)%ymean_rel
                    prof_loc%delta = sbg(is)%delta + disp(i, k)
                    prof_loc%mean = (prof_loc%delta)*C_05_R
                    if (sbg(is)%delta > 0) prof_loc%thick = prof_loc%delta/sbg(is)%delta*sbg(is)%thick

                    do j = 1, jmax
                        s(i, j, k) = PROFILES(prof_loc, ycenter, g(2)%nodes(j))
                    end do

                end do
            end do

        end select

        return
    end subroutine SCAL_FLUCTUATION_PLANE

! ###################################################################
    subroutine SCAL_NORMALIZE(is, s)

        integer(ci) is
        real(cp), dimension(imax, jmax, kmax), intent(inout) :: s

        ! -------------------------------------------------------------------
        real(cp) AVG1V2D, dummy, amplify
        external AVG1V2D

        ! ###################################################################
        amplify = C_0_R                                      ! Maximum across the layer
        do j = 1, jmax
            dummy = AVG1V2D(imax, jmax, kmax, j, 2, s)
            amplify = MAX(dummy, amplify)
        end do

        amplify = norm_ini_s(is)/SQRT(amplify)           ! Scaling factor to normalize to maximum rms

        s = s*amplify

        return
    end subroutine SCAL_NORMALIZE

end module SCAL_LOCAL
