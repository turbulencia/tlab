#include "dns_const.h"
#include "dns_error.h"

module SCAL_LOCAL
    use TLab_Constants, only: wfile,efile, lfile, wp, wi, pi_wp, big_wp
    use TLab_Types, only: profiles_dt, discrete_dt
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_scal, MAX_VARS
    use TLAB_VARS, only: g, sbg
    use TLAB_VARS, only: rtime ! rtime is overwritten in io_read_fields
    use TLAB_ARRAYS, only: wrk1d
    use TLAB_POINTERS_3D, only: p_wrk2d, p_wrk3d
    use TLab_WorkFlow
    use IO_FIELDS
    use AVGS, only: AVG1V2D
    use Profiles
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    implicit none
    save
    private

    ! -------------------------------------------------------------------
    integer(wi) :: flag_s, flag_mixture               ! Type of perturbation in scalar fields
    integer, parameter :: PERT_NONE = 0
    integer, parameter :: PERT_LAYER_BROADBAND = 1
    integer, parameter :: PERT_LAYER_DISCRETE = 2
    integer, parameter :: PERT_PLANE_BROADBAND = 4
    integer, parameter :: PERT_PLANE_DISCRETE = 5
    integer, parameter :: PERT_DELTA_BROADBAND = 6
    integer, parameter :: PERT_DELTA_DISCRETE = 7
    integer, parameter :: PERT_FLUX_BROADBAND = 8
    integer, parameter :: PERT_FLUX_DISCRETE = 9

    type(profiles_dt) :: IniS(MAX_VARS)                          ! Geometry of perturbation of initial boundary condition
    type(profiles_dt) :: prof_loc
    real(wp) :: norm_ini_s(MAX_VARS), norm_ini_radiation         ! Scaling of perturbation
    type(discrete_dt) :: fp                                     ! Discrete perturbation

    public :: flag_s, flag_mixture, IniS, norm_ini_radiation
    public :: PERT_LAYER_BROADBAND, PERT_LAYER_DISCRETE, PERT_PLANE_BROADBAND, PERT_PLANE_DISCRETE
    public :: PERT_DELTA_BROADBAND, PERT_DELTA_DISCRETE, PERT_FLUX_BROADBAND, PERT_FLUX_DISCRETE
    public :: SCAL_READ_LOCAL
    public :: SCAL_FLUCTUATION_PLANE, SCAL_FLUCTUATION_VOLUME

    ! -------------------------------------------------------------------
    integer(wi) i, j, k
    integer(wi) im, idsp, kdsp
    real(wp) wx, wz, wx_1, wz_1
    real(wp), dimension(:), pointer :: xn, zn

contains

    ! ###################################################################
    subroutine SCAL_READ_LOCAL(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        integer(wi) idummy, is
        character*512 sRes
        character*32 bakfile
        character*64 lstr

        integer :: IniSvalid(5) = [PROFILE_NONE, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_ANTISYM, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_SURFACE]

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        call TLAB_WRITE_ASCII(lfile, 'Reading local input data')

        ! ###################################################################
        call TLAB_WRITE_ASCII(bakfile, '#')
        call TLAB_WRITE_ASCII(bakfile, '#[IniFields]')
        call TLAB_WRITE_ASCII(bakfile, '#Scalar=<option>')
        call TLAB_WRITE_ASCII(bakfile, '#NormalizeS=<values>')
        call TLAB_WRITE_ASCII(bakfile, '#Mixture=<string>')

        call SCANINICHAR(bakfile, inifile, 'IniFields', 'Scalar', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; flag_s = 0
        else if (trim(adjustl(sRes)) == 'layerbroadband') then; flag_s = PERT_LAYER_BROADBAND
        else if (trim(adjustl(sRes)) == 'layerdiscrete')  then; flag_s = PERT_LAYER_DISCRETE
        else if (trim(adjustl(sRes)) == 'planebroadband') then; flag_s = PERT_PLANE_BROADBAND
        else if (trim(adjustl(sRes)) == 'planediscrete')  then; flag_s = PERT_PLANE_DISCRETE
        else if (trim(adjustl(sRes)) == 'deltabroadband') then; flag_s = PERT_DELTA_BROADBAND
        else if (trim(adjustl(sRes)) == 'deltadiscrete')  then; flag_s = PERT_DELTA_DISCRETE
        else if (trim(adjustl(sRes)) == 'fluxbroadband')  then; flag_s = PERT_FLUX_BROADBAND
        else if (trim(adjustl(sRes)) == 'fluxdiscrete')   then; flag_s = PERT_FLUX_DISCRETE
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Scalar forcing type unknown')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        call Profiles_ReadBlock(bakfile, inifile, 'IniFields', 'IniS', IniS(1), 'gaussiansurface')      ! if IniS, valid for all
        if (trim(adjustl(sRes)) /= 'none') then
            IniS(2:) = IniS(1)
        else                                                                                            ! if not, read separately
            do is = 1, inb_scal
                write (lstr, *) is
                call Profiles_ReadBlock(bakfile, inifile, 'IniFields', 'IniS'//trim(adjustl(lstr)), IniS(is), 'gaussiansurface')
            end do
        end if
        do is = 1, inb_scal
            if (.not. any(IniSvalid == IniS(is)%type)) then
                call TLAB_WRITE_ASCII(efile, __FILE__//'. Undeveloped IniS type.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end do
        IniS(:)%delta = 1.0_wp
        IniS(:)%mean = 0.0_wp

        call SCANINICHAR(bakfile, inifile, 'IniFields', 'NormalizeS', '-1.0', sRes)
        norm_ini_s(:) = 0.0_wp; idummy = inb_scal
        call LIST_REAL(sRes, idummy, norm_ini_s)
        if (idummy /= inb_scal) then            ! Consistency check
            if (idummy == 1) then
                norm_ini_s(2:) = norm_ini_s(1)
                call TLAB_WRITE_ASCII(wfile, __FILE__//'. Using NormalizeS(1) for all scalars.')
            else
                call TLAB_WRITE_ASCII(efile, __FILE__//'. NormalizeS size does not match number of scalars.')
                call TLAB_STOP(DNS_ERROR_OPTION)
            end if
        end if

        call SCANINIREAL(bakfile, inifile, 'IniFields', 'NormalizeR', '0.0', norm_ini_radiation) ! Radiation field

        ! Additional parameters
        call SCANINICHAR(bakfile, inifile, 'IniFields', 'Mixture', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; flag_mixture = 0
        else if (trim(adjustl(sRes)) == 'equilibrium') then; flag_mixture = 1
        else if (trim(adjustl(sRes)) == 'loadfields') then; flag_mixture = 2
        end if

! Discrete Forcing
        call DISCRETE_READBLOCK(bakfile, inifile, 'Discrete', fp) ! Modulation type in fp%type
!   specific for this tool
        call SCANINIREAL(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', fp%parameters(1))
        call SCANINIREAL(bakfile, inifile, 'Discrete', 'ThickStep', '-1.0', fp%parameters(2))

        return
    end subroutine SCAL_READ_LOCAL

! ###################################################################
    subroutine SCAL_SHAPE(is, prof)
        integer(wi) is
        real(wp), dimension(jmax, 1), intent(out) :: prof

        ! -------------------------------------------------------------------
        real(wp) yr

        real(wp), dimension(:), pointer :: yn

        ! ###################################################################
        yn => g(2)%nodes

        prof_loc = IniS(is)
        prof_loc%delta = 1.0_wp
        prof_loc%mean = 0.0_wp
        do j = 1, jmax
            prof(j, 1) = Profiles_Calculate(prof_loc, yn(j))
        end do

        select case (IniS(is)%type)
        case (PROFILE_GAUSSIAN_SURFACE) ! set perturbation and its normal derivative to zero at the boundaries
            do j = 1, jmax
                yr = 0.5_wp*(yn(j) - yn(1))/IniS(is)%thick
                prof(j, 1) = prof(j, 1)*tanh(yr)**2
                yr = -0.5_wp*(yn(j) - yn(jmax))/IniS(is)%thick
                prof(j, 1) = prof(j, 1)*tanh(yr)**2
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
        case (PERT_LAYER_BROADBAND)
            dummy = rtime   ! rtime is overwritten in io_read_fields
            call IO_READ_FIELDS('scal.rand', IO_SCAL, imax, jmax, kmax, inb_scal, is, tmp)
            rtime = dummy

            amplify = 0.0_wp
            do j = 1, jmax
                dummy = AVG1V2D(imax, jmax, kmax, j, 1, tmp)       ! Calculate mean
                p_wrk3d(:, j, :) = (tmp(:, j, :) - dummy)*wrk1d(j, 1)  ! Remove mean and apply shape function
            end do

        case (PERT_LAYER_DISCRETE)
            wx_1 = 2.0_wp*pi_wp/g(1)%scale ! Fundamental wavelengths
            wz_1 = 2.0_wp*pi_wp/g(3)%scale

            p_wrk2d = 0.0_wp
            do im = 1, fp%size
                wx = real(fp%modex(im), wp)*wx_1
                wz = real(fp%modez(im), wp)*wz_1
                do k = 1, kmax
                    p_wrk2d(:, k, 1) = p_wrk2d(:, k, 1) &
                                       + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*cos(wz*zn(kdsp + k) + fp%phasez(im))
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
        case (PERT_PLANE_BROADBAND, PERT_DELTA_BROADBAND, PERT_FLUX_BROADBAND)
            dummy = rtime   ! rtime is overwritten in io_read_fields
            call IO_READ_FIELDS('scal.rand', IO_SCAL, imax, 1, kmax, 1, 1, disp(:, :))
            rtime = dummy
            dummy = AVG1V2D(imax, 1, kmax, 1, 1, disp(:, :))     ! remove mean
            disp(:, :) = disp(:, :) - dummy

        case (PERT_PLANE_DISCRETE, PERT_DELTA_DISCRETE, PERT_FLUX_DISCRETE)
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
                                  + fp%amplitude(im)*tanh(dummy*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*cos(wz*zn(kdsp + k) + fp%phasez(im)))
                    end do

                else
                    do k = 1, kmax
                        disp(:, k) = disp(:, k) &
                                     + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im))*cos(wz*zn(kdsp + k) + fp%phasez(im))
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
                    rcenter = sqrt(xcenter**2 + zcenter**2)
                    amplify = exp(-0.5_wp*(rcenter/fp%parameters(1))**2)
                    disp(i, k) = disp(i, k)*amplify
                end do
            end do
        end if

        ! ###################################################################
        select case (flag_s)
        case (PERT_PLANE_BROADBAND, PERT_PLANE_DISCRETE)    ! Perturbation in the centerplane
            do k = 1, kmax
                do i = 1, imax
                    do j = 1, jmax
                        s(i, j, k) = Profiles_Calculate(sbg(is), g(2)%nodes(j) - disp(i, k))
                    end do
                end do
            end do

        case (PERT_DELTA_BROADBAND, PERT_DELTA_DISCRETE)    ! Perturbation in the thickness
            prof_loc = sbg(is)

            do k = 1, kmax
                do i = 1, imax
                    prof_loc%thick = sbg(is)%thick + disp(i, k)

                    do j = 1, jmax
                        s(i, j, k) = Profiles_Calculate(prof_loc, g(2)%nodes(j))
                    end do

                end do
            end do

        case (PERT_FLUX_BROADBAND, PERT_FLUX_DISCRETE)      ! Perturbation in the magnitude (constant derivative)
            prof_loc = sbg(is)

            do k = 1, kmax
                do i = 1, imax
                    prof_loc%delta = sbg(is)%delta + disp(i, k)
                    prof_loc%mean = (prof_loc%delta)*0.5_wp
                    if (sbg(is)%delta > 0) prof_loc%thick = prof_loc%delta/sbg(is)%delta*sbg(is)%thick

                    do j = 1, jmax
                        s(i, j, k) = Profiles_Calculate(prof_loc, g(2)%nodes(j))
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
            amplify = max(dummy, amplify)
        end do

        amplify = norm_ini_s(is)/sqrt(amplify)           ! Scaling factor to normalize to maximum rms

        s = s*amplify

        return
    end subroutine SCAL_NORMALIZE

end module SCAL_LOCAL
