#include "dns_const.h"
#include "dns_error.h"

module FLOW_LOCAL
    use TLab_Constants, only: efile, lfile, wfile
    use TLab_Constants, only: wp, wi, pi_wp, BCS_DD, BCS_DN, BCS_ND, BCS_NN
    use Discrete, only: discrete_dt
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Memory, only: inb_wrk2d, inb_txc
    use TLab_Grid, only: x, y, z
    use TLab_Time, only: itime, rtime
    use TLab_WorkFlow, only: stagger_on
    use FDM, only: g
    use Tlab_Background, only: tbg, hbg
    use TLab_Pointers_3D, only: p_wrk1d, p_wrk2d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use THERMO_THERMAL
    use THERMO_AIRWATER
    use Thermo_Anelastic
    use IO_Fields
    use Averages, only: AVG1V2D
    use Profiles, only: profiles_dt, Profiles_ReadBlock, Profiles_Calculate
    use Profiles, only: PROFILE_NONE, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_ANTISYM, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_SURFACE, PROFILE_PARABOLIC_SURFACE
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_offset_i, ims_offset_k
#endif
    use FI_VECTORCALCULUS
    use OPR_Partial
    use OPR_Elliptic
    use OPR_Filters, only: PressureFilter
    use Discrete, only: Discrete_ReadBlock
    implicit none
    save
    private

    ! -------------------------------------------------------------------
    integer(wi) :: flag_u, flag_t               ! Type of perturbation in velocity and thermodynamic fields
    integer, parameter :: PERT_NONE = 0
    integer, parameter :: PERT_DISCRETE = 1
    integer, parameter :: PERT_BROADBAND = 2
    integer, parameter :: PERT_BROADBAND_VORTICITY = 3
    integer, parameter :: PERT_BROADBAND_POTENTIAL = 4

    integer(wi) :: flag_wall                    ! BCs at j1/jmax:   0, No-Slip/No-Slip
    !                                                               1, Free-Slip/No-Slip
    !                                                               2, No-Slip/Free-Slip
    !                                                               3, Free-Slip/Free-Slip
    logical :: RemoveDilatation

    type(profiles_dt) :: IniK                   ! Geometry of perturbation of initial boundary condition
    real(wp) :: norm_ini_u, norm_ini_p          ! Scaling of perturbation
    type(discrete_dt) :: fp                     ! Discrete perturbation

    public :: flag_u, flag_t, IniK, norm_ini_p
    public :: PERT_DISCRETE, PERT_BROADBAND, PERT_BROADBAND_POTENTIAL, PERT_BROADBAND_VORTICITY
    public :: Iniflow_Initialize_Parameters
    public :: VELOCITY_BROADBAND, VELOCITY_DISCRETE
    public :: DENSITY_FLUCTUATION, PRESSURE_FLUCTUATION ! Only in compressible formulation

    ! -------------------------------------------------------------------
    integer(wi) i, j, k
    integer(wi) im, idsp, kdsp
    real(wp) wx, wz, wx_1, wz_1
    real(wp), dimension(:), pointer :: xn, zn

contains

    ! ###################################################################
    subroutine Iniflow_Initialize_Parameters(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes

        integer(wi) :: bcs_flow_jmin, bcs_flow_jmax
 integer :: IniKvalid(6) = [PROFILE_NONE, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_ANTISYM, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_SURFACE, PROFILE_PARABOLIC_SURFACE]

        ! ###################################################################
        call TLab_Write_ASCII(lfile, 'Reading local input data')

        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Inifields'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Velocity=<VelocityDiscrete/VelocityBroadband/PotentialBroadband/VorticityBroadband>')
        call TLab_Write_ASCII(bakfile, '#Temperature=<option>')
        call TLab_Write_ASCII(bakfile, '#ForceDilatation=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#NormalizeK=<value>')
        call TLab_Write_ASCII(bakfile, '#NormalizeP=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Velocity', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; flag_u = PERT_NONE
        else if (trim(adjustl(sRes)) == 'velocitydiscrete') then; flag_u = PERT_DISCRETE
        else if (trim(adjustl(sRes)) == 'velocitybroadband') then; flag_u = PERT_BROADBAND
        else if (trim(adjustl(sRes)) == 'vorticitybroadband') then; flag_u = PERT_BROADBAND_VORTICITY
        else if (trim(adjustl(sRes)) == 'potentialbroadband') then; flag_u = PERT_BROADBAND_POTENTIAL
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Velocity forcing type unknown')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        RemoveDilatation = .true.
        call ScanFile_Char(bakfile, inifile, block, 'ForceDilatation', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'no') RemoveDilatation = .false.

        call Profiles_ReadBlock(bakfile, inifile, block, 'IniK', IniK)
        if (.not. any(IniKvalid == IniK%type)) then
            call TLab_Write_ASCII(efile, __FILE__//'. Undeveloped IniK type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if
        IniK%delta = 1.0_wp
        IniK%mean = 0.0_wp

        call ScanFile_Real(bakfile, inifile, block, 'NormalizeK', '-1.0', norm_ini_u)
        call ScanFile_Real(bakfile, inifile, block, 'NormalizeP', '-1.0', norm_ini_p)

        ! Compressible formulation
        call ScanFile_Char(bakfile, inifile, block, 'Temperature', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; flag_t = 0
        else if (trim(adjustl(sRes)) == 'planediscrete') then; flag_t = PERT_DISCRETE
        else if (trim(adjustl(sRes)) == 'planebroadband') then; flag_t = PERT_BROADBAND
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Temperature forcing type unknown')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! Boundary conditions for the perturbation
        flag_wall = 0
        call ScanFile_Char(bakfile, inifile, 'BoundaryConditions', 'VelocityJmin', 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then; bcs_flow_jmin = DNS_BCS_NONE
        else if (trim(adjustl(sRes)) == 'noslip') then; bcs_flow_jmin = DNS_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; bcs_flow_jmin = DNS_BCS_NEUMANN; flag_wall = flag_wall + 1
        else
            call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.VelocityJmin.')
            call TLab_Stop(DNS_ERROR_IBC)
        end if
        call ScanFile_Char(bakfile, inifile, 'BoundaryConditions', 'VelocityJmax', 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then; bcs_flow_jmax = DNS_BCS_NONE
        else if (trim(adjustl(sRes)) == 'noslip') then; bcs_flow_jmax = DNS_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; bcs_flow_jmax = DNS_BCS_NEUMANN; flag_wall = flag_wall + 2
        else
            call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.VelocityJmax.')
            call TLab_Stop(DNS_ERROR_IBC)
        end if

        ! Discrete Forcing
        call Discrete_ReadBlock(bakfile, inifile, 'Discrete', fp) ! Modulation type in fp%type
        ! specific for this tool
        call ScanFile_Real(bakfile, inifile, 'Discrete', 'Broadening', '-1.0', fp%parameters(1))
        call ScanFile_Real(bakfile, inifile, 'Discrete', 'ThickStep', '-1.0', fp%parameters(2))

        ! ###################################################################
        ! Consistency check
        ! ###################################################################
        ! Staggering of the pressure grid not implemented here
        if (stagger_on) then
            call TLab_Write_ASCII(wfile, __FILE__//'. Staggering of the pressure grid not yet implemented.')
            stagger_on = .false. ! turn staggering off for OPR_Poisson_FourierXZ_Factorize(...)
        end if
        if (any(PressureFilter%type /= DNS_FILTER_NONE)) then
            call TLab_Write_ASCII(wfile, __FILE__//'. Pressure and dpdy Filter not implemented here.')
        end if

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        inb_wrk2d = max(inb_wrk2d, 3)

        inb_txc = 2
        if (flag_u /= PERT_NONE) inb_txc = max(inb_txc, 8)

        return
    end subroutine Iniflow_Initialize_Parameters

    ! ###################################################################
    subroutine VELOCITY_DISCRETE(u, v, w)
        real(wp), dimension(imax, jmax, kmax), intent(OUT) :: u, v, w

        ! -------------------------------------------------------------------
        real(wp) factorx, factorz

        ! ###################################################################
#ifdef USE_MPI
        idsp = ims_offset_i; kdsp = ims_offset_k
#else
        idsp = 0; kdsp = 0
#endif

        xn => x
        zn => z

        call FLOW_SHAPE(p_wrk1d)

        wx_1 = 2.0_wp*pi_wp/g(1)%scale ! Fundamental wavelengths
        wz_1 = 2.0_wp*pi_wp/g(3)%scale

        p_wrk2d = 0.0_wp
        do im = 1, fp%size
            wx = real(fp%modex(im), wp)*wx_1
            wz = real(fp%modez(im), wp)*wz_1

            ! Factor to impose solenoidal constraint
            if (fp%modex(im) == 0 .and. fp%modez(im) == 0) then; exit
            elseif (fp%modez(im) == 0) then; factorx = 1.0_wp/wx; factorz = 0.0_wp
            elseif (fp%modex(im) == 0) then; factorx = 0.0_wp; factorz = 1.0_wp/wz
            else; factorx = 0.5_wp/wx; factorz = 0.5_wp/wz
            end if

            do k = 1, kmax
                p_wrk2d(:, k, 2) = p_wrk2d(:, k, 2) + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im)) &
                                   *cos(wz*zn(kdsp + k) + fp%phasez(im))
                p_wrk2d(:, k, 1) = p_wrk2d(:, k, 1) + fp%amplitude(im)*sin(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im)) &
                                   *cos(wz*zn(kdsp + k) + fp%phasez(im))*factorx
                p_wrk2d(:, k, 3) = p_wrk2d(:, k, 3) + fp%amplitude(im)*cos(wx*xn(idsp + 1:idsp + imax) + fp%phasex(im)) &
                                   *sin(wz*zn(kdsp + k) + fp%phasez(im))*factorz
            end do

        end do

        do k = 1, kmax
            do j = 1, jmax
                u(:, j, k) = p_wrk2d(:, k, 1)*p_wrk1d(j, 2)
                v(:, j, k) = p_wrk2d(:, k, 2)*p_wrk1d(j, 1)
                w(:, j, k) = p_wrk2d(:, k, 3)*p_wrk1d(j, 2)
            end do
        end do

        if (norm_ini_u >= 0.0_wp) call FLOW_NORMALIZE(u, v, w)

        return
    end subroutine VELOCITY_DISCRETE

    ! ###################################################################
    subroutine VELOCITY_BROADBAND(u, v, w, ax, ay, az, tmp4, tmp5)
        use FI_VECTORCALCULUS

        real(wp), dimension(imax, jmax, kmax), intent(OUT) :: u, v, w
        real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: ax, ay, az, tmp4, tmp5

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2), bcs2(2, 2)
        real(wp) dummy, params(0)

        real(wp), allocatable :: bcs_hb(:), bcs_ht(:)

        ! ###################################################################
        bcs = 0

        call IO_Read_Fields('flow.rand', imax, jmax, kmax, itime, 3, 1, u, params)
        call IO_Read_Fields('flow.rand', imax, jmax, kmax, itime, 3, 2, v, params)
        call IO_Read_Fields('flow.rand', imax, jmax, kmax, itime, 3, 3, w, params)

        do j = 1, jmax   ! Remove mean
            dummy = AVG1V2D(imax, jmax, kmax, j, 1, u)
            u(:, j, :) = u(:, j, :) - dummy
            dummy = AVG1V2D(imax, jmax, kmax, j, 1, v)
            v(:, j, :) = v(:, j, :) - dummy
            dummy = AVG1V2D(imax, jmax, kmax, j, 1, w)
            w(:, j, :) = w(:, j, :) - dummy
        end do

        call FLOW_SHAPE(p_wrk1d)

        ! ###################################################################
        select case (flag_u)
        case (PERT_BROADBAND)                           ! Velocity u given
            do j = 1, jmax
                u(:, j, :) = u(:, j, :)*p_wrk1d(j, 2)
                v(:, j, :) = v(:, j, :)*p_wrk1d(j, 1)
                w(:, j, :) = w(:, j, :)*p_wrk1d(j, 2)
            end do

        case (PERT_BROADBAND_POTENTIAL)                 ! Velocity potential u given, calculate u = rot(u)
            do j = 1, jmax
                ax(:, j, :) = u(:, j, :)*p_wrk1d(j, 1)  ! Horizontal components of vector potential give vertical velocity
                ay(:, j, :) = v(:, j, :)*p_wrk1d(j, 2)
                az(:, j, :) = w(:, j, :)*p_wrk1d(j, 1)
            end do

            ! Cannot use fi_curl. I need to impose BCs to zero to get zero velocity there
            bcs2 = 0
            if (any([BCS_DD, BCS_DN] == flag_wall)) bcs2(1, 1) = 1 ! no-slip at ymin
            if (any([BCS_DD, BCS_ND] == flag_wall)) bcs2(2, 1) = 1 ! no-slip at ymax
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs2, g(2), az, u)
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), ay, tmp4)
            u = u - tmp4
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), ax, v)
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), az, tmp4)
            v = v - tmp4
            if (g(3)%size > 1) then
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), ay, w)
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs2, g(2), ax, tmp4)
                w = w - tmp4
            end if

        case (PERT_BROADBAND_VORTICITY)                 ! Vorticity given, solve lap(u) = - rot(vort), vort = rot(u)
            allocate (bcs_hb(imax*kmax), bcs_ht(imax*kmax))

            call FI_CURL(imax, jmax, kmax, u, v, w, ax, ay, az, tmp4)
            do j = 1, jmax
                ax(:, j, :) = -ax(:, j, :)*p_wrk1d(j, 2)
                ay(:, j, :) = -ay(:, j, :)*p_wrk1d(j, 1)
                az(:, j, :) = -az(:, j, :)*p_wrk1d(j, 2)
            end do
            call FI_CURL(imax, jmax, kmax, ax, ay, az, u, v, w, tmp4)

            ! Solve lap(u) = - (rot(vort))_x
            if (g(1)%periodic .and. g(3)%periodic) then
                bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
                call OPR_Poisson(imax, jmax, kmax, flag_wall, u, tmp4, tmp5, bcs_hb, bcs_ht)
            else                                        ! General treatment; undevelop
            end if

            ! Solve lap(v) = - (rot(vort))_y with no penetration bcs
            if (g(1)%periodic .and. g(3)%periodic) then
                bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
                call OPR_Poisson(imax, jmax, kmax, BCS_DD, v, tmp4, tmp5, bcs_hb, bcs_ht)
            else                                        ! General treatment; undevelop
            end if

            ! Solve lap(w) = - (rot(vort))_z
            if (g(3)%size > 1) then
                if (g(1)%periodic .and. g(3)%periodic) then
                    bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
                    call OPR_Poisson(imax, jmax, kmax, flag_wall, w, tmp4, tmp5, bcs_hb, bcs_ht)
                else                                    ! General treatment; undevelop
                end if
            end if

        end select

        ! ###################################################################
        if (RemoveDilatation) then              ! Remove dilatation (vort might not be a vorticity field because it was not solenoidal)
            call FI_SOLENOIDAL(imax, jmax, kmax, u, v, w, ax, ay, az)
        end if

        if (g(3)%size == 1) w = 0.0_wp          ! Impose zero spanwise velocity in 2D case

        if (norm_ini_u >= 0.0_wp) call FLOW_NORMALIZE(u, v, w)

        return
    end subroutine VELOCITY_BROADBAND

    ! ###################################################################
    subroutine FLOW_SHAPE(profs)
        real(wp), dimension(jmax, 5), intent(inout) :: profs

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)
        real(wp) yr

        ! ###################################################################
        bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

        do j = 1, jmax                                              ! Wall-normal velocity
            profs(j, 1) = Profiles_Calculate(IniK, y(j))
        end do
        call OPR_Partial_Y(OPR_P1, 1, jmax, 1, bcs, g(2), profs(1, 1), profs(1, 2))
        profs(:, 2) = -profs(:, 2)                                  ! Negative of the derivative of f, wall-parallel velocity

        select case (IniK%type)
        case (PROFILE_PARABOLIC_SURFACE)
            ! Zero wall-parallel velocity for no-slip condition, multiply by parabolic again, f=f*f
            profs(:, 2) = 2.0_wp*profs(:, 2)*profs(:, 1)            ! Wall-parallel velocity
            profs(:, 1) = profs(:, 1)**2.0_wp                       ! Wall-normal velocity

        case (PROFILE_GAUSSIAN_SURFACE)
            ! Zero wall-normal derivative of wall-parallel velocity for free-slip and potentialvelocity mode, f=f*tanh
            if (any([BCS_DD, BCS_DN] == flag_wall)) then            ! no-slip at jmin
                do j = 1, jmax
                    yr = 0.5_wp*(y(j) - y(1))/IniK%thick
                    profs(j, 2) = profs(j, 2)*tanh(yr)**2 - &       ! Wall-parallel velocity
                                  profs(j, 1)*tanh(yr)/cosh(yr)**2/IniK%thick
                    profs(j, 1) = profs(j, 1)*tanh(yr)**2           ! Wall-normal velocity
                end do
            end if

            if (any([BCS_DD, BCS_ND] == flag_wall)) then            ! no-slip at jmax
                do j = 1, jmax
                    yr = 0.5_wp*(y(jmax) - y(j))/IniK%thick
                    profs(j, 2) = profs(j, 2)*tanh(yr)**2 + &       ! Wall-parallel velocity
                                  profs(j, 1)*tanh(yr)/cosh(yr)**2/IniK%thick
                    profs(j, 1) = profs(j, 1)*tanh(yr)**2           ! Wall-normal velocity
                end do
            end if

        end select

        return
    end subroutine FLOW_SHAPE

    ! ###################################################################
    subroutine FLOW_NORMALIZE(u, v, w)
        real(wp), dimension(imax, jmax, kmax) :: u, v, w

        ! -------------------------------------------------------------------
        real(wp) dummy, amplify

        ! ###################################################################
        amplify = 0.0_wp                                      ! Maximum across the layer
        do j = 1, jmax
            dummy = AVG1V2D(imax, jmax, kmax, j, 2, u) + AVG1V2D(imax, jmax, kmax, j, 2, v) + AVG1V2D(imax, jmax, kmax, j, 2, w)
            amplify = max(dummy, amplify)
        end do
        amplify = 0.5_wp*amplify

        amplify = sqrt(norm_ini_u/amplify)           ! Scaling factor to normalize to maximum TKE

        u = u*amplify
        v = v*amplify
        w = w*amplify

        return
    end subroutine FLOW_NORMALIZE

    ! ###################################################################
    !# Perturbation of the thermodynamic fields by a displacement of the reference center plane.
    !# Array s enters with the scalar total field, including fluctuations.
    !# Only used in compressible formulation
    !# Together discrete and broadband in one procedure
    !########################################################################
    subroutine DENSITY_FLUCTUATION(s, p, rho, T, h)
        use Thermodynamics, only: imixture

        real(wp), dimension(imax, jmax, kmax) :: T, h, rho, p
        real(wp), dimension(imax, jmax, kmax, *) :: s

        ! -------------------------------------------------------------------
        real(wp) dummy
        real(wp) xcenter, amplify, params(0)
        type(profiles_dt) :: prof_loc

        ! ###################################################################

#ifdef USE_MPI
        idsp = ims_offset_i; kdsp = ims_offset_k
#else
        idsp = 0; kdsp = 0
#endif

        ! ###################################################################
        ! Center plane displacement
        ! ###################################################################
#define disp(i,k) p_wrk2d(i,k,1)

        disp(:, :) = 0.0_wp

        select case (flag_t)
        case (PERT_BROADBAND)
            call IO_Read_Fields('scal.rand', imax, 1, kmax, itime, 1, 1, disp(:, :), params)
            dummy = AVG1V2D(imax, 1, kmax, 1, 1, disp(:, :))     ! remove mean
            disp(:, :) = disp(:, :) - dummy

        case (PERT_DISCRETE)
            wx_1 = 2.0_wp*pi_wp/g(1)%scale ! Fundamental wavelengths
            wz_1 = 2.0_wp*pi_wp/g(3)%scale

            do im = 1, fp%size
                wx = real(fp%modex(im), wp)*wx_1
                wz = real(fp%modez(im), wp)*wz_1

                do k = 1, kmax
                    disp(:, k) = disp(:, k) + fp%amplitude(im)*cos(wx*x(idsp + 1:idsp + imax) + fp%phasex(im)) &
                                 *cos(wz*z(kdsp + k) + fp%phasez(im))
                end do

            end do

        end select

        ! -------------------------------------------------------------------
        ! Modulation
        ! -------------------------------------------------------------------
        if (fp%type == PROFILE_GAUSSIAN .and. fp%parameters(1) > 0.0_wp) then
            do k = 1, kmax
                do i = 1, imax
                    xcenter = x(i) - g(1)%scale*0.5_wp - x(1)
                    amplify = exp(-0.5_wp*(xcenter/fp%parameters(1))**2)
                    disp(i, k) = disp(i, k)*amplify
                end do
            end do
        end if

        ! ###################################################################
        ! Perturbation in the thermodynamic fields
        ! ###################################################################
        if (tbg%type /= PROFILE_NONE) then
            prof_loc = tbg
            do k = 1, kmax
                do i = 1, imax
                    prof_loc%ymean = tbg%ymean + disp(i, k)
                    prof_loc%delta = tbg%delta + (tbg%uslope - tbg%lslope)*disp(i, k)*g(2)%scale
                    prof_loc%mean = tbg%mean + 0.5_wp*(tbg%uslope + tbg%lslope)*disp(i, k)*g(2)%scale
                    do j = 1, jmax
                        T(i, j, k) = Profiles_Calculate(prof_loc, y(j))
                    end do
                end do
            end do

            if (imixture == MIXT_TYPE_AIRWATER) then
                call THERMO_AIRWATER_PT(imax*jmax*kmax, s, p, T)
            end if

            call THERMO_THERMAL_DENSITY(imax*jmax*kmax, s, p, T, rho)

        end if

        if (hbg%type /= PROFILE_NONE) then
            prof_loc = hbg
            do k = 1, kmax
                do i = 1, imax
                    prof_loc%ymean = hbg%ymean + disp(i, k)
                    prof_loc%delta = hbg%delta + (hbg%uslope - hbg%lslope)*disp(i, k)*g(2)%scale
                    prof_loc%mean = hbg%mean + 0.5_wp*(hbg%uslope + hbg%lslope)*disp(i, k)*g(2)%scale
                    do j = 1, jmax
                        h(i, j, k) = Profiles_Calculate(prof_loc, y(j))
                    end do
                end do
            end do

            if (imixture == MIXT_TYPE_AIRWATER) then
                call THERMO_AIRWATER_PH_RE(imax*jmax*kmax, s, p, h, T)
            end if

            call THERMO_THERMAL_DENSITY(imax*jmax*kmax, s, p, T, rho)

        end if

        return
    end subroutine DENSITY_FLUCTUATION

    ! ###################################################################
    !# solve Poisson equation for p', nabla^2 p' = d/dx_i d/dx_j (rho_0 u_i u_j),
    !# assuming p/rho^\gamma0 constant(Homentropic conditions)
    ! ###################################################################
    subroutine PRESSURE_FLUCTUATION(u, v, w, rho, p, pprime, txc1, txc2, txc3, txc4)
        use Thermodynamics, only: gama0

        real(wp), dimension(imax, jmax, kmax), intent(in) :: u, v, w
        real(wp), dimension(imax, jmax, kmax), intent(inout) :: rho, p, pprime
        real(wp), dimension(imax, jmax, kmax), intent(inout) :: txc1, txc2, txc3, txc4

        ! -------------------------------------------------------------------
        integer(wi) bcs(2, 2)

        real(wp), allocatable :: bcs_hb(:), bcs_ht(:)

        ! ###################################################################
        ! Calculate RHS d/dx_i d/dx_j (u_i u_j), stored in txc4

        ! terms with u
        txc1 = rho*u*u; txc2 = rho*u*v; txc3 = rho*u*w
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc3, txc4)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc2, txc3)
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc1, txc2)
        txc2 = 2.0_wp*(txc4 + txc3) + txc2

        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), txc2, txc4)

        ! terms with v
        txc1 = rho*v*v; txc2 = rho*v*w
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc2, txc3)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc1, txc2)
        txc2 = txc2 + 2.0_wp*txc3

        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc2, txc1)
        txc4 = txc4 + txc1

        ! terms with w
        txc1 = rho*w*w
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc1, txc2)

        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), txc2, txc1)
        txc4 = txc4 + txc1

        ! Solve Poisson equation; pprime contains fluctuating p' (BCs are equal to zero!)
        if (g(1)%periodic .and. g(3)%periodic) then ! Doubly periodic in xOz
            pprime = -txc4          ! change of forcing term sign

            allocate (bcs_hb(imax*kmax), bcs_ht(imax*kmax))
            bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
            call OPR_Poisson(imax, jmax, kmax, 0, pprime, txc1, txc2, bcs_hb, bcs_ht)
        else                                      ! General treatment
            ! Undevelop
        end if

        ! An amplification factor norm_ini_p is allowed as in previous versions
        rho = (norm_ini_p*pprime/p/gama0 + 1.0_wp)*rho  ! isentropic relation p'/p = \gamma \rho'/rho
        p = norm_ini_p*pprime + p

        return
    end subroutine PRESSURE_FLUCTUATION

end module FLOW_LOCAL
