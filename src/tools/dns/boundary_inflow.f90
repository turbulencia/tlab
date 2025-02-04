#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!#
!# Calculating RHS forcings at the inflow plane in spatially evolving cases
!#
!########################################################################
module BOUNDARY_INFLOW
    use FDM, only: grid_dt
    use TLab_Constants, only: efile, lfile, wp, wi
#ifdef TRACE_ON
    use TLab_Constants, only: tfile
#endif
    use TLAB_VARS, only: imax, jmax, kmax, inb_flow, inb_scal, inb_flow_array, inb_scal_array
    use TLab_WorkFlow, only: flow_on, scal_on
    use NavierStokes, only: nse_eqns
    use FDM, only: g, FDM_Initialize
    use Timer, only: rtime, itime
    use NavierStokes, only: visc 
    use TLab_Background, only: qbg
    use TLab_Arrays, only: wrk1d, wrk2d, wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thermodynamics
    use Discrete, only: discrete_dt
    use THERMO_THERMAL
    use THERMO_CALORIC
    use THERMO_AIRWATER
    use THERMO_ANELASTIC
    use IO_FIELDS
    use OPR_FILTERS
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_VARS, only: ims_offset_k
    use TLabMPI_PROCS, only: TLabMPI_Panic
    use TLabMPI_Transpose, only: TLabMPI_Trp_TypeK_Create, ims_trp_plan_k
#endif
    use OPR_PARTIAL

    implicit none
    save
    private

    type(grid_dt), public :: g_inf(3)
    real(wp), allocatable :: x_inf(:, :), y_inf(:, :), z_inf(:, :)
    real(wp), allocatable :: q_inf(:, :, :, :), s_inf(:, :, :, :)

    integer(wi), public :: inflow_mode, inflow_ifield
    real(wp), public :: inflow_adapt
    type(filter_dt), public :: FilterInflow(3)
    ! integer(wi) :: FilterInflowStep
    type(discrete_dt), public :: fp ! Discrete forcing

    target :: x_inf, y_inf, z_inf

    public :: BOUNDARY_INFLOW_INITIALIZE
    public :: BOUNDARY_INFLOW_BROADBAND
    public :: BOUNDARY_INFLOW_DISCRETE
    public :: BOUNDARY_INFLOW_FILTER

contains
    !########################################################################
    !########################################################################
    !# Initializing inflow fields for broadband forcing case.
    subroutine BOUNDARY_INFLOW_INITIALIZE(etime, txc)
        real(wp) etime
        real(wp), intent(INOUT) :: txc(g_inf(1)%size, g_inf(2)%size, g_inf(3)%size)

        ! -------------------------------------------------------------------
        integer(wi) is, bcs(2, 2)
        integer(wi) joffset, jglobal, j, iwrk_size
        real(wp) tolerance, dy
        real(wp) params(0)
        character*32 fname, sname, str
        character*128 line

#ifdef USE_MPI
        integer(wi) isize_loc !, id
#endif

        ! ###################################################################
#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_INIT')
#endif

#ifdef USE_MPI
        ! I/O routines not yet developed for this particular case
        if (ims_npro_i > 1) then
            call TLab_Write_ASCII(efile, 'BOUNDARY_INIT. I/O routines undeveloped.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
#endif

        g_inf(:)%periodic = g(:)%periodic
        g_inf(:)%uniform = g(:)%uniform
        if (inflow_mode == 2) then
            g_inf(1)%periodic = .true.
            g_inf(1)%uniform = .true.
        end if
        if (g_inf(1)%size > 1 .and. .not. allocated(x_inf)) then ! Grid set only when entering the first time
            call IO_READ_GRID('grid.inf', g_inf(1)%size, g_inf(2)%size, g_inf(3)%size, &
                              g_inf(1)%scale, g_inf(2)%scale, g_inf(3)%scale, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3))
            call FDM_Initialize(x_inf, g_inf(1), wrk1d(:, 1), wrk1d(:, 4))
            call FDM_Initialize(y_inf, g_inf(2), wrk1d(:, 2), wrk1d(:, 4))
            call FDM_Initialize(z_inf, g_inf(3), wrk1d(:, 3), wrk1d(:, 4))
        end if

        if (.not. allocated(q_inf)) allocate (q_inf(g_inf(1)%size, g_inf(2)%size, g_inf(3)%size, inb_flow_array))
        if (.not. allocated(s_inf)) allocate (s_inf(g_inf(1)%size, g_inf(2)%size, g_inf(3)%size, inb_scal_array))

        ! #######################################################################
        ! Definining types for parallel mode
        ! #######################################################################
#ifdef USE_MPI
        if (FilterInflow(1)%type /= DNS_FILTER_NONE) then !  Required for inflow explicit filter
            ! call TLab_Write_ASCII(lfile, 'Initialize MPI types for inflow filter.')
            ! id = TLAB_MPI_TRP_K_INFLOW
            isize_loc = FilterInflow(1)%size*FilterInflow(2)%size
            ! call TLabMPI_TypeK_Create(ims_npro_k, kmax, isize_loc, 1, 1, 1, 1,  id)
            ims_trp_plan_k(TLAB_MPI_TRP_K_INFLOW) = TLabMPI_Trp_TypeK_Create(kmax, isize_loc, 1, 1, 1, 1, 'inflow filter.')
            FilterInflow(3)%mpitype = TLAB_MPI_TRP_K_INFLOW
        end if
#endif

        iwrk_size = g_inf(1)%size*g_inf(2)%size*kmax
        if (imax*jmax*kmax < iwrk_size) then
            call TLab_Write_ASCII(efile, 'BOUNDARY_INFLOW_INIT. Not enough space in array txc.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        ! ###################################################################
        if (inflow_mode == 2 .or. inflow_mode == 3 .or. inflow_mode == 4) then

            ! Checking the matching; we could move this outside...
            tolerance = 1.0e-10_wp
            joffset = (jmax - g_inf(2)%size)/2
            do j = 1, g_inf(2)%size
                jglobal = joffset + j
                dy = abs(g(2)%nodes(jglobal) - g_inf(2)%nodes(j))
                if (dy > tolerance) then
                    call TLab_Write_ASCII(efile, 'BOUNDARY_INFLOW. Inflow domain does not match.')
                    call TLab_Stop(DNS_ERROR_INFLOWDOMAIN)
                end if
            end do

            ! Reading fields
            fname = 'flow.inf'
            sname = 'scal.inf'
            inflow_ifield = int(qbg(1)%mean*etime/g_inf(1)%scale) + 1
            if (inflow_mode == 3) then
                write (str, *) inflow_ifield
                fname = trim(adjustl(fname))//trim(adjustl(str))
                sname = trim(adjustl(sname))//trim(adjustl(str))
                line = 'Reading InflowFile '//trim(adjustl(str))
                call TLab_Write_ASCII(lfile, line)
            end if

            call IO_READ_FIELDS(fname, g_inf(1)%size, g_inf(2)%size, kmax, itime, inb_flow, 0, q_inf, params)
            call IO_READ_FIELDS(sname, g_inf(1)%size, g_inf(2)%size, kmax, itime, inb_scal, 0, s_inf, params)

            ! array p contains the internal energy. Now we put in the pressure
            call THERMO_CALORIC_TEMPERATURE(g_inf(1)%size*g_inf(2)%size*kmax, s_inf, q_inf(1, 1, 1, 4), q_inf(1, 1, 1, 5), txc, wrk3d)
            call THERMO_THERMAL_PRESSURE(g_inf(1)%size*g_inf(2)%size*kmax, s_inf, q_inf(1, 1, 1, 5), txc, q_inf(1, 1, 1, 4))

            ! ###################################################################
            ! Performing the derivatives
            !
            ! Note that a field f_0(x) is convected with a velocity U along OX, i.e.
            ! f(x,t) = f_0(x-Ut), and therefore \partial f/\partial t = -U df_0/dx.
            ! ###################################################################
            bcs = 0

            if (flow_on) then
                do is = 1, inb_flow
                    call OPR_PARTIAL_X(OPR_P1, g_inf(1)%size, g_inf(2)%size, kmax, bcs, g_inf(1), q_inf(1, 1, 1, is), txc)
                    q_inf(:, :, :, is) = -txc(:, :, :)*qbg(1)%mean
                end do
            end if

            if (scal_on) then
                do is = 1, inb_scal
                    call OPR_PARTIAL_X(OPR_P1, g_inf(1)%size, g_inf(2)%size, kmax, bcs, g_inf(1), s_inf(1, 1, 1, is), txc)
                    s_inf(:, :, :, is) = -txc(:, :, :)*qbg(1)%mean
                end do
            end if

        end if

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_INIT')
#endif

        return
    end subroutine BOUNDARY_INFLOW_INITIALIZE

    !########################################################################
    !########################################################################
    subroutine BOUNDARY_INFLOW_BROADBAND(etime, inf_rhs, txc)
        real(wp) etime
        real(wp), intent(OUT) :: inf_rhs(jmax, kmax, inb_flow + inb_scal)
        real(wp), intent(INOUT) :: txc(*)

        ! -------------------------------------------------------------------
        real(wp) xaux, dx_loc, vmult
        integer(wi) joffset, jglobal, ileft, iright, j, k, is, ip
        real(wp) BSPLINES3P, BSPLINES3

        ! ###################################################################
#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_BROADBAND')
#endif

        ! Transient factor
        if (inflow_adapt > 0.0_wp .and. etime <= inflow_adapt) then
            vmult = etime/inflow_adapt
        else
            vmult = 1.0_wp
        end if

        ! check if we need to read again inflow data
        if (inflow_mode == 3 .and. int(qbg(1)%mean*etime/g_inf(1)%scale) + 1 /= inflow_ifield) then
            call BOUNDARY_INFLOW_INITIALIZE(etime, txc)
        end if

        ! ###################################################################
        ! Getting the position
        ! ###################################################################
        joffset = (jmax - g_inf(2)%size)/2

        xaux = qbg(1)%mean*etime
        ! Remove integral length scales of box
        xaux = xaux - int(xaux/g_inf(1)%scale)*g_inf(1)%scale
        ! Set distance from box initial length
        xaux = g_inf(1)%scale - xaux

        dx_loc = g_inf(1)%nodes(2) - g_inf(1)%nodes(1)
        ! Get left index
        ileft = int(xaux/dx_loc) + 1
        ! Check bounds
        if (ileft > g_inf(1)%size) then
            ileft = 1
        end if
        ! Set right index
        iright = ileft + 1
        ! Check bounds
        if (iright > g_inf(1)%size) then
            iright = 1
        end if
        ! Get relative distance from left point
        xaux = (xaux - (g_inf(1)%nodes(ileft) - g_inf(1)%nodes(1)))/dx_loc

        ! ###################################################################
        ! Sampling the information
        ! ###################################################################
        ! -------------------------------------------------------------------
        ! Periodic
        ! -------------------------------------------------------------------
        if (inflow_mode == 2) then
            do k = 1, kmax
                do j = 1, g_inf(2)%size
                    jglobal = joffset + j
                    do is = 1, inb_scal
                        inf_rhs(jglobal, k, is) = inf_rhs(jglobal, k, is) + vmult*BSPLINES3P(q_inf(1, j, k, is), g_inf(1)%size, ileft, xaux)
                    end do

                    if (scal_on) then
                        do is = 1, inb_scal
                            ip = inb_flow + is
                            inf_rhs(jglobal, k, ip) = inf_rhs(jglobal, k, ip) + vmult*BSPLINES3P(s_inf(1, j, k, is), g_inf(1)%size, ileft, xaux)
                        end do
                    end if

                end do
            end do

            ! -------------------------------------------------------------------
            ! Sequential
            ! -------------------------------------------------------------------
        else
            do k = 1, kmax
                do j = 1, g_inf(2)%size
                    jglobal = joffset + j
                    do is = 1, inb_flow
                        inf_rhs(jglobal, k, is) = inf_rhs(jglobal, k, is) + vmult*BSPLINES3(q_inf(1, j, k, is), g_inf(1)%size, ileft, xaux)
                    end do

                    if (scal_on) then
                        do is = 1, inb_scal
                            ip = inb_flow + is
                            inf_rhs(jglobal, k, ip) = inf_rhs(jglobal, k, ip) + vmult*BSPLINES3(s_inf(1, j, k, is), g_inf(1)%size, ileft, xaux)
                        end do
                    end if

                end do
            end do

        end if

        ! ###################################################################
        ! Filling the rest
        ! ###################################################################
        do j = 1, joffset
            inf_rhs(j, :, :) = inf_rhs(j, :, :) + 0.0_wp
        end do
        do j = jmax - joffset + 1, jmax
            inf_rhs(j, :, :) = inf_rhs(j, :, :) + 0.0_wp
        end do

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_BROADBAND')
#endif
        return
    end subroutine BOUNDARY_INFLOW_BROADBAND

    !########################################################################
    !########################################################################
    subroutine BOUNDARY_INFLOW_DISCRETE(etime, inf_rhs)
        use Profiles, only: profiles_dt, Profiles_Calculate, PROFILE_GAUSSIAN, PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_ANTISYM
        use TLab_Constants, only: pi_wp

        real(wp) etime
        real(wp), intent(OUT) :: inf_rhs(jmax, kmax, inb_flow + inb_scal)

        ! -------------------------------------------------------------------
        integer(wi) j, k, im, kdsp
        real(wp) wx, wz, wx_1, wz_1, xaux, vmult, factorx, factorz, dummy

        real(wp) yr
        type(profiles_dt) prof_loc

        real(wp), dimension(:), pointer :: y, z

        ! ###################################################################
#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_DISCRETE')
#endif

        ! Define pointers
        y => g(2)%nodes
        z => g(3)%nodes

#ifdef USE_MPI
        kdsp = ims_offset_k
#else
        kdsp = 0
#endif

        xaux = -qbg(1)%mean*etime

        prof_loc = profiles_dt(type=PROFILE_GAUSSIAN)
        prof_loc%thick = fp%parameters(1)

        ! ###################################################################
        ! Shape function
        ! ###################################################################
        select case (fp%type)
        case (PROFILE_GAUSSIAN)
            prof_loc%ymean = qbg(1)%ymean
            do j = 1, jmax
                yr = y(j) - prof_loc%ymean
                wrk1d(j, 1) = Profiles_Calculate(prof_loc, y(j))
                wrk1d(j, 2) = yr/(prof_loc%thick**2)*wrk1d(j, 1) ! Derivative of f
            end do

        case (PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_ANTISYM)
            prof_loc%ymean = qbg(1)%ymean - 0.5_wp*qbg(1)%diam
            do j = 1, jmax
                yr = y(j) - prof_loc%ymean
                wrk1d(j, 1) = Profiles_Calculate(prof_loc, y(j))
                wrk1d(j, 2) = -yr/(prof_loc%thick**2)*wrk1d(j, 1)
            end do

            prof_loc%ymean = qbg(1)%ymean + 0.5_wp*qbg(1)%diam
            if (fp%type == PROFILE_GAUSSIAN_ANTISYM) then; factorx = -1.0_wp ! varicose
            elseif (fp%type == PROFILE_GAUSSIAN_SYM) then; factorx = 1.0_wp ! Sinuous
            end if
            do j = 1, jmax
                yr = y(j) - prof_loc%ymean
                dummy = Profiles_Calculate(prof_loc, y(j))
                wrk1d(j, 1) = wrk1d(j, 1) + dummy
                wrk1d(j, 2) = wrk1d(j, 2) + yr/(prof_loc%thick**2)*dummy
            end do

        end select

        ! ###################################################################
        ! Fourier series
        ! ###################################################################
        wx_1 = 2.0_wp*pi_wp/fp%parameters(2) ! Fundamental wavelengths
        wz_1 = 2.0_wp*pi_wp/g(3)%scale

        wrk2d = 0.0_wp
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
                wrk2d(k, 2) = wrk2d(k, 2) + fp%amplitude(im)*cos(wx*xaux + fp%phasex(im))*cos(wz*z(kdsp + k) + fp%phasez(im))
                wrk2d(k, 1) = wrk2d(k, 1) + fp%amplitude(im)*sin(wx*xaux + fp%phasex(im))*cos(wz*z(kdsp + k) + fp%phasez(im))*factorx
                wrk2d(k, 3) = wrk2d(k, 3) + fp%amplitude(im)*cos(wx*xaux + fp%phasex(im))*sin(wz*z(kdsp + k) + fp%phasez(im))*factorz
            end do

        end do

        ! ###################################################################
        ! Forcing
        ! ###################################################################
        ! Transient factor
        if (inflow_adapt > 0.0_wp .and. etime <= inflow_adapt) then
            vmult = etime/inflow_adapt
        else
            vmult = 1.0_wp
        end if

        do k = 1, kmax
            do j = 1, jmax
                inf_rhs(j, k, 2) = inf_rhs(j, k, 2) - vmult*qbg(1)%mean*wrk2d(k, 1)*wrk1d(j, 2) ! u
                inf_rhs(j, k, 3) = inf_rhs(j, k, 3) - vmult*qbg(1)%mean*wrk2d(k, 2)*wrk1d(j, 1) ! v
                inf_rhs(j, k, 4) = inf_rhs(j, k, 4) - vmult*qbg(1)%mean*wrk2d(k, 3)*wrk1d(j, 2) ! w
            end do
        end do

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_DISCRETE')
#endif

        return
    end subroutine BOUNDARY_INFLOW_DISCRETE

    !########################################################################
    !########################################################################
    ! Filter
    ! This should be integrated into the inflow buffer, as the filter contribution
    ! BufferFilter should then be a block in tlab.ini as [Filter], which is read in TLab_Initialize_Parameters.

    subroutine BOUNDARY_INFLOW_FILTER(bcs_vi, bcs_vi_scal, q, s, txc)
        real(wp), dimension(imax, jmax, kmax, *), intent(INOUT) :: q, s
        real(wp), dimension(jmax, kmax, *), intent(IN) :: bcs_vi, bcs_vi_scal
        real(wp), dimension(imax*jmax*kmax, 2), intent(INOUT) :: txc

        target q

        ! -----------------------------------------------------------------------
        integer(wi) i, j, k, ip, iq, iq_loc(inb_flow), is
        integer(wi) j1, imx, jmx, ifltmx, jfltmx

        ! Pointers to existing allocated space
        real(wp), dimension(:, :, :), pointer :: e, rho, p, T, vis

        ! ###################################################################
        ! #######################################################################
        call TLab_Write_ASCII(efile, 'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
        ! FilterInflow needs to be initiliazed
        call TLab_Stop(DNS_ERROR_UNDEVELOP)

        ! Define pointers
        if (nse_eqns == DNS_EQNS_TOTAL .or. nse_eqns == DNS_EQNS_INTERNAL) then
            e => q(:, :, :, 4)
            rho => q(:, :, :, 5)
            p => q(:, :, :, 6)
            T => q(:, :, :, 7)

            if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) vis => q(:, :, :, 8)

        end if

        ! Define counters
        imx = FilterInflow(1)%size
        j1 = (jmax - FilterInflow(2)%size)/2 + 1
        jmx = (jmax + FilterInflow(2)%size)/2
        j1 = min(max(j1, 1), jmax)
        jmx = min(max(jmx, 1), jmax)

        ifltmx = imx - 1 + 1
        jfltmx = jmx - j1 + 1

        if (any([DNS_EQNS_TOTAL, DNS_EQNS_INTERNAL] == nse_eqns)) then
            iq_loc = (/5, 1, 2, 3, 6/) ! Filtered variables: rho, u,v,w, p
        else
            iq_loc = (/1, 2, 3/)
        end if

        ! #######################################################################
        do iq = 1, inb_flow

            ! -----------------------------------------------------------------------
            ! Remove mean field
            ! -----------------------------------------------------------------------
            ip = 1
            do k = 1, kmax
                do j = j1, jmx
                    do i = 1, imx
                        txc(ip, 1) = q(i, j, k, iq_loc(iq)) - bcs_vi(j, k, iq_loc(iq))
                        ip = ip + 1
                    end do
                end do
            end do

            ! -----------------------------------------------------------------------
            call OPR_FILTER(ifltmx, jfltmx, kmax, FilterInflow, txc(:, 1), txc(:, 2))

            ! -----------------------------------------------------------------------
            ! Add mean field
            ! -----------------------------------------------------------------------
            ip = 1
            do k = 1, kmax
                do j = j1, jmx
                    do i = 1, imx
                        q(i, j, k, iq_loc(iq)) = txc(ip, 1) + bcs_vi(j, k, iq_loc(iq))
                        ip = ip + 1
                    end do
                end do
            end do

        end do

        ! #######################################################################
        do is = 1, inb_scal

            ! -----------------------------------------------------------------------
            ! Remove mean field
            ! -----------------------------------------------------------------------
            ip = 1
            do k = 1, kmax
                do j = j1, jmx
                    do i = 1, imx
                        txc(ip, 1) = s(i, j, k, is) - bcs_vi_scal(j, k, is)
                        ip = ip + 1
                    end do
                end do
            end do

            ! -----------------------------------------------------------------------
            call OPR_FILTER(ifltmx, jfltmx, kmax, FilterInflow, txc(:, 1), txc(:, 2))

            ! -----------------------------------------------------------------------
            ! Add mean field
            ! -----------------------------------------------------------------------
            ip = 1
            do k = 1, kmax
                do j = j1, jmx
                    do i = 1, imx
                        s(i, j, k, is) = txc(ip, 1) + bcs_vi_scal(j, k, is)
                        ip = ip + 1
                    end do
                end do
            end do

        end do

        ! #######################################################################
        ! recalculation of diagnostic variables
        if (inb_flow_array > inb_flow .or. inb_scal_array > inb_scal) then
            select case (imode_thermo)
            case (THERMO_TYPE_ANELASTIC)
                if (imixture == MIXT_TYPE_AIRWATER) then                            ! Calculate liquid content q_l
                    call THERMO_ANELASTIC_PH(imax, jmax, kmax, s(1, 1, 1, 2), s(1, 1, 1, 1))

                end if

            case (THERMO_TYPE_LINEAR)
                if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then                     ! Calculate liquid content q_l
                    call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, 1, 1, inb_scal_array))

                end if

            case (THERMO_TYPE_COMPRESSIBLE)
                if (imixture == MIXT_TYPE_AIRWATER) then
                    call THERMO_AIRWATER_RP(imax*jmax*kmax, s, p, rho, T, wrk3d)
                else
                    call THERMO_THERMAL_TEMPERATURE(imax*jmax*kmax, s, p, rho, T)
                end if
                call THERMO_CALORIC_ENERGY(imax*jmax*kmax, s, T, e)

                ! This recalculation of T and p is made to make sure that the same numbers are
                ! obtained in statistics postprocessing as in the simulation; avg* files
                ! can then be compared with diff command.
                if (imixture == MIXT_TYPE_AIRWATER) then
                    call THERMO_CALORIC_TEMPERATURE(imax*jmax*kmax, s, e, rho, T, wrk3d)
                    call THERMO_THERMAL_PRESSURE(imax*jmax*kmax, s, rho, T, p)
                end if

                if (itransport == EQNS_TRANS_SUTHERLAND .or. itransport == EQNS_TRANS_POWERLAW) call THERMO_VISCOSITY(imax*jmax*kmax, T, vis)
                
            end select

        end if

        return
    end subroutine BOUNDARY_INFLOW_FILTER

end module BOUNDARY_INFLOW
