#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!#
!# Calculating RHS forcings at the inflow plane in spatially evolving cases
!#
!########################################################################
module BOUNDARY_INFLOW
    use TLAB_TYPES, only: filter_dt, grid_dt, discrete_dt
    use TLAB_CONSTANTS, only: efile, lfile, wp, wi
#ifdef TRACE_ON
    use TLAB_CONSTANTS, only: tfile
#endif
    use TLAB_VARS, only: imax, jmax, kmax, inb_flow, inb_scal, inb_flow_array, inb_scal_array, flow_on, scal_on
    use TLAB_VARS, only: imode_eqns, itransport
    use TLAB_VARS, only: g, qbg
    use TLAB_VARS, only: rtime, itime
    use TLAB_VARS, only: visc, damkohler
    use TLAB_ARRAYS, only: wrk1d, wrk2d, wrk3d
    use TLAB_PROCS
    use Thermodynamics, only: imixture
    use THERMO_THERMAL
    use THERMO_CALORIC
    use THERMO_AIRWATER
    use THERMO_ANELASTIC
    use IO_FIELDS
    use OPR_FILTERS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_i, ims_npro_k
    use TLAB_MPI_VARS, only: ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLAB_MPI_VARS, only: ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLAB_MPI_VARS, only: ims_offset_k
    use TLAB_MPI_PROCS
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
        integer(wi) is, itimetmp, bcs(2, 2)
        integer(wi) joffset, jglobal, j, iwrk_size
        real(wp) tolerance, dy
        real(wp) visctmp, rtimetmp
        character*32 fname, sname, str
        character*128 line

#ifdef USE_MPI
        integer(wi) isize_loc, id
#endif

        ! ###################################################################
#ifdef TRACE_ON
        call TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_INIT')
#endif

#ifdef USE_MPI
        ! I/O routines not yet developed for this particular case
        if (ims_npro_i > 1) then
            call TLAB_WRITE_ASCII(efile, 'BOUNDARY_INIT. I/O routines undeveloped.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
#endif

        if (.not. allocated(x_inf)) allocate (x_inf(g_inf(1)%size, g_inf(1)%inb_grid))
        if (.not. allocated(y_inf)) allocate (y_inf(g_inf(2)%size, g_inf(2)%inb_grid))
        if (.not. allocated(z_inf)) allocate (z_inf(g_inf(3)%size, g_inf(3)%inb_grid))
        if (.not. allocated(q_inf)) allocate (q_inf(g_inf(1)%size, g_inf(2)%size, g_inf(3)%size, inb_flow_array))
        if (.not. allocated(s_inf)) allocate (s_inf(g_inf(1)%size, g_inf(2)%size, g_inf(3)%size, inb_scal_array))

        if (g_inf(1)%size > 1) then ! Inflow fields for spatial simulations
            if (.not. associated(g_inf(1)%nodes)) &
                call IO_READ_GRID('grid.inf', g_inf(1)%size, g_inf(2)%size, g_inf(3)%size, &
                                  g_inf(1)%scale, g_inf(2)%scale, g_inf(3)%scale, x_inf, y_inf, z_inf)
            call FDM_INITIALIZE(x_inf, g_inf(1), wrk1d)
            if (.not. associated(g_inf(2)%nodes)) g_inf(2)%nodes => y_inf(:, 1)
            if (.not. associated(g_inf(3)%nodes)) g_inf(3)%nodes => z_inf(:, 1)
        end if

        ! #######################################################################
        ! Definining types for parallel mode
        ! #######################################################################
#ifdef USE_MPI
        if (FilterInflow(1)%type /= DNS_FILTER_NONE) then !  Required for inflow explicit filter
            call TLAB_WRITE_ASCII(lfile, 'Initialize MPI types for inflow filter.')
            id = TLAB_MPI_K_INFLOW
            isize_loc = FilterInflow(1)%size*FilterInflow(2)%size
            call TLAB_MPI_TYPE_K(ims_npro_k, kmax, isize_loc, 1, 1, 1, 1, &
                                 ims_size_k(id), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
            FilterInflow(3)%mpitype = id
        end if
#endif

        iwrk_size = g_inf(1)%size*g_inf(2)%size*kmax
        if (imax*jmax*kmax < iwrk_size) then
            call TLAB_WRITE_ASCII(efile, 'BOUNDARY_INFLOW_INIT. Not enough space in array txc.')
            call TLAB_STOP(DNS_ERROR_WRKSIZE)
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
                    call TLAB_WRITE_ASCII(efile, 'BOUNDARY_INFLOW. Inflow domain does not match.')
                    call TLAB_STOP(DNS_ERROR_INFLOWDOMAIN)
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
                call TLAB_WRITE_ASCII(lfile, line)
            end if

            rtimetmp = rtime
            itimetmp = itime
            visctmp = visc
            call IO_READ_FIELDS(fname, IO_FLOW, g_inf(1)%size, g_inf(2)%size, kmax, inb_flow, 0, q_inf)
            call IO_READ_FIELDS(sname, IO_SCAL, g_inf(1)%size, g_inf(2)%size, kmax, inb_scal, 0, s_inf)
            rtime = rtimetmp
            itime = itimetmp
            visc = visctmp

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
        call TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_INIT')
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
        call TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_BROADBAND')
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
        call TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_BROADBAND')
#endif
        return
    end subroutine BOUNDARY_INFLOW_BROADBAND

    !########################################################################
    !########################################################################
    subroutine BOUNDARY_INFLOW_DISCRETE(etime, inf_rhs)
        use TLAB_TYPES, only: profiles_dt
        use TLAB_CONSTANTS, only: pi_wp
        use PROFILES

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
        call TLAB_WRITE_ASCII(tfile, 'ENTERING BOUNDARY_INFLOW_DISCRETE')
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

        prof_loc%type = PROFILE_GAUSSIAN
        prof_loc%thick = fp%parameters(1)
        prof_loc%delta = 1.0_wp
        prof_loc%mean = 0.0_wp
        prof_loc%lslope = 0.0_wp
        prof_loc%uslope = 0.0_wp
        prof_loc%parameters = 0.0_wp

        ! ###################################################################
        ! Shape function
        ! ###################################################################
        select case (fp%type)
        case (PROFILE_GAUSSIAN)
            prof_loc%ymean = qbg(1)%ymean
            do j = 1, jmax
                yr = y(j) - prof_loc%ymean
                wrk1d(j, 1) = PROFILES_CALCULATE(prof_loc, y(j))
                wrk1d(j, 2) = yr/(prof_loc%thick**2)*wrk1d(j, 1) ! Derivative of f
            end do

        case (PROFILE_GAUSSIAN_SYM, PROFILE_GAUSSIAN_ANTISYM)
            prof_loc%ymean = qbg(1)%ymean - 0.5_wp*qbg(1)%diam
            do j = 1, jmax
                yr = y(j) - prof_loc%ymean
                wrk1d(j, 1) = PROFILES_CALCULATE(prof_loc, y(j))
                wrk1d(j, 2) = -yr/(prof_loc%thick**2)*wrk1d(j, 1)
            end do

            prof_loc%ymean = qbg(1)%ymean + 0.5_wp*qbg(1)%diam
            if (fp%type == PROFILE_GAUSSIAN_ANTISYM) then; factorx = -1.0_wp ! varicose
            elseif (fp%type == PROFILE_GAUSSIAN_SYM) then; factorx = 1.0_wp ! Sinuous
            end if
            do j = 1, jmax
                yr = y(j) - prof_loc%ymean
                dummy = PROFILES_CALCULATE(prof_loc, y(j))
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
        call TLAB_WRITE_ASCII(tfile, 'LEAVING BOUNDARY_INFLOW_DISCRETE')
#endif

        return
    end subroutine BOUNDARY_INFLOW_DISCRETE

    !########################################################################
    !########################################################################
    ! Filter
    ! This should be integrated into the inflow buffer, as the filter contribution
    ! BufferFilter should then be a block in tlab.ini as [Filter], which is read in io_read_global.

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
        call TLAB_WRITE_ASCII(efile, 'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
        ! FilterInflow needs to be initiliazed
        call TLAB_STOP(DNS_ERROR_UNDEVELOP)

        ! Define pointers
        if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
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

        if (imode_eqns == DNS_EQNS_TOTAL .or. imode_eqns == DNS_EQNS_INTERNAL) then
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
        if (any([DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC] == imode_eqns)) then
            if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
                call THERMO_ANELASTIC_PH(imax, jmax, kmax, s(1, 1, 1, 2), s(1, 1, 1, 1))

            else if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                call THERMO_AIRWATER_LINEAR(imax*jmax*kmax, s, s(1, 1, 1, inb_scal_array))

            end if

        else
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

        end if

        return
    end subroutine BOUNDARY_INFLOW_FILTER

end module BOUNDARY_INFLOW
