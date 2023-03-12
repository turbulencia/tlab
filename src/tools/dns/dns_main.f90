#include "dns_const.h"

program DNS

    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_PROCS
#endif
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS
    use DNS_LOCAL
    use DNS_ARRAYS
    use TIME
    use DNS_TOWER
    use IBM_VARS
    use PLANES
    use BOUNDARY_INFLOW
    use BOUNDARY_BUFFER
    use BOUNDARY_BCS
    use STATISTICS
    use PARTICLE_TRAJECTORIES
    use AVG_SCAL_ZT
    use IO_FIELDS
    use OPR_FILTERS
    use OPR_FOURIER
    implicit none
    save

    ! -------------------------------------------------------------------
    character(len=32) fname, str
    integer ig
    integer, parameter :: i0 = 0, i1 = 1

    ! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call PARTICLE_READ_GLOBAL(ifile)
    call DNS_READ_LOCAL(ifile)
    if (imode_ibm == 1) then
        call IBM_READ_INI(ifile)
        call IBM_READ_CONSISTENCY_CHECK(imode_rhs, BcsFlowJmin%type(:),    BcsFlowJmax%type(:), &
                                                   BcsScalJmin%type(:),    BcsScalJmax%type(:), &
                                                   BcsScalJmin%SfcType(:), BcsScalJmax%SfcType(:))
    end if
#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#ifdef USE_PSFFT
    if (imode_rhs == EQNS_RHS_NONBLOCKING) call DNS_NB3DFFT_INITIALIZE
#endif
#endif

    ! #######################################################################
    ! Initialize memory space and grid data
    ! #######################################################################
    call TLAB_ALLOCATE(__FILE__)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE()

    call TLAB_ALLOCATE_ARRAY_DOUBLE(__FILE__, hq, [isize_field, inb_flow], 'flow-rhs')
    call TLAB_ALLOCATE_ARRAY_DOUBLE(__FILE__, hs, [isize_field, inb_scal], 'scal-rhs')

    call PARTICLE_ALLOCATE(__FILE__)
    call TLAB_ALLOCATE_ARRAY_DOUBLE(__FILE__, l_hq, [isize_part, inb_part], 'part-rhs')

    call STATISTICS_INITIALIZE()

    call PLANES_INITIALIZE()

    if (use_tower) then
        call DNS_TOWER_INITIALIZE(tower_stride)
    end if

    if (imode_ibm == 1) then
        call IBM_ALLOCATE(__FILE__)
    end if

    ! ###################################################################
    ! Initialize operators and reference data
    ! ###################################################################
    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig))
        call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
        call OPR_FILTER_INITIALIZE(g(ig), PressureFilter(ig))
    end do

    if (fourier_on) then
        call OPR_FOURIER_INITIALIZE()
    end if

    call OPR_CHECK()

    ! ###################################################################
    ! Initialize fields
    ! ###################################################################
    itime = nitera_first

    visc_stop = visc ! Value read in ifile

    if (scal_on) then
        write (fname, *) nitera_first; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, 0, s)
    end if

    write (fname, *) nitera_first; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
    call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, inb_flow, 0, q)

    call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)  ! Initialize diagnostic thermodynamic quantities

    if (part%type /= PART_TYPE_NONE) then
        write (fname, *) nitera_first; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
        call IO_READ_PARTICLE(fname, l_g, l_q)
        call PARTICLE_INITIALIZE()

        if (imode_traj /= TRAJ_TYPE_NONE) then
            call PARTICLE_TRAJECTORIES_INITIALIZE()
        end if

    end if

    if (imode_sim == DNS_MODE_SPATIAL .and. nitera_stats_spa > 0) then
        write (fname, *) nitera_first; fname = 'st'//trim(adjustl(fname))
        call IO_READ_AVG_SPATIAL(fname, mean_flow, mean_scal)
    end if

    ! ###################################################################
    ! Initialize change in viscosity
    ! ###################################################################
    flag_viscosity = .false.
    if (visc /= visc_stop) then
        write (str, *) visc
        call TLAB_WRITE_ASCII(lfile, 'Changing original viscosity '//trim(adjustl(str))//' to new value.')
        if (visc_time > 0.0_wp) then
            visc_rate = (visc_stop - visc)/visc_time
            visc_time = rtime + visc_time                 ! Stop when this time is reached
            flag_viscosity = .true.
        else
            visc = visc_stop
        end if
    end if

    ! ###################################################################
    ! Initialize data for boundary conditions
    ! ###################################################################
    call BOUNDARY_BUFFER_INITIALIZE(q, s, txc)

    call BOUNDARY_BCS_INITIALIZE()

    if (imode_sim == DNS_MODE_SPATIAL) then
        call BOUNDARY_INFLOW_INITIALIZE(rtime, txc)
    end if

    ! ###################################################################
    ! Initialize IBM
    ! ###################################################################
    if (imode_ibm == 1) then
        call IBM_INITIALIZE_GEOMETRY(txc, wrk3d)
        call IBM_BCS_FIELD_COMBINED(i0, q)
        if (scal_on) call IBM_INITIALIZE_SCAL(i1, s)
    end if

    ! ###################################################################
    ! Check
    ! ###################################################################
    logs_data(1) = 0 ! Status
    call DNS_BOUNDS_CONTROL()
    call DNS_BOUNDS_LIMIT()

    ! ###################################################################
    ! Initialize time marching scheme
    ! ###################################################################
    call TIME_INITIALIZE()
    call TIME_COURANT()

    ! ###################################################################
    ! Check-pointing: Initialize logfiles
    ! ###################################################################
    call DNS_LOGS_INITIALIZE() ! headers
    call DNS_LOGS() ! first line

    ! ###################################################################
    ! Do simulation: Integrate equations
    ! ###################################################################
    itime = nitera_first

    write (str, *) itime
    call TLAB_WRITE_ASCII(lfile, 'Starting time integration at It'//trim(adjustl(str))//'.')

    do
        if (itime >= nitera_last) exit
        if (int(logs_data(1)) /= 0) exit

        call TIME_RUNGEKUTTA()

        itime = itime + 1
        rtime = rtime + dtime

        if (mod(itime - nitera_first, nitera_filter) == 0) then
            call DNS_FILTER()
            if (imode_ibm == 1) then
                call IBM_BCS_FIELD_COMBINED(i0, q) ! apply IBM BCs
                if (scal_on) call IBM_INITIALIZE_SCAL(i0, s)
            end if
        end if

        if (flag_viscosity) then                ! Change viscosity if necessary
            visc = visc + visc_rate*dtime
            if (rtime > visc_time) then
                visc = visc_stop                ! Fix new value without any roundoff
                flag_viscosity = .false.
            end if
        end if

        call TIME_COURANT()

        ! -------------------------------------------------------------------
        ! The rest: Logging, postprocessing and check-pointing
        ! -------------------------------------------------------------------
        call DNS_BOUNDS_CONTROL()
        if (mod(itime - nitera_first, nitera_log) == 0 .or. int(logs_data(1)) /= 0) then
            call DNS_LOGS()
        end if

        if (use_tower) then
            call DNS_TOWER_ACCUMULATE(q, 1, wrk1d)
            call DNS_TOWER_ACCUMULATE(s, 2, wrk1d)
        end if

        if (imode_traj /= TRAJ_TYPE_NONE) then
            call PARTICLE_TRAJECTORIES_ACCUMULATE()
        end if

        if (mod(itime - nitera_first, nitera_stats_spa) == 0) then  ! Accumulate statistics in spatially evolving cases
            if (flow_on) call AVG_FLOW_ZT_REDUCE(q, hq, txc, mean_flow)
            if (scal_on) call AVG_SCAL_ZT_REDUCE(q, s, hq, txc, mean_scal)
        end if

        if (mod(itime - nitera_first, nitera_stats) == 0) then      ! Calculate statistics
            if (imode_sim == DNS_MODE_TEMPORAL) call STATISTICS_TEMPORAL()
            if (imode_sim == DNS_MODE_SPATIAL) call STATISTICS_SPATIAL()
        end if

        if (mod(itime - nitera_first, nitera_save) == 0 .or. &      ! Check-pointing: Save restart files
            itime == nitera_last .or. int(logs_data(1)) /= 0 .or. & ! Secure that one restart file is saved 
            wall_time > nruntime_sec) then                          ! If max runtime of the code is reached 

            if (flow_on) then
                write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
                call IO_WRITE_FIELDS(fname, IO_FLOW, imax, jmax, kmax, inb_flow, q)
            end if

            if (scal_on) then
                write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
                call IO_WRITE_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, s)
            end if

            if (use_tower) then
                call DNS_TOWER_WRITE(wrk3d)
            end if

            if (part%type /= PART_TYPE_NONE) then
                write (fname, *) itime; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
                call IO_WRITE_PARTICLE(fname, l_g, l_q)
                if (imode_traj /= TRAJ_TYPE_NONE) then
                    write (fname, *) itime; fname = trim(adjustl(tag_traj))//trim(adjustl(fname))
                    call PARTICLE_TRAJECTORIES_WRITE(fname)
                end if
            end if

            if (imode_sim == DNS_MODE_SPATIAL .and. nitera_stats_spa > 0) then ! Spatial; running averages
                write (fname, *) itime; fname = 'st'//trim(adjustl(fname))
                call IO_WRITE_AVG_SPATIAL(fname, mean_flow, mean_scal)
            end if

        end if

        if (mod(itime - nitera_first, nitera_pln) == 0) then
            call PLANES_SAVE()
        end if

        if (wall_time > nruntime_sec) then 
            write (str, *) wall_time
            ! write to efile so that job is not resubmitted
            call TLAB_WRITE_ASCII(efile, 'Maximum walltime of '//trim(adjustl(str))//' seconds is reached.')
            exit
        end if
        
    end do

    ! ###################################################################
    call TLAB_STOP(int(logs_data(1)))

contains
!########################################################################
! Create headers or dns.out file
!
!# logs_data01 State (>0 if error)
!#
!# logs_data02 Maximum CFL number
!# logs_data03 Maximum diffusion number
!# logs_data04 Maximum source number
!#
!# logs_data05 Minimum pressure
!# logs_data06 Maximum pressure
!# logs_data07 Minimum density
!# logs_data08 Maximum density
!# logs_data09 NEWTONRAPHSON_ERROR
!#
!# logs_data10 Minimum dilatation
!# logs_data11 Maximum dilatation
!########################################################################
    subroutine DNS_LOGS_INITIALIZE()
        use THERMO_VARS, only: imixture

        integer ip
        character(len=256) line1

        line1 = '#'; ip = 1
        line1 = line1(1:ip)//' '//' Itn.'; ip = ip + 1 + 7
        line1 = line1(1:ip)//' '//' time'; ip = ip + 1 + 13
        line1 = line1(1:ip)//' '//' dt'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' CFL#'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' D#'; ip = ip + 1 + 10
        line1 = line1(1:ip)//' '//' visc'; ip = ip + 1 + 10

        select case (imode_eqns)
        case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
            line1 = line1(1:ip)//' '//' DilMin'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' DilMax'; ip = ip + 1 + 13

        case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
            line1 = line1(1:ip)//' '//' PMin'; ip = ip + 1 + 10
            line1 = line1(1:ip)//' '//' PMax'; ip = ip + 1 + 10
            line1 = line1(1:ip)//' '//' RMin'; ip = ip + 1 + 10
            line1 = line1(1:ip)//' '//' RMax'; ip = ip + 1 + 10

        end select

        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
            line1 = line1(1:ip)//' '//' NewtonRs'; ip = ip + 1 + 10
        end if

        line1 = line1(1:ip - 1)//'#'
        call TLAB_WRITE_ASCII(ofile, repeat('#', len_trim(line1)))
        call TLAB_WRITE_ASCII(ofile, trim(adjustl(line1)))
        call TLAB_WRITE_ASCII(ofile, repeat('#', len_trim(line1)))

    end subroutine DNS_LOGS_INITIALIZE

!########################################################################
!########################################################################
    subroutine DNS_LOGS()
        use THERMO_VARS, only: imixture, NEWTONRAPHSON_ERROR
#ifdef USE_MPI
        use MPI
        use TLAB_MPI_VARS, only: ims_err
        real(wp) dummy
#endif

        integer ip
        character(len=256) line1, line2

        write (line1, 100) int(logs_data(1)), itime, rtime, dtime, (logs_data(ip), ip=2, 3), visc
100     format((1x, I1), (1x, I7), (1x, E13.6), 4(1x, E10.3))

        select case (imode_eqns)
        case (DNS_EQNS_INCOMPRESSIBLE, DNS_EQNS_ANELASTIC)
            write (line2, 200) logs_data(10), logs_data(11)
200         format(2(1x, E13.6))
            line1 = trim(line1)//trim(line2)

        case (DNS_EQNS_INTERNAL, DNS_EQNS_TOTAL)
            write (line2, 300) (logs_data(ip), ip=5, 8)
300         format(4(1x, E10.3))
            line1 = trim(line1)//trim(line2)

        end select

        if (imixture == MIXT_TYPE_AIRWATER .and. damkohler(3) <= 0.0_wp) then
#ifdef USE_MPI
            call MPI_ALLREDUCE(NEWTONRAPHSON_ERROR, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
            NEWTONRAPHSON_ERROR = dummy
#endif
            write (line2, 400) NEWTONRAPHSON_ERROR
400         format(1(1x, E10.3))
            line1 = trim(line1)//trim(line2)
        end if

        call TLAB_WRITE_ASCII(ofile, trim(adjustl(line1)))

    end subroutine DNS_LOGS

end program DNS
