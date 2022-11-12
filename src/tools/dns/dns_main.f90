#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "DNS"

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
    implicit none
    save

    ! -------------------------------------------------------------------
    character*32 fname, str
    integer ig
    integer, parameter :: i0 = 0, i1 = 1, i2 = 2

    ! ###################################################################
    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call PARTICLE_READ_GLOBAL(ifile)
    call DNS_READ_LOCAL(ifile)
    if (imode_ibm == 1) then
        call IBM_READ_INI(ifile)
        call IBM_READ_CONSISTENCY_CHECK(imode_rhs, BcsFlowJmin%type(:), &
                                        BcsScalJmin%type(:), BcsScalJmax%type(:), &
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
    call TLAB_ALLOCATE(C_FILE_LOC)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE(wrk1d)

    call PARTICLE_ALLOCATE(C_FILE_LOC)

    call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, hq, inb_flow, isize_field, 'flow-rhs')
    call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, hs, inb_scal, isize_field, 'scal-rhs')
    if (imode_part /= PART_TYPE_NONE) then
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, l_hq, inb_part, isize_part, 'part-rhs')
        call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, l_comm, isize_l_comm, 'l_comm')
    end if

    call STATISTICS_INITIALIZE()

    call PLANES_INITIALIZE()

    if (use_tower) then
        call DNS_TOWER_INITIALIZE(tower_stride)
    end if

    if (imode_ibm == 1) then
        call IBM_ALLOCATE(C_FILE_LOC)
    end if

    ! ###################################################################
    ! Initialize operators and reference data
    ! ###################################################################
    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig), wrk1d)
        call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig), wrk1d)
    end do

    if (ifourier == 1) then
        call OPR_FOURIER_INITIALIZE(txc, wrk1d, wrk2d, wrk3d)
    end if

    call OPR_CHECK(imax, jmax, kmax, q, txc, wrk2d, wrk3d)

    ! ###################################################################
    ! Initialize fields
    ! ###################################################################
    itime = nitera_first

    visc_stop = visc ! Value read in ifile

    if (icalc_scal == 1) then
        write (fname, *) nitera_first; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
        call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, 0, s, wrk3d)
    end if

    write (fname, *) nitera_first; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
    call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, inb_flow, 0, q, wrk3d)

    call FI_DIAGNOSTIC(imax, jmax, kmax, q, s, wrk3d)  ! Initialize diagnostic thermodynamic quantities

    if (imode_part /= PART_TYPE_NONE) then
        write (fname, *) nitera_first; fname = trim(adjustl(tag_part))//trim(adjustl(fname))
        call IO_READ_PARTICLE(fname, l_g, l_q)
        call PARTICLE_INITIALIZE()

        if (imode_traj /= TRAJ_TYPE_NONE) then
            call PARTICLE_TRAJECTORIES_INITIALIZE(nitera_save, nitera_last)
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
    call BOUNDARY_BUFFER_INITIALIZE(q, s, txc, wrk3d)

    call BOUNDARY_BCS_INITIALIZE(wrk3d)

    if (imode_sim == DNS_MODE_SPATIAL) then
        call BOUNDARY_INFLOW_INITIALIZE(rtime, txc, wrk1d, wrk2d, wrk3d)
    end if

    ! ###################################################################
    ! Initialize IBM
    ! ###################################################################
    if (imode_ibm == 1) then
        call IBM_INITIALIZE_GEOMETRY(txc, wrk3d)
        call IBM_BCS_FIELD_COMBINED(i0, q)
        if (icalc_scal == 1) call IBM_INITIALIZE_SCAL(i1, s)
    end if

    ! ###################################################################
    ! Check
    ! ###################################################################
    logs_data(1) = 0 ! Status
    call DNS_CONTROL(0)

    ! ###################################################################
    ! Initialize time marching scheme
    ! ###################################################################
    call TIME_INITIALIZE()
    call TIME_COURANT(q, wrk3d)

    ! ###################################################################
    ! Check-pointing: Initialize logfiles
    ! ###################################################################
    call DNS_LOGS(i1) ! headers
    call DNS_LOGS(i2) ! first line

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
                if (icalc_scal == 1) call IBM_INITIALIZE_SCAL(i0, s)
            end if
        end if

        if (flag_viscosity) then                ! Change viscosity if necessary
            visc = visc + visc_rate*dtime
            if (rtime > visc_time) then
                visc = visc_stop                ! Fix new value without any roundoff
                flag_viscosity = .false.
            end if
        end if

        call TIME_COURANT(q, wrk3d)

        ! -------------------------------------------------------------------
        ! The rest: Logging, postprocessing and saving
        ! -------------------------------------------------------------------
        if (mod(itime - nitera_first, nitera_log) == 0 .or. int(logs_data(1)) /= 0) then
            call DNS_LOGS(i2)
        end if

        if (use_tower) then
            call DNS_TOWER_ACCUMULATE(q, 1, wrk1d)
            call DNS_TOWER_ACCUMULATE(s, 2, wrk1d)
        end if

        if (imode_traj /= TRAJ_TYPE_NONE) then
            call PARTICLE_TRAJECTORIES_ACCUMULATE(q, s, txc, l_g, l_q, l_hq, l_txc, l_comm, wrk2d, wrk3d)
        end if

        if (mod(itime - nitera_first, nitera_stats_spa) == 0) then  ! Accumulate statistics in spatially evolving cases
            if (icalc_flow == 1) call AVG_FLOW_ZT_REDUCE(q, hq, txc, mean_flow, wrk2d, wrk3d)
            if (icalc_scal == 1) call AVG_SCAL_ZT_REDUCE(q, s, hq, txc, mean_scal, wrk2d, wrk3d)
        end if

        if (mod(itime - nitera_first, nitera_stats) == 0) then      ! Calculate statistics
            if (imode_sim == DNS_MODE_TEMPORAL) call STATISTICS_TEMPORAL()
            if (imode_sim == DNS_MODE_SPATIAL) call STATISTICS_SPATIAL()
        end if

        if (mod(itime - nitera_first, nitera_save) == 0 .or. &      ! Save restart files
            itime == nitera_last .or. int(logs_data(1)) /= 0) then  ! Secure that one restart file is saved

            if (icalc_flow == 1) then
                write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
                call IO_WRITE_FIELDS(fname, IO_FLOW, imax, jmax, kmax, inb_flow, q, wrk3d)
            end if

            if (icalc_scal == 1) then
                write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
                call IO_WRITE_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, s, wrk3d)
            end if

            if (use_tower) then
                call DNS_TOWER_WRITE(wrk3d)
            end if

            if (imode_part /= PART_TYPE_NONE) then
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
            call PLANES_SAVE(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk1d, wrk2d, wrk3d)
        end if

    end do

    ! ###################################################################
    call TLAB_STOP(int(logs_data(1)))

end program DNS
