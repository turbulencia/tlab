#include "dns_const.h"
#include "dns_error.h"

program DNS

    use TLab_Constants, only: ifile, efile, wfile, lfile, gfile, tag_flow, tag_scal, tag_part, tag_traj
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
    use TLAB_VARS, only: imode_sim, fourier_on, scal_on, flow_on
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: itime
    use TLAB_VARS, only: damkohler
    use TLAB_VARS, only: FilterDomain, Dealiasing, PressureFilter
    use Tlab_Background, only: TLab_Initialize_Background, pbg, rbg
    use FDM, only: g, FDM_Initialize
    use TLab_Arrays
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Transpose_Initialize
#endif
    use Thermodynamics, only: Thermodynamics_Initialize_Parameters
    use Gravity, only: Gravity_Initialize
    use Radiation, only: Radiation_Initialize
    use Microphysics, only: Microphysics_Initialize
    use Chemistry, only: Chemistry_Initialize
    use SpecialForcing, only: SpecialForcing_Initialize
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
    use DNS_STATISTICS, only: DNS_STATISTICS_INITIALIZE, DNS_STATISTICS_SPATIAL, DNS_STATISTICS_TEMPORAL, mean_flow, mean_scal
    use ParticleTrajectories
    use AVG_SCAL_ZT
    use IO_FIELDS
    use OPR_ELLIPTIC
    use OPR_FILTERS
    use OPR_FOURIER
    use AVG_PHASE
    use Avg_Spatial, only: IO_READ_AVG_SPATIAL, IO_WRITE_AVG_SPATIAL
    implicit none
    save

    ! -------------------------------------------------------------------
    character(len=32) fname, str
    integer ig
    integer, parameter :: i0 = 0, i1 = 1

    ! ###################################################################
    call SYSTEM_CLOCK(start_clock)
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Transpose_Initialize(ifile)
#endif
    call Particle_Initialize_Parameters(ifile)
    call IBM_READ_INI(ifile)
    if (imode_ibm == 1) then
        call IBM_READ_CONSISTENCY_CHECK()
    end if

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call Chemistry_Initialize(ifile)

    call TLab_Consistency_Check()

    call DNS_READ_LOCAL(ifile)
#ifdef USE_PSFFT
    if (imode_rhs == EQNS_RHS_NONBLOCKING) call DNS_NB3DFFT_INITIALIZE
#endif

    ! #######################################################################
    ! Initialize memory space and grid data
    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3))
    call FDM_Initialize(x, g(1), wrk1d(:, 1), wrk1d(:, 4))
    call FDM_Initialize(y, g(2), wrk1d(:, 2), wrk1d(:, 4))
    call FDM_Initialize(z, g(3), wrk1d(:, 3), wrk1d(:, 4))

    call SpecialForcing_Initialize(ifile)

    call TLab_Initialize_Background(ifile)

    call TLab_Allocate_Real(__FILE__, hq, [isize_field, inb_flow], 'flow-rhs')
    call TLab_Allocate_Real(__FILE__, hs, [isize_field, inb_scal], 'scal-rhs')

    call ParticleTrajectories_Initialize(ifile)
    call Particle_Initialize_Memory(__FILE__)
    call TLab_Allocate_Real(__FILE__, l_hq, [isize_part, inb_part], 'part-rhs')

    call DNS_STATISTICS_INITIALIZE()

    call PLANES_INITIALIZE()

    if (PhAvg%active) then
        call AvgPhaseInitializeMemory(__FILE__, nitera_save)
    end if

    if (use_tower) then
        call DNS_TOWER_INITIALIZE(tower_stride)
    end if

    if (imode_ibm == 1) then
        call IBM_ALLOCATE(__FILE__)
    end if

    ! ###################################################################
    ! Initialize operators and reference data
    ! ###################################################################
    call OPR_Elliptic_Initialize(ifile)

    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig))
        call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
        call OPR_FILTER_INITIALIZE(g(ig), PressureFilter(ig))
    end do

    if (fourier_on) call OPR_FOURIER_INITIALIZE()

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
        call Particle_Initialize_Fields()
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
        call TLab_Write_ASCII(lfile, 'Changing original viscosity '//trim(adjustl(str))//' to new value.')
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
    if (bound_p%min < 0.0_wp) bound_p%min = pbg%mean*1.0e-6_wp
    if (bound_p%max < 0.0_wp) bound_p%max = pbg%mean*1.0e6_wp
    if (bound_r%min < 0.0_wp) bound_r%min = rbg%mean*1.0e-6_wp
    if (bound_r%max < 0.0_wp) bound_r%max = rbg%mean*1.0e6_wp

    logs_data(1) = 0; obs_data(1) = 0 ! Status
    call DNS_BOUNDS_CONTROL()
    call DNS_OBS_CONTROL()
    call DNS_BOUNDS_LIMIT()

    ! ###################################################################
    ! Initialize time marching scheme
    ! ###################################################################
    call TIME_INITIALIZE()
    call TIME_COURANT()

    ! ###################################################################
    ! Check-pointing: Initialize logfiles, write header & first line
    ! ###################################################################
    call DNS_LOGS_PATH_INITIALIZE()
    call DNS_LOGS_INITIALIZE()
    call DNS_LOGS()
    if (dns_obs_log /= OBS_TYPE_NONE) then
        call DNS_OBS_INITIALIZE()
        call DNS_OBS()
    end if

    ! ###################################################################
    ! Do simulation: Integrate equations
    ! ###################################################################
    itime = nitera_first

    write (str, *) itime
    call TLab_Write_ASCII(lfile, 'Starting time integration at It'//trim(adjustl(str))//'.')

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
        call DNS_OBS_CONTROL()
        if (mod(itime - nitera_first, nitera_log) == 0 .or. int(logs_data(1)) /= 0) then
            call DNS_LOGS()
            if (dns_obs_log /= OBS_TYPE_NONE) then
                call DNS_OBS()
            end if
        end if

        if (PhAvg%active) then
            if (mod(itime, PhAvg%stride) == 0) then
                call AvgPhaseSpace(wrk2d, inb_flow, itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, 1)
                call AvgPhaseSpace(wrk2d, inb_scal, itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, 2)
                ! Pressure is taken from the RHS subroutine
                ! call AvgPhaseSpace(wrk2d, 6       , itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride, 8)
                call AvgPhaseStress(q, itime/PhAvg%stride, nitera_first, nitera_save/PhAvg%stride)
                if (mod(itime - nitera_first, nitera_save) == 0) then
                    call IO_Write_AvgPhase(avg_planes, inb_flow, IO_FLOW, nitera_save, PhAvg%stride, avgu_name, 1, avg_flow)
                    call IO_Write_AvgPhase(avg_planes, inb_scal, IO_SCAL, nitera_save, PhAvg%stride, avgs_name, 2, avg_scal)
                    call IO_Write_AvgPhase(avg_planes, 1, IO_SCAL, nitera_save, PhAvg%stride, avgp_name, 4, avg_p)
                    call IO_Write_AvgPhase(avg_planes, 6, IO_FLOW, nitera_save, PhAvg%stride, avgstr_name, 8, avg_stress)

                    call AvgPhaseResetVariable()
                end if
            end if
        end if

        if (use_tower) then
            call DNS_TOWER_ACCUMULATE(q, 1, wrk1d)
            call DNS_TOWER_ACCUMULATE(s, 2, wrk1d)
        end if
        if (imode_traj /= TRAJ_TYPE_NONE) then
            call ParticleTrajectories_Accumulate()
        end if
        if (mod(itime - nitera_first, nitera_stats_spa) == 0) then  ! Accumulate statistics in spatially evolving cases
            if (flow_on) call AVG_FLOW_ZT_REDUCE(q, hq, txc, mean_flow)
            if (scal_on) call AVG_SCAL_ZT_REDUCE(q, s, hq, txc, mean_scal)
        end if
        if (mod(itime - nitera_first, nitera_stats) == 0) then      ! Calculate statistics
            if (imode_sim == DNS_MODE_TEMPORAL) call DNS_STATISTICS_TEMPORAL()
            if (imode_sim == DNS_MODE_SPATIAL) call DNS_STATISTICS_SPATIAL()
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
                    call ParticleTrajectories_Write(fname)
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
            call TLab_Write_ASCII(efile, 'Maximum walltime of '//trim(adjustl(str))//' seconds is reached.')
            exit
        end if
    end do

    ! ###################################################################
    call TLab_Stop(int(logs_data(1)))

contains

!########################################################################
!# Initialize path to write dns.out & tlab.logs
!########################################################################
    subroutine DNS_LOGS_PATH_INITIALIZE()

        integer env_status, path_len

        call get_environment_variable("DNS_LOGGER_PATH", logger_path, path_len, env_status, .true.)

        select case (env_status)
        case (-1)
            call TLab_Write_ASCII(efile, "DNS_START. The environment variable  is too long and cannot be handled in the foreseen array.")
            call TLab_Stop(DNS_ERROR_OPTION)
        case (0)
            if (.not. logger_path(path_len:path_len) == '/') then
                logger_path = trim(adjustl(logger_path))//'/'
            end if

        case (1:)
            logger_path = trim(adjustl(''))

        end select

    end subroutine DNS_LOGS_PATH_INITIALIZE

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
        use Thermodynamics, only: imixture

        integer ip
        character(len=256) line1

        ofile = trim(adjustl(logger_path))//trim(adjustl(ofile_base))
        ofile = trim(adjustl(ofile))

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
        call TLab_Write_ASCII(ofile, repeat('#', len_trim(line1)))
        call TLab_Write_ASCII(ofile, trim(adjustl(line1)))
        call TLab_Write_ASCII(ofile, repeat('#', len_trim(line1)))

    end subroutine DNS_LOGS_INITIALIZE

!########################################################################

    subroutine DNS_LOGS()
        use Thermodynamics, only: imixture, NEWTONRAPHSON_ERROR
#ifdef USE_MPI
        use MPI
        use TLabMPI_VARS, only: ims_err
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

        call TLab_Write_ASCII(ofile, trim(adjustl(line1)))

    end subroutine DNS_LOGS

!########################################################################
!# Create headers or dns.obs file
!########################################################################
    subroutine DNS_OBS_INITIALIZE()

        implicit none

        integer(wi) :: ip, is
        character(len=256) :: line1

        vfile = trim(adjustl(logger_path))//trim(adjustl(vfile_base))
        vfile = trim(adjustl(vfile))

        line1 = '#'; ip = 1
        line1 = line1(1:ip)//' '//' Itn.'; ip = ip + 1 + 7
        line1 = line1(1:ip)//' '//' time'; ip = ip + 1 + 13

        select case (dns_obs_log)
        case (OBS_TYPE_EKMAN)
            line1 = line1(1:ip)//' '//' u_bulk'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' w_bulk'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' u_y(1)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' w_y(1)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' alpha(1)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' alpha(ny)'; ip = ip + 1 + 13
            line1 = line1(1:ip)//' '//' entstrophy'; ip = ip + 1 + 13
            if (scal_on) then
                do is = 1, inb_scal
                    write (str, *) is
                    line1 = line1(1:ip)//' '//' s'//trim(adjustl(str))//'_y(1)'; ip = ip + 1 + 13
                end do
            end if
        end select

        line1 = line1(1:ip - 1)//'#'
        call TLab_Write_ASCII(vfile, repeat('#', len_trim(line1)))
        call TLab_Write_ASCII(vfile, trim(adjustl(line1)))
        call TLab_Write_ASCII(vfile, repeat('#', len_trim(line1)))

    end subroutine DNS_OBS_INITIALIZE

!########################################################################

    subroutine DNS_OBS()

        implicit none

        integer(wi) :: ip, is
        character(len=256) :: line1, line2

        write (line1, 100) int(obs_data(1)), itime, rtime
100     format((1x, I1), (1x, I7), (1x, E13.6))

        select case (dns_obs_log)
        case (OBS_TYPE_EKMAN)
            write (line2, 200) (obs_data(ip), ip=2, 8)
200         format(7(1x, E13.6))
            line1 = trim(line1)//trim(line2)
            if (scal_on) then
                do is = 1, inb_scal
                    write (line2, 300) obs_data(8 + is)
300                 format(1x, E13.6)
                    line1 = trim(line1)//trim(line2)
                end do
            end if
        end select

        call TLab_Write_ASCII(vfile, trim(adjustl(line1)))

    end subroutine DNS_OBS

end program DNS
