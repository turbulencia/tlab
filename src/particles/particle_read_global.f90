#include "dns_error.h"
#include "dns_const.h"

subroutine PARTICLE_READ_GLOBAL(inifile)
    use TLAB_TYPES, only: cp, ci
    use TLAB_CONSTANTS, only: efile, lfile
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use TLAB_VARS, only: imax, jmax, kmax, isize_wrk2d
    use PARTICLE_VARS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro
#endif

    implicit none

    character*(*) inifile

! -------------------------------------------------------------------
    character*512 sRes
    character*32 bakfile, block
    integer(ci) idummy
    real(cp) memory_factor

! ###################################################################
    bakfile = TRIM(ADJUSTL(inifile))//'.bak'
    block = 'Particles'
    call TLAB_WRITE_ASCII(lfile, 'Reading '//trim(adjustl(block))//' input data.')

    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Number=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#MemoryFactor=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#IniType=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#IniYMean=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#IniThick=<value>')

    call SCANINICHAR(bakfile, inifile, block, 'Type', 'None', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; imode_part = PART_TYPE_NONE
    else if (TRIM(ADJUSTL(sRes)) == 'tracer') then; imode_part = PART_TYPE_TRACER
    else if (TRIM(ADJUSTL(sRes)) == 'simplesettling') then; imode_part = PART_TYPE_SIMPLE_SETT
    else if (TRIM(ADJUSTL(sRes)) == 'bilinearcloudthree') then; imode_part = PART_TYPE_BIL_CLOUD_3
    else if (TRIM(ADJUSTL(sRes)) == 'bilinearcloudfour') then; imode_part = PART_TYPE_BIL_CLOUD_4
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Particles.Type.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

! -------------------------------------------------------------------
    if (imode_part /= PART_TYPE_NONE) then
        call SCANINILONGINT(bakfile, inifile, block, 'Number', '0', isize_part_total)
        call SCANINIREAL(bakfile, inifile, block, 'MemoryFactor', '2.0', memory_factor)

! -------------------------------------------------------------------
        call SCANINIINT(bakfile, inifile, block, 'IniType', '1', part_ini_mode)
        call SCANINIREAL(bakfile, inifile, block, 'IniYMean', '0.5', part_ini_ymean)
        call SCANINIREAL(bakfile, inifile, block, 'IniThick', '1.0', part_ini_thick)

! -------------------------------------------------------------------
        call SCANINICHAR(bakfile, inifile, block, 'CalculatePdf', 'no', sRes)
        if (TRIM(ADJUSTL(sRes)) == 'yes') then; icalc_part_pdf = 1
        elseif (TRIM(ADJUSTL(sRes)) == 'no') then; icalc_part_pdf = 0
        end if
        call SCANINICHAR(bakfile, inifile, block, 'PdfSubdomain', '-1', sRes)
        particle_pdf_subdomain = 0.0_cp; idummy = 6
        call LIST_REAL(sRes, idummy, particle_pdf_subdomain)
        call SCANINIREAL(bakfile, inifile, block, 'PdfMax', '10', particle_pdf_max)
        call SCANINIREAL(bakfile, inifile, block, 'PdfInterval', '0.5', particle_pdf_interval)

        call SCANINICHAR(bakfile, inifile, block, 'ResidenceReset', 'yes', sRes)
        if (TRIM(ADJUSTL(sRes)) == 'yes') then; residence_reset = 1
        elseif (TRIM(ADJUSTL(sRes)) == 'no') then; residence_reset = 0
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. ResidenceReset must be yes or no')
            call TLAB_STOP(DNS_ERROR_RESIDENCERESET)
        end if

        call SCANINICHAR(bakfile, inifile, block, 'Parameters', '0.0', sRes)
        idummy = MAX_LAGPARAM
        call LIST_REAL(sRes, idummy, PARTICLE_param)

! -------------------------------------------------------------------
        inb_traj = inb_flow_array + inb_scal_array
        call SCANINICHAR(bakfile, inifile, block, 'TrajectoryType', 'first', sRes)
        if (TRIM(ADJUSTL(sRes)) == 'first') then; imode_traj = TRAJ_TYPE_FIRST
        elseif (TRIM(ADJUSTL(sRes)) == 'vorticity') then; imode_traj = TRAJ_TYPE_VORTICITY
            inb_traj = inb_traj + 3 ! + vorticity
        elseif (TRIM(ADJUSTL(sRes)) == 'largest') then; imode_traj = TRAJ_TYPE_LARGEST
        elseif (TRIM(ADJUSTL(sRes)) == 'none') then; imode_traj = TRAJ_TYPE_NONE
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Invalid option in TrajectoryType')
            call TLAB_STOP(DNS_ERROR_CALCTRAJECTORIES)
        end if

        call SCANINIINT(bakfile, inifile, block, 'TrajectoryNumber', '0', isize_traj)
        if (isize_traj > isize_part_total) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of trajectories must be less or equal than number of particles.')
            call TLAB_STOP(DNS_ERROR_CALCTRAJECTORIES)
        end if
        if (isize_traj <= 0) imode_traj = TRAJ_TYPE_NONE
        if (imode_part == PART_TYPE_NONE) imode_traj = TRAJ_TYPE_NONE

! ###################################################################
! Initializing size of particle arrays
! ###################################################################
! -------------------------------------------------------------------
! default
        inb_part_array = 3          ! # of particle properties in array
        inb_part = 3                ! # of particle properties in Runge-Kutta
        inb_part_txc = 1            ! # of particle auxiliary properties for intermediate calculations
        inb_part_interp = 3

        if (imode_part == PART_TYPE_BIL_CLOUD_3) then
            inb_part_array = 5
            inb_part = 5
            inb_part_txc = 4
            inb_part_interp = inb_part_interp + 4
            PARTICLE_SPNAME(1) = 'droplet_diff_3'
            PARTICLE_SPNAME(2) = 'droplet_nodiff_3'

        elseif (imode_part == PART_TYPE_BIL_CLOUD_4) then
            inb_part_array = 5 + 2 ! Space for residence time pdf
            inb_part = 5
            inb_part_txc = 4
            inb_part_interp = inb_part_interp + 4
            PARTICLE_SPNAME(1) = 'droplet_diff_3'
            PARTICLE_SPNAME(2) = 'droplet_nodiff_3'
            PARTICLE_SPNAME(3) = 'residence_part'

        end if

#ifdef USE_MPI
        isize_part = INT(isize_part_total/INT(ims_npro, KIND=8))
        if (MOD(isize_part_total, INT(ims_npro, KIND=8)) /= 0) then ! All PEs with equal memory
            isize_part = isize_part + 1
        end if
        isize_part = isize_part*INT(memory_factor*100)/100          ! extra memory space to allow particles concentrating into a few processors
#else
        isize_part = INT(isize_part_total)
#endif

        if (imode_traj /= TRAJ_TYPE_NONE) then
            inb_part_txc = MAX(inb_part_txc, inb_flow_array + inb_scal_array - 3)
            inb_part_interp = MAX(inb_part_interp, inb_traj)
        end if

        isize_pbuffer = int(isize_part/4*(inb_part_array*2 + 1)) !same size for both buffers
        isize_l_comm = 2*jmax*kmax*inb_part_interp &
                       + imax*jmax*2*inb_part_interp &
                       + 2*jmax*2*inb_part_interp
        isize_l_comm = MAX(isize_l_comm, 2*isize_pbuffer)

        idummy = MAX((imax + 1)*jmax, MAX((imax + 1)*kmax, jmax*(kmax + 1)))
        isize_wrk2d = MAX(isize_wrk2d, idummy)

    end if

    return
end subroutine PARTICLE_READ_GLOBAL
