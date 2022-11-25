#include "dns_error.h"
#include "dns_const.h"

subroutine PARTICLE_READ_GLOBAL(inifile)
    use TLAB_CONSTANTS, only: wp, wi, longi, efile, lfile
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use TLAB_VARS, only: imax, jmax, kmax, isize_wrk2d
    use PARTICLE_VARS
    use TLAB_PROCS
    use PROFILES
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro
#endif

    implicit none

    character(len=*) inifile

! -------------------------------------------------------------------
    character*512 sRes
    character*32 bakfile, block
    integer(wi) idummy
#ifdef USE_MPI
    real(wp) memory_factor
#endif

! ###################################################################
    bakfile = trim(adjustl(inifile))//'.bak'
    block = 'Particles'
    call TLAB_WRITE_ASCII(lfile, 'Reading '//trim(adjustl(block))//' input data.')

    call TLAB_WRITE_ASCII(bakfile, '#')
    call TLAB_WRITE_ASCII(bakfile, '#['//trim(adjustl(block))//']')
    call TLAB_WRITE_ASCII(bakfile, '#Type=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#Number=<value>')
    call TLAB_WRITE_ASCII(bakfile, '#MemoryFactor=<value>')

    call SCANINICHAR(bakfile, inifile, block, 'Type', 'None', sRes)
    if (trim(adjustl(sRes)) == 'none') then; part%type = PART_TYPE_NONE
    else if (trim(adjustl(sRes)) == 'tracer') then; part%type = PART_TYPE_TRACER
    else if (trim(adjustl(sRes)) == 'simplesettling') then; part%type = PART_TYPE_SIMPLE_SETT
    else if (trim(adjustl(sRes)) == 'bilinearcloudthree') then; part%type = PART_TYPE_BIL_CLOUD_3
    else if (trim(adjustl(sRes)) == 'bilinearcloudfour') then; part%type = PART_TYPE_BIL_CLOUD_4
    else
        call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Particles.Type.')
        call TLAB_STOP(DNS_ERROR_OPTION)
    end if

    isize_part = 0
    if (part%type /= PART_TYPE_NONE) then
        call SCANINICHAR(bakfile, inifile, block, 'Parameters', '0.0', sRes)
        idummy = MAX_PARS
        call LIST_REAL(sRes, idummy, part%parameters)

        call SCANINILONGINT(bakfile, inifile, block, 'Number', '0', isize_part_total)

! -------------------------------------------------------------------
        call PROFILES_READBLOCK(bakfile, inifile, block, 'IniP', IniP, 'gaussian') ! using gaussian as dummy to read rest of profile information
        call SCANINICHAR(bakfile, inifile, block, 'ProfileIniP', 'None', sRes)
        if (trim(adjustl(sRes)) == 'scalar') IniP%type = PART_INITYPE_SCALAR

! -------------------------------------------------------------------
        particle_pdf_calc = .false.
        call SCANINICHAR(bakfile, inifile, block, 'CalculatePdf', 'no', sRes)
        if (trim(adjustl(sRes)) == 'yes') particle_pdf_calc = .true.

        call SCANINICHAR(bakfile, inifile, block, 'PdfSubdomain', '-1', sRes)
        particle_pdf_subdomain = 0.0_wp; idummy = 6
        call LIST_REAL(sRes, idummy, particle_pdf_subdomain)
        call SCANINIREAL(bakfile, inifile, block, 'PdfMax', '10', particle_pdf_max)
        call SCANINIREAL(bakfile, inifile, block, 'PdfInterval', '0.5', particle_pdf_interval)

        call SCANINICHAR(bakfile, inifile, block, 'ResidenceReset', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'yes') then; residence_reset = 1
        elseif (trim(adjustl(sRes)) == 'no') then; residence_reset = 0
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. ResidenceReset must be yes or no')
            call TLAB_STOP(DNS_ERROR_RESIDENCERESET)
        end if

! -------------------------------------------------------------------
        call SCANINIINT(bakfile, inifile, block, 'TrajNumber', '0', isize_traj)
        if (isize_traj > isize_part_total) then
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Number of trajectories must be less or equal than number of particles.')
            call TLAB_STOP(DNS_ERROR_CALCTRAJECTORIES)
        end if
        if (isize_traj <= 0) imode_traj = TRAJ_TYPE_NONE
        if (part%type == PART_TYPE_NONE) imode_traj = TRAJ_TYPE_NONE

        call SCANINICHAR(bakfile, inifile, block, 'TrajFileName', 'void', traj_filename)

        inb_traj = inb_flow_array + inb_scal_array
        call SCANINICHAR(bakfile, inifile, block, 'TrajType', 'basic', sRes)
        if (trim(adjustl(sRes)) == 'basic') then; imode_traj = TRAJ_TYPE_BASIC
        elseif (trim(adjustl(sRes)) == 'vorticity') then; imode_traj = TRAJ_TYPE_VORTICITY
            inb_traj = inb_traj + 3 ! + vorticity
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Invalid option in TrajectoryType')
            call TLAB_STOP(DNS_ERROR_CALCTRAJECTORIES)
        end if

! ###################################################################
! Initializing size of particle arrays
! ###################################################################
! -------------------------------------------------------------------
! default
        inb_part_array = 3          ! # of particle properties in array
        inb_part = 3                ! # of particle properties in Runge-Kutta
        inb_part_txc = 1            ! # of particle auxiliary properties for intermediate calculations
        inb_part_interp = 3         ! # of interpolated fields into lagrangian framework

        if (part%type == PART_TYPE_BIL_CLOUD_3) then
            inb_part_array = 5
            inb_part = 5
            inb_part_txc = 4
            inb_part_interp = inb_part_interp + 4
            part_spname(1) = 'droplet_diff_3'
            part_spname(2) = 'droplet_nodiff_3'

        elseif (part%type == PART_TYPE_BIL_CLOUD_4) then
            inb_part_array = 5 + 2 ! Space for residence time pdf
            inb_part = 5
            inb_part_txc = 4
            inb_part_interp = inb_part_interp + 4
            part_spname(1) = 'droplet_diff_3'
            part_spname(2) = 'droplet_nodiff_3'
            part_spname(3) = 'residence_part'

        end if

#ifdef USE_MPI
        call SCANINIREAL(bakfile, inifile, block, 'MemoryFactor', '2.0', memory_factor)
        isize_part = int(isize_part_total/int(ims_npro, longi))
        if (mod(isize_part_total, int(ims_npro, longi)) /= 0) then ! All PEs with equal memory
            isize_part = isize_part + 1
        end if
        isize_part = isize_part*int(memory_factor*100)/100          ! extra memory space to allow particles concentrating into a few processors
#else
        isize_part = int(isize_part_total)
#endif

        if (imode_traj /= TRAJ_TYPE_NONE) then
            inb_part_txc = max(inb_part_txc, inb_flow_array + inb_scal_array - 3)
            inb_part_interp = max(inb_part_interp, inb_traj)
        end if

        isize_l_comm = (2*jmax*kmax + imax*jmax*2 + 2*jmax*2)*inb_part_interp
#ifdef USE_MPI
        isize_pbuffer = int(isize_part/4*(inb_part_array*2 + 1)) ! same size for both buffers
        isize_l_comm = max(isize_l_comm, 2*isize_pbuffer)
#endif

        ! I do not think this is needed; ony used in particle_to_field.f90
        idummy = max((imax + 1)*jmax, max((imax + 1)*kmax, jmax*(kmax + 1)))
        isize_wrk2d = max(isize_wrk2d, idummy)

    end if

    return
end subroutine PARTICLE_READ_GLOBAL
