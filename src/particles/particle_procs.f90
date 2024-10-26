#include "dns_const.h"
#include "dns_error.h"

module PARTICLE_PROCS
    use TLab_Constants, only: wp, wi, longi, efile, lfile, MAX_PARS
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use TLAB_VARS, only: imax, jmax, kmax, isize_wrk3d
    use TLab_WorkFlow
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro
#endif
    implicit none
    save
    private

    public :: Particle_Initialize_Parameters
    public :: Particle_Initialize_Memory
    public :: Particle_Initialize_Fields

contains
    ! ###################################################################
    ! ###################################################################
    subroutine Particle_Initialize_Parameters(inifile)
        ! use PARTICLE_TINIA

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
        else if (trim(adjustl(sRes)) == 'inertia') then; part%type = PART_TYPE_INERTIA
        else if (trim(adjustl(sRes)) == 'bilinearcloudthree') then; part%type = PART_TYPE_BIL_CLOUD_3
        else if (trim(adjustl(sRes)) == 'bilinearcloudfour') then; part%type = PART_TYPE_BIL_CLOUD_4
        else if (trim(adjustl(sRes)) == 'tiniaone') then; part%type = PART_TYPE_TINIA_1
        else
            call TLAB_WRITE_ASCII(efile, __FILE__//'. Wrong Particles.Type.')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

        isize_part = 0
        if (part%type == PART_TYPE_NONE) return

! -------------------------------------------------------------------
        part_bcs = PART_BCS_NONE
        if (part%type == PART_TYPE_INERTIA) part_bcs = PART_BCS_SPECULAR
        if (part%type == PART_TYPE_TINIA_1) part_bcs = PART_BCS_STICK
        call SCANINICHAR(bakfile, inifile, block, 'BoundaryCondition', 'Void', sRes)
        if (trim(adjustl(sRes)) == 'none') then; part_bcs = PART_BCS_NONE
        else if (trim(adjustl(sRes)) == 'specular') then; part_bcs = PART_BCS_SPECULAR
        else if (trim(adjustl(sRes)) == 'stick') then; part_bcs = PART_BCS_STICK
        end if

        call SCANINICHAR(bakfile, inifile, block, 'Parameters', '0.0', sRes)
        idummy = MAX_PARS
        call LIST_REAL(sRes, idummy, part%parameters)

        call SCANINILONGINT(bakfile, inifile, block, 'Number', '0', isize_part_total)

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

! ###################################################################
! Initializing size of particle arrays
! ###################################################################
        inb_part = 3                ! # of particle properties in Runge-Kutta (prognostic properties): at least, the position
        inb_part_array = inb_part   ! # of particle properties in array (prognostic plus diagnostic): at leat, prognostic variables
        inb_part_txc = 1            ! # of particle auxiliary properties for intermediate calculations
        inb_part_interp = 3         ! # of interpolated fields into lagrangian framework

        select case (part%type)

        case (PART_TYPE_INERTIA)
            inb_part = inb_part + 3         ! add particle velocity to particle position
            inb_part_array = inb_part
            inb_part_txc = 3

        case (PART_TYPE_BIL_CLOUD_3)
            inb_part = 5
            inb_part_array = inb_part
            inb_part_txc = 4
            inb_part_interp = inb_part_interp + 4
            part_spname(1) = 'droplet_diff_3'
            part_spname(2) = 'droplet_nodiff_3'

        case (PART_TYPE_BIL_CLOUD_4)
            inb_part = 5
            inb_part_array = inb_part + 2   ! add for residence time pdf
            inb_part_txc = 4
            inb_part_interp = inb_part_interp + 4
            part_spname(1) = 'droplet_diff_3'
            part_spname(2) = 'droplet_nodiff_3'
            part_spname(3) = 'residence_part'

            ! case (PART_TYPE_NEW_CASES)
        case (PART_TYPE_TINIA_1)
            ! call PARTICLE_TINIA_READBLOCK(bakfile, inifile, block)

        end select

#ifdef USE_MPI
        call SCANINIREAL(bakfile, inifile, block, 'MemoryFactor', '2.0', memory_factor)
        isize_part = int(isize_part_total/int(ims_npro, longi))
        if (mod(isize_part_total, int(ims_npro, longi)) /= 0) then ! All PEs with equal memory
            isize_part = isize_part + 1
        end if
        isize_part = isize_part*int(memory_factor*100)/100         ! extra memory space to allow particles concentrating into a few processors
#else
        isize_part = int(isize_part_total)
#endif

        isize_wrk3d = max(isize_wrk3d, (imax + 1)*jmax*(kmax + 1))
        isize_wrk3d = max(isize_wrk3d, (jmax*(kmax + 1)*inb_part_interp*2))
        isize_wrk3d = max(isize_wrk3d, (jmax*(imax + 1)*inb_part_interp*2))

        return
    end subroutine Particle_Initialize_Parameters

    ! ###################################################################
    ! ###################################################################
    subroutine Particle_Initialize_Memory(C_FILE_LOC)
        use TLab_Memory

        character(len=*) C_FILE_LOC

        ! -------------------------------------------------------------------
        integer(wi) :: isize_l_work, ip, idummy
        integer iv

        ! ###################################################################
        call TLab_Allocate_LONG_INT(C_FILE_LOC, l_g%tags, [isize_part], 'l_tags')
        call TLab_Allocate_INT(C_FILE_LOC, l_g%nodes, [isize_part], 'l_g')
#ifdef USE_MPI
        allocate (ims_np_all(ims_npro))
#endif

        call TLab_Allocate_DOUBLE(C_FILE_LOC, l_q, [isize_part, inb_part_array], 'l_q')
        call TLab_Allocate_DOUBLE(C_FILE_LOC, l_txc, [isize_part, inb_part_txc], 'l_txc')

        ! Work array
        isize_l_work = (2*jmax*(kmax + 1) + (imax + 1)*jmax*2)*inb_part_interp  ! halos for particle_interpolate
#ifdef USE_MPI
        idummy = (2*jmax*(kmax + 1) + (imax + 1)*jmax*2)*inb_part_interp        ! transfer halos between MPI tasks
        isize_l_work = isize_l_work + idummy
        idummy = int(isize_part/4*(inb_part_array*2 + 1))                       ! transfer particles between MPI tasks
        isize_l_work = max(isize_l_work, 2*idummy)
#endif
        call TLab_Allocate_DOUBLE(C_FILE_LOC, l_work, [isize_l_work], 'l_work')

#ifdef USE_MPI
        p_buffer_1 => l_work(1:idummy)
        p_buffer_2 => l_work(idummy + 1:idummy*2)
#endif

        ip = 1; idummy = 2*jmax*kmax*inb_part_interp
        halo_i(1:2, 1:jmax, 1:kmax, 1:inb_part_interp) => l_work(ip:ip + idummy - 1)

        ip = ip + idummy; idummy = imax*jmax*2*inb_part_interp
        halo_k(1:imax, 1:jmax, 1:2, 1:inb_part_interp) => l_work(ip:ip + idummy - 1)

        ip = ip + idummy; idummy = 2*jmax*2*inb_part_interp
        halo_ik(1:2, 1:jmax, 1:2, 1:inb_part_interp) => l_work(ip:ip + idummy - 1)

#ifdef USE_MPI
        ip = ip + idummy; idummy = imax*jmax*inb_part_interp
        halo_mpi_send_k(1:imax, 1:jmax, 1:inb_part_interp) => l_work(ip:ip + idummy - 1)

        ip = ip + idummy; idummy = imax*jmax*inb_part_interp
        halo_mpi_recv_k(1:imax, 1:jmax, 1:inb_part_interp) => l_work(ip:ip + idummy - 1)

        ip = ip + idummy; idummy = jmax*(kmax + 1)*inb_part_interp
        halo_mpi_send_i(1:jmax, 1:kmax + 1, 1:inb_part_interp) => l_work(ip:ip + idummy - 1)

        ip = ip + idummy; idummy = jmax*(kmax + 1)*inb_part_interp
        halo_mpi_recv_i(1:jmax, 1:kmax + 1, 1:inb_part_interp) => l_work(ip:ip + idummy - 1)
#endif

        ! arrays of pointers that is used in partcile_interpolate
        allocate (p_halo_i(inb_part_interp))
        allocate (p_halo_k(inb_part_interp))
        allocate (p_halo_ik(inb_part_interp))

        do iv = 1, inb_part_interp
            p_halo_i(iv)%field => halo_i(:, :, :, iv)
            p_halo_k(iv)%field => halo_k(:, :, :, iv)
            p_halo_ik(iv)%field => halo_ik(:, :, :, iv)
        end do

        return
    end subroutine Particle_Initialize_Memory

! ###################################################################
! ###################################################################
    subroutine Particle_Initialize_Fields()

        ! set boundarys for residence time pdf
        if (part%type == PART_TYPE_BIL_CLOUD_4) then
            ! to be rewritten
            ! l_y_lambda = (g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(1)%ymean_rel - 2.0_wp
            ! l_y_base = ((g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(1)%ymean_rel - &
            !             (g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(3)%ymean_rel)/2.0_wp &
            !            + (g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(3)%ymean_rel
            if (residence_reset == 1) then
                l_q(:, 6:7) = 0.0_wp
            end if
        end if

        return
    end subroutine Particle_Initialize_Fields

end module PARTICLE_PROCS
