#include "dns_const.h"

#define C_FILE_LOC "INIPART"

program INIPART
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal, tag_part
    use TLab_Time, only: itime
    use TLab_Memory, only: inb_scal, inb_txc, inb_scal_array, inb_flow_array
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use FDM, only: FDM_Initialize
    use Thermodynamics
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize
    use Rotation, only: Rotation_Initialize
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use PARTICLE_PROCS
    use Profiles, only: profiles_dt, Profiles_ReadBlock
    use TLab_Grid

    implicit none

    type(profiles_dt) :: IniP                           ! Information about the initialization
    integer, parameter :: PART_INITYPE_HARDCODED = 101  ! Special type of particle initialization for testing
    integer, parameter :: PART_INITYPE_SCALAR = 102     ! Special type of particle initialization not included in default profile data

    character*512 sRes
    character*32 bakfile

    !########################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif
    call Particle_Initialize_Parameters(ifile)

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
    call Rotation_Initialize(ifile)

    call TLab_Consistency_Check()

    if (part%type /= PART_TYPE_NONE) then

        ! -------------------------------------------------------------------
        ! Read partcile parameters
        ! -------------------------------------------------------------------
        bakfile = trim(adjustl(ifile))//'.bak'

        call Profiles_ReadBlock(bakfile, ifile, 'Particles', 'IniP', IniP, 'gaussian') ! using gaussian as dummy to read rest of profile information
        call ScanFile_Char(bakfile, ifile, 'Particles', 'ProfileIniP', 'None', sRes)
        if (trim(adjustl(sRes)) == 'scalar') IniP%type = PART_INITYPE_SCALAR
        if (trim(adjustl(sRes)) == 'hardcoded') IniP%type = PART_INITYPE_HARDCODED

        ! -------------------------------------------------------------------
        ! Allocating memory space
        ! -------------------------------------------------------------------
        inb_flow_array = 0
        inb_scal_array = 0
        inb_txc = inb_scal

        call TLab_Initialize_Memory(C_FILE_LOC)

        call Particle_Initialize_Memory(C_FILE_LOC)

        ! problem if I enter with inb_scal_array = 0
        inb_scal_array = inb_scal
        call TLab_Initialize_Background(ifile)
        if (IniP%relative) IniP%ymean = y%nodes(1) + y%scale*IniP%ymean_rel

        call Particle_Initialize_Fields()

        ! -------------------------------------------------------------------
        ! Initialize particle information
        ! -------------------------------------------------------------------
        call Particle_Initialize_Variables(l_g, l_q, l_txc, txc)

        call IO_WRITE_PARTICLE(trim(adjustl(tag_part))//'ics', l_g, l_q)

    end if

    call TLab_Stop(0)

contains

    ! ###################################################################
    subroutine Particle_Initialize_Variables(l_g, l_q, l_txc, txc)
        use TLab_Constants, only: wp, wi, longi
        use TLab_Pointers, only: pointers_dt
        use TLab_Pointers_3D, only: pointers3d_dt
        use TLab_Memory, only: imax, jmax, kmax, inb_scal
        use TLab_Grid, only: x, y, z
        use Tlab_Background, only: sbg
        use PARTICLE_TYPES, only: particle_dt
        use PARTICLE_VARS
        use PARTICLE_INTERPOLATE
        use Thermodynamics, only: imixture
        use THERMO_AIRWATER
#ifdef USE_MPI
        use mpi_f08
        use TLabMPI_VARS
        use PARTICLE_ARRAYS, only: ims_np_all
#endif
        ! use PARTICLE_TINIA
        use IO_Fields
        implicit none

        type(particle_dt), intent(inout) :: l_g
        real(wp), intent(inout) :: l_q(isize_part, inb_part_array)
        real(wp), intent(inout), target :: l_txc(isize_part, 2)
        real(wp), intent(inout), target :: txc(imax, jmax, kmax, inb_scal)

! -------------------------------------------------------------------
        integer(wi) i, is
        integer(wi), allocatable :: x_seed(:)
        integer(wi) size_seed

        real(wp) xref, zref, xscale, zscale, dy_loc, dummy, dy_frac
        real(wp) y_limits(2)
        integer(wi) j_limits(2)

        real(wp) rnd_number(3), rnd_number_second
        integer(wi) rnd_scal(3)

        integer(wi) nvar
        type(pointers3d_dt), dimension(2) :: data
        type(pointers_dt), dimension(2) :: data_out

        integer(longi) count

        real(wp) params(0)

!########################################################################
#ifdef USE_MPI
        l_g%np = int(isize_part_total/int(ims_npro, KIND=8))
        if (ims_pro < int(mod(isize_part_total, int(ims_npro, KIND=8)))) then
            l_g%np = l_g%np + 1
        end if
        call MPI_ALLGATHER(l_g%np, 1, MPI_INTEGER4, ims_np_all, 1, MPI_INTEGER4, MPI_COMM_WORLD, ims_err)

#else
        l_g%np = int(isize_part_total)

#endif

! Create tags
        count = 0
#ifdef USE_MPI
        do i = 1, ims_pro
            count = count + int(ims_np_all(i), KIND=8)
        end do
#endif
        do i = 1, l_g%np
            l_g%tags(i) = int(i, KIND=8) + count
        end do

! Generate seed - different seed for each processor
        call random_seed(SIZE=size_seed)
        allocate (x_seed(size_seed))
#ifdef USE_MPI
        x_seed = [(i, i=1 + ims_pro, size_seed + ims_pro)]
#else
        x_seed = [(i, i=1, size_seed)]
#endif
        call random_seed(PUT=x_seed)

#ifdef USE_MPI
        xref = x%nodes(ims_offset_i + 1)
        xscale = x%scale/real(ims_npro_i, wp)
        zref = z%nodes(ims_offset_k + 1)
        zscale = z%scale/real(ims_npro_k, wp)
#else
        xref = x%nodes(1)
        xscale = x%scale
        zref = z%nodes(1)
        zscale = z%scale
#endif
        if (z%size == 1) zscale = 0.0_wp ! 2D case

!########################################################################
! Particle position
!########################################################################
        select case (IniP%type)

        case default

            call random_number(l_q(1:l_g%np, 1))
            l_q(1:l_g%np, 1) = xref + l_q(1:l_g%np, 1)*xscale

            call random_number(l_q(1:l_g%np, 3))
            l_q(1:l_g%np, 3) = zref + l_q(1:l_g%np, 3)*zscale

            call random_number(l_q(1:l_g%np, 2))
            l_q(1:l_g%np, 2) = IniP%ymean + (l_q(1:l_g%np, 2) - 0.5_wp)*IniP%diam

        case (PART_INITYPE_HARDCODED)       ! For testing
            l_q(1:l_g%np, 1) = xref         ! First point in the local domain, to ensure that the particles are distributed over all MPI tasks
            l_q(1:l_g%np, 3) = zref

            l_q(1:l_g%np, 2) = IniP%ymean

        case (PART_INITYPE_SCALAR)          ! Use the scalar field to create the particle distribution
            call IO_Read_Fields('scal.ics', imax, jmax, kmax, itime, inb_scal, 0, txc, params)
            is = 1 ! Reference scalar

            y_limits(1) = IniP%ymean - 0.5_wp*IniP%diam
            y_limits(2) = IniP%ymean + 0.5_wp*IniP%diam
            call LOCATE_Y(2, y_limits, j_limits, y%size, y%nodes)
            dy_loc = y%nodes(j_limits(2)) - y%nodes(j_limits(1))

            i = 1
            do while (i <= l_g%np)
                call random_number(rnd_number(1:3))

                rnd_scal(1) = 1 + floor(rnd_number(1)*imax)
                rnd_scal(3) = 1 + floor(rnd_number(3)*kmax)
                rnd_scal(2) = j_limits(1) + floor(rnd_number(2)*(j_limits(2) - j_limits(1) + 1))
                dy_frac = rnd_number(2)*(j_limits(2) - j_limits(1) + 1) &
                          - floor(rnd_number(2)*(j_limits(2) - j_limits(1) + 1))

                dummy = (txc(rnd_scal(1), rnd_scal(2), rnd_scal(3), is) - sbg(is)%mean)/sbg(is)%delta
                dummy = abs(dummy + 0.5_wp)

                call random_number(rnd_number_second)

                if (rnd_number_second <= dummy) then

                    l_q(i, 1) = xref + rnd_number(1)*xscale
                    l_q(i, 3) = zref + rnd_number(3)*zscale
                    l_q(i, 2) = y%nodes(rnd_scal(2)) + dy_frac*dy_loc

                    i = i + 1

                end if

            end do

        end select

! Calculating closest node below in Y direction
        call LOCATE_Y(l_g%np, l_q(1, 2), l_g%nodes, y%size, y%nodes)

!########################################################################
! Remaining particle properties
!########################################################################
        select case (part%type)

        case (PART_TYPE_INERTIA)                            ! velocity
            l_q(:, 4:6) = 0.0_wp

        case (PART_TYPE_BIL_CLOUD_3, PART_TYPE_BIL_CLOUD_4) ! scalar fields
            call IO_Read_Fields('scal.ics', imax, jmax, kmax, itime, inb_scal, 0, txc, params)

            if (imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                nvar = 0
                nvar = nvar + 1; data(nvar)%field => txc(:, :, :, 1); data_out(nvar)%field => l_txc(:, 1)
                nvar = nvar + 1; data(nvar)%field => txc(:, :, :, 2); data_out(nvar)%field => l_txc(:, 2)
                l_txc(:, 1:2) = 0.0_wp
                call FIELD_TO_PARTICLE(data(1:nvar), data_out(1:nvar), l_g, l_q)

                l_q(:, 4) = 0.0_wp
                call THERMO_AIRWATER_LINEAR(l_g%np, l_txc(1, 1), l_q(1, 4))

                l_q(:, 5) = l_q(:, 4) ! l_q(:,6) for bil_cloud_4 is set =0 in dns_main at initialization

            end if

            ! case (PART_TYPE_NEW_CASES)
        case (PART_TYPE_TINIA_1)
            ! call PARTICLE_TINIA_INITIALIZE()

        end select

        return
    end subroutine Particle_Initialize_Variables

end program INIPART
