#include "dns_error.h"
#include "dns_const.h"

#define USE_ACCESS_STREAM

module ParticleTrajectories

    use TLab_Constants, only: efile, lfile, wp, sp, wi, longi, sizeofint
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use TLab_WorkFlow
    use TLab_Memory
    use PARTICLE_VARS, only: inb_part, isize_part_total, inb_part_interp, inb_part_txc, part, PART_TYPE_NONE
    use DNS_LOCAL, only: nitera_save
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_pro, ims_err
#endif
    use FI_VECTORCALCULUS
    implicit none
    save
    private

    real(sp), allocatable :: l_traj(:, :, :)
    integer(longi), allocatable :: l_traj_tags(:)
    integer(wi) :: counter
#ifdef USE_MPI
    real(sp), allocatable :: mpi_tmp(:, :)
#endif

    ! Posible values of imode_traj
    integer, parameter :: TRAJ_TYPE_NONE = 0
    integer, parameter :: TRAJ_TYPE_BASIC = 1           ! save particle prognostic properties
    integer, parameter :: TRAJ_TYPE_EULERIAN = 2        ! add the Eulerian prognostic properties
    integer, parameter :: TRAJ_TYPE_VORTICITY = 3       ! add the Eulerian vorticity

    integer(wi) :: imode_traj = TRAJ_TYPE_NONE          ! Type of trajectory information that is saved
    integer(wi) :: isize_traj                           ! # of saved trajectories
    integer(wi) :: inb_traj                             ! # of properties saved along trajectories
    character(len=32) :: traj_filename                  ! file with the particle tags to be tracked; if void, then the first isize_traj particles are used

    public :: imode_traj, TRAJ_TYPE_NONE
    public :: ParticleTrajectories_Initialize
    public :: ParticleTrajectories_Accumulate
    public :: ParticleTrajectories_Write

contains
    !#######################################################################
    !#######################################################################
    subroutine ParticleTrajectories_Initialize(inifile)

        character(len=*) inifile

        ! -------------------------------------------------------------------
        character*512 sRes
        character*32 bakfile, block

        character(len=32) name
        integer ims_npro_loc
        integer(wi) j
        integer(longi) stride

        ! -------------------------------------------------------------------
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Particles'

        call SCANINICHAR(bakfile, inifile, block, 'TrajType', 'basic', sRes)
        if (trim(adjustl(sRes)) == 'basic') then; imode_traj = TRAJ_TYPE_BASIC
        elseif (trim(adjustl(sRes)) == 'eulerian') then; imode_traj = TRAJ_TYPE_EULERIAN
        elseif (trim(adjustl(sRes)) == 'vorticity') then; imode_traj = TRAJ_TYPE_VORTICITY
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Invalid option in TrajectoryType')
            call TLab_Stop(DNS_ERROR_CALCTRAJECTORIES)
        end if

        call SCANINIINT(bakfile, inifile, block, 'TrajNumber', '0', isize_traj)
        if (isize_traj > isize_part_total) then
            call TLab_Write_ASCII(efile, __FILE__//'. Number of trajectories must be less or equal than number of particles.')
            call TLab_Stop(DNS_ERROR_CALCTRAJECTORIES)
        end if
        if (isize_traj <= 0) imode_traj = TRAJ_TYPE_NONE
        if (part%type == PART_TYPE_NONE) imode_traj = TRAJ_TYPE_NONE

        if (imode_traj == TRAJ_TYPE_NONE) return

        call SCANINICHAR(bakfile, inifile, block, 'TrajFileName', 'void', traj_filename)

        ! -------------------------------------------------------------------
        inb_traj = inb_part             ! save particle prognostic properties

        select case (imode_traj)
        case (TRAJ_TYPE_BASIC)          ! save particle prognostic properties

        case (TRAJ_TYPE_EULERIAN)       ! add the Eulerian prognostic properties
            inb_traj = inb_traj + inb_flow_array + inb_scal_array

        case (TRAJ_TYPE_VORTICITY)      ! add the Eulerian vorticity
            inb_traj = inb_traj + 3

        end select

        inb_part_txc = max(inb_part_txc, inb_traj)
        inb_part_interp = max(inb_part_interp, inb_traj)

        !#######################################################################
        call TLab_Allocate_SINGLE(__FILE__, l_traj, [isize_traj + 1, nitera_save, inb_traj], 'l_traj')
        call TLab_Allocate_LONG_INT(__FILE__, l_traj_tags, [isize_traj], 'l_traj_tags')
#ifdef USE_MPI
        call TLab_Allocate_SINGLE(__FILE__, mpi_tmp, [isize_traj + 1, nitera_save], 'mpi_tmp')
#endif

        !#######################################################################
        ! set the particle tags to be tracked
        select case (trim(adjustl(traj_filename)))
        case ('void')               ! track isize_traj particles unformly distributed over the population
            stride = isize_part_total/isize_traj
            do j = 1, isize_traj
                l_traj_tags(j) = 1 + (j - 1)*stride
            end do
            if (l_traj_tags(isize_traj) > isize_part_total) then
                call TLab_Write_ASCII(efile, __FILE__//'. Tags of trajectories out of range.')
                call TLab_Stop(DNS_ERROR_CALCTRAJECTORIES)
            end if

        case default                ! track the ones given in a file
            ! write (name, *) nitera_last; name = trim(adjustl(traj_filename))//trim(adjustl(name))
            name = trim(adjustl(traj_filename))
#ifdef USE_MPI
            if (ims_pro == 0) then
#endif
#define LOC_UNIT_ID 117
#define LOC_STATUS 'old'
#include "dns_open_file.h"
                read (LOC_UNIT_ID) ims_npro_loc
                !        READ(LOC_UNIT_ID) ims_np_all(1:ims_npro_loc)
                read (LOC_UNIT_ID, POS=SIZEOFINT*(ims_npro_loc + 1) + 1) l_traj_tags
                close (LOC_UNIT_ID)
#undef LOC_UNIT_ID
#undef LOC_STATUS
#ifdef USE_MPI
            end if
            call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
            call MPI_BCAST(l_traj_tags, isize_traj, MPI_INTEGER8, 0, MPI_COMM_WORLD, ims_err)
#endif

        end select

        ! initialize values
        l_traj = 0.0_sp
        counter = 0

        return
    end subroutine ParticleTrajectories_Initialize

    !#######################################################################
    !#######################################################################
    subroutine ParticleTrajectories_Accumulate()
        use TLab_Types, only: pointers_dt, pointers3d_dt
        use TLAB_VARS, only: inb_flow_array, inb_scal_array
        use TLAB_VARS, only: imax, jmax, kmax
        use TLAB_VARS, only: rtime
        use TLab_Arrays
        use FI_VECTORCALCULUS
        use DNS_ARRAYS
        use PARTICLE_ARRAYS
        use PARTICLE_INTERPOLATE

        ! -------------------------------------------------------------------
        integer(wi) i, j
        integer iv, nvar
        type(pointers3d_dt), dimension(inb_traj) :: data_in
        type(pointers_dt), dimension(inb_traj) :: data

        !#######################################################################
        counter = counter + 1

        ! -------------------------------------------------------------------
        nvar = 0
        select case (imode_traj)

        case (TRAJ_TYPE_EULERIAN)           ! add the Eulerian prognostic properties
            do iv = 1, inb_flow_array
                nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, iv); data(nvar)%field => l_txc(:, nvar)
                l_txc(:, nvar) = 0.0_wp     ! Field_to_particle is additive
            end do

            do iv = 1, inb_scal_array
                nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => s(:, iv); data(nvar)%field => l_txc(:, nvar)
                l_txc(:, nvar) = 0.0_wp     ! Field_to_particle is additive
            end do

        case (TRAJ_TYPE_VORTICITY)          ! add the Eulerian vorticity
            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 1); data(nvar)%field => l_hq(:, 1)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 2); data(nvar)%field => l_hq(:, 2)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 3); data(nvar)%field => l_hq(:, 3)
            l_hq(:, 1:3) = 0.0_wp           ! Field_to_particle is additive
        end select

        ! Interpolation
        if (nvar > 0) then
            call FIELD_TO_PARTICLE(data_in(1:nvar), data(1:nvar), l_g, l_q)
        end if

        ! -------------------------------------------------------------------
        ! Accumulate time
#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
            l_traj(1, counter, 1:inb_part + nvar) = SNGL(rtime)
#ifdef USE_MPI
        end if
#endif

        ! Accumulating the data
        do i = 1, l_g%np
            do j = 1, isize_traj
                if (l_g%tags(i) == l_traj_tags(j)) then
                    do iv = 1, inb_part     ! save the particle prognostic variables
                        l_traj(1 + j, counter, iv) = SNGL(l_q(i, iv))
                    end do
                    do iv = 1, nvar         ! add additional variables
                        l_traj(1 + j, counter, inb_part + iv) = SNGL(data(iv)%field(i))
                    end do
                end if
            end do
        end do

        return
    end subroutine ParticleTrajectories_Accumulate

    !#######################################################################
    !#######################################################################
    subroutine ParticleTrajectories_Write(fname)
        character(len=*) fname

        ! -------------------------------------------------------------------
        character(len=32) name
        integer iv

        !#######################################################################
        do iv = 1, inb_traj
#ifdef USE_MPI
            mpi_tmp = 0.0_sp
            call MPI_REDUCE(l_traj(1, 1, iv), mpi_tmp, (1 + isize_traj)*nitera_save, MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
            call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
            if (ims_pro == 0) then
#endif

                write (name, *) iv; name = trim(adjustl(fname))//'.'//trim(adjustl(name))
#define LOC_UNIT_ID 115
#define LOC_STATUS 'unknown'
#include "dns_open_file.h"
                rewind (LOC_UNIT_ID)
#ifdef USE_MPI
                write (LOC_UNIT_ID) mpi_tmp
#else
                write (LOC_UNIT_ID) l_traj(:, :, iv)
#endif
                close (LOC_UNIT_ID)
#ifdef USE_MPI
            end if
#endif
        end do

        ! reinitialize
        l_traj = 0.0_sp
        counter = 0

        return
    end subroutine ParticleTrajectories_Write

end module ParticleTrajectories
