#include "dns_error.h"
#include "dns_const.h"
#include "types.h"

module PARTICLE_TRAJECTORIES

    use TLAB_CONSTANTS, only: efile, lfile, wp, sp, wi, longi
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use PARTICLE_VARS, only: isize_part
    use TLAB_PROCS
    use PARTICLE_VARS
    use DNS_LOCAL, only: nitera_save, nitera_last
#ifdef USE_MPI
    use MPI
    use TLAB_MPI_VARS, only: ims_pro, ims_err
#endif
    implicit none
    save
    private

    real(sp), allocatable :: l_trajectories(:, :, :)
    integer(longi), allocatable :: l_trajectories_tags(:)
    integer(wi) :: counter
#ifdef USE_MPI
    real(sp), allocatable :: mpi_tmp(:, :)
#endif

    public :: PARTICLE_TRAJECTORIES_INITIALIZE
    public :: PARTICLE_TRAJECTORIES_ACCUMULATE
    public :: PARTICLE_TRAJECTORIES_WRITE

contains
!#######################################################################
!#######################################################################
    subroutine PARTICLE_TRAJECTORIES_INITIALIZE()

! -------------------------------------------------------------------
        character(len=32) name
        integer ims_npro_loc
        integer(wi) j

!#######################################################################
        call TLAB_ALLOCATE_ARRAY_SINGLE(__FILE__, l_trajectories, [isize_traj + 1, nitera_save, inb_traj], 'l_traj')
        call TLAB_ALLOCATE_ARRAY1_LONG_INT(__FILE__, l_trajectories_tags, isize_traj, 'l_trajectory_tags')
#ifdef USE_MPI
        call TLAB_ALLOCATE_ARRAY_SINGLE(__FILE__, mpi_tmp, [isize_traj + 1, nitera_save], 'mpi_tmp')
#endif

!#######################################################################
        ! set the particle tags to be tracked
        select case (imode_traj)
        case default                ! track only the first isize_traj particles
            do j = 1, isize_traj
                l_trajectories_tags(j) = int(j, KIND=longi)
            end do

        case (TRAJ_TYPE_LARGEST)    ! track the ones given in a file
            write (name, *) nitera_last; name = 'largest_particle.'//trim(adjustl(name))
#ifdef USE_MPI
            if (ims_pro == 0) then
#endif
#define LOC_UNIT_ID 117
#define LOC_STATUS 'old'
#include "dns_open_file.h"
                read (LOC_UNIT_ID) ims_npro_loc
!        READ(LOC_UNIT_ID) ims_np_all(1:ims_npro_loc)
                read (LOC_UNIT_ID, POS=SIZEOFINT*(ims_npro_loc + 1) + 1) l_trajectories_tags
                close (LOC_UNIT_ID)
#undef LOC_UNIT_ID
#undef LOC_STATUS
#ifdef USE_MPI
            end if
            call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
            call MPI_BCAST(l_trajectories_tags, isize_traj, MPI_INTEGER8, 0, MPI_COMM_WORLD, ims_err)
#endif

        end select

        ! initialize values
        l_trajectories = 0.0_sp
        counter = 0

        return
    end subroutine PARTICLE_TRAJECTORIES_INITIALIZE

!#######################################################################
!#######################################################################
    subroutine PARTICLE_TRAJECTORIES_ACCUMULATE()
        use TLAB_TYPES, only: pointers_dt, pointers3d_dt
        use TLAB_VARS, only: imax, jmax, kmax
        use TLAB_VARS, only: rtime
        use TLAB_ARRAYS
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
! Setting pointers to position; we always save the position
        nvar = 0
        nvar = nvar + 1; data(nvar)%field => l_q(:, 1)
        nvar = nvar + 1; data(nvar)%field => l_q(:, 2)
        nvar = nvar + 1; data(nvar)%field => l_q(:, 3)

! Primitive variables
        do iv = 4, inb_flow_array
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, iv); data(nvar)%field => l_txc(:, iv - 3)
            l_txc(:, iv - 3) = 0.0_wp ! Field to particle is additive
        end do

        do iv = 1, inb_scal_array
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => s(:, iv); data(nvar)%field => l_txc(:, iv - 3 + inb_flow_array)
            l_txc(:, iv - 3 + inb_flow_array) = 0.0_wp ! Field to particle is additive
        end do

! -------------------------------------------------------------------
! Additional information
        if (imode_traj == TRAJ_TYPE_VORTICITY) then
            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk2d, wrk3d)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 1); data(nvar)%field => l_hq(:, 1)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 2); data(nvar)%field => l_hq(:, 2)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 3); data(nvar)%field => l_hq(:, 3)
            l_hq(:, 1:3) = 0.0_wp ! Field to particle is additive
        end if

! -------------------------------------------------------------------
! Interpolation
        if (nvar > 3) then
            call FIELD_TO_PARTICLE(nvar - 3, data_in(4:), data(4:), l_g, l_q, l_comm, wrk3d)
        end if

! -------------------------------------------------------------------
! Accumulate time
#ifdef USE_MPI
        if (ims_pro == 0) then
#endif
            l_trajectories(1, counter, 1:nvar) = SNGL(rtime)
#ifdef USE_MPI
        end if
#endif

! Accumulating the data
        do i = 1, l_g%np
            do j = 1, isize_traj
                if (l_g%tags(i) == l_trajectories_tags(j)) then
                    do iv = 1, nvar
                        l_trajectories(1 + j, counter, iv) = SNGL(data(iv)%field(i))
                    end do
                end if
            end do
        end do

        return
    end subroutine PARTICLE_TRAJECTORIES_ACCUMULATE

!#######################################################################
!#######################################################################
    subroutine PARTICLE_TRAJECTORIES_WRITE(fname)
        character*(*) fname

! -------------------------------------------------------------------
        character(len=32) name
        integer iv

!#######################################################################
        do iv = 1, inb_traj
#ifdef USE_MPI
            mpi_tmp = 0.0_sp
    call MPI_REDUCE(l_trajectories(1, 1, iv), mpi_tmp, (1 + isize_traj)*nitera_save, MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
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
                write (LOC_UNIT_ID) l_trajectories(:, :, iv)
#endif
                close (LOC_UNIT_ID)
#ifdef USE_MPI
            end if
#endif
        end do

        l_trajectories = 0.0_sp
        counter = 0

        return
    end subroutine PARTICLE_TRAJECTORIES_WRITE

end module PARTICLE_TRAJECTORIES
