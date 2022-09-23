#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

module PARTICLE_TRAJECTORIES

    use TLAB_CONSTANTS, only: efile, lfile
    use TLAB_VARS, only: inb_flow_array, inb_scal_array
    use PARTICLE_VARS, only: isize_part
    use TLAB_PROCS
    use PARTICLE_VARS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_pro, ims_err
#endif

    implicit none
    save

    real(4), dimension(:, :, :), allocatable :: l_trajectories
    integer(8), dimension(:), allocatable :: l_trajectories_tags
    TINTEGER :: counter, isize_time
#ifdef USE_MPI
    real(4), dimension(:, :), allocatable :: mpi_tmp
#endif

contains

!#######################################################################
!#######################################################################
    subroutine PARTICLE_TRAJECTORIES_INITIALIZE(nitera_save, nitera_last)
#ifdef USE_MPI
        use MPI
#endif

        implicit none

        TINTEGER nitera_save, nitera_last

! -------------------------------------------------------------------
        character*32 name
        character*128 str, line
        TINTEGER ims_npro_loc

        TINTEGER j, ierr

!#######################################################################
        isize_time = nitera_save

! Adding space for saving the time
        write (str, *) (isize_traj + 1)*isize_time; line = 'Allocating array l_trajectories of size '//TRIM(ADJUSTL(str))//'x'
        write (str, *) inb_traj; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (l_trajectories(isize_traj + 1, isize_time, inb_traj), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, 'DNS. Not enough memory for l_trajectories.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

        write (str, *) isize_traj; line = 'Allocating array l_trajectory_tags of size '//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (l_trajectories_tags(isize_traj), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, 'DNS. Not enough memory for l_trajectories_tags.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

#ifdef USE_MPI
        write (str, *) (isize_traj + 1)*isize_time; line = 'Allocating array mpi_tmp of size '//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (mpi_tmp(isize_traj + 1, isize_time), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, 'DNS. Not enough memory for l_trajectories.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if
#endif

! Initialize
        if (imode_traj == TRAJ_TYPE_LARGEST) then ! Read file with tags of largest particles, the ones to track
            write (name, *) nitera_last; name = 'largest_particle.'//TRIM(ADJUSTL(name))
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

        else ! default track only isize_traj particles
            do j = 1, isize_traj
                l_trajectories_tags(j) = INT(j, KIND=8)
            end do

        end if

        l_trajectories = C_0_R
        counter = 0

        return
    end subroutine PARTICLE_TRAJECTORIES_INITIALIZE

!#######################################################################
!#######################################################################
    subroutine PARTICLE_TRAJECTORIES_ACCUMULATE(q, s, txc, l_g, l_q, l_hq, l_txc, l_comm, wrk2d, wrk3d)

        use TLAB_TYPES, only: pointers_dt, pointers3d_dt
        use TLAB_VARS, only: isize_field, imax, jmax, kmax
        use TLAB_VARS, only: rtime
        use PARTICLE_INTERPOLATE

        implicit none

        TREAL, dimension(isize_field, *), target :: q, s, txc
        type(particle_dt) :: l_g
        TREAL, dimension(isize_part, *), target :: l_q, l_hq, l_txc ! l_hq as aux array
        TREAL, dimension(isize_l_comm) :: l_comm
        TREAL, dimension(*) :: wrk2d, wrk3d

! -------------------------------------------------------------------
        TINTEGER i, j
        TINTEGER iv, nvar
        type(pointers3d_dt), dimension(inb_traj) :: data_in
        type(pointers_dt), dimension(inb_traj) :: data

!#######################################################################
        counter = counter + 1

! -------------------------------------------------------------------
! Setting pointers to position
        nvar = 0
        nvar = nvar + 1; data(nvar)%field => l_q(:, 1)
        nvar = nvar + 1; data(nvar)%field => l_q(:, 2)
        nvar = nvar + 1; data(nvar)%field => l_q(:, 3)

! Primitive variables
        do iv = 4, inb_flow_array
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => q(:, iv); data(nvar)%field => l_txc(:, iv - 3)
            l_txc(:, iv - 3) = C_0_R ! Field to particle is additive
        end do

        do iv = 1, inb_scal_array
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => s(:, iv); data(nvar)%field => l_txc(:, iv - 3 + inb_flow_array)
            l_txc(:, iv - 3 + inb_flow_array) = C_0_R ! Field to particle is additive
        end do

! -------------------------------------------------------------------
! Additional information
        if (imode_traj == TRAJ_TYPE_VORTICITY) then
            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), wrk2d, wrk3d)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 1); data(nvar)%field => l_hq(:, 1)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 2); data(nvar)%field => l_hq(:, 2)
            nvar = nvar + 1; data_in(nvar)%field(1:imax, 1:jmax, 1:kmax) => txc(:, 3); data(nvar)%field => l_hq(:, 3)
            l_hq(:, 1:3) = C_0_R ! Field to particle is additive
        end if

! -------------------------------------------------------------------
! Interpolation
        if (nvar - 3 > 0) then
            iv = nvar - 3
            call FIELD_TO_PARTICLE(iv, data_in(4), data(4), l_g, l_q, l_comm, wrk3d)
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
#ifdef USE_MPI
        use MPI
#endif

        implicit none

        character*(*) fname

! -------------------------------------------------------------------
        character(len=32) name
        TINTEGER iv

!#######################################################################
        do iv = 1, inb_traj
#ifdef USE_MPI
            mpi_tmp = C_0_R
     call MPI_REDUCE(l_trajectories(1, 1, iv), mpi_tmp, (1 + isize_traj)*isize_time, MPI_REAL4, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
            call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
            if (ims_pro == 0) then
#endif

                write (name, *) iv; name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(name))
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

        l_trajectories = C_0_R
        counter = 0

        return
    end subroutine PARTICLE_TRAJECTORIES_WRITE

end module PARTICLE_TRAJECTORIES
