#include "dns_const.h"

module PARTICLE_PROCS
    use TLAB_CONSTANTS, only: wp, wi, longi
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use TLAB_VARS, only: g, imax, jmax, kmax
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro
#endif
    implicit none
    save
    private

    public :: PARTICLE_ALLOCATE
    public :: PARTICLE_INITIALIZE

contains
    ! ###################################################################
    ! ###################################################################
    subroutine PARTICLE_ALLOCATE(C_FILE_LOC)

        character(len=*) C_FILE_LOC

        ! -------------------------------------------------------------------
        integer(wi) :: isize_l_work, ip, idummy
        integer iv

        ! ###################################################################
        call TLAB_ALLOCATE_ARRAY_LONG_INT(C_FILE_LOC, l_g%tags, [isize_part], 'l_tags')
        call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, l_g%nodes, [isize_part], 'l_g')
#ifdef USE_MPI
        allocate (ims_np_all(ims_npro))
#endif

        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, l_q, [isize_part, inb_part_array], 'l_q')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, l_txc, [isize_part, inb_part_txc], 'l_txc')

        ! Work array
        isize_l_work = (2*jmax*(kmax + 1) + (imax + 1)*jmax*2)*inb_part_interp  ! halos for particle_interpolate
#ifdef USE_MPI
        idummy = (2*jmax*(kmax + 1) + (imax + 1)*jmax*2)*inb_part_interp        ! transfer halos between MPI tasks
        isize_l_work = isize_l_work + idummy
        idummy = int(isize_part/4*(inb_part_array*2 + 1))                       ! transfer particles between MPI tasks
        isize_l_work = max(isize_l_work, 2*idummy)
#endif
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, l_work, [isize_l_work], 'l_work')

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
    end subroutine PARTICLE_ALLOCATE

! ###################################################################
! ###################################################################
    subroutine PARTICLE_INITIALIZE()

        ! Define profiles to initialize particles
        if (IniP%relative) IniP%ymean = g(2)%nodes(1) + g(2)%scale*IniP%ymean_rel

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
    end subroutine PARTICLE_INITIALIZE

end module PARTICLE_PROCS
