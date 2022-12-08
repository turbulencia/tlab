#include "dns_const.h"

subroutine PARTICLE_ALLOCATE(C_FILE_LOC)
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use TLAB_VARS, only: imax, jmax, kmax
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro
#endif
    implicit none

    character(len=*) C_FILE_LOC

    ! -------------------------------------------------------------------
    integer(wi) :: isize_l_comm                 !

    ! ###################################################################
    call TLAB_ALLOCATE_ARRAY_LONG_INT(C_FILE_LOC, l_g%tags, [isize_part], 'l_tags')
    call TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, l_g%nodes, [isize_part], 'l_g')
#ifdef USE_MPI
    allocate (ims_np_all(ims_npro))
#endif

    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, l_q, [isize_part, inb_part_array], 'l_q')
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, l_txc, [isize_part, inb_part_txc], 'l_txc')

    ! memory space for the halo regions for particle_interpolate
    isize_l_comm = (2*jmax*(kmax + 1) + (imax + 1)*jmax*2)*inb_part_interp
#ifdef USE_MPI
    isize_l_comm = isize_l_comm + (2*jmax*(kmax + 1) + (imax + 1)*jmax*2)*inb_part_interp
    isize_pbuffer = int(isize_part/4*(inb_part_array*2 + 1)) ! same size for both buffers
    isize_l_comm = max(isize_l_comm, 2*isize_pbuffer)
#endif
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, l_comm, [isize_l_comm], 'l_comm')

    allocate (p_halo_i(inb_part_interp))
    allocate (p_halo_k(inb_part_interp))
    allocate (p_halo_ik(inb_part_interp))

    return
end subroutine PARTICLE_ALLOCATE

! ###################################################################
! ###################################################################
subroutine PARTICLE_INITIALIZE()
    use TLAB_VARS, only: g, imax, jmax, kmax
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    implicit none

    integer(wi) ip_i, ip_k, ip_ik, np_i, np_k, np_ik
    integer iv

    ! Define pointers to l_work array to create halo regions
    ip_i = 1;                           np_i = 2*jmax*kmax
    ip_k = ip_i + np_i*inb_part_interp; np_k = imax*jmax*2
    ip_ik= ip_k + np_k*inb_part_interp; np_ik= 2*jmax*2

    halo_i(1:2, 1:jmax, 1:kmax, 1:inb_part_interp) => l_comm(ip_i:ip_i + np_i*inb_part_interp - 1)
    halo_k(1:imax, 1:jmax, 1:2, 1:inb_part_interp) => l_comm(ip_k:ip_k + np_k*inb_part_interp - 1)
    halo_ik(1:2, 1:jmax, 1:2, 1:inb_part_interp) => l_comm(ip_ik:ip_ik + np_ik*inb_part_interp - 1)

    do iv = 1, inb_part_interp
        p_halo_i(iv)%field => halo_i(:,:,:, iv)
        p_halo_k(iv)%field => halo_k(:,:,:, iv)
        p_halo_ik(iv)%field => halo_ik(:,:,:, iv)
    end do

#ifdef USE_MPI
    allocate(halo_mpi_send_k(imax, jmax, inb_part_interp), halo_mpi_recv_k(imax, jmax, inb_part_interp))
    allocate(halo_mpi_send_i(jmax, kmax + 1, inb_part_interp), halo_mpi_recv_i(jmax, kmax + 1, inb_part_interp))
#endif

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
