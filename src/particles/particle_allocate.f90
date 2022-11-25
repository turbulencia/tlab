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
    isize_l_comm = (2*jmax*kmax + imax*jmax*2 + 2*jmax*2)*inb_part_interp
#ifdef USE_MPI
    isize_pbuffer = int(isize_part/4*(inb_part_array*2 + 1)) ! same size for both buffers
    isize_l_comm = max(isize_l_comm, 2*isize_pbuffer)
#endif
    call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, l_comm, [isize_l_comm], 'l_comm')

    return
end subroutine PARTICLE_ALLOCATE

! ###################################################################
! ###################################################################
subroutine PARTICLE_INITIALIZE()
    use TLAB_VARS, only: g !, sbg
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    implicit none

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
