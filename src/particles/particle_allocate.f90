#include "dns_const.h"
#include "dns_error.h"

subroutine PARTICLE_ALLOCATE(C_FILE_LOC)
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro
#endif
    implicit none

    character(LEN=*) C_FILE_LOC

    ! -------------------------------------------------------------------

    ! ###################################################################
    call TLAB_ALLOCATE_ARRAY1_LONG_INT(C_FILE_LOC, l_g%tags, isize_part, 'l_tags')
    call TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, l_g%nodes, isize_part, 'l_g')
#ifdef USE_MPI
    allocate (ims_np_all(ims_npro))
#endif

    call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, l_q, inb_part_array, isize_part, 'l_q')
    call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, l_txc, inb_part_txc, isize_part, 'l_txc')

    call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, l_comm, isize_l_comm, 'l_comm')

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
    if (imode_part == PART_TYPE_BIL_CLOUD_4) then
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
