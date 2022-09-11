#include "dns_const.h"
#include "dns_error.h"

subroutine PARTICLE_ALLOCATE(C_FILE_LOC)
    use TLAB_TYPES, only: ci
    use TLAB_CONSTANTS, only: lfile, efile
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro
#endif

    implicit none

    character(LEN=*) C_FILE_LOC

    ! -------------------------------------------------------------------
    character(LEN=128) str, line
    integer(ci) ierr

    ! ###################################################################
    write (str, *) isize_part; line = 'Allocating array l_g.tags of size '//TRIM(ADJUSTL(str))
    call TLAB_WRITE_ASCII(lfile, line)
    allocate (l_g%tags(isize_part), stat=ierr)
    if (ierr /= 0) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_tags.')
        call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    write (str, *) isize_part; line = 'Allocating array l_g.nodes of size '//TRIM(ADJUSTL(str))
    call TLAB_WRITE_ASCII(lfile, line)
    allocate (l_g%nodes(isize_part), stat=ierr)
    if (ierr /= 0) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_g.')
        call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

#ifdef USE_MPI
    allocate (ims_size_p(ims_npro))
#endif

    write (str, *) isize_part; line = 'Allocating array l_q of size '//TRIM(ADJUSTL(str))//'x'
    write (str, *) inb_part_array; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    call TLAB_WRITE_ASCII(lfile, line)
    allocate (l_q(isize_part, inb_part_array), stat=ierr)
    if (ierr /= 0) then
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_q.')
        call TLAB_STOP(DNS_ERROR_ALLOC)
    end if

    if (inb_part_txc > 0) then
        write (str, *) isize_part; line = 'Allocating array l_txc of size '//TRIM(ADJUSTL(str))//'x'
        write (str, *) inb_part_txc; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (l_txc(isize_part, inb_part_txc), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for l_txc.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if
    end if

    return
end subroutine PARTICLE_ALLOCATE

! ###################################################################
! ###################################################################
subroutine PARTICLE_INITIALIZE()
    use TLAB_VARS, only: g, sbg
    use PARTICLE_VARS
    use PARTICLE_ARRAYS
    implicit none

    ! set boundarys for residence time pdf
    if (imode_part == PART_TYPE_BIL_CLOUD_4) then
        l_y_lambda = (g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(1)%ymean - 2.0_cp
        l_y_base = ((g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(1)%ymean - &
                    (g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(3)%ymean)/2.0_cp &
                   + (g(2)%nodes(g(2)%size) - g(2)%nodes(1))*sbg(3)%ymean
        if (residence_reset == 1) then
            l_q(:, 6:7) = 0.0_cp
        end if
    end if

    return
end subroutine PARTICLE_INITIALIZE
