#include "dns_error.h"

module TLabMPI_PROCS
    use MPI
    use TLab_Constants, only: wp, dp, sp, wi, lfile, efile
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLabMPI_VARS

    public :: TLabMPI_Initialize
    public :: TLabMPI_Panic

contains
    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) id, ip, npage
        integer(wi) dims(2), coord(2)
        logical period(2), remain_dims(2), reorder

        character(len=32) bakfile, block
        character(len=512) sRes, line
        character*64 lstr

        ! #######################################################################
        call TLab_Write_ASCII(lfile, 'Creating MPI communicators.')

        ! the first index in the grid corresponds to k, the second to i
        dims(1) = ims_npro_k; dims(2) = ims_npro_i; period = .true.; reorder = .false.
        ! dims(1) = ims_npro_i; dims(2) = ims_npro_k; period = .true.; reorder = .false.
        call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xz, ims_err)

        call MPI_CART_COORDS(ims_comm_xz, ims_pro, 2, coord, ims_err)
        ims_pro_k = coord(1); ims_pro_i = coord(2)      ! starting at 0
        ! ims_pro_k = coord(2); ims_pro_i = coord(1)
        !
        ! equivalent to:
        ! ims_pro_i = mod(ims_pro, ims_npro_i) ! Starting at 0
        ! ims_pro_k = ims_pro/ims_npro_i  ! Starting at 0
        ! to revert them:

        remain_dims(1) = .false.; remain_dims(2) = .true.
        call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)
        ! call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)

        remain_dims(1) = .true.; remain_dims(2) = .false.
        call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_z, ims_err)
        ! call MPI_CART_SUB(ims_comm_xz, remain_dims, ims_comm_x, ims_err)

        ! ip = ims_pro
        ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_x, ims_err)
        ! print*, 'P:', ims_pro, 'Sum along X', id
        ! CALL MPI_ALLREDUCE(ip, id, 1, MPI_INTEGER4, MPI_SUM, ims_comm_z, ims_err)
        ! print*, 'P:', ims_pro, 'Sum along Z', id

        ims_offset_i = ims_pro_i*imax       ! local offset in grid points
        ims_offset_j = 0
        ims_offset_k = ims_pro_k*kmax

        ! -----------------------------------------------------------------------
        ! Control of MPI type
        select case (wp)
        case (dp)
            TLAB_MPI_REAL_TYPE = MPI_REAL8
        case (sp)
            TLAB_MPI_REAL_TYPE = MPI_REAL4
        end select

    end subroutine TLabMPI_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Panic(location, mpi_error_code)
        character(len=*), intent(in) :: location
        integer, intent(in) :: mpi_error_code

        !##############################
        character error_string*1024
        integer error_local, error_len

        call MPI_Error_String(mpi_error_code, error_string, error_len, error_local)
        call TLab_Write_ASCII(efile, 'MPI-ERROR: Source file'//trim(adjustl(LOCATION)), .true.)
        call TLab_Write_ASCII(efile, error_string, .true.)

        call TLab_Stop(mpi_error_code)
        ! Not supposed to return from this subroutine

    end subroutine TLabMPI_Panic

    !########################################################################
    !# Moving plane information between adjacent PEs, circulant version.
    !# npl is smaller than 2*kmax
    !# The number of plane to move is given by npl
    !########################################################################
    subroutine TLabMPI_COPYPLN_1(ijmax, kmax, npl, a, bl, br)

        integer(wi) ijmax, kmax, npl
        real(wp) a(ijmax, *)
        real(wp) bl(ijmax, *)
        real(wp) br(ijmax, *)

        ! -----------------------------------------------------------------------
        integer status(MPI_STATUS_SIZE, 4)
        integer mpireq(4)
        integer ims_pro_l, ims_pro_r
        integer icount

        ! #######################################################################
        if (ims_npro > 1) then

            ! Careful in case only 2 PEs
            if (ims_npro == 2) then
                call TLab_Write_ASCII(efile, 'TLabMPI_COPYPLN_1. Undeveloped for 2 PEs.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            ! left and right PEs
            ims_pro_l = mod(ims_pro - 1 + ims_npro, ims_npro)
            ims_pro_r = mod(ims_pro + 1 + ims_npro, ims_npro)

            icount = ijmax*npl

            call MPI_IRECV(bl(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
            call MPI_IRECV(br(1, 1), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
            call MPI_ISEND(a(1, kmax + 1 - npl), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

            call MPI_WAITALL(4, mpireq, status, ims_err)

        end if

        return
    end subroutine TLabMPI_COPYPLN_1

    !########################################################################
    !########################################################################
    subroutine TLabMPI_COPYPLN_2(ijmax, kmax, npl, a, bl, br)

        integer(wi) ijmax, kmax, npl
        real(wp) a(ijmax, kmax)
        real(wp) bl(ijmax, npl)
        real(wp) br(ijmax, npl)

        ! -----------------------------------------------------------------------
        integer(wi) npl_loc
        integer status(MPI_STATUS_SIZE, 8)
        integer mpireq(8)
        integer ims_pro_l, ims_pro_r
        integer icount

        ! #######################################################################
        if (ims_npro > 1) then

            ! Careful in case only 2 PEs
            if (ims_npro == 2) then
                call TLab_Write_ASCII(efile, 'TLabMPI_COPYPLN_2. Undeveloped for 2 PEs.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            ! -----------------------------------------------------------------------
            ! left and right PEs. Same as in routine TLabMPI_COPYPLN_1
            ! -----------------------------------------------------------------------
            npl_loc = kmax

            ims_pro_l = mod(ims_pro - 1 + ims_npro, ims_npro)
            ims_pro_r = mod(ims_pro + 1 + ims_npro, ims_npro)

            icount = ijmax*npl_loc

            call MPI_IRECV(bl(1, npl - kmax + 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
            call MPI_IRECV(br(1, 1), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

            call MPI_WAITALL(4, mpireq, status, ims_err)

            ! -----------------------------------------------------------------------
            ! second-left and second-right PEs.
            ! -----------------------------------------------------------------------
            npl_loc = npl - kmax

            ims_pro_l = mod(ims_pro - 2 + ims_npro, ims_npro)
            ims_pro_r = mod(ims_pro + 2 + ims_npro, ims_npro)

            icount = ijmax*npl_loc

            call MPI_IRECV(bl(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(7), ims_err)
            call MPI_IRECV(br(1, 1 + kmax), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(8), ims_err)

            call MPI_ISEND(a(1, 1), icount, MPI_REAL8, ims_pro_l, &
                           ims_tag, MPI_COMM_WORLD, mpireq(5), ims_err)
            call MPI_ISEND(a(1, kmax + 1 - npl_loc), icount, MPI_REAL8, ims_pro_r, &
                           ims_tag, MPI_COMM_WORLD, mpireq(6), ims_err)

            call MPI_WAITALL(4, mpireq(5:), status, ims_err)

        end if

        return
    end subroutine TLabMPI_COPYPLN_2

    ! ######################################################################
    ! Initialization of PSFFT Library for nonblocking communication
    ! ######################################################################
#ifdef USE_PSFFT
#include "nb3dfft_defines.inc"
#endif

    subroutine DNS_NB3DFFT_INITIALIZE
#ifdef USE_PSFFT
        use NB3DFFT, only: nb3dfft_test_setup, nb3dfft_setup, get_dims
#endif

        implicit none

        ! #######################################################################
#ifdef USE_PSFFT
        call TLab_Write_ASCII(lfile, 'Initialize nonblocking communication.')

        ims_nb_proc_grid = (/ims_npro_i, ims_npro_k/)
        call NB3DFFT_SETUP(ims_nb_proc_grid, g(1)%size, g(2)%size, g(3)%size, &
                           ims_nb_msize)

        call GET_DIMS(ims_nb_xsrt, ims_nb_xend, ims_nb_xsiz, 1, 1)
        call GET_DIMS(ims_nb_ysrt, ims_nb_yend, ims_nb_ysiz, 1, 2)
        call GET_DIMS(ims_nb_zsrt, ims_nb_zend, ims_nb_zsiz, 1, 3)

        if (ims_nb_xsrt(1) == 1 .and. ims_nb_xend(1) == g(1)%size &
            ! .and. ims_nb_xsiz(2)*ims_nb_xsiz(3) == ims_size_i(TLAB_MPI_TRP_I_PARTIAL)) then
            .and. ims_nb_xsiz(2)*ims_nb_xsiz(3) == ims_plan_dx%nlines) then
            ! Decomp standing in X okay
        else
            call TLab_Write_ASCII(efile, 'Decomp standing in X-BAD')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        if (ims_nb_ysrt(1) == ims_offset_i + 1 &
            .and. ims_nb_ysrt(2) == ims_offset_j + 1 &
            .and. ims_nb_ysrt(3) == ims_offset_k + 1 &
            .and. ims_nb_ysiz(1) == imax &
            .and. ims_nb_ysiz(2) == jmax &
            .and. ims_nb_ysiz(3) == kmax) then
        else
            call TLab_Write_ASCII(efile, 'Decomp standing in Y--BAD')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        if (ims_nb_zsrt(3) == 1 .and. ims_nb_zend(3) == g(3)%size &
            ! .and. ims_nb_zsiz(1)*ims_nb_zsiz(2) == ims_size_k(TLAB_MPI_TRP_K_PARTIAL)) then
            .and. ims_nb_zsiz(1)*ims_nb_zsiz(2) == ims_plan_dz%nlines) then
            ! Decomp standing in Z okay
        else
            call TLab_Write_ASCII(efile, 'Decomp standing in Z--BAD')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        call TLab_Write_ASCII(lfile, 'Checking that NB3DFFT and DNS domain decompositions agree.')

        call nb3dfft_test_setup()

#else
        call TLab_Write_ASCII(efile, 'Compiler flag USE_PSFFT needs to be used.')
        call TLab_Stop(DNS_ERROR_PARPARTITION)

#endif

        return
    end subroutine DNS_NB3DFFT_INITIALIZE

end module TLabMPI_PROCS
