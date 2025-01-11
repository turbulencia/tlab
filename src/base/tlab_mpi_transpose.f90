#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

! Circular transposition across directional communicators
module TLabMPI_Transpose
    use MPI
    use TLab_Constants, only: lfile, efile, wp, dp, sp, wi, sizeofreal
    use TLAB_VARS, only: imax, jmax, kmax, isize_wrk3d, isize_txc_dimx, isize_txc_dimz
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use TLabMPI_VARS
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    public :: TLabMPI_Transpose_Initialize
    public :: TLabMPI_Trp_TypeI_Create, TLabMPI_Trp_TypeK_Create
    public :: TLabMPI_TransposeK_Forward, TLabMPI_TransposeK_Backward
    public :: TLabMPI_TransposeI_Forward, TLabMPI_TransposeI_Backward

    type, public :: tmpi_transpose_dt
        sequence
        integer :: type_s, type_r                           ! send/recv types
        integer(wi) :: nlines                               !
        integer(wi) :: size3d                               !
        integer(wi), allocatable :: disp_s(:), disp_r(:)    ! send/recv displacements
    end type tmpi_transpose_dt
    type(tmpi_transpose_dt), public :: ims_trp_plan_i(TLAB_MPI_TRP_I_MAXTYPES)
    type(tmpi_transpose_dt), public :: ims_trp_plan_k(TLAB_MPI_TRP_K_MAXTYPES)

    integer :: ims_trp_mode_i, ims_trp_mode_k               ! Mode of transposition
    integer, parameter :: TLAB_MPI_TRP_NONE = 0
    integer, parameter :: TLAB_MPI_TRP_ASYNCHRONOUS = 1
    integer, parameter :: TLAB_MPI_TRP_SENDRECV = 2
    integer, parameter :: TLAB_MPI_TRP_ALLTOALL = 3

    integer(wi) :: ims_sizBlock_i, ims_sizBlock_k                   ! group sizes of rend/recv messagest to use explicit sed/recv
    integer(wi), allocatable :: maps_send_i(:), maps_recv_i(:)      ! PE maps to use explicit sed/recv
    integer(wi), allocatable :: maps_send_k(:), maps_recv_k(:)

    integer, allocatable :: counts(:), types_send(:), types_recv(:) ! to use alltoallw

    integer, allocatable :: ims_status(:, :)
    integer, allocatable :: ims_request(:)

    integer :: ims_trp_type_i, ims_trp_type_k               ! Tranposition in double or single precission

    real(wp), allocatable, target :: wrk_mpi(:)             ! 3D work array for MPI; maybe in tlab_memory
    real(sp), pointer :: a_wrk(:) => null(), b_wrk(:) => null()

contains

    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Transpose_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) ip, npage

        character(len=32) bakfile, block
        character(len=512) sRes, line
        character*64 lstr

        ! #######################################################################
        ! Read data
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Parallel'

        call ScanFile_Char(bakfile, inifile, block, 'TransposeModeI', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'ComModeITranspose', 'asynchronous', sRes)
        if (trim(adjustl(sRes)) == 'none') then; ims_trp_mode_i = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; ims_trp_mode_i = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; ims_trp_mode_i = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; ims_trp_mode_i = TLAB_MPI_TRP_ALLTOALL
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeModeI option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeModeK', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'ComModeKTranspose', 'asynchronous', sRes)
        if (trim(adjustl(sRes)) == 'none') then; ims_trp_mode_k = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; ims_trp_mode_k = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; ims_trp_mode_k = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; ims_trp_mode_k = TLAB_MPI_TRP_ALLTOALL
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeModeK option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeK', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; ims_trp_type_k = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; ims_trp_type_k = MPI_REAL4
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeTypeK.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeI', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; ims_trp_type_i = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; ims_trp_type_i = MPI_REAL4
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeTypeI.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! #######################################################################
        ! Initialize

        ! Size of communication in explicit send/recv
#ifdef HLRS_HAWK
        ! On hawk, we tested that 192 yields optimum performace;
        ! Blocking will thus only take effect in very large cases
        ims_sizBlock_k = 192
        ims_sizBlock_i = 384
#else
        ! We assume that this will help to release some of the very heavy
        ! network load in transpositions on most systems
        ims_sizBlock_k = 64
        ims_sizBlock_i = 128
        ! ims_sizBlock_k=1e5   -- would essentially switch off the blocking
#endif

        if (ims_npro_i > ims_sizBlock_i) then
            write (line, *) ims_sizBlock_i
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_I'
            call TLab_Write_ASCII(lfile, line)
        end if

        if (ims_npro_k > ims_sizBlock_k) then
            write (line, *) ims_sizBlock_k
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_K'
            call TLab_Write_ASCII(lfile, line)
        end if

        allocate (ims_status(MPI_STATUS_SIZE, 2*max(ims_sizBlock_i, ims_sizBlock_k, ims_npro_i, ims_npro_k)))
        allocate (ims_request(2*max(ims_sizBlock_i, ims_sizBlock_k, ims_npro_i, ims_npro_k)))

        ! -----------------------------------------------------------------------
        ! local PE mappings for explicit send/recv
        allocate (maps_send_i(ims_npro_i))
        allocate (maps_recv_i(ims_npro_i))
        do ip = 0, ims_npro_i - 1
            maps_send_i(ip + 1) = ip
            maps_recv_i(ip + 1) = mod(ims_npro_i - ip, ims_npro_i)
        end do
        maps_send_i = cshift(maps_send_i, ims_pro_i)
        maps_recv_i = cshift(maps_recv_i, -ims_pro_i)

        allocate (maps_send_k(ims_npro_k))
        allocate (maps_recv_k(ims_npro_k))
        do ip = 0, ims_npro_k - 1
            maps_send_k(ip + 1) = ip
            maps_recv_k(ip + 1) = mod(ims_npro_k - ip, ims_npro_k)
        end do
        maps_send_k = cshift(maps_send_k, ims_pro_k)
        maps_recv_k = cshift(maps_recv_k, -ims_pro_k)

        ! -----------------------------------------------------------------------
        ! to use alltoallw
        allocate (counts(max(ims_npro_i, ims_npro_j, ims_npro_k)))
        allocate (types_send(max(ims_npro_i, ims_npro_j, ims_npro_k)))
        allocate (types_recv(max(ims_npro_i, ims_npro_j, ims_npro_k)))
        counts(:) = 1

        ! -----------------------------------------------------------------------
        ! to use single transposition when running in double precission
        call TLab_Allocate_Real(__FILE__, wrk_mpi, [isize_wrk3d], 'wrk-mpi')

        ! -----------------------------------------------------------------------
        ! Create basic transposition plans used for partial X and partial Z; could be in another module...
        if (ims_npro_i > 1) then
            npage = kmax*jmax
            ims_trp_plan_i(TLAB_MPI_TRP_I_PARTIAL) = TLabMPI_Trp_TypeI_Create(imax, npage, 1, 1, 1, 1, 'Ox derivatives.')
        end if

        if (ims_npro_k > 1) then
            npage = imax*jmax
            ims_trp_plan_k(TLAB_MPI_TRP_K_PARTIAL) = TLabMPI_Trp_TypeK_Create(kmax, npage, 1, 1, 1, 1, 'Oz derivatives.')
        end if

        return
    end subroutine TLabMPI_Transpose_Initialize

    ! ######################################################################
    ! ######################################################################
    ! Pointers and types for transposition across processors
    function TLabMPI_Trp_TypeI_Create(nmax, npage, nd, md, n1, n2, message) result(trp_plan)
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in) :: nd, md, n1, n2
        character(len=*), intent(in), optional :: message
        type(tmpi_transpose_dt) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) i
        integer ims_tmp1, ims_tmp2, ims_tmp3
        integer ims_ss, ims_rs
        character*64 str, line

        ! #######################################################################
        if (present(message)) &
            call TLab_Write_ASCII(lfile, 'Creating derived MPI types for '//trim(adjustl(message)))

        if (mod(npage, ims_npro_i) == 0) then
            trp_plan%nlines = npage/ims_npro_i
            allocate (trp_plan%disp_s(ims_npro_i), trp_plan%disp_r(ims_npro_i))
            trp_plan%size3d = npage*nmax
        else
            call TLab_Write_ASCII(efile, 'TLabMPI_TypeI_Create. Ratio npage/npro not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        ! Calculate array displacements in Forward Send/Receive
        trp_plan%disp_s(1) = 0
        trp_plan%disp_r(1) = 0
        do i = 2, ims_npro_i
            trp_plan%disp_s(i) = trp_plan%disp_s(i - 1) + nmax*nd*trp_plan%nlines
            trp_plan%disp_r(i) = trp_plan%disp_r(i - 1) + nmax*md
        end do

        ! #######################################################################
        ims_tmp1 = trp_plan%nlines*n1 ! count
        ims_tmp2 = nmax*n2 ! block
        ims_tmp3 = ims_tmp2  ! stride = block because things are together
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, ims_trp_type_i, trp_plan%type_s, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_s, ims_err)

        ims_tmp1 = trp_plan%nlines*n1 ! count
        ims_tmp2 = nmax*n2 ! block
        ims_tmp3 = nmax*ims_npro_i*n2 ! stride is a multiple of nmax_total=nmax*ims_npro_i
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, ims_trp_type_i, trp_plan%type_r, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_r, ims_err)

        call MPI_TYPE_SIZE(trp_plan%type_s, ims_ss, ims_err)
        call MPI_TYPE_SIZE(trp_plan%type_r, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end function TLabMPI_Trp_TypeI_Create

    !########################################################################
    !########################################################################
    function TLabMPI_Trp_TypeK_Create(nmax, npage, nd, md, n1, n2, message) result(trp_plan)
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in) :: nd, md, n1, n2
        character(len=*), intent(in), optional :: message
        type(tmpi_transpose_dt) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) i
        integer ims_tmp1, ims_tmp2, ims_tmp3
        integer ims_ss, ims_rs
        character*64 str, line

        ! #######################################################################
        if (present(message)) &
            call TLab_Write_ASCII(lfile, 'Creating derived MPI types for '//trim(adjustl(message)))

        if (mod(npage, ims_npro_k) == 0) then
            trp_plan%nlines = npage/ims_npro_k
            allocate (trp_plan%disp_s(ims_npro_k), trp_plan%disp_r(ims_npro_k))
            trp_plan%size3d = npage*nmax
        else
            call TLab_Write_ASCII(efile, 'TLabMPI_TypeI_Create. Ratio npage/npro not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        ! Calculate array displacements in Forward Send/Receive
        trp_plan%disp_s(1) = 0
        trp_plan%disp_r(1) = 0
        do i = 2, ims_npro_k
            trp_plan%disp_s(i) = trp_plan%disp_s(i - 1) + trp_plan%nlines*nd
            trp_plan%disp_r(i) = trp_plan%disp_r(i - 1) + trp_plan%nlines*md*nmax
        end do

        ! #######################################################################
        ims_tmp1 = nmax*n1                  ! count
        ims_tmp2 = trp_plan%nlines*n2       ! block
        ims_tmp3 = npage*n2                 ! stride
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, ims_trp_type_k, trp_plan%type_s, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_s, ims_err)

        ims_tmp1 = nmax*n1                  ! count
        ims_tmp2 = trp_plan%nlines*n2       ! block
        ims_tmp3 = ims_tmp2                 ! stride = block to put things together
        call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, ims_trp_type_k, trp_plan%type_r, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_r, ims_err)

        call MPI_TYPE_SIZE(trp_plan%type_s, ims_ss, ims_err)
        call MPI_TYPE_SIZE(trp_plan%type_r, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end function TLabMPI_Trp_TypeK_Create

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TransposeK_Forward(a, b, trp_plan)
        ! use, intrinsic :: iso_c_binding, only: c_int

        ! interface
        !     function DNS_USLEEP(useconds) bind(C, name="usleep")
        !         import
        !         integer(c_int) :: nb3dfft_nbc_usleep
        !         integer(c_int), intent(in), value :: useconds
        !     end function DNS_USLEEP
        ! end interface

        real(wp), dimension(*), intent(in) :: a
        real(wp), dimension(*), intent(out) :: b
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        target b

        ! -----------------------------------------------------------------------
        integer(wi) size

#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif

        if (ims_trp_type_k == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call Transpose_Kernel_Single(a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         ims_comm_z, ims_sizBlock_k, ims_trp_mode_k)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                  b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                  ims_comm_z, ims_sizBlock_k, ims_trp_mode_k)
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_TransposeK_Forward

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TransposeK_Backward(b, a, trp_plan)
        real(wp), dimension(*), intent(in) :: b
        real(wp), dimension(*), intent(out) :: a
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        target a

        ! -----------------------------------------------------------------------
        integer(wi) size
#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif

        if (ims_trp_type_k == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(a), b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), a_wrk, shape=[size])
            b_wrk(1:size) = real(b(1:size), sp)
            call Transpose_Kernel_Single(b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         ims_comm_z, ims_sizBlock_k, ims_trp_mode_k)
            a(1:size) = real(a_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                  a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                  ims_comm_z, ims_sizBlock_k, ims_trp_mode_k)
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_TransposeK_Backward

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TransposeI_Forward(a, b, trp_plan)
        real(wp), dimension(*), intent(in) :: a
        real(wp), dimension(*), intent(out) :: b
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        target b

        ! -----------------------------------------------------------------------
        integer(wi) size

        ! #######################################################################
        if (ims_trp_type_i == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call Transpose_Kernel_Single(a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         ims_comm_x, ims_sizBlock_i, ims_trp_mode_i)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                  b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                  ims_comm_x, ims_sizBlock_i, ims_trp_mode_i)
        end if

        return
    end subroutine TLabMPI_TransposeI_Forward

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TransposeI_Backward(b, a, trp_plan)
        real(wp), dimension(*), intent(in) :: b
        real(wp), dimension(*), intent(out) :: a
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        target a

        ! -----------------------------------------------------------------------
        integer(wi) size

        ! #######################################################################
        if (ims_trp_type_i == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(a), b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), a_wrk, shape=[size])
            b_wrk(1:size) = real(b(1:size), sp)
            call Transpose_Kernel_Single(b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         ims_comm_x, ims_sizBlock_i, ims_trp_mode_i)
            a(1:size) = real(a_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                  a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                  ims_comm_x, ims_sizBlock_i, ims_trp_mode_i)
        end if

        return
    end subroutine TLabMPI_TransposeI_Backward

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        real(wp), intent(in) :: a(*)
        real(wp), intent(out) :: b(*)

        integer, intent(in) :: comm                         ! communicator
        integer, intent(in) :: tsend, trecv                 ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)       ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)       ! maps send/receive
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer(wi) npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(dsend(:))

        select case (mode)
        case (TLAB_MPI_TRP_ASYNCHRONOUS)
            do j = 1, npro, step
                l = 0
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_ISEND(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, comm, ims_request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, ims_request(l), ims_err)
                end do
                call MPI_WAITALL(l, ims_request, ims_status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, ims_status(:, 1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = tsend
            types_recv(1:npro) = trecv
            call MPI_ALLTOALLW(a, counts, dsend*int(sizeof(1.0_wp)), types_send, &
                               b, counts, drecv*int(sizeof(1.0_wp)), types_recv, comm, ims_err)
            ! call MPI_ALLTOALLW(a, spread(1, 1, npro), dsend*int(sizeof(1.0_wp)), spread(tsend, 1, npro), &
            !                    b, spread(1, 1, npro), drecv*int(sizeof(1.0_wp)), spread(trecv, 1, npro), comm, ims_err)
        end select

        return
    end subroutine Transpose_Kernel

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Single(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        real(sp), intent(in) :: a(*)
        real(sp), intent(out) :: b(*)

        integer, intent(in) :: comm                         ! communicator
        integer, intent(in) :: tsend, trecv                 ! types send/receive
        integer(wi), intent(in) :: dsend(:), drecv(:)       ! displacements send/receive
        integer(wi), intent(in) :: msend(:), mrecv(:)       ! maps send/receive
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer(wi) npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(dsend(:))

        select case (mode)
        case (TLAB_MPI_TRP_ASYNCHRONOUS)
            do j = 1, npro, step
                l = 0
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_ISEND(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, comm, ims_request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, ims_request(l), ims_err)
                end do
                call MPI_WAITALL(l, ims_request, ims_status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, ims_status(:, 1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = tsend
            types_recv(1:npro) = trecv
            call MPI_ALLTOALLW(a, counts, dsend*int(sizeof(1.0_sp)), types_send, &
                               b, counts, drecv*int(sizeof(1.0_sp)), types_recv, comm, ims_err)
            ! call MPI_ALLTOALLW(a, spread(1, 1, npro), dsend*int(sizeof(1.0_wp)), spread(tsend, 1, npro), &
            !                    b, spread(1, 1, npro), drecv*int(sizeof(1.0_wp)), spread(trecv, 1, npro), comm, ims_err)

        end select

        return
    end subroutine Transpose_Kernel_Single

end module TLabMPI_Transpose
