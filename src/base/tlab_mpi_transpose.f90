#include "dns_error.h"

! Circular transposition within directional communicators
module TLabMPI_Transpose
    use TLab_Constants, only: lfile, efile, wp, dp, sp, wi, sizeofreal
    use TLab_Memory, only: imax, jmax, kmax, isize_wrk3d
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    use TLabMPI_VARS
    implicit none
    private

    public :: TLabMPI_Trp_Initialize
    public :: TLabMPI_Trp_PlanI, TLabMPI_Trp_PlanK
    public :: TLabMPI_Trp_ExecK_Forward, TLabMPI_Trp_ExecK_Backward
    public :: TLabMPI_Trp_ExecI_Forward, TLabMPI_Trp_ExecI_Backward

    type, public :: tmpi_transpose_dt
        ! sequence
        type(MPI_Datatype) :: type_s, type_r                        ! derived send/recv types
        integer(wi) :: nlines
        integer(wi) :: size3d
        integer(wi), allocatable :: disp_s(:), disp_r(:)            ! send/recv displacements
    end type tmpi_transpose_dt
    type(tmpi_transpose_dt), public :: tmpi_plan_dx                 ! general plans used in derivatives and other operators
    type(tmpi_transpose_dt), public :: tmpi_plan_dz

    ! -----------------------------------------------------------------------
    integer :: trp_mode_i, trp_mode_k                               ! Mode of transposition
    integer, parameter :: TLAB_MPI_TRP_NONE = 0
    integer, parameter :: TLAB_MPI_TRP_ASYNCHRONOUS = 1
    integer, parameter :: TLAB_MPI_TRP_SENDRECV = 2
    integer, parameter :: TLAB_MPI_TRP_ALLTOALL = 3

    integer(wi) :: trp_sizBlock_i, trp_sizBlock_k                   ! explicit sed/recv: group sizes of rend/recv messages
    integer(wi), allocatable :: maps_send_i(:), maps_recv_i(:)      ! PE maps to use explicit sed/recv
    integer(wi), allocatable :: maps_send_k(:), maps_recv_k(:)
    type(MPI_Datatype), allocatable :: types_send(:), types_recv(:) ! alltoallw
    integer, allocatable :: counts(:)

    type(MPI_Datatype) :: trp_datatype_i, trp_datatype_k            ! Transposition in double or single precision

    real(wp), allocatable, target :: wrk_mpi(:)                     ! 3D work array
    real(sp), pointer :: a_wrk(:) => null(), b_wrk(:) => null()

    type(MPI_Status), allocatable :: status(:)
    type(MPI_Request), allocatable :: request(:)

    interface TLabMPI_Trp_ExecK_Forward
        module procedure TLabMPI_Trp_ExecK_Forward_Real, TLabMPI_Trp_ExecK_Forward_Complex
    end interface TLabMPI_Trp_ExecK_Forward
    interface TLabMPI_Trp_ExecK_Backward
        module procedure TLabMPI_Trp_ExecK_Backward_Real, TLabMPI_Trp_ExecK_Backward_Complex
    end interface TLabMPI_Trp_ExecK_Backward

    interface TLabMPI_Trp_ExecI_Forward
        module procedure TLabMPI_Trp_ExecI_Forward_Real, TLabMPI_Trp_ExecI_Forward_Complex
    end interface TLabMPI_Trp_ExecI_Forward
    interface TLabMPI_Trp_ExecI_Backward
        module procedure TLabMPI_Trp_ExecI_Backward_Real, TLabMPI_Trp_ExecI_Backward_Complex
    end interface TLabMPI_Trp_ExecI_Backward

contains

    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Trp_Initialize(inifile)
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
        if (trim(adjustl(sRes)) == 'none') then; trp_mode_i = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; trp_mode_i = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; trp_mode_i = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; trp_mode_i = TLAB_MPI_TRP_ALLTOALL
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeModeI option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeModeK', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'ComModeKTranspose', 'asynchronous', sRes)
        if (trim(adjustl(sRes)) == 'none') then; trp_mode_k = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; trp_mode_k = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; trp_mode_k = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; trp_mode_k = TLAB_MPI_TRP_ALLTOALL
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeModeK option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeK', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; trp_datatype_k = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; trp_datatype_k = MPI_REAL4
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeTypeK.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeI', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; trp_datatype_i = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; trp_datatype_i = MPI_REAL4
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
        trp_sizBlock_k = 192
        trp_sizBlock_i = 384
#else
        ! We assume that this will help to release some of the very heavy
        ! network load in transpositions on most systems
        trp_sizBlock_k = 64
        trp_sizBlock_i = 128
        ! trp_sizBlock_k=1e5   -- would essentially switch off the blocking
#endif

        if (ims_npro_i > trp_sizBlock_i) then
            write (line, *) trp_sizBlock_i
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_I'
            call TLab_Write_ASCII(lfile, line)
        end if

        if (ims_npro_k > trp_sizBlock_k) then
            write (line, *) trp_sizBlock_k
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_K'
            call TLab_Write_ASCII(lfile, line)
        end if

        allocate (status(2*max(trp_sizBlock_i, trp_sizBlock_k, ims_npro_i, ims_npro_k)))
        allocate (request(2*max(trp_sizBlock_i, trp_sizBlock_k, ims_npro_i, ims_npro_k)))

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
            tmpi_plan_dx = TLabMPI_Trp_PlanI(imax, npage, message='Ox derivatives.')
        end if

        if (ims_npro_k > 1) then
            npage = imax*jmax
            tmpi_plan_dz = TLabMPI_Trp_PlanK(kmax, npage, message='Oz derivatives.')
        end if

        return
    end subroutine TLabMPI_Trp_Initialize

    ! ######################################################################
    ! ######################################################################
    ! Pointers and types for transposition across processors
    function TLabMPI_Trp_PlanI(nmax, npage, locStride, locType, message) result(trp_plan)
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in), optional :: locStride
        type(MPI_Datatype), intent(in), optional :: locType
        character(len=*), intent(in), optional :: message
        type(tmpi_transpose_dt) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) i
        type(MPI_Datatype) :: datatype
        integer block_count, block_length, stride
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

        block_count = trp_plan%nlines
        block_length = nmax

        ! Calculate array displacements in Forward Send/Receive
        trp_plan%disp_s(1) = 0
        trp_plan%disp_r(1) = 0
        do i = 2, ims_npro_i
            trp_plan%disp_s(i) = trp_plan%disp_s(i - 1) + block_length*block_count
            trp_plan%disp_r(i) = trp_plan%disp_r(i - 1) + block_length
        end do

        ! #######################################################################
        if (present(locType)) then
            datatype = locType
        else
            datatype = trp_datatype_i
        end if

        stride = block_length       ! stride = block_length because things are together
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_s, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_s, ims_err)

        stride = nmax*ims_npro_i    ! stride is a multiple of nmax_total=nmax*ims_npro_i
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_r, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_r, ims_err)

        ! -----------------------------------------------------------------------
        call MPI_TYPE_SIZE(trp_plan%type_s, ims_ss, ims_err)
        call MPI_TYPE_SIZE(trp_plan%type_r, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end function TLabMPI_Trp_PlanI

    !########################################################################
    !########################################################################
    function TLabMPI_Trp_PlanK(nmax, npage, locStride, locType, message) result(trp_plan)
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in), optional :: locStride
        type(MPI_Datatype), intent(in), optional :: locType
        character(len=*), intent(in), optional :: message
        type(tmpi_transpose_dt) :: trp_plan

        ! -----------------------------------------------------------------------
        integer(wi) i
        type(MPI_Datatype) :: datatype
        integer block_count, block_length, stride
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

        block_count = nmax
        block_length = trp_plan%nlines

        ! Calculate array displacements in Forward Send/Receive
        trp_plan%disp_s(1) = 0
        trp_plan%disp_r(1) = 0
        do i = 2, ims_npro_k
            trp_plan%disp_s(i) = trp_plan%disp_s(i - 1) + block_length
            trp_plan%disp_r(i) = trp_plan%disp_r(i - 1) + block_length*block_count
        end do

        ! #######################################################################
        if (present(locType)) then
            datatype = locType
        else
            datatype = trp_datatype_i
        end if

        stride = npage
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_s, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_s, ims_err)

        stride = block_length       ! stride = block_length to put things together
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, trp_plan%type_r, ims_err)
        call MPI_TYPE_COMMIT(trp_plan%type_r, ims_err)

        ! -----------------------------------------------------------------------
        call MPI_TYPE_SIZE(trp_plan%type_s, ims_ss, ims_err)
        call MPI_TYPE_SIZE(trp_plan%type_r, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        return
    end function TLabMPI_Trp_PlanK

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Forward_Real(a, b, trp_plan)
        real(wp), intent(in) :: a(*)
        real(wp), intent(out) :: b(*)
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

        if (trp_datatype_k == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call Transpose_Kernel_Single(a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         ims_comm_z, trp_sizBlock_k, trp_mode_k)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel_Double(a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         ims_comm_z, trp_sizBlock_k, trp_mode_k)
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_Trp_ExecK_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Forward_Complex(a, b, trp_plan)
        complex(wp), intent(in) :: a(*)
        complex(wp), intent(out) :: b(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! #######################################################################
        call Transpose_Kernel_Complex(a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                      b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                      ims_comm_z, trp_sizBlock_k, trp_mode_k)

        return
    end subroutine TLabMPI_Trp_ExecK_Forward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Backward_Real(b, a, trp_plan)
        real(wp), intent(in) :: b(*)
        real(wp), intent(out) :: a(*)
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

        if (trp_datatype_k == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(a), b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), a_wrk, shape=[size])
            b_wrk(1:size) = real(b(1:size), sp)
            call Transpose_Kernel_Single(b_wrk, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         a_wrk, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         ims_comm_z, trp_sizBlock_k, trp_mode_k)
            a(1:size) = real(a_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel_Double(b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         ims_comm_z, trp_sizBlock_k, trp_mode_k)
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine TLabMPI_Trp_ExecK_Backward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecK_Backward_Complex(b, a, trp_plan)
        complex(wp), intent(in) :: b(*)
        complex(wp), intent(out) :: a(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! #######################################################################
        call Transpose_Kernel_Complex(b, maps_recv_k(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                      a, maps_send_k(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                      ims_comm_z, trp_sizBlock_k, trp_mode_k)

        return
    end subroutine TLabMPI_Trp_ExecK_Backward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Forward_Real(a, b, trp_plan)
        real(wp), dimension(*), intent(in) :: a
        real(wp), dimension(*), intent(out) :: b
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        target b

        ! -----------------------------------------------------------------------
        integer(wi) size

        ! #######################################################################
        if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call Transpose_Kernel_Single(a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         ims_comm_x, trp_sizBlock_i, trp_mode_i)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel_Double(a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         ims_comm_x, trp_sizBlock_i, trp_mode_i)
        end if

        return
    end subroutine TLabMPI_Trp_ExecI_Forward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Forward_Complex(a, b, trp_plan)
        complex(wp), intent(in) :: a(*)
        complex(wp), intent(out) :: b(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! #######################################################################
        call Transpose_Kernel_Complex(a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                      b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                      ims_comm_x, trp_sizBlock_i, trp_mode_i)

        return
    end subroutine TLabMPI_Trp_ExecI_Forward_Complex

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Backward_Real(b, a, trp_plan)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
        real(wp), intent(in) :: b(*)
        real(wp), intent(out) :: a(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        target a

        ! -----------------------------------------------------------------------
        integer(wi) size

        ! #######################################################################
        if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
            size = trp_plan%size3d
            call c_f_pointer(c_loc(a), b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), a_wrk, shape=[size])
            b_wrk(1:size) = real(b(1:size), sp)
            call Transpose_Kernel_Single(b_wrk, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         a_wrk, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         ims_comm_x, trp_sizBlock_i, trp_mode_i)
            a(1:size) = real(a_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel_Double(b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                         a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                         ims_comm_x, trp_sizBlock_i, trp_mode_i)
        end if

        return
    end subroutine TLabMPI_Trp_ExecI_Backward_Real

    !########################################################################
    !########################################################################
    subroutine TLabMPI_Trp_ExecI_Backward_Complex(b, a, trp_plan)
        complex(wp), intent(in) :: b(*)
        complex(wp), intent(out) :: a(*)
        type(tmpi_transpose_dt), intent(in) :: trp_plan

        ! #######################################################################
        call Transpose_Kernel_Complex(b, maps_recv_i(:), trp_plan%disp_r(:), trp_plan%type_r, &
                                      a, maps_send_i(:), trp_plan%disp_s(:), trp_plan%type_s, &
                                      ims_comm_x, trp_sizBlock_i, trp_mode_i)

        return
    end subroutine TLabMPI_Trp_ExecI_Backward_Complex

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Double(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        real(wp), intent(in) :: a(*)
        real(wp), intent(out) :: b(*)

        type(MPI_Comm), intent(in) :: comm                         ! communicator
        type(MPI_Datatype), intent(in) :: tsend, trecv                 ! types send/receive
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
                    call MPI_ISEND(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, comm, request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, status(1), ims_err)
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
    end subroutine Transpose_Kernel_Double

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Single(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        real(sp), intent(in) :: a(*)
        real(sp), intent(out) :: b(*)

        type(MPI_Comm), intent(in) :: comm                         ! communicator
        type(MPI_Datatype), intent(in) :: tsend, trecv                 ! types send/receive
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
                    call MPI_ISEND(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, comm, request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, status(1), ims_err)
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

    !########################################################################
    !########################################################################
    subroutine Transpose_Kernel_Complex(a, msend, dsend, tsend, b, mrecv, drecv, trecv, comm, step, mode)
        complex(wp), intent(in) :: a(*)
        complex(wp), intent(out) :: b(*)

        type(MPI_Comm), intent(in) :: comm                         ! communicator
        type(MPI_Datatype), intent(in) :: tsend, trecv                 ! types send/receive
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
                    call MPI_ISEND(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, comm, request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, status(1), ims_err)
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
    end subroutine Transpose_Kernel_Complex

end module TLabMPI_Transpose
