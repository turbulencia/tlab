#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

module TLabMPI_PROCS
    use MPI
    use TLab_Constants, only: lfile, efile, wp, dp, sp, wi
    use TLAB_VARS, only: imax, jmax, kmax, isize_wrk3d, isize_txc_dimx, isize_txc_dimz
    use TLAB_VARS, only: fourier_on
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use TLabMPI_VARS
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
    implicit none
    private

    public :: TLabMPI_Transpose_Initialize
    public :: TLabMPI_TypeK_Create, TLabMPI_TypeI_Create
    public :: TLabMPI_TransposeK_Forward, TLabMPI_TransposeK_Backward
    public :: TLabMPI_TransposeI_Forward, TLabMPI_TransposeI_Backward

    integer(wi), allocatable, public :: ims_size_i(:), ims_size_k(:)    ! Maybe public in module for TLabMPI_Transpose

    ! type, public :: mpi_transpose_dt
    !     sequence
    !     integer :: type_s, type_r                           ! send/recv types
    !     integer(wi), allocatable :: disp_s(:), disp_r(:)    ! send/recv displacements
    !     integer(wi) :: size
    ! end type mpi_transpose_dt

    integer(wi) :: ims_sizBlock_i, ims_sizBlock_k

    integer :: ims_trp_type_i, ims_trp_type_k
    integer, parameter :: TLAB_MPI_TRP_SINGLE = 1
    integer, parameter :: TLAB_MPI_TRP_DOUBLE = 2

    integer :: ims_trp_mode_i, ims_trp_mode_k
    integer, parameter :: TLAB_MPI_TRP_NONE = 0
    integer, parameter :: TLAB_MPI_TRP_ASYNCHRONOUS = 1
    integer, parameter :: TLAB_MPI_TRP_SENDRECV = 2
    integer, parameter :: TLAB_MPI_TRP_ALLTOALL = 3

    integer, dimension(:), allocatable :: ims_ts_i, ims_tr_i
    integer(wi), dimension(:, :), allocatable :: ims_ds_i, ims_dr_i
    integer(wi), allocatable :: ims_plan_trps_i(:), ims_plan_trpr_i(:)

    integer, dimension(:), allocatable :: ims_ts_k, ims_tr_k
    integer(wi), dimension(:, :), allocatable :: ims_ds_k, ims_dr_k
    integer(wi), allocatable :: ims_plan_trps_k(:), ims_plan_trpr_k(:)

    integer, allocatable :: ims_status(:, :)
    integer, allocatable :: ims_request(:)

    real(wp), allocatable, target :: wrk_mpi(:)      ! 3D work array for MPI
    real(sp), pointer :: a_wrk(:) => null(), b_wrk(:) => null()

contains

    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Transpose_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) id, ip, npage

        character(len=32) bakfile, block
        character(len=512) sRes, line
        character*64 lstr

        ! #######################################################################
        ! Circular transposition
        ! Maybe this should be a separate module
        ! #######################################################################
        ! -----------------------------------------------------------------------
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
        if (trim(adjustl(sRes)) == 'double') then; ims_trp_type_k = TLAB_MPI_TRP_DOUBLE
        elseif (trim(adjustl(sRes)) == 'single') then; ims_trp_type_k = TLAB_MPI_TRP_SINGLE
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeTypeK.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeI', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; ims_trp_type_i = TLAB_MPI_TRP_DOUBLE
        elseif (trim(adjustl(sRes)) == 'single') then; ims_trp_type_i = TLAB_MPI_TRP_SINGLE
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TransposeTypeI.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! -----------------------------------------------------------------------
        ! Allocation
        allocate (ims_ts_i(TLAB_MPI_TRP_I_MAXTYPES))         ! derived MPI types for send/recv
        allocate (ims_tr_i(TLAB_MPI_TRP_I_MAXTYPES))
        allocate (ims_size_i(TLAB_MPI_TRP_I_MAXTYPES))       ! metadata inside/to-calculate MPI types
        allocate (ims_ds_i(ims_npro_i, TLAB_MPI_TRP_I_MAXTYPES))
        allocate (ims_dr_i(ims_npro_i, TLAB_MPI_TRP_I_MAXTYPES))
        allocate (ims_plan_trps_i(ims_npro_i))          ! mappings for explicit send/recv
        allocate (ims_plan_trpr_i(ims_npro_i))

        allocate (ims_ts_k(TLAB_MPI_TRP_K_MAXTYPES))         ! derived MPI types for send/recv
        allocate (ims_tr_k(TLAB_MPI_TRP_K_MAXTYPES))
        allocate (ims_size_k(TLAB_MPI_TRP_K_MAXTYPES))       ! metadata inside/to-calculate MPI types
        allocate (ims_ds_k(ims_npro_k, TLAB_MPI_TRP_K_MAXTYPES))
        allocate (ims_dr_k(ims_npro_k, TLAB_MPI_TRP_K_MAXTYPES))
        allocate (ims_plan_trps_k(ims_npro_k))          ! mappings for explicit send/recv
        allocate (ims_plan_trpr_k(ims_npro_k))

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
            line = 'Using blocking of '//TRIM(ADJUSTL(line))//' in TLabMPI_TRP<F,B>_I'
            call TLab_Write_ASCII(lfile, line)
        end if

        if (ims_npro_k > ims_sizBlock_k) then
            write (line, *) ims_sizBlock_k
            line = 'Using blocking of '//TRIM(ADJUSTL(line))//' in TLabMPI_TRP<F,B>_K'
            call TLab_Write_ASCII(lfile, line)
        end if

        allocate (ims_status(MPI_STATUS_SIZE, 2*max(ims_sizBlock_i, ims_sizBlock_k, ims_npro_i, ims_npro_k)))
        allocate (ims_request(2*max(ims_sizBlock_i, ims_sizBlock_k, ims_npro_i, ims_npro_k)))

        call TLab_Allocate_Real(__FILE__, wrk_mpi, [isize_wrk3d], 'wrk-mpi')

        ! -----------------------------------------------------------------------
        ! Initialize
        if (ims_npro_i > 1) then
            call TLab_Write_ASCII(lfile, 'Creating MPI types for Ox derivatives.')
            id = TLAB_MPI_TRP_I_PARTIAL
            npage = kmax*jmax
            call TLabMPI_TypeI_Create(ims_npro_i, imax, npage, 1, 1, 1, 1, id)
        end if

        if (ims_npro_k > 1) then
            call TLab_Write_ASCII(lfile, 'Creating MPI types for Oz derivatives.')
            id = TLAB_MPI_TRP_K_PARTIAL
            npage = imax*jmax
            call TLabMPI_TypeK_Create(ims_npro_k, kmax, npage, 1, 1, 1, 1, id)
        end if

        ! -----------------------------------------------------------------------
        if (ims_npro_i > 1 .and. fourier_on) then
            call TLab_Write_ASCII(lfile, 'Creating MPI types for Ox FFTW in Poisson solver.')
            id = TLAB_MPI_TRP_I_POISSON1
            npage = isize_txc_dimx ! isize_txc_field/imax
            call TLabMPI_TypeI_Create(ims_npro_i, imax, npage, 1, 1, 1, 1, id)

            call TLab_Write_ASCII(lfile, 'Creating MPI types for Ox FFTW in Poisson solver.')
            id = TLAB_MPI_TRP_I_POISSON2 ! isize_txc_field/(imax+2)
            npage = isize_txc_dimx
            call TLabMPI_TypeI_Create(ims_npro_i, imax + 2, npage, 1, 1, 1, 1, id)

        end if

        if (ims_npro_k > 1 .and. fourier_on) then
            call TLab_Write_ASCII(lfile, 'Creating MPI types for Oz FFTW in Poisson solver.')
            id = TLAB_MPI_TRP_K_POISSON
            npage = isize_txc_dimz ! isize_txc_field/kmax
            call TLabMPI_TypeK_Create(ims_npro_k, kmax, npage, 1, 1, 1, 1, id)
        end if

        ! -----------------------------------------------------------------------
        ! local PE mappings for explicit send/recv
        do ip = 0, ims_npro_i - 1
            ims_plan_trps_i(ip + 1) = ip
            ims_plan_trpr_i(ip + 1) = mod(ims_npro_i - ip, ims_npro_i)
        end do
        ims_plan_trps_i = cshift(ims_plan_trps_i, ims_pro_i)
        ims_plan_trpr_i = cshift(ims_plan_trpr_i, -ims_pro_i)

        do ip = 0, ims_npro_k - 1
            ims_plan_trps_k(ip + 1) = ip
            ims_plan_trpr_k(ip + 1) = mod(ims_npro_k - ip, ims_npro_k)
        end do
        ims_plan_trps_k = cshift(ims_plan_trps_k, ims_pro_k)
        ims_plan_trpr_k = cshift(ims_plan_trpr_k, -ims_pro_k)

        ! do ip = 0, ims_npro_i - 1
        !     if (ims_pro == ip) then
        !         write (*, *) ims_pro, ims_pro_i, 'SEND:', ims_plan_trps_i
        !         write (*, *) ims_pro, ims_pro_i, 'RECV:', ims_plan_trpr_i
        !     end if
        !     call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
        ! end do

        call TLabMPI_TAGRESET

        return
    end subroutine TLabMPI_Transpose_Initialize

    ! ######################################################################
    ! Pointers and types for transposition across ims_npro processors
    ! ######################################################################
    subroutine TLabMPI_TypeI_Create(npro_i, nmax, npage, nd, md, n1, n2, id)
        integer, intent(in) :: npro_i
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in) :: nd, md, n1, n2
        integer, intent(in) :: id

        ! -----------------------------------------------------------------------
        integer(wi) i
        integer ims_ss, ims_rs
        integer ims_tmp1, ims_tmp2, ims_tmp3
        character*64 str, line

#define nsize       ims_size_i(id)
#define sdisp(j)    ims_ds_i(j, id)
#define rdisp(j)    ims_dr_i(j, id)
#define stype       ims_ts_i(id)
#define rtype       ims_tr_i(id)

        ! #######################################################################
        if (mod(npage, npro_i) == 0) then
            nsize = npage/npro_i
        else
            call TLab_Write_ASCII(efile, 'TLabMPI_TypeI_Create. Ratio npage/npro_i not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        ! Calculate array displacements in Forward Send/Receive
        sdisp(1) = 0
        rdisp(1) = 0
        do i = 2, npro_i
            sdisp(i) = sdisp(i - 1) + nmax*nd*nsize
            rdisp(i) = rdisp(i - 1) + nmax*md
        end do

        ! #######################################################################
        ims_tmp1 = nsize*n1 ! count
        ims_tmp2 = nmax*n2 ! block
        ims_tmp3 = ims_tmp2  ! stride = block because things are together
        select case (ims_trp_type_i)
        case (TLAB_MPI_TRP_DOUBLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype, ims_err)
        case (TLAB_MPI_TRP_SINGLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL4, stype, ims_err)
        end select
        call MPI_TYPE_COMMIT(stype, ims_err)

        ims_tmp1 = nsize*n1 ! count
        ims_tmp2 = nmax*n2 ! block
        ims_tmp3 = nmax*npro_i*n2 ! stride is a multiple of nmax_total=nmax*npro_i
        select case (ims_trp_type_i)
        case (TLAB_MPI_TRP_DOUBLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype, ims_err)
        case (TLAB_MPI_TRP_SINGLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL4, rtype, ims_err)
        end select
        call MPI_TYPE_COMMIT(rtype, ims_err)

        call MPI_TYPE_SIZE(stype, ims_ss, ims_err)
        call MPI_TYPE_SIZE(rtype, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            write (str, *) 1  ! i
            line = trim(adjustl(line))//' in message '//trim(adjustl(str))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

#undef nsize
#undef sdisp
#undef rdisp
#undef stype
#undef rtype

        return
    end subroutine TLabMPI_TypeI_Create

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TypeK_Create(npro_k, nmax, npage, nd, md, n1, n2, id)
        integer, intent(in) :: npro_k
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in) :: nd, md, n1, n2
        integer, intent(in) :: id

        ! -----------------------------------------------------------------------
        integer(wi) i
        integer ims_ss, ims_rs
        integer ims_tmp1, ims_tmp2, ims_tmp3

#define nsize       ims_size_k(id)
#define sdisp(j)    ims_ds_k(j, id)
#define rdisp(j)    ims_dr_k(j, id)
#define stype       ims_ts_k(id)
#define rtype       ims_tr_k(id)

! #######################################################################
        if (mod(npage, npro_k) == 0) then
            nsize = npage/npro_k
        else
            call TLab_Write_ASCII(efile, 'TLabMPI_TypeK_Create. Ratio npage/npro_k not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        ! Calculate array displacements in Forward Send/Receive
        sdisp(1) = 0
        rdisp(1) = 0
        do i = 2, npro_k
            sdisp(i) = sdisp(i - 1) + nsize*nd
            rdisp(i) = rdisp(i - 1) + nsize*md*nmax
        end do

        ! #######################################################################
        ims_tmp1 = nmax*n1 ! count
        ims_tmp2 = nsize*n2 ! block
        ims_tmp3 = npage*n2 ! stride
        select case (ims_trp_type_k)
        case (TLAB_MPI_TRP_DOUBLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, stype, ims_err)
        case (TLAB_MPI_TRP_SINGLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL4, stype, ims_err)
        end select
        call MPI_TYPE_COMMIT(stype, ims_err)

        ims_tmp1 = nmax*n1 ! count
        ims_tmp2 = nsize*n2 ! block
        ims_tmp3 = ims_tmp2  ! stride = block to put things together
        select case (ims_trp_type_k)
        case (TLAB_MPI_TRP_DOUBLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL8, rtype, ims_err)
        case (TLAB_MPI_TRP_SINGLE)
            call MPI_TYPE_VECTOR(ims_tmp1, ims_tmp2, ims_tmp3, MPI_REAL4, rtype, ims_err)
        end select
        call MPI_TYPE_COMMIT(rtype, ims_err)

        call MPI_TYPE_SIZE(stype, ims_ss, ims_err)
        call MPI_TYPE_SIZE(rtype, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            print *, 'Message   : ', 1, ' size is wrong' ! i
            print *, 'Send size : ', ims_ss
            print *, 'Recv size : ', ims_rs
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

#undef nsize
#undef sdisp
#undef rdisp
#undef stype
#undef rtype

        return
    end subroutine TLabMPI_TypeK_Create

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TransposeK_Forward(a, b, id)
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
        integer, intent(in) :: id

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

        if (ims_trp_type_k == TLAB_MPI_TRP_SINGLE .and. wp == dp) then
            call MPI_TYPE_SIZE(ims_ts_k(id), size, ims_err)
            size = size/sizeof(1.0_sp)
            size = size*ims_npro_k
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call Transpose_Kernel_Single(a_wrk, ims_plan_trps_k(:), ims_ds_k(:, id), ims_ts_k(id), &
                                         b_wrk, ims_plan_trpr_k(:), ims_dr_k(:, id), ims_tr_k(id), &
                                         ims_comm_z, ims_sizBlock_k, ims_trp_mode_k)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(a, ims_plan_trps_k(:), ims_ds_k(:, id), ims_ts_k(id), &
                                  b, ims_plan_trpr_k(:), ims_dr_k(:, id), ims_tr_k(id), &
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
    subroutine TLabMPI_TransposeK_Backward(b, a, id)
        real(wp), dimension(*), intent(in) :: b
        real(wp), dimension(*), intent(out) :: a
        integer, intent(in) :: id

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

        if (ims_trp_type_k == TLAB_MPI_TRP_SINGLE .and. wp == dp) then
            call MPI_TYPE_SIZE(ims_ts_k(id), size, ims_err)
            size = size/sizeof(1.0_sp)
            size = size*ims_npro_k
            call c_f_pointer(c_loc(a), b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), a_wrk, shape=[size])
            b_wrk(1:size) = real(b(1:size), sp)
            call Transpose_Kernel_Single(b_wrk, ims_plan_trpr_k(:), ims_dr_k(:, id), ims_tr_k(id), &
                                         a_wrk, ims_plan_trps_k(:), ims_ds_k(:, id), ims_ts_k(id), &
                                         ims_comm_z, ims_sizBlock_k, ims_trp_mode_k)
            a(1:size) = real(a_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(b, ims_plan_trpr_k(:), ims_dr_k(:, id), ims_tr_k(id), &
                                  a, ims_plan_trps_k(:), ims_ds_k(:, id), ims_ts_k(id), &
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
    subroutine TLabMPI_TransposeI_Forward(a, b, id)
        real(wp), dimension(*), intent(in) :: a
        real(wp), dimension(*), intent(out) :: b
        integer, intent(in) :: id

        target b

        ! -----------------------------------------------------------------------
        integer(wi) size

        ! #######################################################################
        if (ims_trp_type_i == TLAB_MPI_TRP_SINGLE .and. wp == dp) then
            call MPI_TYPE_SIZE(ims_ts_i(id), size, ims_err)
            size = size/sizeof(1.0_sp)
            size = size*ims_npro_i
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call Transpose_Kernel_Single(a_wrk, ims_plan_trps_i(:), ims_ds_i(:, id), ims_ts_i(id), &
                                         b_wrk, ims_plan_trpr_i(:), ims_dr_i(:, id), ims_tr_i(id), &
                                         ims_comm_x, ims_sizBlock_i, ims_trp_mode_i)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(a, ims_plan_trps_i(:), ims_ds_i(:, id), ims_ts_i(id), &
                                  b, ims_plan_trpr_i(:), ims_dr_i(:, id), ims_tr_i(id), &
                                  ims_comm_x, ims_sizBlock_i, ims_trp_mode_i)
        end if

        return
    end subroutine TLabMPI_TransposeI_Forward

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TransposeI_Backward(b, a, id)
        real(wp), dimension(*), intent(in) :: b
        real(wp), dimension(*), intent(out) :: a
        integer, intent(in) :: id

        target a

        ! -----------------------------------------------------------------------
        integer(wi) size

        ! #######################################################################
        if (ims_trp_type_i == TLAB_MPI_TRP_SINGLE .and. wp == dp) then
            call MPI_TYPE_SIZE(ims_ts_i(id), size, ims_err)
            size = size/sizeof(1.0_sp)
            size = size*ims_npro_i
            call c_f_pointer(c_loc(a), b_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), a_wrk, shape=[size])
            b_wrk(1:size) = real(b(1:size), sp)
            call Transpose_Kernel_Single(b_wrk, ims_plan_trpr_i(:), ims_dr_i(:, id), ims_tr_i(id), &
                                         a_wrk, ims_plan_trps_i(:), ims_ds_i(:, id), ims_ts_i(id), &
                                         ims_comm_x, ims_sizBlock_i, ims_trp_mode_i)
            a(1:size) = real(a_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call Transpose_Kernel(b, ims_plan_trpr_i(:), ims_dr_i(:, id), ims_tr_i(id), &
                                  a, ims_plan_trps_i(:), ims_ds_i(:, id), ims_ts_i(id), &
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

                call TLabMPI_TAGUPDT

            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, ims_status(:, 1), ims_err)
                end do

                call TLabMPI_TAGUPDT

            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            ! call MPI_ALLTOALL(a, 1, tsend, &
            !                   b, 1, trecv, comm, ims_err)

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

                call TLabMPI_TAGUPDT

            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = msend(m) + 1; ips = ns - 1
                    nr = mrecv(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(a(dsend(ns) + 1), 1, tsend, ips, ims_tag, &
                                      b(drecv(nr) + 1), 1, trecv, ipr, ims_tag, comm, ims_status(:, 1), ims_err)
                end do

                call TLabMPI_TAGUPDT

            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            ! call MPI_ALLTOALL(a, 1, tsend, &
            !                   b, 1, trecv, comm, ims_err)

        end select

        return
    end subroutine Transpose_Kernel_Single

end module TLabMPI_PROCS
