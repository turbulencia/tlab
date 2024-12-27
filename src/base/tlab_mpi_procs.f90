#include "dns_const.h"
#include "dns_const_mpi.h"
#include "dns_error.h"

! To be changed into TLabMPI (includes TLabMPI_VARS) and TLabMPI_Transpose

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

    ! Global variables and procedures
    public :: TLabMPI_Initialize
    public :: TLabMPI_TypeK_Create, TLabMPI_TypeI_Create
    public :: TLabMPI_TransposeK_Forward, TLabMPI_TransposeK_Backward
    public :: TLabMPI_TransposeI_Forward, TLabMPI_TransposeI_Backward
    public :: TLabMPI_PANIC
    public :: TLabMPI_WRITE_PE0_SINGLE

    integer, public :: TLAB_MPI_REAL_TYPE

    ! Local variables and procedures; mainly for the transposition--should be a separate module TlabMPI_Transpose

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

    integer :: ims_tag
    integer, allocatable :: ims_status(:, :)
    integer, allocatable :: ims_request(:)

    real(wp), allocatable, target :: wrk_mpi(:)      ! 3D work array for MPI
    real(sp), pointer :: a_wrk(:) => null(), b_wrk(:) => null()

#ifdef USE_PSFFT
    integer :: ims_nb_thrsupp_provided
    integer, dimension(2) :: ims_nb_proc_grid
    integer, dimension(3) :: ims_nb_msize
    integer, dimension(3) :: ims_nb_xsrt, ims_nb_xend, ims_nb_xsiz
    integer, dimension(3) :: ims_nb_ysrt, ims_nb_yend, ims_nb_ysiz
    integer, dimension(3) :: ims_nb_zsrt, ims_nb_zend, ims_nb_zsiz
#endif

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
    end subroutine TLabMPI_Initialize

    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_TAGUPDT

        ims_tag = ims_tag + 1
        if (ims_tag > 32000) then
            call TLabMPI_TAGRESET
        end if

        return
    end subroutine TLabMPI_TAGUPDT

    !########################################################################
    !########################################################################
    subroutine TLabMPI_TAGRESET

        ims_tag = 0

        return
    end subroutine TLabMPI_TAGRESET

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

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_PANIC(location, mpi_error_code)
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

    end subroutine TLabMPI_PANIC

    !########################################################################
    !########################################################################
    subroutine TLabMPI_WRITE_PE0_SINGLE(iunit, nx, ny, nz, subdomain, u, tmp1, tmp2)

        integer(wi) iunit, nx, ny, nz, subdomain(6)
        real(wp), dimension(nx*ny*nz), target :: u, tmp1
        real(wp), dimension(nx*ims_npro_i, *), target :: tmp2

        ! -------------------------------------------------------------------
        integer(wi) nx_total, ny_total, nz_total
        integer(wi) nx_min, nx_max, ny_min, ny_max, nz_min, nz_max
        integer(wi) nyz

        integer(wi) ip_i, ip_k, joffset_loc, koffset_loc, id
        integer(wi) i, jk, j_loc, k_loc
        integer mpio_size, mpio_ip
        integer status(MPI_STATUS_SIZE)

        real(wp), dimension(:), pointer :: p_org

        ! ###################################################################
        nx_total = nx*ims_npro_i
        ny_total = ny
        nz_total = nz*ims_npro_k

        nx_min = subdomain(1); nx_max = subdomain(2)
        ny_min = subdomain(3); ny_max = subdomain(4)
        nz_min = subdomain(5); nz_max = subdomain(6)

        koffset_loc = 0
        joffset_loc = 0

        id = TLAB_MPI_TRP_I_PARTIAL

        ! -------------------------------------------------------------------
        ! Transposing along Ox
        ! -------------------------------------------------------------------
        if (ims_npro_i > 1) then
            call TLabMPI_TransposeI_Forward(u, tmp1, id)
            p_org => tmp1
            nyz = ims_size_i(id)
        else
            p_org => u
            nyz = ny*nz
        end if
        mpio_size = nyz*nx_total

        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

        ! -------------------------------------------------------------------
        ! Passing all data through PE#0
        ! -------------------------------------------------------------------
        if (ims_pro == 0) then

            do ip_k = 1, ims_npro_k
                koffset_loc = nz*(ip_k - 1)

                do ip_i = 1, ims_npro_i
                    joffset_loc = nyz*(ip_i - 1) ! Remember that data is Ox-transposed

                    mpio_ip = ims_npro_i*(ip_k - 1) + ip_i - 1
                    if (mpio_ip == 0) then
                        tmp2(1:mpio_size, 1) = p_org(1:mpio_size)
                    else
                        call MPI_RECV(tmp2, mpio_size, MPI_REAL8, mpio_ip, ims_tag, MPI_COMM_WORLD, status, ims_err)
                    end if

                    do jk = 1, nyz
                        j_loc = mod((jk - 1 + joffset_loc), ny_total) + 1
                        k_loc = ((jk - 1 + joffset_loc)/ny_total) + 1 + koffset_loc

                        if ((j_loc >= ny_min) .and. (j_loc <= ny_max) .and. &
                            (k_loc >= nz_min) .and. (k_loc <= nz_max)) then
                            write (iunit) (SNGL(tmp2(i, jk)), i=nx_min, nx_max)
                        end if

                    end do

                end do
            end do

        else
            call MPI_SEND(p_org, mpio_size, MPI_REAL8, 0, ims_tag, MPI_COMM_WORLD, ims_err)
        end if

        return
    end subroutine TLabMPI_WRITE_PE0_SINGLE

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

            call TLabMPI_TAGUPDT

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

            call TLabMPI_TAGUPDT

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

            call TLabMPI_TAGUPDT

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
            .and. ims_nb_xsiz(2)*ims_nb_xsiz(3) == ims_size_i(TLAB_MPI_TRP_I_PARTIAL)) then
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
            .and. ims_nb_zsiz(1)*ims_nb_zsiz(2) == ims_size_k(TLAB_MPI_TRP_K_PARTIAL)) then
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
