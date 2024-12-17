! debug code for levante where we are experiencing some problems possibly due to stripping of disks.
program vmpi_io_levante
    use MPI
    implicit none

    integer mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE), ims_Err
    integer(KIND=MPI_OFFSET_KIND) mpio_disp
    integer subarray

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: wp = dp               ! working precision
    integer, parameter :: wi = selected_int_kind(9)

    ! using 1536 cores, 12 nodes. Decomposition along X and Z in 48x32 pencils
    integer ims_npro                            ! number of tasks
    ! integer, parameter :: ims_npro_i = 48       ! number of tasks in X
    ! integer, parameter :: ims_npro_k = 32       ! number of tasks in Z
    integer, parameter :: ims_npro_i = 2       ! number of tasks in X
    integer, parameter :: ims_npro_k = 2       ! number of tasks in Z
    integer ims_pro                             ! local task in global communicator
    integer ims_pro_i, ims_pro_k                ! local task in each directional communicator; here used only as offsets

    ! grid points
    integer nx, ny, nz

    ! number of variables
    integer, parameter :: nv = 3
    integer iv

    ! array
    real(sp), allocatable :: a(:, :)
    ! real(dp), allocatable :: a(:, :)

    ! file name
    character(*), parameter :: name = 'test-io.'
    character(len=32) str

    ! ###############################################################
    call MPI_INIT(ims_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ims_npro, ims_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)

    ! initialize grid
    ny = 384
    ! nx = 1536
    ! nz = 1536
    nx = 16
    nz = 24
    nx = nx/ims_npro_i                      ! task-local number of grid points along X
    nz = nz/ims_npro_k                      ! task-local number of grid points along Z
    
    ims_pro_i = mod(ims_pro, ims_npro_i)    ! MPI offset
    ims_pro_k = ims_pro/ims_npro_i          ! MPI offset

    allocate (a(nx*ny*nz, nv))
    a = 0.0_wp

    ! single precission
    subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL4)
    ! double precission
    ! subarray = IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, MPI_REAL8)

    mpio_disp = 0
    mpio_locsize = nx*ny*nz
    do iv = 1, nv
        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

        write (str, *) iv; str = trim(adjustl(name))//trim(adjustl(str))
        call MPI_FILE_OPEN(MPI_COMM_WORLD, str, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, mpio_fh, ims_err)

        ! single precission
        call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
        call MPI_File_write_all(mpio_fh, a(:, iv), mpio_locsize, MPI_REAL4, status, ims_err)
        ! double precission
        ! call MPI_File_set_view(mpio_fh, mpio_disp, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
        ! call MPI_File_write_all(mpio_fh, a(:, :, :, iv), mpio_locsize, MPI_REAL8, status, ims_err)

        call MPI_File_close(mpio_fh, ims_err)
    end do

    call MPI_FINALIZE(ims_err)

    stop

contains
    function IO_CREATE_SUBARRAY_XOZ(nx, ny, nz, mpi_type) result(subarray)
        integer(wi), intent(in) :: nx, ny, nz
        integer, intent(in) :: mpi_type

        integer :: subarray
        integer, parameter :: ndims = 3
        integer(wi) :: sizes(ndims), locsize(ndims), offset(ndims)

        sizes = [nx*ims_npro_i, ny, nz*ims_npro_k]
        locsize = [nx, ny, nz]
        offset = [nx*ims_pro_i, 0, nz*ims_pro_k]

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
                                      MPI_ORDER_FORTRAN, mpi_type, subarray, ims_err)
        call MPI_Type_commit(subarray, ims_err)

    end function IO_CREATE_SUBARRAY_XOZ

end program
