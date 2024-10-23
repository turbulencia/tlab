#include "dns_const.h"
#include "dns_error.h"

module TLAB_PROCS
    use TLAB_CONSTANTS, only: sp, wp, wi, longi, lfile, efile
    use TLAB_VARS
#ifdef USE_OPENMP
    use OMP_LIB
#endif
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_pro, ims_npro, ims_time_max, ims_time_min, ims_time_trans, ims_err
#endif
    implicit none

    save
    private

    character*128 :: str, line
    integer :: ierr

    public :: TLAB_START
    public :: TLAB_STOP
    public :: TLAB_WRITE_ASCII
    public :: TLab_Initialize_Memory

    interface  TLAB_ALLOCATE_LOG
        module procedure TLAB_ALLOCATE_LOG_SHORT, TLAB_ALLOCATE_LOG_LONG
    end interface TLAB_ALLOCATE_LOG

#ifdef NO_ASSUMED_RANKS
    interface TLAB_ALLOCATE_ARRAY_SINGLE
        module procedure TLAB_ALLOCATE_ARRAY_SINGLE1, TLAB_ALLOCATE_ARRAY_SINGLE2, TLAB_ALLOCATE_ARRAY_SINGLE3, TLAB_ALLOCATE_ARRAY_SINGLE4
    end interface TLAB_ALLOCATE_ARRAY_SINGLE

    interface TLAB_ALLOCATE_ARRAY_DOUBLE
        module procedure TLAB_ALLOCATE_ARRAY_DOUBLE1, TLAB_ALLOCATE_ARRAY_DOUBLE2, TLAB_ALLOCATE_ARRAY_DOUBLE3, TLAB_ALLOCATE_ARRAY_DOUBLE4
    end interface TLAB_ALLOCATE_ARRAY_DOUBLE

    interface TLAB_ALLOCATE_ARRAY_INT
        module procedure TLAB_ALLOCATE_ARRAY_INT1, TLAB_ALLOCATE_ARRAY_INT2, TLAB_ALLOCATE_ARRAY_INT3, TLAB_ALLOCATE_ARRAY_INT4
    end interface TLAB_ALLOCATE_ARRAY_INT

    interface TLAB_ALLOCATE_ARRAY_LONG_INT
        module procedure TLAB_ALLOCATE_ARRAY_LONG_INT1, TLAB_ALLOCATE_ARRAY_LONG_INT2, TLAB_ALLOCATE_ARRAY_LONG_INT3, TLAB_ALLOCATE_ARRAY_LONG_INT4
    end interface TLAB_ALLOCATE_ARRAY_LONG_INT
#endif
    public :: TLAB_ALLOCATE_ARRAY_SINGLE
    public :: TLAB_ALLOCATE_ARRAY_DOUBLE
    public :: TLAB_ALLOCATE_ARRAY_INT
    public :: TLAB_ALLOCATE_ARRAY_LONG_INT
    public :: TLAB_ALLOCATE_ARRAY_DOUBLE1_LONG
contains

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Initialize_Memory(C_FILE_LOC)
        use TLAB_ARRAYS

        character(len=*), intent(in) :: C_FILE_LOC

! loop counters over the whole domain are integer*4
        if (isize_txc_field > huge(imax)) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Integer model of 4 bytes is not big enough.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
    
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, x, [g(1)%size, g(1)%inb_grid], g(1)%name)
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, y, [g(2)%size, g(2)%inb_grid], g(2)%name)
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, z, [g(3)%size, g(3)%inb_grid], g(3)%name)

        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, q, [isize_field, inb_flow_array], 'flow')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, s, [isize_field, inb_scal_array], 'scal')

        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, txc, [isize_txc_field, inb_txc], 'txc')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, wrk1d, [isize_wrk1d, inb_wrk1d], 'wrk1d')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, wrk2d, [isize_wrk2d, inb_wrk2d], 'wrk2d')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, wrk3d, [isize_wrk3d], 'wrk3d')

        if (any(Dealiasing(:)%type /= DNS_FILTER_NONE)) then
            call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, wrkdea, [isize_field, 2], 'wrk-dealiasing')
        end if

        call TLAB_DEFINE_POINTERS()

        call TLAB_DEFINE_POINTERS_3D()

        call TLAB_DEFINE_POINTERS_C()

        return
    end subroutine TLab_Initialize_Memory

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_DEFINE_POINTERS()
        use TLAB_ARRAYS
        use TLAB_POINTERS

        integer(wi) idummy(2)

        idummy = shape(q)
        if (idummy(2) >= 1) u(1:isize_field) => q(1:isize_field, 1)
        if (idummy(2) >= 2) v(1:isize_field) => q(1:isize_field, 2)
        if (idummy(2) >= 3) w(1:isize_field) => q(1:isize_field, 3)
        ! compressible flows variables
        if (idummy(2) >= 4) e(1:isize_field) => q(1:isize_field, 4)
        if (idummy(2) >= 5) rho(1:isize_field) => q(1:isize_field, 5)
        if (idummy(2) >= 6) p(1:isize_field) => q(1:isize_field, 6)
        if (idummy(2) >= 7) T(1:isize_field) => q(1:isize_field, 7)
        if (idummy(2) >= 8) vis(1:isize_field) => q(1:isize_field, 8)

        idummy = shape(txc)
        if (idummy(2) >= 1) tmp1(1:isize_field) => txc(1:isize_field, 1)
        if (idummy(2) >= 2) tmp2(1:isize_field) => txc(1:isize_field, 2)
        if (idummy(2) >= 3) tmp3(1:isize_field) => txc(1:isize_field, 3)
        if (idummy(2) >= 4) tmp4(1:isize_field) => txc(1:isize_field, 4)
        if (idummy(2) >= 5) tmp5(1:isize_field) => txc(1:isize_field, 5)
        if (idummy(2) >= 6) tmp6(1:isize_field) => txc(1:isize_field, 6)
        if (idummy(2) >= 7) tmp7(1:isize_field) => txc(1:isize_field, 7)
        if (idummy(2) >= 8) tmp8(1:isize_field) => txc(1:isize_field, 8)
        if (idummy(2) >= 9) tmp9(1:isize_field) => txc(1:isize_field, 9)

        return
    end subroutine TLAB_DEFINE_POINTERS

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_DEFINE_POINTERS_3D()
        use TLAB_ARRAYS
        use TLAB_POINTERS_3D

        integer(wi) idummy(2)

        idummy = shape(q)
        if (idummy(2) >= 1) u(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 1)
        if (idummy(2) >= 2) v(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 2)
        if (idummy(2) >= 3) w(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 3)
        ! compressible flows variables
        if (idummy(2) >= 4) e(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 4)
        if (idummy(2) >= 5) rho(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 5)
        if (idummy(2) >= 6) p(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 6)
        if (idummy(2) >= 7) T(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 7)
        if (idummy(2) >= 8) vis(1:imax, 1:jmax, 1:kmax) => q(1:isize_field, 8)

        if (allocated(q)) p_q(1:imax, 1:jmax, 1:kmax, 1:inb_flow_array) => q(1:isize_field*inb_flow_array, 1)
        if (allocated(s)) p_s(1:imax, 1:jmax, 1:kmax, 1:inb_scal_array) => s(1:isize_field*inb_scal_array, 1)
        if (allocated(wrk3d)) p_wrk3d(1:imax, 1:jmax, 1:kmax) => wrk3d(1:isize_field)
        if (allocated(wrk2d)) p_wrk2d(1:imax, 1:kmax, 1:inb_wrk2d) => wrk2d(1:imax*kmax*inb_wrk2d, 1)    ! this is the most common wrk2d dimensions
        if (allocated(wrk1d)) p_wrk1d(1:jmax, 1:inb_wrk1d) => wrk1d(1:jmax*inb_wrk1d, 1)                 ! this is the most common wrk1d dimensions

        idummy = shape(txc)
        if (idummy(2) >= 1) tmp1(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 1)
        if (idummy(2) >= 2) tmp2(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 2)
        if (idummy(2) >= 3) tmp3(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 3)
        if (idummy(2) >= 4) tmp4(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 4)
        if (idummy(2) >= 5) tmp5(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 5)
        if (idummy(2) >= 6) tmp6(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 6)
        if (idummy(2) >= 7) tmp7(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 7)
        if (idummy(2) >= 8) tmp8(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 8)
        if (idummy(2) >= 9) tmp9(1:imax, 1:jmax, 1:kmax) => txc(1:isize_field, 9)

        return
    end subroutine TLAB_DEFINE_POINTERS_3D

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_DEFINE_POINTERS_C()
        use TLAB_ARRAYS
        use TLAB_POINTERS_C
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

        if (allocated(wrk1d)) call c_f_pointer(c_loc(wrk1d), c_wrk1d, shape=[jmax, inb_wrk1d/2])
        if (allocated(wrk3d)) call c_f_pointer(c_loc(wrk3d), c_wrk3d, shape=[isize_txc_dimz/2, kmax])

        return
    end subroutine TLAB_DEFINE_POINTERS_C

    ! ######################################################################
    ! ######################################################################
#ifndef NO_ASSUMED_RANKS
    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, a, dims, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(wp), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        ! if (any(dims <= 0)) return; better allocate to zero than not allocate; error in supermuc
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        rank (4)
            allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        rank default
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Rank too large while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE

! ### DOUBLE ALLOCATION ROUTINES FOR LARGE 1D ARRAYS
    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE1_LONG(C_FILE_LOC, a, dims, s)
      character(len=*), intent(in) :: C_FILE_LOC,s
      real(8), allocatable, intent(inout) :: a(:)
      integer(8), intent(in) :: dims(1)
      integer id
      !#####################################################################
      call TLAB_ALLOCATE_LOG_LONG(lfile,dims,s)
      allocate (a(dims(1)), stat=ierr)
      call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE1_LONG

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY_SINGLE(C_FILE_LOC, a, dims, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(sp), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLAB_ALLOCATE_ARRAY_SINGLE

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC
        integer(wi), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLAB_ALLOCATE_ARRAY_INT

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY_LONG_INT(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC
        integer(longi), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)

    end subroutine TLAB_ALLOCATE_ARRAY_LONG_INT
#endif
    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_START()

        character*10 clock(2)

        !#####################################################################
        ! Inititalize MPI parallel mode
#ifdef USE_MPI
#ifdef USE_PSFFT
        call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, ims_nb_thrsupp_provided, ims_err)
        if (ims_nb_thrsupp_provided < MPI_THREAD_SERIALIZED) then
            call TLAB_WRITE_ASCII(efile, "DNS_START. MPI Library is missing the needed level of thread support for nonblocking communication.")
            call TLAB_WRITE_ASCII(efile, "DNS_START. Try MPI_THREAD_FUNNELED after reading the documentation.")
            call TLAB_STOP(DNS_ERROR_PSFFT)
        end if
#else
        call MPI_INIT(ims_err)
#endif

        call MPI_COMM_SIZE(MPI_COMM_WORLD, ims_npro, ims_err)
        call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)

        ims_time_min = MPI_WTIME()

        ims_time_trans = 0.0_wp

#endif

        !########################################################################
        ! First output
        call date_and_time(clock(1), clock(2))
        line = 'Starting on '//trim(adjustl(clock(1) (1:8)))//' at '//trim(adjustl(clock(2)))
        call TLAB_WRITE_ASCII(lfile, line)

        line = 'Git-hash '//GITHASH
        call TLAB_WRITE_ASCII(lfile, line)

        line = 'Git-branch '//GITBRANCH
        call TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
        write (line, *) ims_npro
        line = 'Number of MPI tasks '//trim(adjustl(line))
        call TLAB_WRITE_ASCII(lfile, line)

        if (ims_npro == 0) then
            call TLAB_WRITE_ASCII(efile, 'DNS_START. Number of processors is zero.')
            call TLAB_STOP(DNS_ERROR_MINPROC)
        end if
#endif

        !#####################################################################
        ! Inititalize OpenMP mode
#ifdef USE_OPENMP
        dns_omp_numThreads = omp_get_max_threads()
        write (line, *) dns_omp_numThreads
        line = 'Number of OMP threads '//trim(adjustl(line))
        call TLAB_WRITE_ASCII(lfile, line)

#else
        dns_omp_numThreads = 1

#endif

        return
    end subroutine TLAB_START

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_STOP(error_code)

        integer, intent(in) :: error_code

        ! ###################################################################
! #ifdef USE_FFTW
!         if (fourier_on) then
!             call dfftw_destroy_plan(fft_plan_fx)
!             call dfftw_destroy_plan(fft_plan_bx)
!             if (g(3)%size > 1) then
!                 call dfftw_destroy_plan(fft_plan_fz)
!                 call dfftw_destroy_plan(fft_plan_bz)
!             end if
!         end if
! #endif

        ! ###################################################################
        if (error_code /= 0) then
            write (line, *) error_code
            line = 'Error code '//trim(adjustl(line))//'.'
            call TLAB_WRITE_ASCII(efile, line)
        end if

        call GETARG(0, line)
        write (line, *) 'Finalizing program '//trim(adjustl(line))
        if (error_code == 0) then
            line = trim(adjustl(line))//' normally.'
        else
            line = trim(adjustl(line))//' abnormally. Check '//trim(adjustl(efile))
        end if
        call TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
        ims_time_max = MPI_WTIME()
        write (line, 1000) ims_time_max - ims_time_min
        line = 'Time elapse ....................: '//trim(adjustl(line))
        call TLAB_WRITE_ASCII(lfile, line)

#ifdef PROFILE_ON
        write (line, 1000) ims_time_trans
        line = 'Time in array transposition ....: '//trim(ADJUST(line))
        call TLAB_WRITE_ASCII(lfile, line)
#endif

1000    format(G_FORMAT_R)

#endif

        call TLAB_WRITE_ASCII(lfile, '########################################')

#ifdef USE_MPI
        if (error_code == 0) then
            call MPI_FINALIZE(ims_err)
        else
            call MPI_Abort(MPI_COMM_WORLD, error_code, ims_err)
        end if
#endif
        return
    end subroutine TLAB_STOP

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_WRITE_ASCII(file, lineloc, flag_all)

        character*(*), intent(in) :: file, lineloc
        logical, intent(in), optional :: flag_all

        ! -----------------------------------------------------------------------
        character*10 clock(2)

        ! #######################################################################
#ifdef USE_MPI
        if (ims_pro == 0 .or. present(flag_all)) then
#endif

            if (imode_verbosity > 0) then

                open (UNIT=22, FILE=file, STATUS='unknown', POSITION='APPEND')
                if (imode_verbosity == 1) then
                    write (22, '(a)') trim(adjustl(lineloc))
                else if (imode_verbosity == 2) then
                    call date_and_time(clock(1), clock(2))
                    write (22, '(a)') '['//trim(adjustr(clock(2)))//'] '//trim(adjustl(lineloc))
                end if
                close (22)

            end if

#ifndef PARALLEL
            if (file == efile) then
                write (*, *) trim(adjustl(lineloc))
            end if
#endif

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine TLAB_WRITE_ASCII

#ifdef NO_ASSUMED_RANKS
    ! #######################################################
    ! ### INSTANCES FOR INTERFACES TO ALLOCATION ROUTINES
    ! ###
    ! ### SINGLE ALLOCATION ROUTINES
    subroutine TLAB_ALLOCATE_ARRAY_SINGLE1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(4), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_SINGLE1

    subroutine TLAB_ALLOCATE_ARRAY_SINGLE2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(4), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_SINGLE2

    subroutine TLAB_ALLOCATE_ARRAY_SINGLE3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(4), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_SINGLE3

    subroutine TLAB_ALLOCATE_ARRAY_SINGLE4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(4), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_SINGLE4

    ! ### DOUBLE ALLOCATION ROUTINES
    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(8), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE1

    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(8), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE2

    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(8), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE3

    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        real(8), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE4

    ! # INTEGER ALLOCATION ROUTINES
    subroutine TLAB_ALLOCATE_ARRAY_INT1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_INT1

    subroutine TLAB_ALLOCATE_ARRAY_INT2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        !   if (any(dims <= 0)) return; better allocate to zero than not allocate; error in supermuc
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_INT2

    subroutine TLAB_ALLOCATE_ARRAY_INT3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_INT3

    subroutine TLAB_ALLOCATE_ARRAY_INT4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(wi), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_INT4

    ! # LONG INTEGER ALLOCATION ROUTINES
    subroutine TLAB_ALLOCATE_ARRAY_LONG_INT1(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:)
        integer(wi), intent(in) :: dims(:)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_LONG_INT1

    subroutine TLAB_ALLOCATE_ARRAY_LONG_INT2(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:, :)
        integer(wi), intent(in) :: dims(2)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_LONG_INT2

    subroutine TLAB_ALLOCATE_ARRAY_LONG_INT3(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:, :, :)
        integer(wi), intent(in) :: dims(3)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_LONG_INT3

    subroutine TLAB_ALLOCATE_ARRAY_LONG_INT4(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC, s
        integer(longi), allocatable, intent(inout) :: a(:, :, :, :)
        integer(wi), intent(in) :: dims(4)
        integer id
        !#####################################################################
        call TLAB_ALLOCATE_LOG(lfile, dims, s)
        allocate (a(dims(1), dims(2), dims(3), dims(4)), stat=ierr)
        call TLAB_ALLOCATE_ERR(C_FILE_LOC, efile, s)
    end subroutine TLAB_ALLOCATE_ARRAY_LONG_INT4

#endif

    subroutine TLAB_ALLOCATE_LOG_SHORT(log_file, dims, s)
        integer(wi), intent(IN) :: dims(:)
        character(len=*), intent(IN) :: log_file, s
        integer id
        !#####################################################################

        if (any(dims < 0)) then
            ierr = DNS_ERROR_ALLOC
            call TLAB_ALLOCATE_ERR('TLAB_ALLOCATE_LOG', efile, s)
        end if

        if (any(dims == 0)) return      ! do not print out lines when allocation a zero-space array

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//' x '//trim(adjustl(str))
        end do
        call TLAB_WRITE_ASCII(log_file, line)
    end subroutine TLAB_ALLOCATE_LOG_SHORT

    subroutine TLAB_ALLOCATE_LOG_LONG(log_file, dims, s)
        integer(8), intent(IN) :: dims(:)
        character(len=*), intent(IN) :: log_file, s
        integer id
        !#####################################################################

        if (any(dims < 0)) then
            ierr = DNS_ERROR_ALLOC
            call TLAB_ALLOCATE_ERR('TLAB_ALLOCATE_LOG', efile, s)
        end if

        if (any(dims == 0)) return      ! do not print out lines when allocation a zero-space array

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//' x '//trim(adjustl(str))
        end do
        call TLAB_WRITE_ASCII(log_file, line)
    end subroutine TLAB_ALLOCATE_LOG_LONG

    subroutine TLAB_ALLOCATE_ERR(C_FILE_LOC, log_file, s)
        character(len=*) :: C_FILE_LOC, log_file, s

        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(log_file, C_FILE_LOC//'. Error while allocating memory space for'//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ERR

end module TLAB_PROCS
