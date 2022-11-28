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
    use TLAB_MPI_VARS
#endif
    implicit none
    save
    private

    character*128 :: str, line
    integer :: ierr

    public :: TLAB_START
    public :: TLAB_STOP
    public :: TLAB_WRITE_ASCII
    public :: TLAB_ALLOCATE
    public :: TLAB_ALLOCATE_ARRAY_SINGLE
    public :: TLAB_ALLOCATE_ARRAY_DOUBLE
    public :: TLAB_ALLOCATE_ARRAY_INT
    public :: TLAB_ALLOCATE_ARRAY_LONG_INT
    ! to be removed
    public :: TLAB_ALLOCATE_ARRAY1, TLAB_ALLOCATE_ARRAY2
    public :: TLAB_ALLOCATE_ARRAY1_INT!, TLAB_ALLOCATE_ARRAY1_LONG_INT
#ifdef USE_MPI
    public :: TLAB_MPI_PANIC
#endif
contains

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_ALLOCATE(C_FILE_LOC)
        use TLAB_ARRAYS

        character(len=*), intent(in) :: C_FILE_LOC

        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, x, [g(1)%size, g(1)%inb_grid], g(1)%name)
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, y, [g(2)%size, g(2)%inb_grid], g(2)%name)
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, z, [g(3)%size, g(3)%inb_grid], g(3)%name)

        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, q, [isize_field, inb_flow_array], 'flow')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, s, [isize_field, inb_scal_array], 'scal')

        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, txc, [isize_txc_field, inb_txc], 'txc')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, wrk1d, [isize_wrk1d, inb_wrk1d], 'wrk1d')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, wrk2d, [isize_wrk2d, inb_wrk2d], 'wrk2d')
        call TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, wrk3d, [isize_wrk3d], 'wrk3d')

        return
    end subroutine TLAB_ALLOCATE

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, a, i1, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(wp), dimension(:), allocatable, intent(inout) :: a
        integer(wi), intent(in) :: i1
        character(len=*), intent(in) :: s

        !#####################################################################
        if (i1 <= 0) return

        write (str, *) i1; line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY1

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC, a, i1, s)

        character(len=*), intent(in) :: C_FILE_LOC
        integer(wi), dimension(:), allocatable, intent(inout) :: a
        integer(wi), intent(in) :: i1
        character(len=*), intent(in) :: s

        !#####################################################################
        if (i1 <= 0) return

        write (str, *) i1; line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY1_INT

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY1_LONG_INT(C_FILE_LOC, a, i1, s)

        character(len=*), intent(in) :: C_FILE_LOC
        integer(longi), dimension(:), allocatable, intent(inout) :: a
        integer(wi), intent(in) :: i1
        character(len=*), intent(in) :: s

        !#####################################################################
        if (i1 <= 0) return

        write (str, *) i1; line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY1_LONG_INT

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, a, i1, i2, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(wp), dimension(:, :), allocatable, intent(inout) :: a
        integer(wi), intent(in) :: i1, i2
        character(len=*), intent(in) :: s

        !#####################################################################
        if (i1 <= 0 .or. i2 <= 0) return

        write (str, *) i2; line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))//'x'
        write (str, *) i1; line = trim(adjustl(line))//trim(adjustl(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i2, i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY2

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY_DOUBLE(C_FILE_LOC, a, dims, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(wp), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        integer id

        !#####################################################################
        if (any(dims <= 0)) return

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//'x'//trim(adjustl(str))
        end do
        call TLAB_WRITE_ASCII(lfile, line)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY_DOUBLE

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY_SINGLE(C_FILE_LOC, a, dims, s)

        character(len=*), intent(in) :: C_FILE_LOC
        real(sp), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        integer id

        !#####################################################################
        if (any(dims <= 0)) return

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//'x'//trim(adjustl(str))
        end do
        call TLAB_WRITE_ASCII(lfile, line)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY_SINGLE

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY_INT(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC
        integer(wi), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        integer id

        !#####################################################################
        if (any(dims <= 0)) return

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//'x'//trim(adjustl(str))
        end do
        call TLAB_WRITE_ASCII(lfile, line)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY_INT

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY_LONG_INT(C_FILE_LOC, a, dims, s)
        character(len=*), intent(in) :: C_FILE_LOC
        integer(longi), allocatable, intent(inout) :: a(..)
        integer(wi), intent(in) :: dims(:)
        character(len=*), intent(in) :: s

        integer id

        !#####################################################################
        if (any(dims <= 0)) return

        write (str, *) dims(1); line = 'Allocating array '//trim(adjustl(s))//' of size '//trim(adjustl(str))
        do id = 2, size(dims)
            write (str, *) dims(id); line = trim(adjustl(line))//'x'//trim(adjustl(str))
        end do
        call TLAB_WRITE_ASCII(lfile, line)
        select rank (a)
        rank (1)
            allocate (a(dims(1)), stat=ierr)
        rank (2)
            allocate (a(dims(1), dims(2)), stat=ierr)
        rank (3)
            allocate (a(dims(1), dims(2), dims(3)), stat=ierr)
        end select
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//trim(adjustl(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY_LONG_INT

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

        line = 'Git-hash   '//GITHASH
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
#ifdef USE_FFTW
        if (ifourier == 1) then
            call dfftw_destroy_plan(fft_plan_fx)
            call dfftw_destroy_plan(fft_plan_bx)
            if (g(3)%size > 1) then
                call dfftw_destroy_plan(fft_plan_fz)
                call dfftw_destroy_plan(fft_plan_bz)
            end if
        end if
        if (ivfilter == 1) then
            call dfftw_destroy_plan(fft_plan_fy1d)
            call dfftw_destroy_plan(fft_plan_by1d)
        end if
#endif

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
        if (ims_err == 0) then
            call MPI_FINALIZE(ims_err)
        else
            call MPI_Abort(MPI_COMM_WORLD, ims_err, ims_err)
        end if
#endif
        return
    end subroutine TLAB_STOP

#ifdef USE_MPI
    subroutine TLAB_MPI_PANIC(location, mpi_error_code)

        implicit none

        character(len=*), intent(in) :: location
        integer, intent(in) :: mpi_error_code

        !##############################
        character error_string*1024, line*512
        integer error_local, error_len

        call MPI_Error_String(mpi_error_code, error_string, error_len, error_local)
        call TLAB_WRITE_ASCII(efile, 'MPI-ERROR: Source file'//trim(adjustl(LOCATION)), .true.)
        call TLAB_WRITE_ASCII(efile, error_string, .true.)

        call TLAB_STOP(mpi_error_code)
        ! Not supposed to return from this subroutine

    end subroutine TLAB_MPI_PANIC
#endif

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

end module TLAB_PROCS
