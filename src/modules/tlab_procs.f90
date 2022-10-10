#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

module TLAB_PROCS
    use TLAB_CONSTANTS, only: wp, wi, longi, lfile, efile
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
    TINTEGER :: ierr

    public :: TLAB_START
    public :: TLAB_STOP
    public :: TLAB_WRITE_ASCII
    public :: TLAB_ALLOCATE
    public :: TLAB_ALLOCATE_ARRAY1, TLAB_ALLOCATE_ARRAY2
    public :: TLAB_ALLOCATE_ARRAY1_INT, TLAB_ALLOCATE_ARRAY1_LONG_INT
#ifdef USE_MPI
    public :: TLAB_MPI_PANIC
#endif
contains

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_ALLOCATE(C_FILE_LOC)
        use TLAB_ARRAYS

        character(len=*), intent(in) :: C_FILE_LOC

        ! -------------------------------------------------------------------
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, x, g(1)%inb_grid, g(1)%size, 'x')
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, y, g(2)%inb_grid, g(2)%size, 'y')
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, z, g(3)%inb_grid, g(3)%size, 'z')
        ! -------------------------------------------------------------------
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, q, inb_flow_array, isize_field, 'flow')
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, s, inb_scal_array, isize_field, 'scal')
        ! -------------------------------------------------------------------
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, txc, inb_txc, isize_txc_field, 'txc')
        call TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, wrk3d, isize_wrk3d, 'wrk3d')
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, wrk1d, inb_wrk1d, isize_wrk1d, 'wrk1d')
        call TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, wrk2d, inb_wrk2d, isize_wrk2d, 'wrk2d')

        return
    end subroutine TLAB_ALLOCATE

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY1(C_FILE_LOC, a, i1, s)

        character(len=*), intent(in) :: C_FILE_LOC
        TREAL, dimension(:), allocatable, intent(inout) :: a
        TINTEGER, intent(in) :: i1
        character(len=*), intent(in) :: s

        !#####################################################################
        if (i1 <= 0) return

        write (str, *) i1; line = 'Allocating array '//TRIM(ADJUSTL(s))//' of size '//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//TRIM(ADJUSTL(s))//'.')
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

        write (str, *) i1; line = 'Allocating array '//TRIM(ADJUSTL(s))//' of size '//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//TRIM(ADJUSTL(s))//'.')
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

        write (str, *) i1; line = 'Allocating array '//TRIM(ADJUSTL(s))//' of size '//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//TRIM(ADJUSTL(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY1_LONG_INT

    ! ######################################################################
    ! ######################################################################
    subroutine TLAB_ALLOCATE_ARRAY2(C_FILE_LOC, a, i1, i2, s)

        character(len=*), intent(in) :: C_FILE_LOC
        TREAL, dimension(:, :), allocatable, intent(inout) :: a
        TINTEGER, intent(in) :: i1, i2
        character(len=*), intent(in) :: s

        !#####################################################################
        if (i1 <= 0 .or. i2 <= 0) return

        write (str, *) i1; line = 'Allocating array '//TRIM(ADJUSTL(s))//' of size '//TRIM(ADJUSTL(str))//'x'
        write (str, *) i2; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
        call TLAB_WRITE_ASCII(lfile, line)
        allocate (a(i2, i1), stat=ierr)
        if (ierr /= 0) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//TRIM(ADJUSTL(s))//'.')
            call TLAB_STOP(DNS_ERROR_ALLOC)
        end if

    end subroutine TLAB_ALLOCATE_ARRAY2

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

        ims_time_trans = C_0_R

#endif

        !########################################################################
        ! First output
        call DATE_AND_TIME(clock(1), clock(2))
        line = 'Starting on '//TRIM(ADJUSTL(clock(1) (1:8)))//' at '//TRIM(ADJUSTL(clock(2)))
        call TLAB_WRITE_ASCII(lfile, line)

        line = 'Git-hash   '//GITHASH
        call TLAB_WRITE_ASCII(lfile, line)

        line = 'Git-branch '//GITBRANCH
        call TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
        write (line, *) ims_npro
        line = 'Number of MPI tasks '//TRIM(ADJUSTL(line))
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
        line = 'Number of OMP threads '//TRIM(ADJUSTL(line))
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
            line = 'Error code '//TRIM(ADJUSTL(line))//'.'
            call TLAB_WRITE_ASCII(efile, line)
        end if

        call GETARG(0, line)
        write (line, *) 'Finalizing program '//TRIM(ADJUSTL(line))
        if (error_code == 0) then
            line = TRIM(ADJUSTL(line))//' normally.'
        else
            line = TRIM(ADJUSTL(line))//' abnormally. Check '//TRIM(ADJUSTL(efile))
        end if
        call TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
        ims_time_max = MPI_WTIME()
        write (line, 1000) ims_time_max - ims_time_min
        line = 'Time elapse ....................: '//TRIM(ADJUSTL(line))
        call TLAB_WRITE_ASCII(lfile, line)

#ifdef PROFILE_ON
        write (line, 1000) ims_time_trans
        line = 'Time in array transposition ....: '//TRIM(ADJUST(line))
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
        call TLAB_WRITE_ASCII(efile, 'MPI-ERROR: Source file'//TRIM(ADJUSTL(LOCATION)), .true.)
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
        if (ims_pro == 0 .or. PRESENT(flag_all)) then
#endif

            if (imode_verbosity > 0) then

                open (UNIT=22, FILE=file, STATUS='unknown', POSITION='APPEND')
                if (imode_verbosity == 1) then
                    write (22, '(a)') TRIM(ADJUSTL(lineloc))
                else if (imode_verbosity == 2) then
                    call DATE_AND_TIME(clock(1), clock(2))
                    write (22, '(a)') '['//TRIM(ADJUSTR(clock(2)))//'] '//TRIM(ADJUSTL(lineloc))
                end if
                close (22)

            end if

#ifndef PARALLEL
            if (file == efile) then
                write (*, *) TRIM(ADJUSTL(lineloc))
            end if
#endif

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine TLAB_WRITE_ASCII

end module TLAB_PROCS
