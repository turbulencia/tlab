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

        call TLAB_DEFINE_POINTERS()

        call TLAB_DEFINE_POINTERS_3D()

        call TLAB_DEFINE_POINTERS_C()

        return
    end subroutine TLAB_ALLOCATE

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
        if (ims_err == 0) then
            call MPI_FINALIZE(ims_err)
        else
            call MPI_Abort(MPI_COMM_WORLD, ims_err, ims_err)
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

end module TLAB_PROCS
