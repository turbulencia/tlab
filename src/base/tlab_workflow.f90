#include "dns_const.h"
#include "dns_error.h"

module TLab_WorkFlow
    use TLab_Constants, only: sp, wp, wi, longi, lfile, efile
    use TLAB_VARS, only: imode_verbosity
#ifdef USE_OPENMP
    use OMP_LIB
#endif
#ifdef USE_MPI
    use MPI
    use TLabMPI_VARS, only: ims_pro, ims_npro, ims_time_max, ims_time_min, ims_time_trans, ims_err
#endif
    implicit none
    private
    save

    character*128 :: line

    public :: TLAB_START
    public :: TLAB_STOP
    public :: TLAB_WRITE_ASCII

contains

    ! ###################################################################
    ! ###################################################################
    subroutine TLAB_START()
        use TLab_OpenMP
        
        character*10 clock(2)

        !#####################################################################
        ! Inititalize MPI parallel mode
#ifdef USE_MPI
#ifdef USE_PSFFT
        call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED, ims_nb_thrsupp_provided, ims_err)
        if (ims_nb_thrsupp_provided < MPI_THREAD_SERIALIZED) then
            call TLAB_WRITE_ASCII(efile, __FILE__//". MPI Library is missing the needed level of thread support for nonblocking communication.")
            call TLAB_WRITE_ASCII(efile, __FILE__//". Try MPI_THREAD_FUNNELED after reading the documentation.")
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
        TLab_OMP_numThreads = omp_get_max_threads()
        write (line, *) TLab_OMP_numThreads
        line = 'Number of OMP threads '//trim(adjustl(line))
        call TLAB_WRITE_ASCII(lfile, line)

#else
        TLab_OMP_numThreads = 1

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

end module TLab_WorkFlow
