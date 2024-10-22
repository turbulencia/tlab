#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

subroutine OPR_CHECK()
    use TLAB_CONSTANTS, only: lfile, wp, wi
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, inb_flow_array, inb_txc
    use TLAB_VARS, only: g
    use TLAB_VARS, only: fourier_on
    use TLab_WorkFlow
    use TLAB_ARRAYS
    use OPR_FOURIER
#ifdef USE_MPI
    use TLAB_VARS, only: itime
    use MPI
    use TLabMPI_VARS, only: ims_err
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_VARS, only: ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
    use TLabMPI_VARS, only: ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
    use TLabMPI_VARS, only: ims_sizBlock_i, ims_sizBlock_k
    use TLabMPI_PROCS
#endif

    implicit none

! -------------------------------------------------------------------
    real(wp) residual
    integer(wi) t_srt, t_end, t_dif, PROC_CYCLES, MAX_CYCLES
    real(wp) norm
    character*64 str
    character*256 line

#ifdef USE_MPI
    real(wp) dummy
    integer(wi) idummy, id
#endif

! ###################################################################
    if ( inb_flow_array < 3 .or. inb_txc < 4 ) return ! not enough memory

! Create random array
    call RANDOM_NUMBER(q(1:isize_field, 1))

! -------------------------------------------------------------------
! Transposition along OX
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_npro_i > 1) then
        id = TLabMPI_I_PARTIAL

        call SYSTEM_CLOCK(t_srt, PROC_CYCLES, MAX_CYCLES)
        call TLabMPI_TRPF_I(q(1, 1), wrk3d, ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        call TLabMPI_TRPB_I(wrk3d, q(1, 2), ims_ds_i(1, id), ims_dr_i(1, id), ims_ts_i(1, id), ims_tr_i(1, id))
        call SYSTEM_CLOCK(t_end, PROC_CYCLES, MAX_CYCLES)

        idummy = t_end - t_srt
        call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        write (str, 100) real(t_dif,wp)/PROC_CYCLES

        dummy = MAXVAL(ABS(q(1:isize_field, 1) - q(1:isize_field, 2)))
        call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

        write (line, 100) residual
        line = 'Checking MPI transposition for Ox derivatives: Residual ' &
               //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
        call TLAB_WRITE_ASCII(lfile, line)

        if (ims_npro_i > ims_sizBlock_i) then
            line = ''
            write (line, *) ims_sizBlock_i
            line = '   using blocking of '//TRIM(ADJUSTL(line))//' in  TLabMPI_TRP<F,B>_I'
            call TLAB_WRITE_ASCII(lfile, line)
        end if
    end if
#endif

! -------------------------------------------------------------------
! Transposition along OZ
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_npro_k > 1) then
        id = TLabMPI_K_PARTIAL

        call SYSTEM_CLOCK(t_srt, PROC_CYCLES, MAX_CYCLES)
        idummy = itime; itime = -1  ! set itime to -1 for this call to trigger interruption
        call TLabMPI_TRPF_K(q(1, 1), wrk3d, ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
        itime = idummy
        call TLabMPI_TRPB_K(wrk3d, q(1, 2), ims_ds_k(1, id), ims_dr_k(1, id), ims_ts_k(1, id), ims_tr_k(1, id))
        call SYSTEM_CLOCK(t_end, PROC_CYCLES, MAX_CYCLES)

        idummy = t_end - t_srt
        call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        write (str, 100) real(t_dif,wp)/PROC_CYCLES

        dummy = MAXVAL(ABS(q(1:isize_field, 1) - q(1:isize_field, 2)))
        call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

        write (line, 100) residual
        line = 'Checking MPI transposition for Oz derivatives: Residual ' &
               //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
        call TLAB_WRITE_ASCII(lfile, line)

        if (ims_npro_k > ims_sizBlock_k) then
            line = ''
            write (line, *) ims_sizBlock_k
            line = '   using blocking of '//TRIM(ADJUSTL(line))//' in  TLabMPI_TRP<F,B>_K'
            call TLAB_WRITE_ASCII(lfile, line)
        end if
    end if
#endif

! -------------------------------------------------------------------
! Poisson FFT
! -------------------------------------------------------------------
    if (fourier_on) then

        wrk2d(:, 1:2) = 0.0_wp
        txc(1:isize_field, 3) = q(1:isize_field, 1)

!     fft_reordering = .true.
        call SYSTEM_CLOCK(t_srt, PROC_CYCLES, MAX_CYCLES)
        call OPR_FOURIER_F(2, imax, jmax, kmax, txc(1, 3), txc(1, 1), txc(1, 2))
        call OPR_FOURIER_B(2, imax, jmax, kmax, txc(1, 1), txc(1, 2))
        call SYSTEM_CLOCK(t_end, PROC_CYCLES, MAX_CYCLES)
!     fft_reordering = .false.

        q(1:isize_field, 2) = txc(1:isize_field, 2)

#ifdef USE_MPI
        idummy = t_end - t_srt
        call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
#else
        t_dif = t_end - t_srt
#endif
        write (str, 100) real(t_dif,wp)/PROC_CYCLES

        norm = 1.0_wp/real(g(1)%size*g(3)%size,wp)

#ifdef USE_MPI
        dummy = MAXVAL(ABS(norm*q(1:isize_field, 2) - q(1:isize_field, 1)))
        call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
#else
        residual = MAXVAL(ABS(norm*q(1:isize_field, 2) - q(1:isize_field, 1)))
#endif

        write (line, 100) residual
        line = 'Checking FFT routines: Residual ' &
               //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
        call TLAB_WRITE_ASCII(lfile, line)

    end if

    return

100 format(G_FORMAT_R)

end subroutine OPR_CHECK
