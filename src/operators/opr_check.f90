#include "dns_const.h"
#ifdef USE_MPI

#endif

subroutine OPR_CHECK()
    use TLab_Constants, only: lfile, wp, wi
    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_flow_array, inb_txc
    use FDM, only: g
    use TLab_WorkFlow, only: fourier_on
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Arrays
    use OPR_FOURIER
#ifdef USE_MPI
    use TLab_Time, only: itime
    use mpi_f08
    use TLabMPI_VARS, only: ims_err
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_Transpose
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
    integer(wi) idummy !, id
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
        ! id = TLAB_MPI_TRP_I_PARTIAL

        call SYSTEM_CLOCK(t_srt, PROC_CYCLES, MAX_CYCLES)
        call TLabMPI_TransposeI_Forward(q(1, 1), wrk3d, tmpi_plan_dx)
        call TLabMPI_TransposeI_Backward(wrk3d, q(1, 2), tmpi_plan_dx)
        call SYSTEM_CLOCK(t_end, PROC_CYCLES, MAX_CYCLES)

        idummy = t_end - t_srt
        call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        write (str, 100) real(t_dif,wp)/PROC_CYCLES

        dummy = MAXVAL(ABS(q(1:isize_field, 1) - q(1:isize_field, 2)))
        call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

        write (line, 100) residual
        line = 'Checking MPI transposition for Ox derivatives: Residual ' &
               //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
        call TLab_Write_ASCII(lfile, line)

    end if
#endif

! -------------------------------------------------------------------
! Transposition along OZ
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_npro_k > 1) then
        ! id = TLAB_MPI_TRP_K_PARTIAL

        call SYSTEM_CLOCK(t_srt, PROC_CYCLES, MAX_CYCLES)
        idummy = itime; itime = -1  ! set itime to -1 for this call to trigger interruption
        call TLabMPI_TransposeK_Forward(q(:, 1), wrk3d, tmpi_plan_dz)
        itime = idummy
        call TLabMPI_TransposeK_Backward(wrk3d, q(:, 2), tmpi_plan_dz)
        call SYSTEM_CLOCK(t_end, PROC_CYCLES, MAX_CYCLES)

        idummy = t_end - t_srt
        call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
        write (str, 100) real(t_dif,wp)/PROC_CYCLES

        dummy = MAXVAL(ABS(q(1:isize_field, 1) - q(1:isize_field, 2)))
        call MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

        write (line, 100) residual
        line = 'Checking MPI transposition for Oz derivatives: Residual ' &
               //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
        call TLab_Write_ASCII(lfile, line)

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
        call TLab_Write_ASCII(lfile, line)

    end if

    return

100 format(G_FORMAT_R)

end subroutine OPR_CHECK
