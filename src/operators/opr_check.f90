#include "types.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

SUBROUTINE OPR_CHECK(nx,ny,nz, a, txc, wrk2d,wrk3d)

  USE TLAB_VARS, ONLY : isize_field,isize_txc_field, isize_wrk2d
  USE TLAB_VARS, ONLY : itime
  USE TLAB_VARS, ONLY : g
  USE TLAB_VARS, ONLY : ifourier !, fft_reordering
  USE TLAB_CONSTANTS, ONLY : lfile
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_err
  USE TLAB_MPI_VARS, ONLY : ims_npro_i, ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_VARS, ONLY : ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_VARS, ONLY : ims_sizBlock_i, ims_sizBlock_k 
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx,ny,nz
  TREAL, DIMENSION(isize_field,2)     :: a
  TREAL, DIMENSION(isize_txc_field,3) :: txc
  TREAL, DIMENSION(isize_txc_field)   :: wrk3d
  TREAL, DIMENSION(isize_wrk2d,2)     :: wrk2d

! -------------------------------------------------------------------
  TREAL residual
  TINTEGER t_srt,t_end,t_dif, PROC_CYCLES, MAX_CYCLES
  TREAL norm
  CHARACTER*64 str
  CHARACTER*256 line

#ifdef USE_MPI
  TREAL dummy
  TINTEGER idummy, id
#endif

! ###################################################################
! Create random array
  CALL RANDOM_NUMBER(a(1:isize_field,1))

! -------------------------------------------------------------------
! Transposition along OX
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     id = TLAB_MPI_I_PARTIAL

     CALL SYSTEM_CLOCK(t_srt,PROC_CYCLES,MAX_CYCLES)
     CALL TLAB_MPI_TRPF_I(a(1,1), wrk3d, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     CALL TLAB_MPI_TRPB_I(wrk3d, a(1,2), ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     CALL SYSTEM_CLOCK(t_end,PROC_CYCLES,MAX_CYCLES)

     idummy = t_end-t_srt
     CALL MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
     WRITE(str,100) M_REAL(t_dif)/PROC_CYCLES

     dummy = MAXVAL(ABS(a(1:isize_field,1)-a(1:isize_field,2)))
     CALL MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

     WRITE(line,100) residual
     line = 'Checking MPI transposition for Ox derivatives: Residual '&
          //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
     CALL TLAB_WRITE_ASCII(lfile,line)

     IF ( ims_npro_i .GT. ims_sizBlock_i ) THEN 
        line=''
        WRITE(line,*) ims_sizBlock_i
        line = '   using blocking of ' // TRIM(ADJUSTL(line)) // ' in  TLAB_MPI_TRP<F,B>_I'
        CALL TLAB_WRITE_ASCII(lfile,line) 
     ENDIF
  ENDIF
#endif

! -------------------------------------------------------------------
! Transposition along OZ
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     id = TLAB_MPI_K_PARTIAL

     CALL SYSTEM_CLOCK(t_srt,PROC_CYCLES,MAX_CYCLES)
     idummy=itime; itime=-1  ! set itime to -1 for this call to trigger interruption
     CALL TLAB_MPI_TRPF_K(a(1,1), wrk3d, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     itime=idummy
     CALL TLAB_MPI_TRPB_K(wrk3d, a(1,2), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     CALL SYSTEM_CLOCK(t_end,PROC_CYCLES,MAX_CYCLES)

     idummy = t_end-t_srt
     CALL MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
     WRITE(str,100) M_REAL(t_dif)/PROC_CYCLES

     dummy = MAXVAL(ABS(a(1:isize_field,1)-a(1:isize_field,2)))
     CALL MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

     WRITE(line,100) residual
     line = 'Checking MPI transposition for Oz derivatives: Residual '&
          //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
     CALL TLAB_WRITE_ASCII(lfile,line)

     IF ( ims_npro_k .GT. ims_sizBlock_k ) THEN 
        line=''
        WRITE(line,*) ims_sizBlock_k
        line = '   using blocking of ' // TRIM(ADJUSTL(line)) // ' in  TLAB_MPI_TRP<F,B>_K'
        CALL TLAB_WRITE_ASCII(lfile,line) 
     ENDIF
  ENDIF
#endif

! -------------------------------------------------------------------
! Poisson FFT
! -------------------------------------------------------------------
  IF ( ifourier .EQ. 1) THEN

     wrk2d(:,1:2) = C_0_R
     txc(1:isize_field,3) = a(1:isize_field,1)

!     fft_reordering = 1
     CALL SYSTEM_CLOCK(t_srt,PROC_CYCLES,MAX_CYCLES)
     CALL OPR_FOURIER_F(i2, nx,ny,nz, txc(1,3),txc(1,1), txc(1,2), wrk2d,wrk3d)
     CALL OPR_FOURIER_B(i2, nx,ny,nz, txc(1,1),txc(1,2), wrk3d)
     CALL SYSTEM_CLOCK(t_end,PROC_CYCLES,MAX_CYCLES)
!     fft_reordering = 0

     a(1:isize_field,2) = txc(1:isize_field,2)

#ifdef USE_MPI
     idummy = t_end-t_srt
     CALL MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
#else
     t_dif  = t_end-t_srt
#endif
     WRITE(str,100) M_REAL(t_dif)/PROC_CYCLES

     norm = C_1_R/M_REAL( g(1)%size *g(3)%size )
!     norm = norm /M_REAL(g(2)%size) ! for large domains we need to do it in two steps !

#ifdef USE_MPI
     dummy = MAXVAL(ABS(norm*a(1:isize_field,2)-a(1:isize_field,1)))
     CALL MPI_REDUCE(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)
#else
     residual = MAXVAL(ABS(norm*a(1:isize_field,2)-a(1:isize_field,1)))
#endif

     WRITE(line,100) residual
     line = 'Checking FFT routines: Residual '&
          //TRIM(ADJUSTL(line))//'. Max. elapsed time '//TRIM(ADJUSTL(str))//' sec.'
     CALL TLAB_WRITE_ASCII(lfile,line)

  ENDIF

  RETURN

100 FORMAT(G_FORMAT_R)

END SUBROUTINE OPR_CHECK
