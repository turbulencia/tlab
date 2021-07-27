#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2012/04/26 - C. Ansorge
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Forward/Backward Fourier transform of the array a extended by two planes
!# in Oy direction (that is, jmax+2).
!#
!# In the case of the pressure solver these planes represent the
!# boundary conditions.
!#
!# The transformed complex field is saved to an array with shape
!# (isize_txc_dimz/2,nz), that is, a stride isize_txc_dimz/2 between z-planes
!# (z-planes need not be contiguous).
!#
!# Each z-plane contains jmax+2 lines of nx/2+1 complex numbers. The first
!# nx/2 Fourier coefficients are contiguous and the element nx/2+1 is the
!# Nyquist frequency; if OX MPI decomposition, this element is just a padding used
!# to homogeneize arrays across PEs with ims_npro_i-PE, which contains then the
!# Nyquist frequency of the corresponding line.
!#
!########################################################################
SUBROUTINE OPR_FOURIER_INITIALIZE(tmp, wrk1d,wrk2d,wrk3d)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL, ONLY : fft_plan_fx_bcs,fft_plan_fx,fft_plan_bx
  USE DNS_GLOBAL, ONLY : fft_plan_fy,fft_plan_by
  USE DNS_GLOBAL, ONLY : fft_plan_fz,fft_plan_bz
  USE DNS_GLOBAL, ONLY : imax,jmax, isize_txc_field
  USE DNS_GLOBAL, ONLY : g
  USE TLAB_PROCS

#ifdef USE_MPI
USE TLAB_MPI_VARS, ONLY : ims_npro_i, ims_npro_k
USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_size_k
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif
#ifdef USE_FFTW
#include "fftw3.f"
#endif

  TREAL, DIMENSION(isize_txc_field), INTENT(INOUT) :: tmp, wrk3d
  TREAL, DIMENSION(g(1)%size+2),     INTENT(INOUT) :: wrk1d, wrk2d

  ! -----------------------------------------------------------------------
  TINTEGER isize_stride, isize_disp, isize_fft_z, isize_fft_y, isize_fft_x

  ! #######################################################################
  ! FFTW library
  ! #######################################################################
#ifdef USE_FFTW

#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
    IF ( ims_size_i(DNS_MPI_I_POISSON1) .NE. ims_size_i(DNS_MPI_I_POISSON2) ) THEN
      CALL TLAB_WRITE_ASCII(efile,'OPR_FOURIER_INITIALIZE. Error in the size in the transposition arrays.')
      CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
    ENDIF
  ENDIF
#endif

  ! -----------------------------------------------------------------------
  ! Oz direction
  ! -----------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
    isize_fft_z = ims_size_k(DNS_MPI_K_POISSON )/2 ! divide by 2 bcs. we work w complex #
  ELSE
#endif
    isize_fft_z = (imax/2+1)*(jmax+2)
#ifdef USE_MPI
  ENDIF
#endif

  isize_stride = isize_fft_z

  IF ( g(3)%size .GT. 1 ) THEN
#ifdef _DEBUG
    CALL dfftw_plan_many_dft(fft_plan_fz, i1, g(3)%size, isize_fft_z, &
        tmp,   g(3)%size, isize_stride, i1, &
        wrk3d, g(3)%size, isize_stride, i1, FFTW_FORWARD, FFTW_ESTIMATE)

    CALL dfftw_plan_many_dft(fft_plan_bz, i1, g(3)%size, isize_fft_z, &
        tmp,   g(3)%size, isize_stride, i1, &
        wrk3d, g(3)%size, isize_stride, i1, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
    CALL dfftw_plan_many_dft(fft_plan_fz, i1, g(3)%size, isize_fft_z, &
        tmp,   g(3)%size, isize_stride, i1, &
        wrk3d, g(3)%size, isize_stride, i1, FFTW_FORWARD, FFTW_MEASURE)

    CALL dfftw_plan_many_dft(fft_plan_bz, i1, g(3)%size, isize_fft_z, &
        tmp,   g(3)%size, isize_stride, i1, &
        wrk3d, g(3)%size, isize_stride, i1, FFTW_BACKWARD, FFTW_MEASURE)
#endif
  ENDIF

  ! -----------------------------------------------------------------------
  ! Ox direction
  ! -----------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
    isize_fft_x = ims_size_i(DNS_MPI_I_POISSON1)
  ELSE
#endif
    isize_fft_x = jmax
#ifdef USE_MPI
  ENDIF
#endif

#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
    isize_disp = (imax/2+1)*ims_npro_i
  ELSE
#endif
    isize_disp = g(1)%size/2+1
#ifdef _DEBUG
    CALL dfftw_plan_dft_r2c_1d(fft_plan_fx_bcs, g(1)%size, wrk1d,wrk2d, FFTW_ESTIMATE)
#else
    CALL dfftw_plan_dft_r2c_1d(fft_plan_fx_bcs, g(1)%size, wrk1d,wrk2d, FFTW_MEASURE)
#endif
#ifdef USE_MPI
  ENDIF
#endif

#ifdef _DEBUG
  CALL dfftw_plan_many_dft_r2c(fft_plan_fx, i1, g(1)%size, isize_fft_x, &
      tmp,   g(1)%size,     i1, g(1)%size, &
      wrk3d, g(1)%size/2+1, i1, isize_disp, FFTW_ESTIMATE)
  CALL dfftw_plan_many_dft_c2r(fft_plan_bx, i1, g(1)%size, isize_fft_x, &
      tmp,   g(1)%size/2+1, i1, isize_disp, &
      wrk3d, g(1)%size,     i1, g(1)%size, FFTW_ESTIMATE)
#else
  CALL dfftw_plan_many_dft_r2c(fft_plan_fx, i1, g(1)%size, isize_fft_x, &
      tmp,   g(1)%size,     i1, g(1)%size, &
      wrk3d, g(1)%size/2+1, i1, isize_disp, FFTW_MEASURE)
  CALL dfftw_plan_many_dft_c2r(fft_plan_bx, i1, g(1)%size, isize_fft_x, &
      tmp,   g(1)%size/2+1, i1, isize_disp, &
      wrk3d, g(1)%size,     i1, g(1)%size, FFTW_MEASURE)
#endif

  ! -----------------------------------------------------------------------
  ! Oy direction
  ! -----------------------------------------------------------------------
  isize_fft_y = imax/2+1

  isize_stride = isize_fft_y

  IF ( g(2)%size .GT. 1 ) THEN
#ifdef _DEBUG
    CALL dfftw_plan_many_dft(fft_plan_fy, i1, g(2)%size, isize_fft_y, &
        tmp,   g(2)%size, isize_stride, i1, &
        wrk3d, g(2)%size, isize_stride, i1, FFTW_FORWARD, FFTW_ESTIMATE)

    CALL dfftw_plan_many_dft(fft_plan_by, i1, g(2)%size, isize_fft_y, &
        tmp,   g(2)%size, isize_stride, i1, &
        wrk3d, g(2)%size, isize_stride, i1, FFTW_BACKWARD, FFTW_ESTIMATE)
#else
    CALL dfftw_plan_many_dft(fft_plan_fy, i1, g(2)%size, isize_fft_y, &
        tmp,   g(2)%size, isize_stride, i1, &
        wrk3d, g(2)%size, isize_stride, i1, FFTW_FORWARD, FFTW_MEASURE)

    CALL dfftw_plan_many_dft(fft_plan_by, i1, g(2)%size, isize_fft_y, &
        tmp,   g(2)%size, isize_stride, i1, &
        wrk3d, g(2)%size, isize_stride, i1, FFTW_BACKWARD, FFTW_MEASURE)
#endif
  ENDIF

#else
  CALL TLAB_WRITE_ASCII(efile,'OPR_FOURIER_INITIALIZE. FFTW needed for POISSON solver.')
  CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)

#endif

  RETURN
END SUBROUTINE OPR_FOURIER_INITIALIZE

!########################################################################
!########################################################################
SUBROUTINE OPR_FOURIER_F_X_EXEC(nx,ny,nz, in,in_bcs_hb,in_bcs_ht, out, wrk1,wrk2)

  USE DNS_GLOBAL, ONLY : fft_plan_fx, fft_plan_fx_bcs
  USE DNS_GLOBAL, ONLY : isize_txc_dimz, isize_txc_field
  USE DNS_GLOBAL, ONLY : g
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
#endif
  USE, INTRINSIC :: ISO_C_binding, ONLY : c_f_pointer, c_loc

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL,    DIMENSION(nx*ny,*)               :: in
  TREAL,    DIMENSION(nx,nz)                 :: in_bcs_hb,in_bcs_ht
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz)   :: out
  TCOMPLEX, DIMENSION( nx/2+1,            *) :: wrk2
#ifdef USE_MPI
  TCOMPLEX, DIMENSION((nx/2+1)*ims_npro_i,*) :: wrk1
#else
  TCOMPLEX, DIMENSION(g(1)%size/2+1,*)       :: wrk1
#endif

  TARGET wrk1, wrk2, out

  ! -----------------------------------------------------------------------
  TINTEGER j, k, ip, isize_page, isize_line
  TCOMPLEX, DIMENSION(:), POINTER :: p_tmp

#ifdef USE_MPI
  TINTEGER i, id, iold, inew
#endif

  TREAL, POINTER :: r_wrk1(:) => NULL(), r_wrk2(:) => NULL()

  ! #######################################################################
  isize_line = nx/2+1; isize_page = isize_line*ny

  ! #######################################################################
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN

    ! Pass memory address from complex array to real array
    CALL c_f_POINTER(c_LOC(wrk1), r_wrk1, shape=[isize_txc_field])
    CALL c_f_POINTER(c_LOC(wrk2), r_wrk2, shape=[isize_txc_field])

    ! Add bcs into array a; there must be space !
    ip = 1
    ip = ip + nx*ny*nz; in(ip:ip+nx*nz-1,1) = in_bcs_hb(1:nx*nz,1)
    ip = ip + nx*nz;    in(ip:ip+nx*nz-1,1) = in_bcs_ht(1:nx*nz,1)

    ! Transpose array a into b
    id = DNS_MPI_I_POISSON1
    CALL DNS_MPI_TRPF_I(in, r_wrk2, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

    ! ims_size_k(id) FFTWs
    CALL dfftw_execute_dft_r2c(fft_plan_fx, wrk2, wrk1)

    ! reorganize wrk1 (FFTW make a stride in wrk1 already before)
    id = DNS_MPI_I_POISSON1
    DO k = 1,ims_size_i(id)
      inew = (nx/2+1)*ims_npro_i
      iold = g(1)%size/2 + 1
      wrk1(inew,k) = wrk1(iold,k)
      DO ip = ims_npro_i,2,-1
        DO i = nx/2,1,-1
          inew = (ip-1)*(nx/2+1) + i
          iold = (ip-1)* nx/2    + i
          wrk1(inew,k) = wrk1(iold,k)
        ENDDO
      ENDDO
    ENDDO

    ! Transpose array back
    id = DNS_MPI_I_POISSON2
    CALL DNS_MPI_TRPB_I(r_wrk1, r_wrk2, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

    ! reorganize wrk2 into b
    DO k = 1,nz
      ip = 1  + ny*(k-1); out(1:isize_page            ,k) = wrk2(1:isize_page,ip)
      ip = k  + ny* nz;   out(isize_page           +1:,k) = wrk2(1:isize_line,ip)
      ip = ip +     nz;   out(isize_page+isize_line+1:,k) = wrk2(1:isize_line,ip)
    ENDDO

    NULLIFY(r_wrk1,r_wrk2)

  ELSE
#endif

    ! #######################################################################
    DO k = 1,nz
      p_tmp => out( 1:,k)
      CALL dfftw_execute_dft_r2c(fft_plan_fx,     in(1,k),        p_tmp)
      j = 1+ny; ip = (j-1)*isize_line + 1; p_tmp => out(ip:,k)
      CALL dfftw_execute_dft_r2c(fft_plan_fx_bcs, in_bcs_hb(1,k), p_tmp)
      j = 2+ny; ip = (j-1)*isize_line + 1; p_tmp => out(ip:,k)
      CALL dfftw_execute_dft_r2c(fft_plan_fx_bcs, in_bcs_ht(1,k), p_tmp)
    ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE OPR_FOURIER_F_X_EXEC

!########################################################################
!########################################################################
SUBROUTINE OPR_FOURIER_B_X_EXEC(nx,ny,nz, in,out, wrk)

  USE DNS_GLOBAL, ONLY : isize_txc_dimz, isize_txc_field
  USE DNS_GLOBAL, ONLY : fft_plan_bx
#ifdef USE_MPI
  USE DNS_GLOBAL, ONLY : g
  USE TLAB_MPI_VARS, ONLY : ims_npro_i
  USE TLAB_MPI_VARS, ONLY : ims_size_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i
  USE TLAB_MPI_PROCS
#endif
  USE, INTRINSIC :: ISO_C_binding, ONLY : c_f_pointer, c_loc

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz),   INTENT(IN) :: in
#ifdef USE_MPI
  TCOMPLEX, DIMENSION((nx/2+1)*ims_npro_i,*), INTENT(OUT)          :: out
#else
  TREAL,    DIMENSION(nx,ny,nz),              INTENT(OUT)          :: out
#endif
  TCOMPLEX, DIMENSION( nx/2+1,            *), INTENT(INOUT)        :: wrk

  TARGET in, out, wrk

  ! -----------------------------------------------------------------------
#ifdef USE_MPI
  TINTEGER i, k, ip, id, iold, inew, isize_page
#else
  TINTEGER k
#endif
  TCOMPLEX, DIMENSION(:), POINTER :: p_tmp

  TREAL, POINTER :: r_wrk(:) => NULL(), r_out(:) => NULL()

  !########################################################################
  ! Ox Parallel Decomposition
  !########################################################################
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN

    ! Pass memory address from complex array to real array
    CALL c_f_POINTER(c_LOC(wrk), r_wrk, shape=[isize_txc_field])
    CALL c_f_POINTER(c_LOC(out), r_out, shape=[isize_txc_field])

    ! reorganize in into wrk
    isize_page =(nx/2+1)*ny
    DO k = 1,nz
      ip = 1  + ny*(k-1); wrk(1:isize_page,ip) = in(1:isize_page,k)
      !     ip = 1  +ny*(nx/2+1)* nz;   wrk(1: nx/2+1,    ip) = in(   ny   *(nx/2+1)+1:,k) !Idonotneed
      !     ip = ip +   (nx/2+1)* nz;   wrk(1: nx/2+1,    ip) = in(  (ny+1)*(nx/2+1)+1:,k) !Idonotneed
    ENDDO

    ! Transpose array
    id  = DNS_MPI_I_POISSON2
    CALL DNS_MPI_TRPF_I(r_wrk, r_out, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

    ! reorganize a (FFTW make a stride in a already before)
    id = DNS_MPI_I_POISSON1
    DO k = 1,ims_size_i(id)
      DO ip = 2,ims_npro_i
        DO i = 1,nx/2
          iold = (ip-1)*(nx/2+1) + i
          inew = (ip-1)* nx/2    + i
          out(inew,k) = out(iold,k)
        ENDDO
      ENDDO
      iold = ims_npro_i*(nx/2+1)
      inew = g(1)%size/2 + 1
      out(inew,k) = out(iold,k)
    ENDDO

    ! ims_size_i(id) FFTWs
    CALL dfftw_execute_dft_c2r(fft_plan_bx, out, wrk)

    ! Transpose array wrk into out
    id  = DNS_MPI_I_POISSON1
    CALL DNS_MPI_TRPB_I(r_wrk, r_out, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

    NULLIFY(r_wrk,r_out)

  ELSE
#endif

    !########################################################################
    ! No Ox Parallel Decomposition
    !########################################################################
    DO k = 1,nz
      p_tmp => in(1:,k)
#ifdef USE_MPI
      ! divide by two b/o array is declared complex
      ip = (nx*ny*(k-1))/2 + 1
      CALL dfftw_execute_dft_c2r(fft_plan_bx, p_tmp, out(ip,1))
#else
      CALL dfftw_execute_dft_c2r(fft_plan_bx, p_tmp, out(1,1,k))
#endif
    ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE OPR_FOURIER_B_X_EXEC

!########################################################################
!########################################################################
SUBROUTINE OPR_FOURIER_F_Z_EXEC(in,out)

  USE DNS_GLOBAL, ONLY : fft_plan_fz, fft_reordering
  USE DNS_GLOBAL, ONLY : isize_txc_dimz, isize_txc_field
  USE DNS_GLOBAL, ONLY : g
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_PROCS
#endif
  USE, INTRINSIC :: ISO_C_binding, ONLY : c_f_pointer, c_loc

  IMPLICIT NONE

#include "integers.h"

#ifdef USE_MPI
  TCOMPLEX, DIMENSION(ims_size_k(DNS_MPI_K_POISSON)/2,g(3)%size), TARGET :: in,out
#else
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,g(3)%size),                TARGET :: in,out
#endif

  ! -----------------------------------------------------------------------
  TCOMPLEX, DIMENSION(:,:), POINTER :: p_org, p_dst

  TINTEGER k, k_old1, k_old2, k_new1, k_new2
#ifdef USE_MPI
  TINTEGER id
#endif

  TREAL, POINTER :: r_in(:) => NULL(), r_out(:) => NULL()

  ! #######################################################################
  ! Forward complex FFT in z
#ifdef USE_MPI
  id = DNS_MPI_K_POISSON

  IF ( ims_npro_k .GT. 1 ) THEN
    ! Pass memory address from complex array to real array
    CALL c_f_POINTER(c_LOC(in), r_in, shape=[isize_txc_field])
    CALL c_f_POINTER(c_LOC(out), r_out, shape=[isize_txc_field])

    CALL DNS_MPI_TRPF_K(r_in, r_out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    p_org => out
    p_dst => in
  ELSE
#endif
    p_org => in
    p_dst => out
#ifdef USE_MPI
  ENDIF
#endif

  CALL dfftw_execute_dft(fft_plan_fz, p_org, p_dst)

  IF ( fft_reordering .EQ. i1 ) THEN ! re-shuffle spectra in z
    DO k = 1,g(3)%size/2
      k_old1 = k + g(3)%size/2
      k_new1 = k
      k_old2 = k
      k_new2 = k + g(3)%size/2

      p_org(:,k_new1) = p_dst(:,k_old1)
      p_org(:,k_new2) = p_dst(:,k_old2)
    ENDDO
    DO k = 1,g(3)%size
      p_dst(:,k) = p_org(:,k)
    ENDDO
  ENDIF

#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
    CALL DNS_MPI_TRPB_K(r_in, r_out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    NULLIFY(r_in,r_out)
  ENDIF
#endif

  NULLIFY(p_org,p_dst)

  RETURN
END SUBROUTINE OPR_FOURIER_F_Z_EXEC

!########################################################################
!########################################################################
SUBROUTINE OPR_FOURIER_B_Z_EXEC(in,out)

  USE DNS_GLOBAL, ONLY : fft_plan_bz, fft_reordering
  USE DNS_GLOBAL, ONLY : isize_txc_dimz, isize_txc_field
  USE DNS_GLOBAL, ONLY : g
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_k
  USE TLAB_MPI_VARS, ONLY : ims_size_k, ims_ds_k, ims_dr_k, ims_ts_k, ims_tr_k
  USE TLAB_MPI_PROCS
#endif
  USE, INTRINSIC :: ISO_C_binding, ONLY : c_f_pointer, c_loc

  IMPLICIT NONE

#include "integers.h"

#ifdef USE_MPI
  TCOMPLEX, DIMENSION(ims_size_k(DNS_MPI_K_POISSON)/2,g(3)%size), TARGET :: in,out
#else
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,g(3)%size),                TARGET :: in,out
#endif

  ! -----------------------------------------------------------------------
  TCOMPLEX, DIMENSION(:,:), POINTER :: p_org, p_dst

  TINTEGER k, k_old1, k_old2, k_new1, k_new2
#ifdef USE_MPI
  TINTEGER id
#endif

  TREAL, POINTER :: r_in(:) => NULL(), r_out(:) => NULL()

  ! #######################################################################
  ! Forward complex FFT in z
#ifdef USE_MPI
  id = DNS_MPI_K_POISSON

  IF ( ims_npro_k .GT. 1 ) THEN
    ! Pass memory address from complex array to real array
    CALL c_f_POINTER(c_LOC(in), r_in, shape=[isize_txc_field])
    CALL c_f_POINTER(c_LOC(out), r_out, shape=[isize_txc_field])

    CALL DNS_MPI_TRPF_K(r_in, r_out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    p_org => out
    p_dst => in
  ELSE
#endif
    p_org => in
    p_dst => out
#ifdef USE_MPI
  ENDIF
#endif

  IF ( fft_reordering .EQ. i1 ) THEN ! re-shuffle spectra in z
    DO k = 1,g(3)%size/2
      k_new1 = k + g(3)%size/2
      k_old1 = k
      k_new2 = k
      k_old2 = k + g(3)%size/2

      p_dst(:,k_new1) = p_org(:,k_old1)
      p_dst(:,k_new2) = p_org(:,k_old2)
    ENDDO
    DO k = 1,g(3)%size
      p_org(:,k) = p_dst(:,k)
    ENDDO
  ENDIF

  CALL dfftw_execute_dft(fft_plan_bz, p_org, p_dst)

#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
    CALL DNS_MPI_TRPB_K(r_in, r_out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
    NULLIFY(r_in,r_out)
  ENDIF
#endif

  NULLIFY(p_org,p_dst)

  RETURN
END SUBROUTINE OPR_FOURIER_B_Z_EXEC

!########################################################################
!########################################################################
! SUBROUTINE OPR_FOURIER_B_Z_EXEC(nx,ny,nz, in,out)

!   USE DNS_GLOBAL, ONLY : fft_plan_bz, isize_txc_dimz
! #ifdef USE_MPI
!   USE TLAB_MPI_VARS,    ONLY : ims_npro_k, ims_err, ims_ts_k, ims_tr_k, ims_ds_k, ims_dr_k
! #endif

!   IMPLICIT NONE

! #include "integers.h"
! #ifdef USE_FFTW
! #include "fftw3.f"
! #endif

!   TINTEGER                                                 :: nx,ny,nz
!   TCOMPLEX, DIMENSION(isize_txc_dimz/2, nz), INTENT(INOUT) :: in, out

! ! -----------------------------------------------------------------------
! #ifdef USE_MPI
!   TINTEGER id
! #endif

! !########################################################################
! ! backwards complex FFT in z
! #ifdef USE_MPI
!   id = DNS_MPI_K_POISSON

!   IF ( ims_npro_k .GT. 1 ) THEN
!      CALL DNS_MPI_TRPF_K(in, out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

!      CALL dfftw_execute_dft(fft_plan_bz,out,in)

!      CALL DNS_MPI_TRPB_K(in, out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

!   ELSE ! ims_npro_k .EQ. 1
! #endif

!      CALL dfftw_execute_dft(fft_plan_bz,in,out)

! #ifdef USE_MPI
!   ENDIF ! ims_npro_k .GT. 1
! #endif

!   RETURN
! END SUBROUTINE OPR_FOURIER_B_Z_EXEC

! #######################################################################
! #######################################################################
SUBROUTINE OPR_FOURIER_SPECTRA_3D(nx,ny,nz, isize_psd, u, psd, wrk1d)

  USE DNS_GLOBAL, ONLY : isize_txc_dimz
  USE DNS_GLOBAL, ONLY : g
#ifdef USE_MPI
  USE TLAB_MPI_VARS,    ONLY : ims_offset_i, ims_offset_k, ims_pro, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER ,                           INTENT(IN)    :: nx,ny,nz, isize_psd
  TREAL, DIMENSION(isize_txc_dimz,nz), INTENT(IN)    :: u
  TREAL, DIMENSION(isize_psd),         INTENT(OUT)   :: psd
  TREAL, DIMENSION(isize_psd),         INTENT(INOUT) :: wrk1d

  ! -----------------------------------------------------------------------
  TINTEGER i,j,k, r, iglobal,kglobal, ip
  TREAL fr, fi,fj,fk

  ! #######################################################################
  psd = C_0_R

  DO k = 1,nz
#ifdef USE_MPI
    kglobal = k + ims_offset_k
#else
    kglobal = k
#endif
    IF ( kglobal .LE. g(3)%size/2+1 ) THEN; fk = M_REAL(kglobal-1)
    ELSE;                                   fk =-M_REAL(g(3)%size+1-kglobal)
    ENDIF

    DO j = 1,ny
      IF ( j       .LE. g(2)%size/2+1 ) THEN; fj = M_REAL(j-1)
      ELSE;                                   fj =-M_REAL(g(2)%size+1-j)
      ENDIF


      DO i = 1,nx/2+1
#ifdef USE_MPI
        iglobal = i + ims_offset_i/2
#else
        iglobal = i
#endif
        IF ( iglobal .LE. g(1)%size/2+1 ) THEN; fi = M_REAL(iglobal-1)
        ELSE;                                   fi =-M_REAL(g(1)%size+1-iglobal)
        ENDIF


        fr = CEILING( SQRT( M_REAL(fi**2 + fj**2 + fk**2) ) )

        ! -----------------------------------------------------------------------
        ! psd
        ! -----------------------------------------------------------------------
        ip = (nx+2)*(j-1) + 2*i

        r = INT(fr)
        IF ( r .EQ. 0 ) THEN ! zero is not written
          !              print*, i,j,k, r
        ELSE
          psd(r) = psd(r) + u(ip-1,k)**2 + u(ip,k)**2
        ENDIF

      ENDDO
    ENDDO
  ENDDO

#ifdef USE_MPI
  wrk1d = psd
  CALL MPI_Reduce(wrk1d, psd, isize_psd, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
  IF ( ims_pro .EQ. 0 ) THEN
    psd = wrk1d
  ENDIF
#endif

  RETURN
END SUBROUTINE OPR_FOURIER_SPECTRA_3D
