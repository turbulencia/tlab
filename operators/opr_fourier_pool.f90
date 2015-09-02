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

  USE DNS_GLOBAL, ONLY : fft_plan_fx_bcs,fft_plan_fx,fft_plan_bx
  USE DNS_GLOBAL, ONLY : fft_plan_fy,fft_plan_by
  USE DNS_GLOBAL, ONLY : fft_plan_fz,fft_plan_bz
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax, imax_total,jmax_total,kmax_total
  USE DNS_GLOBAL, ONLY : isize_txc_field, isize_txc_dimz
  USE DNS_CONSTANTS, ONLY : efile

#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif
#ifdef USE_FFTW
#include "fftw3.f"
#endif

!  TREAL, DIMENSION(*) :: tmp, wrk1d, wrk2d, wrk3d
  TREAL, DIMENSION(isize_txc_field), INTENT(INOUT) :: tmp, wrk3d
  TREAL, DIMENSION(imax_total+2),    INTENT(INOUT) :: wrk1d, wrk2d

! -----------------------------------------------------------------------
  TINTEGER isize_stride, isize_disp, isize_fft_z, isize_fft_y, isize_fft_x

! #######################################################################
! FFTW library
! #######################################################################
#ifdef USE_FFTW

#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     IF ( ims_size_i(DNS_MPI_I_POISSON1) .NE. ims_size_i(DNS_MPI_I_POISSON2) ) THEN
        CALL IO_WRITE_ASCII(efile,'OPR_FOURIER_INITIALIZE. Error in the size in the transposition arrays.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
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

  IF ( kmax_total .GT. 1 ) THEN
     CALL dfftw_plan_many_dft(fft_plan_fz, i1, kmax_total, isize_fft_z, &
          tmp,   kmax_total, isize_stride, i1, &
          wrk3d, kmax_total, isize_stride, i1, FFTW_FORWARD, FFTW_MEASURE)
     
     CALL dfftw_plan_many_dft(fft_plan_bz, i1, kmax_total, isize_fft_z, &
          tmp,   kmax_total, isize_stride, i1, &
          wrk3d, kmax_total, isize_stride, i1, FFTW_BACKWARD, FFTW_MEASURE)
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
     isize_disp = imax_total/2+1
     CALL dfftw_plan_dft_r2c_1d(fft_plan_fx_bcs, imax_total, wrk1d,wrk2d, FFTW_MEASURE)
#ifdef USE_MPI
  ENDIF
#endif

  CALL dfftw_plan_many_dft_r2c(fft_plan_fx, i1, imax_total, isize_fft_x, &
       tmp,   imax_total,     i1, imax_total, &
       wrk3d, imax_total/2+1, i1, isize_disp, FFTW_MEASURE)
  CALL dfftw_plan_many_dft_c2r(fft_plan_bx, i1, imax_total, isize_fft_x, &
       tmp,   imax_total/2+1, i1, isize_disp, &
       wrk3d, imax_total,     i1, imax_total, FFTW_MEASURE)
  
! -----------------------------------------------------------------------
! Oy direction
! -----------------------------------------------------------------------
  isize_fft_y = imax/2+1

  isize_stride = isize_fft_y

  IF ( jmax_total .GT. 1 ) THEN
     CALL dfftw_plan_many_dft(fft_plan_fy, i1, jmax_total, isize_fft_y, &
          tmp,   jmax_total, isize_stride, i1, &
          wrk3d, jmax_total, isize_stride, i1, FFTW_FORWARD, FFTW_MEASURE)
     
     CALL dfftw_plan_many_dft(fft_plan_by, i1, jmax_total, isize_fft_y, &
          tmp,   jmax_total, isize_stride, i1, &
          wrk3d, jmax_total, isize_stride, i1, FFTW_BACKWARD, FFTW_MEASURE)
  ENDIF

#else
  CALL IO_WRITE_ASCII(efile,'OPR_FOURIER_INITIALIZE. FFTW needed for POISSON solver.')
  CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     
#endif
  
  RETURN
END SUBROUTINE OPR_FOURIER_INITIALIZE

!########################################################################
!########################################################################
SUBROUTINE OPR_FOURIER_F_X_EXEC(nx,ny,nz, in,in_bcs_hb,in_bcs_ht, out, wrk1,wrk2)

  USE DNS_GLOBAL, ONLY : imax_total, isize_txc_dimz
  USE DNS_GLOBAL, ONLY : fft_plan_fx, fft_plan_fx_bcs
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TREAL,    DIMENSION(nx*ny,*)                     :: in
  TREAL,    DIMENSION(nx,nz)                       :: in_bcs_hb,in_bcs_ht  
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz), TARGET :: out
  TCOMPLEX, DIMENSION( nx/2+1,            *)       :: wrk2
#ifdef USE_MPI
  TCOMPLEX, DIMENSION((nx/2+1)*ims_npro_i,*)       :: wrk1
#else
  TCOMPLEX, DIMENSION(imax_total/2+1,*)            :: wrk1
#endif

! -----------------------------------------------------------------------
  TINTEGER j, k, ip, isize_page, isize_line
  TCOMPLEX, DIMENSION(:), POINTER :: p_tmp

#ifdef USE_MPI
  TINTEGER i, id, iold, inew
#endif

! #######################################################################
  isize_line = nx/2+1; isize_page = isize_line*ny

! #######################################################################
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN

! Add bcs into array a; there must be space !
     ip = 1
     ip = ip + nx*ny*nz; in(ip:ip+nx*nz-1,1) = in_bcs_hb(1:nx*nz,1)
     ip = ip + nx*nz;    in(ip:ip+nx*nz-1,1) = in_bcs_ht(1:nx*nz,1)
     
! Transpose array a into b
     id = DNS_MPI_I_POISSON1
     CALL DNS_MPI_TRPF_I(in, wrk2, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     
! ims_size_k(id) FFTWs
     CALL dfftw_execute_dft_r2c(fft_plan_fx, wrk2, wrk1)

! reorganize wrk1 (FFTW make a stride in wrk1 already before)
     id = DNS_MPI_I_POISSON1
     DO k = 1,ims_size_i(id)
        inew = (nx/2+1)*ims_npro_i
        iold = imax_total/2 + 1
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
     CALL DNS_MPI_TRPB_I(wrk1, wrk2, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))
     
! reorganize wrk2 into b
     DO k = 1,nz
        ip = 1  + ny*(k-1); out(1:isize_page            ,k) = wrk2(1:isize_page,ip)
        ip = k  + ny* nz;   out(isize_page           +1:,k) = wrk2(1:isize_line,ip)
        ip = ip +     nz;   out(isize_page+isize_line+1:,k) = wrk2(1:isize_line,ip)
     ENDDO

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

  USE DNS_GLOBAL, ONLY : imax_total, isize_txc_dimz
  USE DNS_GLOBAL, ONLY : fft_plan_bx
#ifdef USE_MPI 
  USE DNS_MPI,    ONLY : ims_npro_i, ims_ds_i, ims_dr_i, ims_ts_i, ims_tr_i, ims_size_i
#endif 

  IMPLICIT NONE

  TINTEGER nx,ny,nz
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz),   INTENT(IN),   TARGET :: in
#ifdef USE_MPI 
  TCOMPLEX, DIMENSION((nx/2+1)*ims_npro_i,*), INTENT(OUT)          :: out
#else
  TREAL,    DIMENSION(nx,ny,nz),              INTENT(OUT)          :: out
#endif
  TCOMPLEX, DIMENSION( nx/2+1,            *), INTENT(INOUT)        :: wrk

! -----------------------------------------------------------------------
#ifdef USE_MPI 
  TINTEGER i, k, ip, id, iold, inew, isize_page
#else
  TINTEGER k
#endif
  TCOMPLEX, DIMENSION(:), POINTER :: p_tmp

!########################################################################
! Ox Parallel Decomposition
!########################################################################
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN

! reorganize in into wrk
  isize_page =(nx/2+1)*ny
  DO k = 1,nz
     ip = 1  + ny*(k-1); wrk(1:isize_page,ip) = in(1:isize_page,k) 
!     ip = 1  +ny*(nx/2+1)* nz;   wrk(1: nx/2+1,    ip) = in(   ny   *(nx/2+1)+1:,k) !Idonotneed
!     ip = ip +   (nx/2+1)* nz;   wrk(1: nx/2+1,    ip) = in(  (ny+1)*(nx/2+1)+1:,k) !Idonotneed
  ENDDO

! Transpose array
  id  = DNS_MPI_I_POISSON2
  CALL DNS_MPI_TRPF_I(wrk, out, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

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
     inew = imax_total/2 + 1
     out(inew,k) = out(iold,k)
  ENDDO
  
! ims_size_i(id) FFTWs
  CALL dfftw_execute_dft_c2r(fft_plan_bx, out, wrk)

! Transpose array wrk into out
  id  = DNS_MPI_I_POISSON1
  CALL DNS_MPI_TRPB_I(wrk, out, ims_ds_i(1,id), ims_dr_i(1,id), ims_ts_i(1,id), ims_tr_i(1,id))

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
  USE DNS_GLOBAL, ONLY : isize_txc_dimz, imax_total, jmax_total, kmax_total 
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_npro_k, ims_ts_k, ims_tr_k, ims_ds_k, ims_dr_k, ims_err,ims_size_k
#endif 

  IMPLICIT NONE 
  
#include "integers.h"

#ifdef USE_MPI
  TCOMPLEX, DIMENSION(ims_size_k(DNS_MPI_K_POISSON)/2,kmax_total), TARGET :: in,out
#else 
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,kmax_total),                TARGET :: in,out
#endif

! -----------------------------------------------------------------------
  TCOMPLEX, DIMENSION(:,:), POINTER :: p_org, p_dst

  TINTEGER k, k_old1, k_old2, k_new1, k_new2
#ifdef USE_MPI
  TINTEGER id
#endif

! #######################################################################
! Forward complex FFT in z
#ifdef USE_MPI
  id = DNS_MPI_K_POISSON

  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(in, out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
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
     DO k = 1,kmax_total/2
        k_old1 = k + kmax_total/2
        k_new1 = k 
        k_old2 = k
        k_new2 = k + kmax_total/2

        p_org(:,k_new1) = p_dst(:,k_old1) 
        p_org(:,k_new2) = p_dst(:,k_old2)  
     ENDDO
     DO k = 1,kmax_total
        p_dst(:,k) = p_org(:,k)
     ENDDO
  ENDIF

#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(in, out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
  ENDIF
#endif

  NULLIFY(p_org,p_dst)

  RETURN
END SUBROUTINE OPR_FOURIER_F_Z_EXEC

!########################################################################
!########################################################################
SUBROUTINE OPR_FOURIER_B_Z_EXEC(in,out) 
  
  USE DNS_GLOBAL, ONLY : fft_plan_bz, fft_reordering 
  USE DNS_GLOBAL, ONLY : isize_txc_dimz, kmax_total 
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_npro_k, ims_ts_k, ims_tr_k, ims_ds_k, ims_dr_k, ims_size_k
#endif 

  IMPLICIT NONE 
  
#include "integers.h"

#ifdef USE_MPI
  TCOMPLEX, DIMENSION(ims_size_k(DNS_MPI_K_POISSON)/2,kmax_total), TARGET :: in,out
#else 
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,kmax_total),                TARGET :: in,out
#endif

! -----------------------------------------------------------------------
  TCOMPLEX, DIMENSION(:,:), POINTER :: p_org, p_dst

  TINTEGER k, k_old1, k_old2, k_new1, k_new2
#ifdef USE_MPI
  TINTEGER id
#endif

! #######################################################################
! Forward complex FFT in z
#ifdef USE_MPI
  id = DNS_MPI_K_POISSON

  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPF_K(in, out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
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
     DO k = 1,kmax_total/2
        k_new1 = k + kmax_total/2
        k_old1 = k 
        k_new2 = k
        k_old2 = k + kmax_total/2

        p_dst(:,k_new1) = p_org(:,k_old1) 
        p_dst(:,k_new2) = p_org(:,k_old2)  
     ENDDO
     DO k = 1,kmax_total
        p_org(:,k) = p_dst(:,k)
     ENDDO
  ENDIF

  CALL dfftw_execute_dft(fft_plan_bz, p_org, p_dst)

#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     CALL DNS_MPI_TRPB_K(in, out, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
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
!   USE DNS_MPI,    ONLY : ims_npro_k, ims_err, ims_ts_k, ims_tr_k, ims_ds_k, ims_dr_k 
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

  USE DNS_GLOBAL, ONLY : imax_total, jmax_total, kmax_total, isize_txc_dimz
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_offset_i, ims_offset_k, ims_pro, ims_err 
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
     IF ( kglobal .LE. kmax_total/2+1 ) THEN; fk = M_REAL(kglobal-1)
     ELSE;                                    fk =-M_REAL(kmax_total+1-kglobal); ENDIF

     DO j = 1,ny
        IF ( j       .LE. jmax_total/2+1 ) THEN; fj = M_REAL(j-1)
        ELSE;                                    fj =-M_REAL(jmax_total+1-j); ENDIF
              

        DO i = 1,nx/2+1
#ifdef USE_MPI
           iglobal = i + ims_offset_i/2
#else
           iglobal = i
#endif
           IF ( iglobal .LE. imax_total/2+1 ) THEN; fi = M_REAL(iglobal-1)
           ELSE;                                    fi =-M_REAL(imax_total+1-iglobal); ENDIF


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

