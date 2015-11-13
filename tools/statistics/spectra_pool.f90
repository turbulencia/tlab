#include "types.h" 
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

!########################################################################
!########################################################################
SUBROUTINE INTEGRATE_SPECTRUM(nx,ny,nz, kr_total, isize_aux, &
     spec_2d, data_x,data_z,spec_r, tmp_x,tmp_z,wrk2d)

  USE DNS_GLOBAL, ONLY : kmax_total
#ifdef USE_MPI 
  USE DNS_MPI
#endif 

  IMPLICIT NONE 

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN)  :: nx,ny,nz, kr_total, isize_aux
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN)  :: spec_2d ! power spectral density

  TREAL, DIMENSION(kr_total,ny), INTENT(OUT) :: spec_r 

  TREAL, DIMENSION(nx,ny),          INTENT(OUT)   :: data_x
  TREAL, DIMENSION(nx,ny,2),        INTENT(INOUT) :: tmp_x
  TREAL, DIMENSION(nz/2,ny),        INTENT(OUT)   :: data_z 
  TREAL, DIMENSION(isize_aux,nz,2), INTENT(INOUT) :: tmp_z ! need space for transpostion

#ifdef USE_MPI 
  TREAL, DIMENSION(isize_aux/ims_npro_k,kmax_total,2), INTENT(INOUT) :: wrk2d
#else
  TREAL, DIMENSION(ny,kmax_total,2), INTENT(INOUT) :: wrk2d
#endif

! -----------------------------------------------------------------------
  TINTEGER :: i,k, kx_global,kz_global,kr_global, flag, ny_local
#ifdef USE_MPI 
  TINTEGER count, id
#endif

! #######################################################################
  tmp_x = 0; tmp_z = 0; wrk2d = 0

  DO k=1,nz  
#ifdef USE_MPI 
     kz_global = k+ ims_offset_k 
#else 
     kz_global = k
#endif
     
     tmp_x(:,:,1) = tmp_x(:,:,1) + spec_2d(:,:,k)

     IF ( kz_global .LE. kmax_total/2 ) THEN; kz_global = kmax_total/2 - kz_global + 2; flag = 1
     ELSE;                                    kz_global = kz_global - kmax_total/2;     flag = 0
     ENDIF
! drop the Nyquist frequency; we add it to the previous mode to keep structure below
     IF ( kz_global .EQ. kmax_total/2+1 ) kz_global = MAX(kmax_total/2,1)
        
     DO i=1,nx
#ifdef USE_MPI 
        kx_global = i+ ims_offset_i/2 
#else 
        kx_global = i
#endif 

        tmp_z(1:ny,k,1) = tmp_z(1:ny,k,1) + spec_2d(i,1:ny,k)

! need to use CEILING here; otherwise the zero-wavenumber 
! for radial spectra gets contaminated 
!        kr_global = CEILING( SQRT( M_REAL(kx_global*kx_global + kz_global*kz_global) ) )
        kr_global = INT( SQRT( M_REAL((kx_global-1)**2 + (kz_global-1)**2) ) ) + 1           
        spec_r(kr_global,:) = spec_r(kr_global,:) + spec_2d(i,:,k)

! correction; half of the modes kx_global=1 have been counted twice 
        IF ( kx_global .EQ. 1 .AND. flag .EQ. 1 ) THEN
           tmp_z(1:ny,k,1) = tmp_z(1:ny,k,1) - spec_2d(i,1:ny,k)
 
           spec_r(kr_global,:) = spec_r(kr_global,:) - spec_2d(i,:,k)
        ENDIF

! correction; half of the modes kz_global=1 need to be counted twice
        IF ( kz_global .EQ. 1 .AND. kx_global .GT. 1 ) THEN
           tmp_z(1:ny,k,1) = tmp_z(1:ny,k,1) + spec_2d(i,1:ny,k)
        ENDIF
        
     ENDDO
  ENDDO

! Finalize Ox spectrum
#ifdef USE_MPI
  count = nx*ny
  CALL MPI_ALLREDUCE(tmp_x(:,:,1), tmp_x(:,:,2), count, MPI_REAL8, MPI_SUM, ims_comm_z, ims_err)
  data_x(:,:) = data_x(:,:) + tmp_x(:,:,2)

#else
  data_x(:,:) = data_x(:,:) + tmp_x(:,:,1)

#endif

! Finalize Oz spectrum
#ifdef USE_MPI
  count = isize_aux*nz
  CALL MPI_ALLREDUCE(tmp_z(:,:,1), tmp_z(:,:,2), count, MPI_REAL8, MPI_SUM, ims_comm_x, ims_err)

  IF ( ims_npro_k .GT. 1 ) THEN
     id = DNS_MPI_K_AUX2
     CALL DNS_MPI_TRPF_K(tmp_z(:,:,2), wrk2d(:,:,1), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))

  ELSE
     wrk2d(1:ny*nz,1,1) = tmp_z(1:ny*nz,1,2)

  ENDIF

  ny_local = isize_aux/ims_npro_k

#else
  wrk2d(1:ny*nz,1,1) = tmp_z(1:ny*nz,1,1)

  ny_local = ny

#endif


  DO k = 1,kmax_total
     kz_global = k
     IF ( kz_global .LE. kmax_total/2 ) THEN; kz_global = kmax_total/2 - kz_global + 2
     ELSE;                                    kz_global = kz_global - kmax_total/2
     ENDIF
! drop the Nyquist frequency; we add it to the previous mode to keep structure below
     IF ( kz_global .EQ. kmax_total/2+1 ) kz_global = MAX(kmax_total/2,1)

     wrk2d(1:ny_local,kz_global,2) = wrk2d(1:ny_local,kz_global,2) + wrk2d(1:ny_local,k,1)
        
  ENDDO

#ifdef USE_MPI
  IF ( ims_npro_k .GT. 1 ) THEN
     count = kmax_total/2 /ims_npro_k ! add strides for the transposition
     DO k = 1,kmax_total/2,count
        wrk2d(1:ny_local*count,(k-1)*2+1,1) =  wrk2d(1:ny_local*count,k,2)
     ENDDO

     CALL DNS_MPI_TRPB_K(wrk2d(:,:,1), tmp_z(:,:,1), ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     
  ELSE
#endif

     tmp_z(1:isize_aux*nz,1,1) = wrk2d(1:ny_local*kmax_total,1,2)

#ifdef USE_MPI
  ENDIF
#endif

  DO k = 1,nz/2
     data_z(k,1:ny) = data_z(k,1:ny) + tmp_z(1:ny,k,1)
  ENDDO

  RETURN
END SUBROUTINE INTEGRATE_SPECTRUM

!########################################################################
!########################################################################
SUBROUTINE REDUCE_SPECTRUM(nx,ny,nz, nblock, in,out, tmp1,variance)

  USE DNS_GLOBAL, ONLY : isize_txc_dimz

! need to know about domain decomposition in x b/o 
! nyquist frequency and zero frequency account different for the variance 
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_offset_i, ims_pro_i, ims_npro_i, ims_err 
#endif

  IMPLICIT NONE 

#include "integers.h" 
#ifdef USE_MPI
#include "mpif.h"  
#endif

  TINTEGER,                                 INTENT(IN)  :: nx,ny,nz, nblock
  TCOMPLEX, DIMENSION(isize_txc_dimz/2,nz), INTENT(IN)  :: in, tmp1
  TREAL,    DIMENSION(nx/2,ny/nblock,2*nz), INTENT(OUT) :: out ! Amplitude (1:nz) and Phase (nz+1:2*nz)
  TREAL,    DIMENSION(ny,2)                             :: variance

! -----------------------------------------------------------------------
  TINTEGER :: kx,y,kz, ipy, ip, kx_global, ny_loc
  TCOMPLEX :: cdummy
  TREAL    :: power

! #######################################################################
! Calculate PSD and phase
! #######################################################################
  variance = C_0_R ! use variance to control result

! Drop the uppermost ny%nblock lines as there would not be 
! nblock levels contributing to the output 
  ny_loc = ny - MOD(ny,nblock) 

  DO kz=1,nz

     DO y=1,ny_loc
        ipy = (y-1)/nblock + 1 

! Drop the Nyquist frequency nx/2+1; 
! keeping it would make writing inhomogeneous across processors 
        DO kx=1,nx/2
#ifdef USE_MPI
     kx_global = kx + ims_offset_i
#else 
     kx_global = kx
#endif  
           ip = (nx/2+1) * (y-1) + kx 

           power = REAL(in(ip,kz))
           out(kx, ipy, kz) = out(kx, ipy, kz) + power

! phase; calculate phase only if there is some power 
           IF ( power .GE. 1.e-16 ) THEN
              cdummy = tmp1(ip,kz)
              out(kx, ipy, kz+nz) = out(kx, ipy, kz+nz) + ATAN2(AIMAG(cdummy),REAL(cdummy))
           ENDIF

! use variance to control result
           IF ( kx_global .EQ. i1 ) THEN; variance(y,1) = variance(y,1) + power 
           ELSE;                          variance(y,1) = variance(y,1) + power*C_2_R
           ENDIF

        ENDDO

! add variance from nyquist frequency for the Parseval identity to hold
#ifdef USE_MPI
        IF ( ims_pro_i .EQ. ims_npro_i - i1 ) THEN 
#endif
           ip = (nx/2+1) * (y-1) + nx/2+1 
           variance(y,1) = variance(y,1) + REAL(in(ip,kz))
#ifdef USE_MPI
        ENDIF
#endif 

     ENDDO
  ENDDO

! use variance to control result
#ifdef USE_MPI 
  variance(:,2) = variance(:,1)
  CALL MPI_AllReduce(variance(1,2),variance(1,1),ny_loc,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ims_err) 
#endif 

  RETURN

END SUBROUTINE REDUCE_SPECTRUM

!########################################################################
!########################################################################
! Retaining just 1/4 of the lags domain
! ARGUMENTS: 
! icalc_radial   - whether to reduce radial correlations or not 
SUBROUTINE REDUCE_CORRELATION(nx,ny,nz, nblock, nr_total, &
     in, data_2d,data_x,data_z,data_r, variance1, variance2,icalc_radial_in)

  USE DNS_GLOBAL, ONLY : isize_wrk1d,imax_total,kmax_total
#ifdef USE_MPI 
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k
#endif 

  IMPLICIT NONE 

  TINTEGER,                             INTENT(IN)  :: nx,ny,nz, nblock, nr_total
  TINTEGER,OPTIONAL,                    INTENT(IN)  :: icalc_radial_in  ! defaults to 1 
  TREAL, DIMENSION(nx,ny,nz),           INTENT(IN)  :: in
  TREAL, DIMENSION(nx,ny/nblock,nz),    INTENT(OUT) :: data_2d 
  TREAL, DIMENSION(nx,ny/nblock),       INTENT(OUT) :: data_x 
  TREAL, DIMENSION(nz,ny/nblock),       INTENT(OUT) :: data_z 
  TREAL, DIMENSION(nr_total,ny/nblock), INTENT(OUT) :: data_r
  TREAL, DIMENSION(isize_wrk1d,2),      INTENT(IN)  :: variance1 ! to normalize
  TREAL, DIMENSION(isize_wrk1d),        INTENT(OUT) :: variance2 ! to validate

! -----------------------------------------------------------------------
  TINTEGER i,j,k, ipy,ny_loc, i_global,k_global, r_global, icalc_radial 
  TREAL norm

  IF ( present(icalc_radial_in) ) THEN 
     icalc_radial = icalc_radial_in 
  ELSE 
     icalc_radial = 1 
  ENDIF

! #######################################################################
  variance2(:) = C_0_R ! use variance to control result

! Drop the uppermost ny%nblock lines as there would not be 
! nblock levels contributing to the output 
  ny_loc = ny - MOD(ny,nblock) 
  
  DO j=1,ny_loc
     ipy = (j-1)/nblock + 1 
     norm = SQRT(variance1(j,1))*SQRT(variance1(j,2))
     IF ( norm .GT. C_0_R ) THEN; norm = C_1_R/norm
     ELSE;                        norm = C_1_R
     ENDIF
     
     DO k=1,nz  
#ifdef USE_MPI 
        k_global = k+ ims_offset_k 
#else 
        k_global = k
#endif
        
        DO i=1,nx
#ifdef USE_MPI 
           i_global = i+ ims_offset_i 
#else 
           i_global = i
#endif 
           
! Reduce 2D data
           data_2d(i, ipy, k) = data_2d(i, ipy, k) + in(i,j,k)*norm
           
! Reduce 1D data
           IF ( k_global .EQ. 1 ) data_x(i,ipy) = data_x(i,ipy) + in(i,j,k)*norm

           IF ( i_global .EQ. 1 ) data_z(k,ipy) = data_z(k,ipy) + in(i,j,k)*norm

           ! Avoid contamination by periodicity of domain 
           ! --> has to be taken into account in RADIAL_SAMPLESIZE as well 
           IF ( icalc_radial .EQ. 1 .AND. & 
                k_global .LT. kmax_total/2 .AND. i_global .LT. imax_total/2 ) THEN 
              
              !r_global = SQRT( M_REAL( (k_global-1)**2 + (i_global-1)**2) )  + 1
              !r_global = CEILING(SQRT( M_REAL( (k_global-1)**2 + (i_global-1)**2) )  + 1)
              !IF ( r_global .LE. nr_total ) data_r(r_global,ipy) = data_r(r_global,ipy) + data_2d(i,ipy,k)
              ! Check: Do we still need the 'IF'? 

              r_global = INT(SQRT( M_REAL( (k_global-1)**2 + (i_global-1)**2)))  + 1
              IF ( r_global .LE. nr_total ) data_r(r_global,ipy) = data_r(r_global,ipy) + in(i,j,k)

           ENDIF
! use variance to control result
           IF ( i_global .EQ. 1 .AND. k_global .EQ. 1 ) variance2(j) = in(i,j,k)
           
        ENDDO
     ENDDO
  ENDDO
     
  RETURN
END SUBROUTINE REDUCE_CORRELATION

!########################################################################
!########################################################################
SUBROUTINE RADIAL_SAMPLESIZE(nx,ny,nz, nr_total, samplesize)

  USE DNS_GLOBAL, ONLY : imax_total, kmax_total 

#ifdef USE_MPI 
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k
#endif 

  IMPLICIT NONE 

  TINTEGER,                   INTENT(IN)  :: nx,ny,nz, nr_total
  TREAL, DIMENSION(nr_total), INTENT(OUT) :: samplesize

! -----------------------------------------------------------------------
  TINTEGER i,k, i_global,k_global, r_global

! #######################################################################
  DO k=1,nz  
#ifdef USE_MPI 
     k_global = k+ ims_offset_k 
#else 
     k_global = k
#endif
        
     DO i=1,nx
#ifdef USE_MPI 
        i_global = i+ ims_offset_i 
#else 
        i_global = i
#endif 
        ! Avoid contamination by periodicity of the domain. 
        IF (i_global .LT. imax_total/2 .AND. k_global.LT. kmax_total/2) THEN  
           !r_global = SQRT( M_REAL( (k_global-1)**2 + (i_global-1)**2) ) + 1
           r_global = INT(SQRT( M_REAL( (k_global-1)**2 + (i_global-1)**2) )) + 1
           IF ( r_global .LE. nr_total ) THEN
              samplesize(r_global) = samplesize(r_global) + C_1_R
           ENDIF
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE RADIAL_SAMPLESIZE

!########################################################################
!########################################################################
SUBROUTINE WRITE_SPECTRUM1D(fname, varname, nxy, nvar, pow)

  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : ims_pro, ims_err
#endif

  IMPLICIT NONE

  CHARACTER*(*),                 INTENT(IN) :: fname
  CHARACTER*32, DIMENSION(nvar), INTENT(IN) :: varname(nvar)
  TINTEGER,                      INTENT(IN) :: nxy, nvar
  TREAL, DIMENSION(nxy,nvar),    INTENT(IN) :: pow

! -----------------------------------------------------------------------
  TINTEGER iv
  CHARACTER*64 name

! #######################################################################
  IF ( nxy .EQ. 1 ) RETURN

#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

     DO iv = 1,nvar
        name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))

#include "dns_open_file.h"
        WRITE(LOC_UNIT_ID) SNGL(pow(1:nxy,iv))
        CLOSE(LOC_UNIT_ID)

     ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

END SUBROUTINE WRITE_SPECTRUM1D

!########################################################################
!########################################################################
#ifdef USE_MPI

SUBROUTINE DEFINE_ARRAY_TYPE(opt_main, nblock, subarray)

  USE DNS_GLOBAL, ONLY : imax_total,jmax_total,kmax_total, imax,jmax,kmax
  USE DNS_MPI,    ONLY : ims_offset_i, ims_offset_j, ims_offset_k, ims_err 

  IMPLICIT NONE 

#include "mpif.h" 

  TINTEGER, INTENT(IN)                :: opt_main, nblock 
  TINTEGER, INTENT(OUT), DIMENSION(3) :: subarray

! -----------------------------------------------------------------------
  TINTEGER                :: ndims, i
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset

! #######################################################################
  IF     ( opt_main .EQ. 1 .OR. opt_main .EQ. 2 ) THEN
     ndims = 2 ! Subarray for the output of the Ox spectrum
     sizes(1)   = imax_total/2;   sizes(2)   = jmax_total/nblock
     locsize(1) = imax/2;         locsize(2) = jmax/nblock
     offset(1)  = ims_offset_i/2; offset(2)  = ims_offset_j/nblock

     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, subarray(1), ims_err)
     CALL MPI_Type_commit(subarray(1), ims_err)

! -----------------------------------------------------------------------
     ndims = 2 ! Subarray for the output of the Oz spectrum
     sizes(1)   = kmax_total/2;   sizes(2)   = jmax_total/nblock
     locsize(1) = kmax/2;         locsize(2) = jmax/nblock
     offset(1)  = ims_offset_k/2; offset(2)  = ims_offset_j/nblock

     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, subarray(2), ims_err)
     CALL MPI_Type_commit(subarray(2), ims_err)

! -----------------------------------------------------------------------
     ndims = 3 ! Subarray for the output of the 2D data
     sizes(1)   = imax_total/2;   sizes(2)   = jmax_total/nblock;   sizes(3)   = kmax_total 
     locsize(1) = imax/2;         locsize(2) = jmax/nblock;         locsize(3) = kmax 
     offset(1)  = ims_offset_i/2; offset(2)  = ims_offset_j/nblock; offset(3)  = ims_offset_k
     
     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, subarray(3), ims_err)
     CALL MPI_Type_commit(subarray(3), ims_err)

! -----------------------------------------------------------------------
! I do not know how to save just have of the data
  ! ELSE IF ( opt_main .EQ. 3 ) THEN
  !    ndims = 2 ! Subarray for the output of the Ox correlation
  !    sizes(1)   = imax_total/2;   sizes(2)   = jmax_total/nblock
  !    locsize(1) = imax;           locsize(2) = jmax/nblock
  !    offset(1)  = ims_offset_i;   offset(2)  = ims_offset_j/nblock

  !    CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
  !         MPI_ORDER_FORTRAN, MPI_REAL4, subarray(1), ims_err)
  !    CALL MPI_Type_commit(subarray(1), ims_err)

  !    ndims = 2 ! Subarray for the output of the Oz correlation
  !    sizes(1)   = kmax_total/2;   sizes(2)   = jmax_total/nblock
  !    locsize(1) = kmax;           locsize(2) = jmax/nblock
  !    offset(1)  = ims_offset_k;   offset(2)  = ims_offset_j/nblock

  !    CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
  !         MPI_ORDER_FORTRAN, MPI_REAL4, subarray(2), ims_err)
  !    CALL MPI_Type_commit(subarray(2), ims_err)

! #######################################################################
  ELSE IF ( opt_main .EQ. 4 .OR. opt_main .EQ. 3 ) THEN
     ndims = 2 ! Subarray for the output of the Ox cross-correlation
     sizes(1)   = imax_total;     sizes(2)   = jmax_total/nblock
     locsize(1) = imax;           locsize(2) = jmax/nblock
     offset(1)  = ims_offset_i;   offset(2)  = ims_offset_j/nblock

     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, subarray(1), ims_err)
     CALL MPI_Type_commit(subarray(1), ims_err)

! -----------------------------------------------------------------------
     ndims = 2 ! Subarray for the output of the Oz cross-correlation
     sizes(1)   = kmax_total;     sizes(2)   = jmax_total/nblock
     locsize(1) = kmax;           locsize(2) = jmax/nblock
     offset(1)  = ims_offset_k;   offset(2)  = ims_offset_j/nblock

     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, subarray(2), ims_err)
     CALL MPI_Type_commit(subarray(2), ims_err)

! -----------------------------------------------------------------------
     ndims = 3 ! Subarray for the output of the 2D data     
     sizes(1)   = imax_total;   sizes(2)   = jmax_total/nblock;   sizes(3)   = kmax_total 
     locsize(1) = imax;         locsize(2) = jmax/nblock;         locsize(3) = kmax 
     offset(1)  = ims_offset_i; offset(2)  = ims_offset_j/nblock; offset(3)  = ims_offset_k

     CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
          MPI_ORDER_FORTRAN, MPI_REAL4, subarray(3), ims_err)
     CALL MPI_Type_commit(subarray(3), ims_err)
     
  ENDIF

  RETURN
END SUBROUTINE DEFINE_ARRAY_TYPE
#endif 
