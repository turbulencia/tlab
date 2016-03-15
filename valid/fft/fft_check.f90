#include "dns_error.h" 


SUBROUTINE FFT_CHECK(check_mode, err_count,case_count, &
     trans,  &
     trans_ref, &
     tmp1 ,  &
     tmp2 ,  & 
     tmp3 ,  &
     tmp4 ,  &
     wrk3d,  &
     wrk2d  ) 

#include "types.h"

USE DNS_GLOBAL, ONLY   :  imax_total, jmax_total, kmax_total, imax, jmax, kmax
USE DNS_CONSTANTS, ONLY : lfile
USE DNS_GLOBAL, ONLY   :  isize_txc_dimx, isize_txc_dimz 

#ifdef USE_MPI 
USE DNS_MPI
#endif 

IMPLICIT NONE 

#include "integers.h" 
#ifdef USE_MPI 
#include "mpif.h" 
#include "dns_const_mpi.h" 
#endif 


TINTEGER, INTENT(IN)                                  :: check_mode
TINTEGER, INTENT(INOUT)                               :: err_count,case_count 
TREAL, DIMENSION((imax+2)*(jmax+2)*kmax) :: tmp1,tmp2,tmp3,tmp4,wrk2d,wrk3d
TREAL, DIMENSION(imax*jmax*kmax)         :: trans, trans_ref 

TREAL                                    :: scalex, scaley, dummy,field_variance 
TREAL, DIMENSION(jmax)                   :: spec_variance, dummy_variance 

TINTEGER                            :: i,j,k,iv,ip,ip_ref,bad_count,good_count,bad,nxz 
TINTEGER                            :: good_planes,bad_planes 
TINTEGER                            :: isize_fft3d, isize_wrk3d, isize_trn3d
TINTEGER                            :: isize_page_ij, isize_line_i 
TREAL                               :: norm, residual, re, im, arg  
CHARACTER                           :: fname*32, line*256,check_name*32,label*32

isize_page_ij = imax*jmax
isize_line_i  = imax

isize_wrk3d = imax*jmax*kmax
isize_fft3d = (imax+2)*(jmax+2)*kmax 
isize_trn3d = (imax/2)* jmax   *(2*kmax)


CALL OPR_FOURIER_INITIALIZE(tmp1,tmp2,wrk2d,wrk3d)


norm = C_1_R / M_REAL(imax_total*kmax_total)

tmp1(:)  = C_0_R
tmp3(:)  = C_0_R
tmp4(:)  = C_0_R
wrk2d(:) = C_0_R
wrk3d(:) = C_0_R; tmp2(:) = C_0_R; trans(:) = C_0_R 

DO k=1,kmax
   DO j=1,jmax 
      DO i=1,imax 
         ip = (k-1)*imax*jmax + (j-1)*imax+i 
#ifdef USE_MPI 
         tmp1(ip) = setup_check(check_mode,i+ims_offset_i,   j, k+ims_offset_k)
#else 
         tmp1(ip) = setup_check(check_mode,i,                j, k)
#endif
      ENDDO
   ENDDO
ENDDO

SELECT CASE ( check_mode ) 
CASE (1)
   label = 'fft_check:    delta'
   nxz = imax_total*kmax_total
   dummy          = M_REAL(nxz) / (M_REAL(nxz)-C_1_R)
   field_variance = C_2_R * ( M_REAL(nxz) - C_1_R) / nxz / nxz 
   field_variance = dummy * field_variance * M_REAL(jmax_total) / C_2_R
CASE (2)
   label = 'fft_check:   cosine'
   field_variance = C_05_R * jmax_total
CASE (3) 
   label = 'fft_check:   random' 
CASE DEFAULT
   CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
END SELECT



#ifdef USE_MPI
DO i=0,ims_npro-1 
   IF (ims_pro .EQ. i ) THEN 
      DO k=1,kmax
         WRITE(*,1011) ims_pro, 'input', tmp1((k-1)*imax*jmax+1: k*imax*jmax)
1011     FORMAT (i3,a,16(g10.3,1x))
      ENDDO
   ENDIF
   CALL MPI_Barrier(MPI_COMM_WORLD,ims_err) 
ENDDO 
#endif 

wrk3d(:) = C_0_R
tmp2(:)  = C_0_R
trans(:) = C_0_R
wrk2d(:) = C_0_R
CALL OPR_FOURIER_F(i2, imax,jmax,kmax, tmp1,tmp2, wrk3d, wrk2d, trans)

tmp4(:) = tmp2(:) 
wrk3d(:) = C_0_R
tmp3(:)  = C_0_R 
CALL OPR_FOURIER_B(i2, imax,jmax,kmax, tmp4,tmp3, wrk3d) 

ip = imax*jmax*kmax 
residual = MAXVAL(ABS(norm*tmp3(1:ip)-tmp1(1:ip)))
#ifdef USE_MPI 
dummy = residual 
CALL MPI_Reduce(dummy,residual,1,MPI_REAL8,MPI_MAX,0,MPI_COMM_WORLD,ims_err)  


DO i=0,ims_npro-1 
   IF ( i.EQ.ims_pro ) THEN 
      WRITE(*,*) 'pro', ims_pro 
      DO k=1,kmax
         WRITE(*,1018) norm*tmp3((k-1)*imax*jmax+1:k*imax*jmax)
1018     FORMAT(16(G10.3,1x))
      ENDDO
   ENDIF
   CALL MPI_Barrier(MPI_COMM_WORLD,ims_err) 
END DO

IF (ims_pro .EQ. 0) THEN 
#endif 
   
   IF (  residual  .GT. 1.e-14 ) THEN 
      WRITE(line,1000) TRIM(ADJUSTL(label)), 'FAILED', residual 
      err_count = err_count + 1 
   ELSE 
      WRITE(line,1000) TRIM(ADJUSTL(label)), 'PASSED', residual 
   ENDIF
1000 FORMAT(a, 1x, a6, ' Transform-check.    Max. Residual: ', G13.6)   
   case_count= case_count+1
   CALL IO_WRITE_ASCII(lfile, line) 

#ifdef USE_MPI 
ENDIF
#endif 

IF ( check_mode .EQ. i3 ) GOTO 999

trans(:) = C_0_R
spec_variance(:) = C_0_R

!DO i=0,ims_npro-1 
!   IF ( i .EQ. ims_pro ) THEN 
!      WRITE(*,*) ims_pro, 'after poisson'
!      DO k=1,kmax
!         WRITE(*,*) tmp2((k-1)*2*(imax/2+1)*(jmax+2)+1 : k*2*(imax/2+1)*(jmax+2))
!      ENDDO
!   ENDIF
!   CALL MPI_Barrier(MPI_COMM_WORLD,ims_err) 
!ENDDO 


CALL KXZ_PSD(imax,jmax,kmax,i1,norm,tmp2,trans,spec_variance) 

#ifdef USE_MPI
dummy_variance(:) = spec_variance(:)
CALL MPI_Reduce(dummy_variance,spec_variance,jmax,MPI_REAL8,MPI_SUM,i0,MPI_COMM_WORLD,ims_err) 
#endif 


DO k=1,kmax
   DO j=1,jmax
      DO i=1,imax/2
         ip = (k-1)*imax/2*jmax + (j-1)*imax/2 + i
#ifdef USE_MPI
         trans_ref(ip) = power_check(check_mode,i+ims_offset_i/2,j,k+ims_offset_k)
#else 
         trans_ref(ip) = power_check(check_mode,i,             j,k)
#endif
      ENDDO
   ENDDO
ENDDO


ip = imax/2*jmax*kmax
residual = MAXVAL(ABS(trans_ref(1:ip) - trans(1:ip)))

#ifdef USE_MPI
dummy = residual 
CALL MPI_Reduce(dummy,residual,1,MPI_REAL8,0,MPI_MAX,MPI_COMM_WORLD,ims_err) 

IF ( ims_pro .EQ. 0 ) THEN 
#endif
   IF ( MAXVAL(ABS(trans_ref(1:ip) - trans(1:ip))) .GT. 1.e-16 ) THEN
      WRITE(line,1017) TRIM(ADJUSTL(label)),'FAILED',residual
      err_count = err_count + 1 
   ELSE 
      WRITE(line,1017) TRIM(ADJUSTL(label)),'PASSED',residual
   ENDIF
   case_count = case_count + 1
   
1017 FORMAT(a,1x,a6,1x, 'PSD check.    Max. Residual:', G13.8)
   CALL IO_WRITE_ASCII(lfile,line) 
#ifdef USE_MPI
ENDIF 
#endif 

#ifdef USE_MPI

DO i=1,ims_npro
   IF ( ims_pro .EQ. i-1 ) THEN 
      DO k=1,kmax
!         WRITE(*,*) 'k=',k
!         WRITE(*,1010) ims_pro, trans(    (k-1)*imax/2*jmax+1 : k*imax/2*jmax) 
!         WRITE(*,1010) ims_pro, trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax) 
!         WRITE(*,1010) ims_pro, trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax) - trans((k-1)*imax/2*jmax+1 : k*imax/2*jmax)  
      ENDDO
   ENDIF
   CALL MPI_Barrier(MPI_COMM_WORLD,ims_err) 
ENDDO
#else 
DO k=1,kmax
!   WRITE(*,*) 'k=',k
!   WRITE(*,1010) trans(    (k-1)*imax/2*jmax+1 : k*imax/2*jmax) 
!   WRITE(*,1010) trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax) 
!   WRITE(*,1010) trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax) - trans((k-1)*imax/2*jmax+1 : k*imax/2*jmax)  
ENDDO
#endif  
1010 FORMAT(16(g12.3,1x))


#ifdef USE_MPI
IF ( ims_pro .EQ. 0 ) THEN  
#endif 

   DO j=2,jmax 
      spec_variance(1) = spec_variance(1) + spec_variance(j) 
   ENDDO

   IF ( ABS(field_variance - spec_variance(1)) .LT. 1.e-07 ) THEN 
      WRITE(line,1001) TRIM(ADJUSTL(label)), 'PASSED', &
           ABS(field_variance - spec_variance(1)), field_variance, spec_variance(1)
   ELSE 
      WRITE(line,1001) TRIM(ADJUSTL(label)), 'FAILED', &
           ABS(field_variance - spec_variance(1)), field_variance, spec_variance(1) 
      err_count = err_count + 1
   ENDIF
   case_count = case_count + 1 
   
   CALL IO_WRITE_ASCII(lfile,line) 
1001 FORMAT(a, 1x,a6, 1x, 'PARSEVAL-Identity check. Residual:', &
          1x, G13.8, ' Field:', G13.8, ' Spectrum: ', G13.8)

#ifdef USE_MPI
ENDIF
#endif 

999 CONTINUE 

CONTAINS

  TREAL FUNCTION SETUP_CHECK(check_mode,i,j,k) 

    IMPLICIT NONE 
    
    TINTEGER, INTENT(IN) :: check_mode,i,j,k 

    SELECT CASE ( check_mode )
    CASE (1) ! delta check
       IF ( i.EQ.i1  ) THEN 
          IF ( k.LE.i2 ) THEN 
             SETUP_CHECK= C_1_R * MOD(j,i2)
          ELSE
             SETUP_CHECK = C_0_R 
          ENDIF
       ELSE 
          SETUP_CHECK = C_0_R
       ENDIF
    CASE (2) ! cosine check
       SETUP_CHECK = COS(C_2_R*C_PI_R * M_REAL(i-1) / M_REAL(imax_total))
    CASE (3) ! random check 
       CALL RANDOM_NUMBER(SETUP_CHECK) 
    CASE default
       CALL DNS_STOP(DNS_ERROR_UNDEVELOP) 
    END SELECT
    
  END FUNCTION SETUP_CHECK



  TREAL FUNCTION POWER_CHECK(check_mode,i,j,k) 

    TINTEGER, INTENT(IN) ::  check_mode,i,j,k

    TREAL arg,re,im

    SELECT CASE ( check_mode ) 
    CASE (1)  ! delta-function check 
       IF  ( MOD(j,i2) .EQ. 1 ) THEN
          arg = C_PI_R*C_2_R*M_REAL(k-1)/M_REAL(kmax_total)  
          re = cos(arg)  + cos(  C_2_R*arg) 
          im = sin(-arg) + sin(- C_2_R*arg) 
          power_check = ( re*re+im*im ) / M_REAL(imax_total*imax_total*kmax_total*kmax_total)
       ELSE 
          power_check = C_0_R
       ENDIF
    CASE (2)  ! cosine check 
       IF ( i.EQ.i2 .AND. k.EQ.i1  ) THEN 
          power_check = C_025_R
       ELSE 
          power_check = C_0_R
       ENDIF
    CASE (3)
       power_check = - C_1_R
    END SELECT

  END FUNCTION POWER_CHECK


END SUBROUTINE FFT_CHECK



