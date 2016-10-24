#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2012/07/20 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculation of FFT using the OPR_FOURIER_ routines
!#
!########################################################################
SUBROUTINE OPR_FOURIER_F(flag_mode, nx,ny,nz, in,out, tmp1,wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : jmax_total, kmax_total
  USE DNS_GLOBAL, ONLY : isize_txc_field, isize_txc_dimz, isize_wrk2d
  USE DNS_GLOBAL, ONLY : fft_plan_fy

  IMPLICIT NONE
  
  TINTEGER                                           :: flag_mode  ! 1D, 2D or 3D
  TINTEGER                                           :: nx,ny,nz
  TREAL, DIMENSION(isize_txc_field),   INTENT(INOUT) :: in         ! extended w/ BCs below
  TREAL, DIMENSION(isize_txc_dimz,nz), INTENT(OUT)   :: out
  TREAL, DIMENSION(isize_wrk2d,2),     INTENT(INOUT) :: wrk2d      ! BCs padding
  TREAL, DIMENSION(isize_txc_dimz,nz), INTENT(INOUT) :: tmp1, wrk3d

! -----------------------------------------------------------------------
  TINTEGER k

! #######################################################################
  wrk2d = C_0_R
  
  IF ( flag_mode .EQ. 3 .AND. jmax_total .GT. 1 ) THEN ! 3D FFT (unless 2D sim)
     IF ( kmax_total .GT. 1 ) THEN
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in,wrk2d(1,1),wrk2d(1,2),out, tmp1,wrk3d) 
        CALL OPR_FOURIER_F_Z_EXEC(out,wrk3d)
     ELSE  
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in,wrk2d(1,1),wrk2d(1,2),wrk3d, out,tmp1) 
     ENDIF

     DO k = 1,nz
        CALL dfftw_execute_dft(fft_plan_fy, wrk3d(1,k), out(1,k))
     ENDDO

  ELSE
     IF ( flag_mode .EQ. 2 .AND. kmax_total .GT. 1 ) THEN ! 2D FFT (unless 1D sim)
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in,wrk2d(1,1),wrk2d(1,2),wrk3d, out,tmp1) 
        CALL OPR_FOURIER_F_Z_EXEC(wrk3d,out)
     ELSE  
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in,wrk2d(1,1),wrk2d(1,2),out, tmp1,wrk3d) 
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE OPR_FOURIER_F

! #######################################################################
! #######################################################################
SUBROUTINE OPR_FOURIER_B(flag_mode, nx,ny,nz, in,out, wrk3d)

  USE DNS_GLOBAL, ONLY : jmax_total, kmax_total
  USE DNS_GLOBAL, ONLY : isize_txc_field, isize_txc_dimz
  USE DNS_GLOBAL, ONLY : fft_plan_by

  IMPLICIT NONE
  
  TINTEGER                                           :: flag_mode  ! 1D, 2D or 3D
  TINTEGER                                           :: nx,ny,nz
  TREAL, DIMENSION(isize_txc_dimz,nz), INTENT(INOUT) :: in, wrk3d
  TREAL, DIMENSION(isize_txc_field),   INTENT(OUT)   :: out

! -----------------------------------------------------------------------
  TINTEGER k

! #######################################################################
  IF ( flag_mode .EQ. 3 .AND. jmax_total .GT. 1 ) THEN ! 3D FFT (unless 2D sim)
     DO k = 1,nz
        CALL dfftw_execute_dft(fft_plan_by, in(1,k), wrk3d(1,k))
     ENDDO

     IF ( kmax_total .GT. 1 ) THEN
        CALL OPR_FOURIER_B_Z_EXEC(wrk3d,in) 
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, in,out, wrk3d) 
     ELSE
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,out, in) 
     ENDIF

  ELSE
     IF ( flag_mode .EQ. 2 .AND. kmax_total .GT. 1 ) THEN ! 2D FFT (unless 1D sim)
        CALL OPR_FOURIER_B_Z_EXEC(in,wrk3d) 
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,out, in) 
     ELSE
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, in,out, wrk3d) 
     ENDIF

  ENDIF

  RETURN
END SUBROUTINE OPR_FOURIER_B

! #######################################################################
! #######################################################################
! Calculate spectrum in array in1
! Calculate correlation in array in2
SUBROUTINE OPR_FOURIER_CONVOLUTION_FXZ(flag1,flag2, nx,ny,nz, in1,in2, tmp1,tmp2, wrk2d,wrk3d)

  USE DNS_GLOBAL, ONLY : kmax_total, fft_reordering
  USE DNS_GLOBAL, ONLY : isize_txc_field, isize_wrk2d

  IMPLICIT NONE
  
#include "integers.h"

  TINTEGER                                              :: flag1,flag2, nx,ny,nz
  TCOMPLEX, DIMENSION(isize_txc_field/2), INTENT(INOUT) :: in1,in2, tmp1,tmp2, wrk3d
  TREAL, DIMENSION(isize_wrk2d,2),        INTENT(INOUT) :: wrk2d ! BCs padding

! -----------------------------------------------------------------------

! #######################################################################
  wrk2d = C_0_R
  fft_reordering = i1
 
  IF ( kmax_total .GT. 1 ) THEN
     CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in1,wrk2d(1,1),wrk2d(1,2),wrk3d, tmp1,tmp2) 
     CALL OPR_FOURIER_F_Z_EXEC(wrk3d,tmp1)
  ELSE  
     CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in1,wrk2d(1,1),wrk2d(1,2),tmp1, tmp2,wrk3d) 
  ENDIF

  IF      ( flag1 .EQ. 1 ) THEN    ! Auto-spectra
     in1 = tmp1*CONJG(tmp1)

  ELSE IF ( flag1 .EQ. 2 ) THEN    ! Cross-spectra
     IF ( kmax_total .GT. 1 ) THEN
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in2,wrk2d(1,1),wrk2d(1,2),wrk3d, in1,tmp2) 
        CALL OPR_FOURIER_F_Z_EXEC(wrk3d,in1)
     ELSE  
        CALL OPR_FOURIER_F_X_EXEC(nx,ny,nz, in2,wrk2d(1,1),wrk2d(1,2),in1, tmp2,wrk3d) 
     ENDIF
     in1 = in1*CONJG(tmp1)

  ENDIF

! -----------------------------------------------------------------------
  IF ( flag2 .EQ. 2 ) THEN         ! Calculate correlation in array in2
     tmp2 = in1                    ! the routines below can overwrite the entry array
     IF ( kmax_total .GT. 1 ) THEN
        CALL OPR_FOURIER_B_Z_EXEC(tmp2,wrk3d) 
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, wrk3d,in2, tmp1) 
     ELSE
        CALL OPR_FOURIER_B_X_EXEC(nx,ny,nz, tmp2,in2, wrk3d) 
     ENDIF
  ENDIF

! -----------------------------------------------------------------------
  fft_reordering = i0

  RETURN
END SUBROUTINE OPR_FOURIER_CONVOLUTION_FXZ
