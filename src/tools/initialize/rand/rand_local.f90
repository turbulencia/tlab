#include "types.h"
#include "dns_error.h"

MODULE RAND_LOCAL

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_VARS, ONLY : imax,jmax,kmax, isize_field, isize_txc_field
  USE TLAB_VARS, ONLY : g
  USE TLAB_PROCS
  
  IMPLICIT NONE
  SAVE

  ! -------------------------------------------------------------------
  TINTEGER :: ispectrum
  TREAL    :: spc_param(5)  ! Fundamental frequency, fmin, fmax, sigma

  TINTEGER :: ipdf
  TREAL    :: ucov(6)

  TINTEGER :: seed          ! Random number generator

  ! -------------------------------------------------------------------
  TINTEGER i

#include "integers.h"

CONTAINS

  ! ###################################################################
  SUBROUTINE RAND_FIELD(variance, a, tmp1,tmp2,tmp3, wrk2d,wrk3d)
    IMPLICIT NONE

    TREAL,                              INTENT(IN)    :: variance
    TREAL, DIMENSION(isize_field),      INTENT(OUT)   :: a
    TREAL, DIMENSION(isize_txc_field),  INTENT(INOUT) :: tmp1, tmp2, tmp3
    TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk2d,wrk3d

    TINTEGER idim
    TREAL RAN0, RANG
    EXTERNAL RAN0, RANG

    ! -------------------------------------------------------------------
    SELECT CASE( ipdf )
    CASE( 1 )     ! Uniform distribution
      DO i = 1,isize_field
        tmp2(i) = RAN0(seed) -C_05_R
      ENDDO

    CASE( 2 )     ! Gaussian distribution
      DO i = 1,isize_field
        tmp2(i) = RANG(C_0_R, C_1_R, seed)
      ENDDO

    END SELECT

    IF ( ispectrum .GT. 0 ) THEN
      IF ( g(2)%size .EQ. 1 ) THEN; idim = 2;           ! 2D Fourier transform
      ELSE;                         idim = 3; ENDIF     ! 3D Fourier transform

      IF ( ipdf .GT. 0 ) CALL OPR_FOURIER_F(idim, imax,jmax,kmax, tmp2,tmp1, tmp3,wrk2d,wrk3d)
      CALL RAND_PSD(imax,jmax,kmax, tmp1)
      CALL OPR_FOURIER_B(idim, imax,jmax,kmax, tmp1, tmp2, wrk3d)

    ENDIF

    CALL RAND_NORMALIZE( variance, tmp2 )
    a(1:isize_field) = tmp2(1:isize_field)

    RETURN
  END SUBROUTINE RAND_FIELD

  !########################################################################
  SUBROUTINE RAND_COVARIANCE(cov, u,v,w)
    IMPLICIT NONE

    TREAL cov(6)
    TREAL, DIMENSION(isize_field), INTENT(OUT)   :: u,v,w

    ! -------------------------------------------------------------------
    TREAL trace, lambda1, lambda2, alpha, calpha, salpha
    TREAL rdummy

#define Rxx cov(1)
#define Ryy cov(2)
#define Rzz cov(3)
#define Rxy cov(4)
#define Rxz cov(5)
#define Ryz cov(6)

    ! ###################################################################
    IF ( g(3)%size .GT. 1 ) THEN
      IF ( Rxz .NE. C_0_R .OR. Ryz .NE. C_0_R ) THEN ! only 2D case developed
        CALL TLAB_WRITE_ASCII(efile,'Terms Rxz and Ryz not developed yet.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
      ENDIF

      CALL RAND_NORMALIZE(Rzz, w)

    ENDIF

    IF ( Rxy .EQ. C_0_R ) THEN  ! Diagonal case
      CALL RAND_NORMALIZE(Rxx, u)
      CALL RAND_NORMALIZE(Ryy, v)

    ELSE                        ! Nondiagonal case
      ! get eigenvalues
      trace = Rxx+Ryy
      lambda1 = C_05_R*(trace + SQRT(trace*trace-C_4_R*(Rxx*Ryy-Rxy*Rxy)))
      lambda2 = trace - lambda1

      ! define fields in rotated uncorrelated frame
      CALL RAND_NORMALIZE(lambda1, u)
      CALL RAND_NORMALIZE(lambda2, v)

      ! rotate to XY correlated frame
      alpha  = ATAN((lambda1-Rxx)/Rxy)
      calpha = COS(alpha)
      salpha = SIN(alpha)

      DO i =1,isize_field
        rdummy = calpha*u(i) - salpha*v(i)
        v(i)   = salpha*u(i) + calpha*v(i)
        u(i)   = rdummy
      ENDDO

    ENDIF

    RETURN

#undef Rxx
#undef Ryy
#undef Rzz
#undef Rxy
#undef Rxz
#undef Ryz

  END SUBROUTINE RAND_COVARIANCE

  ! ###################################################################
  SUBROUTINE RAND_NORMALIZE(variance, a)
    IMPLICIT NONE

    TREAL,                            INTENT(IN)    :: variance
    TREAL, DIMENSION(imax,jmax,kmax), INTENT(INOUT) :: a

    ! -------------------------------------------------------------------
    TREAL AVG1V2D, dummy
    EXTERNAL AVG1V2D

    ! ###################################################################
    dummy = AVG1V2D(imax*jmax,i1,kmax, i1, i1, a) ! 3D average
    a = a -dummy

    dummy = AVG1V2D(imax*jmax,i1,kmax, i1, i2, a) ! 3D average
    IF ( dummy .GT. C_0_R ) THEN
      dummy = SQRT(variance/dummy)
      a  = a *dummy
    ENDIF

    RETURN
  END SUBROUTINE RAND_NORMALIZE

END MODULE RAND_LOCAL
