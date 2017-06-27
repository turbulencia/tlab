#include "types.h"
#include "dns_error.h"
  
!########################################################################
!# DESCRIPTION
!# 
!# Transforming a random velocity field to have a given covariance matrix
!#
!########################################################################
SUBROUTINE RAND_COVARIANCE(imax,jmax,kmax, u,v,w, cov)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL, ONLY    : g
  
  IMPLICIT NONE
  
  TINTEGER imax, jmax, kmax
  TREAL, DIMENSION(imax,jmax,kmax) :: u, v, w
  TREAL cov(6)
  
! -------------------------------------------------------------------
  TINTEGER i, j, k
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
     
! only 2D case developed
     IF ( Rxz .NE. C_0_R .OR. Ryz .NE. C_0_R ) THEN
        CALL IO_WRITE_ASCII(efile,'Terms Rxz and Ryz not developed yet.')
        CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF
     
     CALL ISORMS(imax, jmax, kmax, Rzz, w)
     
  ENDIF
      
! -------------------------------------------------------------------
! Diagonal case
! -------------------------------------------------------------------
  IF ( Rxy .EQ. C_0_R ) THEN
     CALL ISORMS(imax, jmax, kmax, Rxx, u)
     CALL ISORMS(imax, jmax, kmax, Ryy, v)
     
  ELSE
! -------------------------------------------------------------------
! Nondiagonal case
! -------------------------------------------------------------------
! get eigenvalues
     trace = Rxx+Ryy
     lambda1 = C_05_R*(trace + SQRT(trace*trace-C_4_R*(Rxx*Ryy-Rxy*Rxy)))
     lambda2 = trace - lambda1
     
! define fields in rotated uncorrelated frame
     CALL ISORMS(imax, jmax, kmax, lambda1, u)
     CALL ISORMS(imax, jmax, kmax, lambda2, v)
     
! rotate to XY correlated frame
     alpha  = ATAN((lambda1-Rxx)/Rxy)
     calpha = COS(alpha)
     salpha = SIN(alpha)
     
     DO k = 1,kmax
        DO j = 1,jmax
           DO i = 1,imax
              rdummy   = calpha*u(i,j,k) - salpha*v(i,j,k)
              v(i,j,k) = salpha*u(i,j,k) + calpha*v(i,j,k)
              u(i,j,k) = rdummy
           ENDDO
        ENDDO
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

!########################################################################
!########################################################################
SUBROUTINE ISORMS(imax,jmax,kmax, var1, u)

  IMPLICIT NONE

  TINTEGER imax, jmax, kmax
  TREAL var1
  TREAL u(imax,jmax,kmax)

! -------------------------------------------------------------------
  TINTEGER i1, i2 
  TREAL AVG1V3D, var0, factor

! ###################################################################
  i1 = 1
  i2 = 2

  var0 = AVG1V3D(imax,jmax,kmax, i1, u)
  u = u - var0

  var0 = AVG1V3D(imax,jmax,kmax, i2, u)
  IF ( var0 .GT. C_0_R ) THEN
     factor = SQRT(var1/var0)
     u = u*factor
  ENDIF
  
  RETURN
END SUBROUTINE ISORMS
