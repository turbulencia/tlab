#include "types.h"

FUNCTION SIMPSON(imax, u, dx)

  IMPLICIT NONE

  TINTEGER,               INTENT(IN) :: imax
  TREAL, DIMENSION(imax), INTENT(IN) :: u, dx
  
  TINTEGER i,nn
  TREAL SIMPSON, slast
  TREAL dx21, dx20, dx10, du20, du10, du21, b, c
  TINTEGER i2

  i2 = 2

  IF ( MOD(imax,i2) .EQ. 0 ) THEN 
     dx21 = dx(imax)
     dx20 = dx(imax)+dx(imax-1)
     dx10 = dx(imax-1)
     du20 = u(imax)-u(imax-2)
     du10 = u(imax-1)-u(imax-2)
     du21 = u(imax)-u(imax-1)
     c = (du21/dx21-du10/dx10)/dx20
     b = (du21/dx21-c*dx21)/C_2_R
     slast = dx21*(u(imax-1)+dx21*(b+c*dx21/C_3_R))
     nn = imax - 1
  ELSE
     slast = C_0_R
     nn = imax
  ENDIF

  SIMPSON = u(1)*dx(1) + u(nn)*dx(nn)

  DO i=2, nn-1,2
     SIMPSON = SIMPSON + C_4_R * u(i) * dx(i)
  ENDDO

  DO i=3, nn-2,2
     SIMPSON = SIMPSON + C_2_R * u(i) * dx(i)
  ENDDO

  SIMPSON = SIMPSON/C_3_R + slast

  RETURN
END FUNCTION SIMPSON

! ###################################################################
! ###################################################################
FUNCTION SIMPSON_NU(imax, u, x)

  IMPLICIT NONE

  TINTEGER,               INTENT(IN) :: imax
  TREAL, DIMENSION(imax), INTENT(IN) :: u, x

  TINTEGER i, nn, i2
  TREAL SIMPSON_NU
  TREAL dx21, dx20, dx10, du20, du10, du21, b, c

  i2 = 2 

! Correct the last element contribution
  IF ( MOD(imax,i2) .EQ. 0 ) THEN 
     dx21 = x(imax)-x(imax-1)
     dx20 = x(imax)-x(imax-2)
     dx10 = x(imax-1)-x(imax-2)
     du20 = u(imax)-u(imax-2)
     du10 = u(imax-1)-u(imax-2)
     du21 = u(imax)-u(imax-1)
     c = (du21/dx21-du10/dx10)/dx20
     b = (du21/dx21-c*dx21)/C_2_R
     SIMPSON_NU = dx21*(u(imax-1)+dx21*(b+c*dx21/C_3_R))
     nn = imax - 1
  ELSE
     SIMPSON_NU = C_0_R
     nn = imax
  ENDIF

  DO i=2, nn-1,2
     dx21 = x(i+1)-x(i)
     dx20 = x(i+1)-x(i-1)
     dx10 = x(i)-x(i-1)
     du20 = u(i+1)-u(i-1)
     du10 = u(i)-u(i-1)
     c = (du20/dx20-du10/dx10)/dx21
     b = (du20/dx20-c*dx20)/C_2_R
     SIMPSON_NU = SIMPSON_NU + dx20*(u(i-1)+dx20*(b+c*dx20/C_3_R))
  ENDDO

  RETURN
END FUNCTION SIMPSON_NU
