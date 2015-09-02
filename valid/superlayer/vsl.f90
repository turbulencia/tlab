! This program just generates a given scalar field in order to check the
! routines of superlayer analysis
PROGRAM VSL

  USE DNS_GLOBAL

  IMPLICIT NONE

#include "types.h"
#include "integers.h"

  TREAL x(:)
  ALLOCATABLE x
  TREAL y(:)
  ALLOCATABLE y
  TREAL z(:)
  ALLOCATABLE z
  TREAL z1(:,:,:)
  ALLOCATABLE z1
  TREAL wrk3d(:)
  ALLOCATABLE wrk3d

  TINTEGER iwrk_size, i, j, k
  TREAL mean, delta, thick, y_center, param, FLOW_SHEAR_TEMPORAL

! ###################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL('dns.ini')

  ALLOCATE(x(imax))
  ALLOCATE(y(jmax))
  ALLOCATE(z(kmax_total))
  ALLOCATE(z1(imax,jmax,kmax))
  ALLOCATE(wrk3d(imax*jmax*kmax))

  CALL IO_READ_GRID(gfile, imax,jmax,kmax_total, scalex,scaley,scalez, x,y,z)

! -------------------------------------------------------------------
! define the field
! -------------------------------------------------------------------
  delta = C_1_R
  mean  = C_05_R
  thick = C_1EM2_R
  DO k = 1,kmax_total
     DO j = 1,jmax
        DO i = 1,imax
           y_center = C_05_R + C_025_R*SIN(2*C_PI_R*x(i))*SIN(2*C_PI_R*z(k))
           z1(i,j,k) = FLOW_SHEAR_TEMPORAL(i3, thick, delta, mean, y_center, param, y(j))
        ENDDO
     ENDDO
  ENDDO

  iwrk_size = imax*jmax*kmax
  CALL DNS_WRITE_FIELDS('scalar_field', i1, imax,jmax,kmax, i1, iwrk_size, z1, wrk3d)

  STOP
END PROGRAM VSL
