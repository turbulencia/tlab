!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

SUBROUTINE BOUNDARY_INFLOW_FILTER(u, txc, mean, wrk1d,wrk2d,wrk3d)
  
  USE DNS_GLOBAL
  USE DNS_LOCAL, ONLY : ifilt_inflow_iwidth, ifilt_inflow_jwidth

  IMPLICIT NONE
  
#include "integers.h"
  
  TREAL u(imax,jmax,kmax)
  TREAL mean(jmax,kmax)
  TREAL, DIMENSION(*) :: txc, wrk1d,wrk2d,wrk3d

! -----------------------------------------------------------------------
  TINTEGER i,j,k,ip
  TINTEGER j1, imx, jmx, ifltmx, jfltmx

  TINTEGER id, ibc_x(4), ibc_y(4), ibc_z(4)

! #######################################################################
  imx = ifilt_inflow_iwidth
  j1  = (jmax-ifilt_inflow_jwidth)/2+1
  jmx = (jmax+ifilt_inflow_jwidth)/2
  j1  = MIN(MAX(j1,i1),jmax)
  jmx = MIN(MAX(jmx,i1),jmax)

! BCs for the filters (see routine FILTER)
  ibc_x(1) = 1; ibc_x(2) = i1bc; ibc_x(3) = 1; ibc_x(4) = 1
  ibc_y(1) = 1; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0 
  ibc_z(1) = 1; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0 

  IF ( imx-i1 .LT. 6 ) THEN
     ibc_x(1) = 0
!     CALL IO_WRITE_ASCII(efile, 'Error: Explicit Filter width too small')
!     CALL DNS_STOP(DNS_ERROR_FLTWIDTH)
  ENDIF

  IF ( jmx-j1 .LT. 6 ) THEN
     ibc_y(1) = 0
!     CALL IO_WRITE_ASCII(efile, 'Error: Explicit Filter width too small')
!     CALL DNS_STOP(DNS_ERROR_FLTWIDTH)
  ENDIF

! -----------------------------------------------------------------------
! Remove mean field
! -----------------------------------------------------------------------
  IF ( ibc_x(1) .EQ. 1 .OR. ibc_y(1) .EQ. 1 .OR. ibc_z(1) .EQ. 1 ) THEN
     ip = 1
     DO k = 1,kmax
        DO j = j1,jmx
           DO i = i1,imx
              txc(ip) = u(i,j,k) - mean(j,k)
              ip = ip + 1
           ENDDO
        ENDDO
     ENDDO
     ifltmx = imx-i1+1
     jfltmx = jmx-j1+1
  ENDIF

! -----------------------------------------------------------------------
  id = DNS_MPI_K_INFLOW
  ! CALL OPR_FILTER(i2, ifltmx,jfltmx,kmax, ibc_x,ibc_y,ibc_z, id, &
  !      txc, wrk3d,wrk3d,wrk3d, wrk1d,wrk2d,wrk3d)
! Needs to be updated
  
! -----------------------------------------------------------------------
! Add mean field
! -----------------------------------------------------------------------
  IF ( ibc_x(1) .EQ. 1 .OR. ibc_y(1) .EQ. 1 .OR. ibc_z(1) .EQ. 1 ) THEN
     ip = 1
     DO k = 1,kmax
        DO j = j1,jmx
           DO i = i1,imx
              u(i,j,k) = txc(ip) + mean(j,k)
              ip = ip + 1
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  RETURN
END SUBROUTINE BOUNDARY_INFLOW_FILTER
