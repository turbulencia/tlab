#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Needs to be updated to the the new filter formulation
!#
!########################################################################
SUBROUTINE BOUNDARY_INFLOW_FILTER(u, txc, mean, wrk1d,wrk2d,wrk3d)
  
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax
  USE DNS_LOCAL,     ONLY : FilterInflow, ifilt_inflow_iwidth, ifilt_inflow_jwidth

  IMPLICIT NONE
  
#include "integers.h"
  
  TREAL, DIMENSION(imax,jmax,kmax),   INTENT(INOUT) :: u ! Filter in place
  TREAL, DIMENSION(     jmax,kmax),   INTENT(IN)    :: mean
  TREAL, DIMENSION(imax*jmax*kmax,2), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! -----------------------------------------------------------------------
  TINTEGER i,j,k,ip
  TINTEGER j1, imx, jmx, ifltmx, jfltmx

! #######################################################################
  CALL IO_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
  CALL DNS_STOP(DNS_ERROR_UNDEVELOP)

  imx = ifilt_inflow_iwidth
  j1  = (jmax-ifilt_inflow_jwidth)/2+1
  jmx = (jmax+ifilt_inflow_jwidth)/2
  j1  = MIN(MAX(j1,i1),jmax)
  jmx = MIN(MAX(jmx,i1),jmax)

! BCs for the filters (see routine FILTER)
  ! ibc_x(1) = 1; ibc_x(2) = i1bc; ibc_x(3) = 1; ibc_x(4) = 1
  ! ibc_y(1) = 1; ibc_y(2) = j1bc; ibc_y(3) = 0; ibc_y(4) = 0 
  ! ibc_z(1) = 1; ibc_z(2) = k1bc; ibc_z(3) = 0; ibc_z(4) = 0 

  ! IF ( imx-i1 .LT. 6 ) THEN ! Turn if off
  !    ibc_x(1) = 0
  ! ENDIF

  ! IF ( jmx-j1 .LT. 6 ) THEN ! Turn it off
  !    ibc_y(1) = 0
  ! ENDIF

! -----------------------------------------------------------------------
! Remove mean field
! -----------------------------------------------------------------------
  ip = 1
  DO k = 1,kmax
     DO j = j1,jmx
        DO i = i1,imx
           wrk3d(ip) = u(i,j,k) - mean(j,k)
           ip = ip + 1
        ENDDO
     ENDDO
  ENDDO
  ifltmx = imx-i1+1
  jfltmx = jmx-j1+1

! -----------------------------------------------------------------------
!  id = DNS_MPI_K_INFLOW
! Need to pass a type flag for the transposition
  CALL OPR_FILTER(ifltmx,jfltmx,kmax, FilterInflow, wrk3d, wrk1d,wrk2d,txc)
  
! -----------------------------------------------------------------------
! Add mean field
! -----------------------------------------------------------------------
  ip = 1
  DO k = 1,kmax
     DO j = j1,jmx
        DO i = i1,imx
           u(i,j,k) = wrk3d(ip) + mean(j,k)
           ip = ip + 1
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE BOUNDARY_INFLOW_FILTER
