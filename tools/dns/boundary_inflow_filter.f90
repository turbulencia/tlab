#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

SUBROUTINE BOUNDARY_INFLOW_FILTER(bcs_vi, q,s, txc, wrk1d,wrk2d,wrk3d)
  
  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_GLOBAL,    ONLY : imax,jmax,kmax, inb_flow,inb_scal, imode_eqns
  USE DNS_LOCAL,     ONLY : FilterInflow

  IMPLICIT NONE
  
#include "integers.h"
  
  TREAL, DIMENSION(imax,jmax,kmax,*), INTENT(INOUT) :: q,s
  TREAL, DIMENSION(jmax,kmax,*),      INTENT(IN)    :: bcs_vi
  TREAL, DIMENSION(imax*jmax*kmax,2), INTENT(INOUT) :: txc
  TREAL, DIMENSION(*),                INTENT(INOUT) :: wrk1d,wrk2d,wrk3d

! -----------------------------------------------------------------------
  TINTEGER i,j,k,ip, iq, iq_loc(inb_flow), is
  TINTEGER j1, imx, jmx, ifltmx, jfltmx

! #######################################################################
  CALL IO_WRITE_ASCII(efile,'BOUNDARY_BUFFER_FILTER. Needs to be updated to new filter routines.')
  ! FilterInflow needs to be initiliazed
  CALL DNS_STOP(DNS_ERROR_UNDEVELOP)

  imx = FilterInflow(1)%size
  j1  = ( jmax -FilterInflow(2)%size )/2 +1
  jmx = ( jmax +FilterInflow(2)%size )/2
  j1  = MIN(MAX(j1,i1),jmax)
  jmx = MIN(MAX(jmx,i1),jmax)

  ifltmx = imx-i1+1
  jfltmx = jmx-j1+1

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

  IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
     iq_loc = (/ 5,1,2,3,6 /) ! Filtered variables: rho, u,v,w, p
  ELSE
     iq_loc = (/ 1,2,3 /)
  ENDIF
  
! #######################################################################
  DO iq = 1,inb_flow
     
! -----------------------------------------------------------------------
! Remove mean field
! -----------------------------------------------------------------------
     ip = 1
     DO k = 1,kmax
        DO j = j1,jmx
           DO i = i1,imx
              wrk3d(ip) = q(i,j,k,iq_loc(iq)) - bcs_vi(j,k,iq_loc(iq))
              ip = ip + 1
           ENDDO
        ENDDO
     ENDDO
     
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
              wrk3d(ip) = q(i,j,k,iq_loc(iq)) + bcs_vi(j,k,iq_loc(iq))
              ip = ip + 1
           ENDDO
        ENDDO
     ENDDO
     
  ENDDO

! #######################################################################
  DO is = 1,inb_scal
     
! -----------------------------------------------------------------------
! Remove mean field
! -----------------------------------------------------------------------
     ip = 1
     DO k = 1,kmax
        DO j = j1,jmx
           DO i = i1,imx
              wrk3d(ip) = s(i,j,k,is) - bcs_vi(j,k,is+inb_flow)
              ip = ip + 1
           ENDDO
        ENDDO
     ENDDO
     
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
              wrk3d(ip) = s(i,j,k,is) + bcs_vi(j,k,is+inb_flow)
              ip = ip + 1
           ENDDO
        ENDDO
     ENDDO
     
  ENDDO

  RETURN
END SUBROUTINE BOUNDARY_INFLOW_FILTER
