#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/11/12 - J.P. Mellado
!#              Extracted to be used in buffer 
!# 2008/11/26 - J.P. Mellado
!#              From FILT6E to use with different filter types
!#
!########################################################################
!# DESCRIPTION
!#
!# This subroutine applies a filter 
!#
!########################################################################
!# ARGUMENTS 
!#
!# impi_id    In    MPI type for the transposition
!# itype      In    Flag indicating type of filter:
!#                  1 Compact 4th order
!#                  2 Explicit 6th order
!#                  3 Explicit 4th order
!#                  4 Explicit ADM
!# txc        In    Auxilar 3D array of size:
!#                  2 if ADM type
!#                  1 otherwise
!# ibc        In    Data about BCs for each direction
!#                  1 Flag active/inactive
!#                  2 Periodic/Nonperiodic
!#                  3 Specific condition on imin
!#                  4 Specific condition on imax
!#
!########################################################################
SUBROUTINE OPR_FILTER(itype, nx,ny,nz, ibc_x,ibc_y,ibc_z, impi_id, u, cx,cy,cz, wrk1d,wrk2d,txc)
        
  USE DNS_GLOBAL, ONLY : imax_total, kmax_total, isize_txc_field
  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

  TINTEGER itype, nx,ny,nz, impi_id
  TINTEGER, DIMENSION(*)                         :: ibc_x, ibc_y, ibc_z
  TREAL,    DIMENSION(nx*ny*nz),          TARGET :: u
  TREAL,    DIMENSION(*)                         :: cx,cy,cz, wrk1d, wrk2d
  TREAL,    DIMENSION(isize_txc_field,*), TARGET :: txc

! ###################################################################
  IF ( ibc_x(1) .EQ. 1 .AND. imax_total .GT. 1 ) THEN
     CALL OPR_FILTER_X(itype, nx,ny,nz, ibc_x(2),ibc_x(3),ibc_x(4), u, cx, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF

  IF ( ibc_y(1) .EQ. 1 .AND. ny .GT. 1 ) THEN
     CALL OPR_FILTER_Y(itype, nx,ny,nz, ibc_y(2),ibc_y(3),ibc_y(4), u, cy, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF

  IF ( ibc_z(1) .EQ. 1 .AND. kmax_total .GT. 1 ) THEN
     CALL OPR_FILTER_Z(itype, nx,ny,nz, ibc_z(2),ibc_z(3),ibc_z(4), u, cz, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF

  RETURN
END SUBROUTINE OPR_FILTER
