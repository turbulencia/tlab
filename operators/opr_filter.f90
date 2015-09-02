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
        
  USE DNS_GLOBAL, ONLY : kmax_total, isize_txc_field
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

! -------------------------------------------------------------------
  TINTEGER nxy

  TREAL, DIMENSION(:), POINTER :: p_a, p_b

! ###################################################################
#ifdef USE_MPI
  IF ( ims_npro_i .GT. 1 ) THEN
     CALL IO_WRITE_ASCII(efile,'OPR_FILTER. Filter in 2D decomposition yet undeveloped.')
     CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  ENDIF
#endif

! -------------------------------------------------------------------
! X direction
! -------------------------------------------------------------------
  IF ( ibc_x(1) .EQ. 1 .AND. nx .GT. 1 ) THEN
     CALL OPR_FILTER_X(itype, nx,ny,nz, ibc_x(2),ibc_x(3),ibc_x(4), u, cx, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF

! -------------------------------------------------------------------
! Y direction
! -------------------------------------------------------------------
  IF ( ibc_y(1) .EQ. 1 .AND. ny .GT. 1 ) THEN
     CALL OPR_FILTER_Y(itype, nx,ny,nz, ibc_y(2),ibc_y(3),ibc_y(4), u, cy, txc(1,2), wrk1d,wrk2d,txc(1,1))
  ENDIF

! -------------------------------------------------------------------
! Z direction
! -------------------------------------------------------------------
  IF ( ibc_z(1) .EQ. 1 .AND. kmax_total .GT. 1 ) THEN

#ifdef USE_MPI
     IF ( ims_npro_k .GT. 1 ) THEN
! Transpose Matrix u -> txc
        CALL DNS_MPI_TRPF_K(u, txc(1,1), ims_ds_k(1,impi_id), ims_dr_k(1,impi_id), &
             ims_ts_k(1,impi_id), ims_tr_k(1,impi_id))
        p_a => txc(:,1)
        p_b => u
        nxy = ims_size_k(impi_id)
     ELSE
#endif
        p_a => u
        p_b => txc(:,1)
        nxy = nx*ny
#ifdef USE_MPI
     ENDIF
#endif

! Filter            
     CALL OPR_FILTER_Z(itype, nxy,kmax_total, ibc_z(2),ibc_z(3),ibc_z(4), p_a, cz, txc(1,2), wrk1d,wrk2d,p_b)

! Transpose back Matrix txc -> u; we know that u=txc(1,1), check OPR_FILTER_Z
#ifdef USE_MPI
     IF ( ims_npro_k .GT. 1 ) THEN
        CALL DNS_MPI_TRPB_K(txc(1,1), u, ims_ds_k(1,impi_id), ims_dr_k(1,impi_id), &
             ims_ts_k(1,impi_id), ims_tr_k(1,impi_id))

     ENDIF
#endif
  NULLIFY(p_a,p_b)

  ENDIF

  RETURN
END SUBROUTINE OPR_FILTER
