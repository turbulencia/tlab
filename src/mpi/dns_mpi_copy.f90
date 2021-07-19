#include "types.h"
#include "dns_error.h"

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
!# Moving plane information between adjacent PEs, circulant version.
!# npl is smaller than 2*kmax
!# The number of plane to move is given by npl
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE DNS_MPI_COPYPLN_1(ijmax, kmax, npl, a, bl, br)

  USE DNS_CONSTANTS, ONLY : efile
  USE TLAB_CORE
  USE DNS_MPI

  IMPLICIT NONE

#include "mpif.h"

  TINTEGER ijmax, kmax, npl
  TREAL a(ijmax,*)
  TREAL bl(ijmax,*)
  TREAL br(ijmax,*)

! -----------------------------------------------------------------------
  INTEGER status(MPI_STATUS_SIZE,4)
  INTEGER mpireq(4)
  INTEGER ims_pro_l, ims_pro_r
  INTEGER icount

! #######################################################################
  IF ( ims_npro .GT. 1 ) THEN

! Careful in case only 2 PEs
     IF ( ims_npro .EQ. 2 ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_COPYPLN_1. Undeveloped for 2 PEs.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

! left and right PEs
     ims_pro_l = MOD(ims_pro-1+ims_npro,ims_npro)
     ims_pro_r = MOD(ims_pro+1+ims_npro,ims_npro)

     icount   = ijmax*npl

     CALL MPI_IRECV(bl(1,1), icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
     CALL MPI_IRECV(br(1,1), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

     CALL MPI_ISEND(a(1,1),          icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
     CALL MPI_ISEND(a(1,kmax+1-npl), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

     CALL MPI_WAITALL(4, mpireq, status, ims_err)

     CALL DNS_MPI_TAGUPDT

  ENDIF

  RETURN
END SUBROUTINE DNS_MPI_COPYPLN_1

!########################################################################
!########################################################################
SUBROUTINE DNS_MPI_COPYPLN_2(ijmax, kmax, npl, a, bl, br)

  USE DNS_CONSTANTS, ONLY : efile
  USE DNS_MPI

  IMPLICIT NONE

#include "mpif.h"

  TINTEGER ijmax, kmax, npl
  TREAL a(ijmax,kmax)
  TREAL bl(ijmax,npl)
  TREAL br(ijmax,npl)

! -----------------------------------------------------------------------
  TINTEGER npl_loc
  INTEGER status(MPI_STATUS_SIZE,8)
  INTEGER mpireq(8)
  INTEGER ims_pro_l, ims_pro_r
  INTEGER icount

! #######################################################################
  IF ( ims_npro .GT. 1 ) THEN

! Careful in case only 2 PEs
     IF ( ims_npro .EQ. 2 ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'DNS_MPI_COPYPLN_2. Undeveloped for 2 PEs.')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
     ENDIF

! -----------------------------------------------------------------------
! left and right PEs. Same as in routine DNS_MPI_COPYPLN_1
! -----------------------------------------------------------------------
     npl_loc = kmax

     ims_pro_l = MOD(ims_pro-1+ims_npro,ims_npro)
     ims_pro_r = MOD(ims_pro+1+ims_npro,ims_npro)

     icount   = ijmax*npl_loc

     CALL MPI_IRECV(bl(1,npl-kmax+1), icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(3), ims_err)
     CALL MPI_IRECV(br(1,1),          icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(4), ims_err)

     CALL MPI_ISEND(a(1,1), icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(1), ims_err)
     CALL MPI_ISEND(a(1,1), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(2), ims_err)

     CALL MPI_WAITALL(4, mpireq, status, ims_err)

     CALL DNS_MPI_TAGUPDT

! -----------------------------------------------------------------------
! second-left and second-right PEs.
! -----------------------------------------------------------------------
     npl_loc = npl-kmax

     ims_pro_l = MOD(ims_pro-2+ims_npro,ims_npro)
     ims_pro_r = MOD(ims_pro+2+ims_npro,ims_npro)

     icount   = ijmax*npl_loc

     CALL MPI_IRECV(bl(1,1),      icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(7), ims_err)
     CALL MPI_IRECV(br(1,1+kmax), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(8), ims_err)

     CALL MPI_ISEND(a(1,1),              icount, MPI_REAL8, ims_pro_l, &
          ims_tag, MPI_COMM_WORLD, mpireq(5), ims_err)
     CALL MPI_ISEND(a(1,kmax+1-npl_loc), icount, MPI_REAL8, ims_pro_r, &
          ims_tag, MPI_COMM_WORLD, mpireq(6), ims_err)

     CALL MPI_WAITALL(4, mpireq(5), status, ims_err)

     CALL DNS_MPI_TAGUPDT

  ENDIF

  RETURN
END SUBROUTINE DNS_MPI_COPYPLN_2
