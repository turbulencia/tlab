!########################################################################
!# Tool/Library PDF
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/07/11 - J.P. Mellado
!#              Cleaned
!# 2008/04/02 - J.P. Mellado
!#              Reformulation of the gate signal
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE INTER1V2D(nx, ny, nz, j, igate, gate, inter)

  IMPLICIT NONE

#include "types.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER nx, ny, nz, j
  TREAL inter
  INTEGER(1) gate(nx,ny,nz), igate

! -------------------------------------------------------------------
  TINTEGER i, k, nz_total

#ifdef USE_MPI
  TREAL inter_p
  INTEGER ims_err
#endif

! ###################################################################
  inter = C_0_R
  DO k = 1,nz
     DO i = 1,nx
        IF ( gate(i,j,k) .EQ. igate ) THEN
           inter = inter + C_1_R
        ENDIF
     ENDDO
  ENDDO

#ifdef USE_MPI
  CALL MPI_ALLREDUCE(inter, inter_p, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  inter = inter_p

  CALL MPI_ALLREDUCE(nz, nz_total, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ims_err)
#else
  nz_total = nz
#endif

  inter = inter/M_REAL(nx*nz_total)

  RETURN
END SUBROUTINE INTER1V2D
