#include "types.h"

!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/03/09 - J.P. Mellado
!#              Reducing number of ops.
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the average of the plane j in array a over nonuniform grids.
!#
!########################################################################
TREAL FUNCTION AVG_IK(nx,ny,nz, j, a, dx,dz, area)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN) :: nx,ny,nz, j
  TREAL, DIMENSION(*),        INTENT(IN) :: dx,dz
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN) :: a
  TREAL,                      INTENT(IN) :: area

! -------------------------------------------------------------------
  TINTEGER i, k, idsp, kdsp
  TREAL sum

! ###################################################################
#ifdef USE_MPI
  idsp = ims_offset_i 
  kdsp = ims_offset_k 
#else 
  idsp = 0
  kdsp = 0
#endif

! number of + ops: nx*nz*1 + nz*1
! number of * ops: nx*nz*1 + nz*1
  AVG_IK = C_0_R
  DO k = 1,nz
     sum = C_0_R
     DO i = 1,nx
        sum = sum + a(i,j,k)*dx(idsp+i)
     ENDDO
     AVG_IK = AVG_IK + sum*dz(kdsp+k)
  ENDDO

#ifdef USE_MPI
  sum = AVG_IK/area
  CALL MPI_ALLREDUCE(sum, AVG_IK, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#else
  AVG_IK = AVG_IK/area
#endif

  RETURN
END FUNCTION AVG_IK

!########################################################################
!########################################################################
! Vector form
SUBROUTINE AVG_IK_V(nx,ny,nz, jm, a, dx,dz, avg, wrk, area)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_offset_i, ims_offset_k, ims_err
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz, jm
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dz
  TREAL, DIMENSION(nx,ny,nz), INTENT(IN)    :: a
  TREAL, DIMENSION(jm),       INTENT(OUT)   :: avg
  TREAL, DIMENSION(jm),       INTENT(INOUT) :: wrk
  TREAL,                      INTENT(IN)    :: area

! -------------------------------------------------------------------
#ifdef USE_MPI
  INTEGER len
#endif

  TINTEGER i, j, k, idsp, kdsp

! ###################################################################
#ifdef USE_MPI
  idsp = ims_offset_i 
  kdsp = ims_offset_k 
#else 
  idsp = 0
  kdsp = 0
#endif

  avg = C_0_R

  DO k = 1,nz
     DO j = 1,jm
        DO i = 1,nx
           avg(j) = avg(j) + a(i,j,k)*dx(i+idsp)*dz(k+kdsp)
        ENDDO
     ENDDO
  ENDDO

#ifdef USE_MPI
  len = jm
  CALL MPI_REDUCE(avg, wrk, len, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ims_err)
  avg = wrk/area
#else
  avg = avg/area
#endif

  RETURN
END SUBROUTINE AVG_IK_V

!########################################################################
! Reynolds fluctuations of array a
!########################################################################
SUBROUTINE REYFLUCT2D(nx,ny,nz, dx,dz, area, a)

  IMPLICIT NONE

  TINTEGER,                   INTENT(IN)    :: nx,ny,nz
  TREAL, DIMENSION(*),        INTENT(IN)    :: dx,dz
  TREAL, DIMENSION(nx,ny,nz), INTENT(INOUT) :: a
  TREAL,                      INTENT(IN)    :: area

! -------------------------------------------------------------------
  TREAL dummy, AVG_IK
  TINTEGER j

! ###################################################################
  DO j = 1,ny
     dummy = AVG_IK(nx,ny,nz, j, a, dx,dz, area)
     a(:,j,:) = a(:,j,:) - dummy
  ENDDO

  RETURN
END SUBROUTINE REYFLUCT2D
