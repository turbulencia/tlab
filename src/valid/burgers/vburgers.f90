#include "types.h"
#include "dns_const.h"

PROGRAM VBURGERS

  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_PROCS
  USE TLAB_MPI_VARS
#endif
  USE IO_FIELDS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
  TREAL error2, dummy2
#else
   TINTEGER, PARAMETER                 :: ims_pro=0
#endif

  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: a, b, c
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: wrk1d, wrk2d
  TREAL, DIMENSION(:),     ALLOCATABLE :: wrk3d, tmp1

  TINTEGER i, j, k,  bcs(2,2)
  TREAL dummy, error

! ###################################################################
  CALL TLAB_START()

  CALL IO_READ_GLOBAL('dns.ini')
#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

  isize_wrk3d = isize_txc_field

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d+1))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(a(imax,jmax,kmax),b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(tmp1(isize_txc_field),wrk3d(isize_wrk3d))

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  CALL FI_BACKGROUND_INITIALIZE(wrk1d)

  bcs = 0

! ###################################################################
! Define forcing term
! ###################################################################
  CALL IO_READ_FIELDS('field.inp', IO_SCAL, imax,jmax,kmax, i1,i0, a, wrk3d)

! ###################################################################
  CALL OPR_PARTIAL_X(OPR_P2_P1, imax,jmax,kmax, bcs, g(1), a,b, c, wrk2d,wrk3d)
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
!           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
           b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  CALL OPR_BURGERS_X(i0,i0, imax,jmax,kmax, bcs, g(1), a,a,a, c, tmp1, wrk2d,wrk3d)

  c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
  error2 = error; dummy2 = dummy
  CALL MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  CALL MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  IF (ims_pro == 0) THEN
     WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
  ENDIF
! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)

! ###################################################################
  CALL OPR_PARTIAL_Y(OPR_P2_P1, imax,jmax,kmax, bcs, g(2), a,b, c, wrk2d,wrk3d)
  DO k = 1,kmax
     DO j = 1,jmax
        DO i = 1,imax
!           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
           b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
        ENDDO
     ENDDO
  ENDDO

  CALL OPR_BURGERS_Y(i0,i0, imax,jmax,kmax, bcs, g(2), a,a,a, c, tmp1, wrk2d,wrk3d)

  c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
  error2 = error; dummy2 = dummy
  CALL MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  CALL MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  IF (ims_pro == 0) THEN
     WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
  ENDIF
! CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)

! ###################################################################
  IF ( g(3)%size .GT. 1 ) THEN

     CALL OPR_PARTIAL_Z(OPR_P2_P1, imax,jmax,kmax, bcs, g(3), a,b, c, wrk2d,wrk3d)
     DO k = 1,kmax
        DO j = 1,jmax
           DO i = 1,imax
   !           b(i,j,k) = b(i,j,k) *visc - a(i,j,k) *c(i,j,k)
              b(i,j,k) = b(i,j,k) *visc *ribackground(j)- a(i,j,k) *c(i,j,k)
           ENDDO
        ENDDO
     ENDDO

     CALL OPR_BURGERS_Z(i0,i0, imax,jmax,kmax, bcs, g(3), a,a,a, c, tmp1, wrk2d,wrk3d)

     c = c - b; error = sum(c**2); dummy = sum(b**2)
#ifdef USE_MPI
     error2 = error; dummy2 = dummy
     CALL MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
     CALL MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
     IF (ims_pro == 0) THEN
        WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
     ENDIF
!    CALL IO_WRITE_FIELDS('field.dif', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)

  END IF

  CALL TLAB_STOP(0)
END PROGRAM VBURGERS
