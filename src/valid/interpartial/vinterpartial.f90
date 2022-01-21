#include "types.h"
#include "dns_const.h"
!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2022/01/21 - J. Kostelecky
!#              Created           
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate compact interpolation schemes.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM VINTERPARTIAL

  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_PROCS
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#else
  TINTEGER, PARAMETER                  :: ims_pro=0, ims_npro=1
#endif
 
  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: a, a_int, a_dif
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: b, c
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: wrk1d, wrk2d
  TREAL, DIMENSION(:),     ALLOCATABLE :: wrk3d, tmp1

  TINTEGER i, j, k,  bcs(2,2)
  TREAL dummy, error
! ###################################################################
  CALL TLAB_START()

  CALL DNS_READ_GLOBAL('dns.ini')
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
  ALLOCATE(a(imax,jmax,kmax),a_int(imax,jmax,kmax),a_dif(imax,jmax,kmax))
  ALLOCATE(b(imax,jmax,kmax),c(imax,jmax,kmax))
  ALLOCATE(tmp1(isize_txc_field),wrk3d(isize_wrk3d))

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  bcs = 0
! ###################################################################
! Define forcing term
! ###################################################################
  CALL DNS_READ_FIELDS('field.inp', i1, imax,jmax,kmax, i1,i0, isize_wrk3d, a, wrk3d)

! ###################################################################
! x-direction: Interpolation + interpolatory 1st derivative
! ###################################################################
! Interpolation: vel. <--> pre. grid
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), a,     a_int, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P0_INT_PV, imax,jmax,kmax, bcs, g(1), a_int, b,     tmp1, wrk2d,wrk3d)
! Difference field and error
  a_dif = a - b; error = C_0_R; dummy = C_0_R
  DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
     error = error + a_dif(i,j,k)*a_dif(i,j,k)
     dummy = dummy + a(i,j,k)*a(i,j,k)
  ENDDO; ENDDO; ENDDO
  error = sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
  dummy = error
  CALL MPI_ALLREDUCE(dummy, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  error = error / ims_npro
#endif
  IF (ims_pro == 0) THEN
     WRITE(*,*) 'x-direction: Interpolation + interp. 1st derivative'
     WRITE(*,*) 'Relative error .............: ', error
  ENDIF
  ! CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, a_dif, wrk3d)
! -------------------------------------------------------------------
! 1st interp. deriv + Interpolation: vel. <--> pre. grid
  CALL OPR_PARTIAL_X(OPR_P1_INT_VP, imax,jmax,kmax, bcs, g(1), a,     a_int, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P0_INT_PV, imax,jmax,kmax, bcs, g(1), a_int, b,     tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1,        imax,jmax,kmax, bcs, g(1), a,     c,     tmp1, wrk2d,wrk3d)
! Difference field and error
  a_dif = c - b; error = C_0_R; dummy = C_0_R
  DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
       error = error + a_dif(i,j,k)*a_dif(i,j,k)
     dummy = dummy + a(i,j,k)*a(i,j,k)
  ENDDO; ENDDO; ENDDO
  error = sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
  dummy = error
  CALL MPI_ALLREDUCE(dummy, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  error = error / ims_npro
#endif
  IF (ims_pro == 0) THEN
     WRITE(*,*) 'Relative error .............: ', error
  ENDIF
  ! CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, a_dif, wrk3d)
! -------------------------------------------------------------------
! 1st interp. deriv + Interpolation: vel. <--> pre. grid
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), a,     a_int, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1_INT_PV, imax,jmax,kmax, bcs, g(1), a_int, b,     tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P1,        imax,jmax,kmax, bcs, g(1), a,     c,     tmp1, wrk2d,wrk3d)
! Difference field and error
  a_dif = c - b; error = C_0_R; dummy = C_0_R
  DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
     error = error + a_dif(i,j,k)*a_dif(i,j,k)
     dummy = dummy + a(i,j,k)*a(i,j,k)
  ENDDO; ENDDO; ENDDO
  error = sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
  dummy = error
  CALL MPI_ALLREDUCE(dummy, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  error = error / ims_npro
#endif
  IF (ims_pro == 0) THEN
     WRITE(*,*) 'Relative error .............: ', error
     WRITE(*,*) '--------------------------------------------------------'
  ENDIF
  ! CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, a_dif, wrk3d)
! ###################################################################
! z-direction: Interpolation + interpolatory 1st derivative
! ###################################################################  
  IF ( g(3)%size .GT. 1 ) THEN
   ! Interpolation: vel. <--> pre. grid
     CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(3), a,     a_int, tmp1, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P0_INT_PV, imax,jmax,kmax, bcs, g(3), a_int, b,     tmp1, wrk2d,wrk3d)
   ! Difference field and error
     a_dif = a - b; error = C_0_R; dummy = C_0_R
     DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
        error = error + a_dif(i,j,k)*a_dif(i,j,k)
        dummy = dummy + a(i,j,k)*a(i,j,k)
     ENDDO; ENDDO; ENDDO
     error = sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
     dummy = error
     CALL MPI_ALLREDUCE(dummy, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
     error = error / ims_npro
#endif
     IF (ims_pro == 0) THEN
        WRITE(*,*) 'z-direction: Interpolation + interp. 1st derivative'
        WRITE(*,*) 'Relative error .............: ', error
     ENDIF
     ! CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, a_dif, wrk3d)
   ! -------------------------------------------------------------------
   ! 1st interp. deriv + Interpolation: vel. <--> pre. grid
     CALL OPR_PARTIAL_Z(OPR_P1_INT_VP, imax,jmax,kmax, bcs, g(3), a,     a_int, tmp1, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P0_INT_PV, imax,jmax,kmax, bcs, g(3), a_int, b,     tmp1, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1,        imax,jmax,kmax, bcs, g(3), a,     c,     tmp1, wrk2d,wrk3d)
   ! Difference field and error
     a_dif = c - b; error = C_0_R; dummy = C_0_R
     DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
        error = error + a_dif(i,j,k)*a_dif(i,j,k)
        dummy = dummy + a(i,j,k)*a(i,j,k)
     ENDDO; ENDDO; ENDDO
     error = sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
     dummy = error
     CALL MPI_ALLREDUCE(dummy, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
     error = error / ims_npro
#endif
     IF (ims_pro == 0) THEN
        WRITE(*,*) 'Relative error .............: ', error
     ENDIF
     ! CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, a_dif, wrk3d)
   ! -------------------------------------------------------------------
   ! 1st interp. deriv + Interpolation: vel. <--> pre. grid
     CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(3), a,     a_int, tmp1, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1_INT_PV, imax,jmax,kmax, bcs, g(3), a_int, b,     tmp1, wrk2d,wrk3d)
     CALL OPR_PARTIAL_Z(OPR_P1,        imax,jmax,kmax, bcs, g(3), a,     c,     tmp1, wrk2d,wrk3d)
   ! Difference field and error
     a_dif = c - b; error = C_0_R; dummy = C_0_R
     DO k = 1,kmax; DO j = 1,jmax; DO i = 1,imax
       error = error + a_dif(i,j,k)*a_dif(i,j,k)
       dummy = dummy + a(i,j,k)*a(i,j,k)
     ENDDO; ENDDO; ENDDO
     error = sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
     dummy = error
     CALL MPI_ALLREDUCE(dummy, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
     error = error / ims_npro
#endif
     IF (ims_pro == 0) THEN
        WRITE(*,*) 'Relative error .............: ', error
        WRITE(*,*) '--------------------------------------------------------'
     ENDIF
     ! CALL DNS_WRITE_FIELDS('field.dif', i1, imax,jmax,kmax, i1, isize_wrk3d, a_dif, wrk3d)
  ENDIF
! ###################################################################
  CALL TLAB_STOP(0)
END PROGRAM VINTERPARTIAL