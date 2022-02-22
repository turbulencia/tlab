#include "types.h"
#include "dns_const.h"
!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2022/02/21 - J. Kostelecky
!#              Created           
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate compact interpolation schemes.
!# Always turn on staggering in dns.ini! Error otherwise...
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM VINTERPDERIV

  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE TLAB_PROCS
  USE IO_FIELDS
#ifdef USE_MPI
  USE TLAB_MPI_PROCS
  USE TLAB_MPI_VARS
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#else
  TINTEGER, PARAMETER                  :: ims_pro=0
#endif
 
  TREAL, DIMENSION(:,:),   ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: a, a_int, a_dif
  TREAL, DIMENSION(:,:,:), ALLOCATABLE :: b, c
  TREAL, DIMENSION(:,:),   ALLOCATABLE :: wrk1d, wrk2d
  TREAL, DIMENSION(:),     ALLOCATABLE :: wrk3d, tmp1, d

  TINTEGER i, j, k,  bcs(2,2), ip_a,ip_b,ip_t
  TREAL dummy, dummy2, error, error2
! ###################################################################
  CALL TLAB_START()

  CALL DNS_READ_GLOBAL('dns.ini')
#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

  isize_wrk3d = isize_txc_field

! -------------------------------------------------------------------
! Check input
! ------------------------------------------------------------------- 
  IF (istagger .EQ. 0) THEN
    CALL TLAB_WRITE_ASCII(efile,'VINTERPARTIAL. Set "StaggerGrid=yes" in dns.ini!')
    CALL TLAB_STOP(i0)
  ENDIF 

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))
  ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d+1))
  ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))
  ALLOCATE(a(imax,jmax,kmax),a_int(imax,jmax,kmax),a_dif(imax,jmax,kmax))
  ALLOCATE(b(imax,jmax,kmax),c(imax,jmax,kmax),d(imax*jmax*kmax))
  ALLOCATE(tmp1(isize_txc_field),wrk3d(isize_wrk3d))

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  bcs = 0
! ###################################################################
! Define forcing term
! ###################################################################
  CALL IO_READ_FIELDS('flow.0', IO_SCAL, imax,jmax,kmax, i1,i0, a, wrk3d)

! ###################################################################
! x-direction: Interpolation + interpolatory 1st derivative
! ###################################################################
! ! Interpolation + derivative
!   CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), a,     b,     tmp1, wrk2d,wrk3d)
!   CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(3), b,     a_int, tmp1, wrk2d,wrk3d)
!   CALL OPR_PARTIAL_X(OPR_P1,        imax,jmax,kmax, bcs, g(1), a_int, b,     tmp1, wrk2d,wrk3d)
! ! Interpolatory derivative
!   CALL OPR_PARTIAL_X(OPR_P1_INT_VP, imax,jmax,kmax, bcs, g(1), a,     a_int, tmp1, wrk2d,wrk3d)
!   CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(3), a_int, c,     tmp1, wrk2d,wrk3d)

! Interpolation + derivative
  ! CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), a,     a_int, tmp1, wrk2d,wrk3d)
  ! CALL OPR_PARTIAL_X(OPR_P1,        imax,jmax,kmax, bcs, g(1), a_int, b,     tmp1, wrk2d,wrk3d)

  CALL OPR_PARTIAL_X(OPR_P1,        imax,jmax,kmax, bcs, g(1), a,     a_int, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), a_int, b,     tmp1, wrk2d,wrk3d)

! Interpolatory derivative
  CALL OPR_PARTIAL_X(OPR_P1_INT_VP, imax,jmax,kmax, bcs, g(1), a,     c,     tmp1, wrk2d,wrk3d)




! Difference field and error
  a_dif = b - c; error = sum(a_dif**2); dummy = sum(c**2)
#ifdef USE_MPI
  error2 = error; dummy2 = dummy
  CALL MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  CALL MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  IF (ims_pro == 0) THEN
    WRITE(*,*) 'x-direction: Interpolation + interp. 1st derivative'
    WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
  ENDIF
  CALL IO_WRITE_FIELDS('Pressure000000', IO_SCAL, imax,jmax,kmax, i1, b, wrk3d)
  CALL IO_WRITE_FIELDS('Pressure000001', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)
  
! ###################################################################
! x-direction: Interpolation + interpolatory 1st derivative
! ###################################################################
! Interpolation + derivative
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), a,     b,     tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(3), b,     a_int, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_Z(OPR_P1,        imax,jmax,kmax, bcs, g(3), a_int, b,     tmp1, wrk2d,wrk3d)
! Interpolatory derivative
  CALL OPR_PARTIAL_Z(OPR_P1_INT_VP, imax,jmax,kmax, bcs, g(3), a,     a_int, tmp1, wrk2d,wrk3d)
  CALL OPR_PARTIAL_X(OPR_P0_INT_VP, imax,jmax,kmax, bcs, g(1), a_int, c,     tmp1, wrk2d,wrk3d)

! Difference field and error
  a_dif = b - c; error = sum(a_dif**2); dummy = sum(c**2)
#ifdef USE_MPI
  error2 = error; dummy2 = dummy
  CALL MPI_ALLREDUCE(dummy2, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
  CALL MPI_ALLREDUCE(error2, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
#endif
  IF (ims_pro == 0) THEN
    WRITE(*,*) 'x-direction: Interpolation + interp. 1st derivative'
    WRITE(*,*) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
  ENDIF
  CALL IO_WRITE_FIELDS('Pressure000002', IO_SCAL, imax,jmax,kmax, i1, b, wrk3d)
  CALL IO_WRITE_FIELDS('Pressure000003', IO_SCAL, imax,jmax,kmax, i1, c, wrk3d)

  CALL TLAB_STOP(0)
END PROGRAM VINTERPDERIV