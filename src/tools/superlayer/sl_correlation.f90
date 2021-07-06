#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# HISTORY
!#
!# 2007/09/04 - J.P. Mellado
!#              Created
!#
!########################################################################
PROGRAM SL_CORRELATION

  USE DNS_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
! Grid and associated arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z

! Flow variables
  TREAL, DIMENSION(:,:), ALLOCATABLE :: q

  TREAL, DIMENSION(:),   ALLOCATABLE :: z1
  TREAL, DIMENSION(:),   ALLOCATABLE :: tmp1, tmp2, tmp3, tmp4, tmp5
  TREAL, DIMENSION(:),   ALLOCATABLE :: wrk1d, wrk2d, wrk3d
  TREAL, DIMENSION(:),   ALLOCATABLE :: profiles

  TARGET q

  TINTEGER  ilog
  CHARACTER*32 fname

  TINTEGER itime_size_max, itime_size, i
  PARAMETER(itime_size_max=128)
  TINTEGER itime_vec(itime_size_max)
  TINTEGER iopt_size_max, iopt_size
  PARAMETER(iopt_size_max=10)
  TREAL opt_vec(iopt_size_max)
  CHARACTER*512 sRes
  CHARACTER*32 line
#ifdef USE_MPI
  INTEGER icount
#endif

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u, v, w

  TREAL, DIMENSION(:,:), POINTER :: dx, dy, dz

! ###################################################################
  CALL DNS_START

  CALL DNS_READ_GLOBAL('dns.ini')

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  isize_wrk3d = imax*jmax*kmax

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(q(imax*jmax*kmax,3))
  ALLOCATE(z1(imax*jmax*kmax))
  ALLOCATE(tmp1(imax*jmax*kmax))
  ALLOCATE(tmp2(imax*jmax*kmax))
  ALLOCATE(tmp3(imax*jmax*kmax))
  ALLOCATE(tmp4(imax*jmax*kmax))
  ALLOCATE(tmp5(imax*jmax*kmax))
  ALLOCATE(profiles(jmax*10))
  ALLOCATE(wrk1d(isize_wrk1d*5))
  ALLOCATE(wrk2d(isize_wrk2d  ))
  ALLOCATE(wrk3d(isize_wrk3d  ))

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CALL SCANINICHAR&
          (lfile, 'dns.ini', 'PostProcessing', 'Files', '-1', sRes)
     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'Integral Iterations ?'
        READ(*,'(A512)') sRes
     ENDIF
     itime_size = itime_size_max
     CALL LIST_INTEGER(sRes, itime_size, itime_vec)
#ifdef USE_MPI
  ENDIF
  CALL MPI_BCAST(itime_size, 1,     MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  icount = itime_size
  CALL MPI_BCAST(itime_vec, icount, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Read local options
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CALL SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'ParamSlCorr', '-1', sRes)
     iopt_size = iopt_size_max
     CALL LIST_REAL(sRes, iopt_size, opt_vec)

     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'Use Log Normal Variable (y/n)?'
        READ(*,'(A1)') line
        IF ( line(1:1) .EQ. 'y' ) THEN; ilog = 1
        ELSE;                           ilog = 0; ENDIF
     ELSE
        ilog = INT(opt_vec(1))
     ENDIF

#ifdef USE_MPI
  ENDIF

  CALL MPI_BCAST(ilog, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
! NONE

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! ###################################################################
! Define pointers
! ###################################################################
  dx => x(:,2:) ! to be removed
  dy => y(:,2:)
  dz => z(:,2:)

  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)

! ###################################################################
! Postprocess given list of files
! ###################################################################
  DO i=1, itime_size

     itime = itime_vec(i)

! read data
     WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
     CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i3,i0, isize_wrk3d, q, wrk3d)

     WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
     CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal,inb_scal, isize_wrk3d, z1, wrk3d)

! do correlations
     CALL SL_CORRELATION_1(ilog, y, dx, dy, dz, u, v, w, z1, profiles, &
          tmp1, tmp2, tmp3, tmp4, tmp5, wrk1d, wrk2d, wrk3d)

  ENDDO

  CALL DNS_STOP(0)
END PROGRAM SL_CORRELATION
