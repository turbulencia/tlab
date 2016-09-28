!########################################################################
!# Tool/Library SUPERLAYER
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

PROGRAM SL_NORMAL_ANALYSIS
  
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
  TREAL, DIMENSION(:,:), POINTER :: q

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u, v, w, p

  TREAL z1(:)
  ALLOCATABLE z1
  TREAL field(:)
  ALLOCATABLE field
  TREAL sl(:)
  ALLOCATABLE sl
  TREAL txc(:)
  ALLOCATABLE txc
  TREAL profiles(:)
  ALLOCATABLE profiles
  TREAL mean(:)
  ALLOCATABLE mean
  TREAL wrk1d(:)
  ALLOCATABLE wrk1d
  TREAL wrk2d(:)
  ALLOCATABLE wrk2d
  TREAL wrk3d(:)
  ALLOCATABLE wrk3d

  TINTEGER iopt, isl, ith, isize_wrk3d, itxc_size, iavg
  TREAL threshold
  TINTEGER ibuffer_npy
  TINTEGER nmax, istep, kstep, nprof_size, nfield
  CHARACTER*32 fname, inifile, bakfile

  TINTEGER itime_size_max, itime_size, i
  PARAMETER(itime_size_max=128)
  TINTEGER itime_vec(itime_size_max)
  TINTEGER iopt_size_max, iopt_size
  PARAMETER(iopt_size_max=10)
  TREAL opt_vec(iopt_size_max)
  CHARACTER*512 sRes
#ifdef USE_MPI
  INTEGER icount
#endif

  TREAL, DIMENSION(:,:), POINTER :: dx, dy, dz
  
! ###################################################################
  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL('dns.ini')

#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  CALL SCANINIINT(bakfile, inifile, 'BufferZone', 'NumPointsY', '0', ibuffer_npy)

  isize_wrk3d = imax*jmax*kmax
  itxc_size = imax*jmax*kmax*7

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(u(imax*jmax*kmax))
  ALLOCATE(v(imax*jmax*kmax))
  ALLOCATE(w(imax*jmax*kmax))
  ALLOCATE(p(imax*jmax*kmax))
  ALLOCATE(z1(imax*jmax*kmax))
  ALLOCATE(field(imax*jmax*kmax))
  ALLOCATE(sl(imax*kmax))
  ALLOCATE(txc(itxc_size))
  ALLOCATE(wrk1d(isize_wrk1d*5))
  ALLOCATE(wrk2d(isize_wrk2d))
  ALLOCATE(wrk3d(isize_wrk3d))

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CALL SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'Files', '-1', sRes)
     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'Integral Iterations ?'
        READ(*,'(A512)') sRes
     ENDIF
     itime_size = itime_size_max
     CALL LIST_INTEGER(sRes, itime_size, itime_vec)
#ifdef USE_MPI
  ENDIF
  CALL MPI_BCAST(itime_size, 1,      MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  icount = itime_size
  CALL MPI_BCAST(itime_vec,  icount, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Read local options
! -------------------------------------------------------------------
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CALL SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'Superlayer', '-1', sRes)
     iopt_size = iopt_size_max
     CALL LIST_REAL(sRes, iopt_size, opt_vec)

     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'Option ?'
        WRITE(*,*) '1. Based on vorticity magnitude w_i w_i'
        WRITE(*,*) '2. Based on scalar gradient G_i G_i'
        WRITE(*,*) '3. Based on strain rate s_ij s_ij'
        READ(*,*) iopt
        WRITE(*,*) 'Upper (1) or lower (2) envelope surface ?'
        READ(*,*) isl
        WRITE(*,*) 'Number of normal points ?'
        READ(*,*) nmax
        WRITE(*,*) 'Averages (1) or instantaneous profiles (2) ?'
        READ(*,*) iavg
        WRITE(*,*) 'Step along OX ?'
        READ(*,*) istep
        WRITE(*,*) 'Step along OZ ?'
        READ(*,*) kstep
        WRITE(*,*) 'Threshold based on maximum (1) or mean (2) ?'
        READ(*,*) ith
        WRITE(*,*) 'Threshold value ?'
        READ(*,*) threshold
     ELSE
        iopt  = DINT(opt_vec(1))
        isl   = DINT(opt_vec(2))
        nmax  = DINT(opt_vec(3))
        iavg  = DINT(opt_vec(4))
        istep = DINT(opt_vec(5))
        kstep = DINT(opt_vec(6))
        ith   = DINT(opt_vec(7))
        threshold = opt_vec(8)
     ENDIF

#ifdef USE_MPI
  ENDIF

  CALL MPI_BCAST(iopt,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(isl,   1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(nmax,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(iavg,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(istep, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(kstep, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(ith,  1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(threshold, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
  nfield = 13
  nprof_size = (imax/istep)*(kmax/kstep)*nmax*nfield

  ALLOCATE(profiles(nprof_size))
  ALLOCATE(mean(nmax*nfield*2))

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

! ###################################################################
! Define pointers
! ###################################################################
  u   => q(:,1)
  v   => q(:,2)
  w   => q(:,3)
  p   => q(:,4)

! ###################################################################
! Postprocess given list of files
! ###################################################################
  DO i=1, itime_size

     itime = itime_vec(i)

! -------------------------------------------------------------------
! Binary data
! -------------------------------------------------------------------
     WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
     CALL DNS_READ_FIELDS(fname, i2, imax,jmax,kmax, i4,i0, isize_wrk3d, q, wrk3d)
     
     WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
     CALL DNS_READ_FIELDS(fname, i1, imax,jmax,kmax, inb_scal,inb_scal, isize_wrk3d, z1, wrk3d)
     
     CALL THERMO_CALORIC_TEMPERATURE(imax, jmax, kmax, z1, p, field, txc, wrk3d)
     CALL THERMO_THERMAL_PRESSURE(imax, jmax, kmax, z1, field, txc, p)

! -------------------------------------------------------------------
! Vorticity analysis
! -------------------------------------------------------------------
     IF ( iopt .EQ. 1 ) THEN
        CALL SL_NORMAL_VORTICITY(isl, ith, iavg, nmax, istep, kstep, nfield, itxc_size,&
             threshold, ibuffer_npy, x, y, z, dx, dy, dz, u, v, w, p, z1, field, sl, profiles, &
             txc, mean, wrk1d, wrk2d, wrk3d)

! -------------------------------------------------------------------
! Scalar gradient analysis
! -------------------------------------------------------------------
!     ELSE IF ( iopt .EQ. 2 ) THEN
!        CALL SL_NORMAL_GRADIENT(isl, nmax, istep, kstep, ibuffer_npy, x, y, z, dx, dy, dz, 
! $           u, v, w, z1, field, sl, profiles, txc, wrk1d, wrk2d, wrk3d)
     ENDIF

  ENDDO

  CALL DNS_END(0)

  STOP
END PROGRAM SL_NORMAL_ANALYSIS
