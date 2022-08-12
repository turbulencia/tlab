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
!# DESCRIPTION
!#
!########################################################################
PROGRAM SL_BOUNDARY

  USE TLAB_VARS
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#include "integers.h"

! -------------------------------------------------------------------
! Grid and associated arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z

! Flow variables
  TREAL, DIMENSION(:,:), ALLOCATABLE, TARGET :: q
  TREAL, DIMENSION(:),   ALLOCATABLE         :: s, field

! Work arrays
  TREAL, DIMENSION(:),   ALLOCATABLE         :: wrk1d, wrk2d, wrk3d

! Surface arrays
  TREAL, DIMENSION(:,:), ALLOCATABLE         :: sl, samples

  TREAL txc(:)
  ALLOCATABLE txc

  TREAL pdf(:)
  ALLOCATABLE pdf

! Pointers to existing allocated space
  TREAL, DIMENSION(:),   POINTER :: u, v, w

  TINTEGER iopt, iint, isl, ith,  itxc_size, nfield, np
  TINTEGER iread_flow, iread_scal, jmin_loc, jmax_loc, idummy
  TREAL threshold, vmin, vmax
  TINTEGER buff_nps_u_jmin, buff_nps_u_jmax
  CHARACTER*64 str
  CHARACTER*32 fname, bakfile

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
  bakfile = TRIM(ADJUSTL(ifile))//'.bak'

  CALL DNS_START

  CALL IO_READ_GLOBAL(ifile)
#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

  CALL SCANINIINT(bakfile, ifile, 'BufferZone', 'PointsUJmin', '0', buff_nps_u_jmin)
  CALL SCANINIINT(bakfile, ifile, 'BufferZone', 'PointsUJmax', '0', buff_nps_u_jmax)

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(sl(imax*kmax,6))
  ALLOCATE(wrk1d(isize_wrk1d* 5))
  ALLOCATE(wrk2d(isize_wrk2d*10))
  ALLOCATE(wrk3d(isize_field   ))

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
     CALL SCANINICHAR(lfile, 'dns.ini', 'PostProcessing', 'ParamSuperlayer', '-1', sRes)
     iopt_size = iopt_size_max
     CALL LIST_REAL(sRes, iopt_size, opt_vec)

     IF ( sRes .EQ. '-1' ) THEN
        WRITE(*,*) 'Option ?'
        WRITE(*,*) '1. Extract scalar gradient envelope surfaces to file'
        WRITE(*,*) '2. PDFs conditioned on the envelope surface'
        WRITE(*,*) '3. RQ JPDF conditioned on the envelope surface'
        WRITE(*,*) '4. WS JPDF conditioned on the envelope surface'
        READ(*,*) iopt

        WRITE(*,*) 'Intermittency conditioning ?'
        WRITE(*,*) ' 1. Based on scalar'
        WRITE(*,*) ' 2. Based on vorticity'
        WRITE(*,*) ' 3. Based on scalar gradient'
        READ(*,*) iint

        WRITE(*,*) 'Threshold based on relative (1) or absolute (2) values?'
        READ(*,*) ith
        WRITE(*,*) 'Threshold value ?'
        READ(*,*) threshold
        IF ( iopt .GT. 1 ) THEN
           WRITE(*,*) 'Upper (1), lower (2) or both (3) envelope surfaces ?'
           READ(*,*) isl
           WRITE(*,*) 'Number of PDF bins ?'
           READ(*,*) np
        ENDIF
     ELSE
        iopt = INT(opt_vec(1))
        iint = INT(opt_vec(2))
        ith  = INT(opt_vec(3))
        threshold = opt_vec(4)
        isl  = INT(opt_vec(5))
        np   = INT(opt_vec(6))
     ENDIF

#ifdef USE_MPI
  ENDIF

  CALL MPI_BCAST(iopt,      1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(iint,      1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(isl,       1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(np,        1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(ith,       1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(threshold, 1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
#endif

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
  IF      ( iopt .EQ. 1 ) THEN; itxc_size = isize_field*2; nfield = 1;
  ELSE IF ( iopt .EQ. 2 ) THEN; itxc_size = isize_field*6; nfield = 5; iread_flow = 1; iread_scal = 1
  ELSE IF ( iopt .GE. 3 ) THEN; itxc_size = isize_field*6; nfield = 4; iread_flow = 1; iread_scal = 0
  ENDIF

  IF      ( iint .EQ. 1 ) THEN; iread_scal = 1
  ELSE IF ( iint .EQ. 2 ) THEN; iread_flow = 1
  ELSE IF ( iint .EQ. 3 ) THEN; iread_scal = 1
  ENDIF

  IF ( iread_flow .EQ. 1 ) ALLOCATE(q(isize_field,3))
  IF ( iread_scal .EQ. 1 ) ALLOCATE(s(isize_field))

  ALLOCATE(pdf(nfield*np))
  ALLOCATE(samples(imax*kmax,nfield*2))
  ALLOCATE(txc(itxc_size))

  IF ( iopt .LE. 2 ) ALLOCATE(field(isize_field))

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
CALL FDM_INITIALIZE(x, g(1), wrk1d)
CALL FDM_INITIALIZE(y, g(2), wrk1d)
CALL FDM_INITIALIZE(z, g(3), wrk1d)

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

     IF ( iread_flow .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_flow))//TRIM(ADJUSTL(fname))
        CALL IO_READ_FIELDS(fname, IO_FLOW, imax,jmax,kmax, i3,i0, q, wrk3d)
     ENDIF

     IF ( iread_scal .EQ. 1 ) THEN
        WRITE(fname,*) itime; fname = TRIM(ADJUSTL(tag_scal))//TRIM(ADJUSTL(fname))
        CALL IO_READ_FIELDS(fname, IO_SCAL, imax,jmax,kmax, inb_scal,inb_scal, s, wrk3d)
     ENDIF

! -------------------------------------------------------------------
! Extract scalar gradient surfaces to file
! -------------------------------------------------------------------
     IF ( iopt .EQ. 1 ) THEN
        jmin_loc = MAX(1,buff_nps_u_jmin)                 ! remove buffers
        jmax_loc = MIN(jmax,jmax - buff_nps_u_jmax +1)

! Based on scalar
        IF ( iint .EQ. 1 ) THEN
           IF      ( ith .EQ. 1 ) THEN ! relative to max
              CALL MINMAX(imax,jmax,kmax, s, vmin,vmax)
              vmin = threshold*threshold*vmax
           ELSE IF ( ith .EQ. 2 ) THEN ! absolute
              vmin = threshold
           ENDIF
           CALL SL_UPPER_BOUNDARY(imax,jmax,kmax, jmax_loc, vmin, y, s, txc, sl(1,1), wrk2d)
           CALL SL_LOWER_BOUNDARY(imax,jmax,kmax, jmin_loc, vmin, y, s, txc, sl(1,2), wrk2d)

! Based on vorticity
        ELSE IF ( iint .EQ. 2 ) THEN
           CALL TLAB_WRITE_ASCII(lfile,'Calculating vorticity...')
           CALL FI_VORTICITY(imax,jmax,kmax, u,v,w, field, txc(1),txc(1+isize_field), wrk2d,wrk3d)
           IF      ( ith .EQ. 1 ) THEN ! relative to max
              CALL MINMAX(imax,jmax,kmax, field, vmin,vmax)
              vmin = threshold*threshold*vmax
           ELSE IF ( ith .EQ. 2 ) THEN ! absolute
              vmin = threshold
           ENDIF

           CALL SL_UPPER_BOUNDARY(imax,jmax,kmax, jmax_loc, vmin, y, field, txc, sl(1,1), wrk2d)
           CALL SL_LOWER_BOUNDARY(imax,jmax,kmax, jmin_loc, vmin, y, field, txc, sl(1,2), wrk2d)

! Based on scalar gradient
        ELSE IF ( iint .EQ. 3 ) THEN
           CALL TLAB_WRITE_ASCII(lfile,'Calculating scalar gradient...')
           CALL FI_GRADIENT(imax,jmax,kmax, s,field, txc, wrk2d,wrk3d)
           CALL MINMAX(imax,jmax,kmax, field, vmin,vmax)
           WRITE(str,'(E22.15E3,E22.15E3)') vmin,vmax; str='Bounds '//TRIM(ADJUSTL(str))
           CALL TLAB_WRITE_ASCII(lfile,str)
           IF      ( ith .EQ. 1 ) THEN ! relative to max
              CALL MINMAX(imax,jmax,kmax, field, vmin,vmax)
              vmin = threshold*threshold*vmax
           ELSE IF ( ith .EQ. 2 ) THEN ! absolute
              vmin = threshold
           ENDIF

           CALL SL_UPPER_BOUNDARY(imax,jmax,kmax, jmax_loc, vmin, y, field, txc, sl(1,1), wrk2d)
           CALL SL_LOWER_BOUNDARY(imax,jmax,kmax, jmin_loc, vmin, y, field, txc, sl(1,2), wrk2d)

        ENDIF

! write threshold
        WRITE(str,'(E22.15E3)') vmin; str='Threshold '//TRIM(ADJUSTL(str))
        CALL TLAB_WRITE_ASCII(lfile,str)

! save surfaces w/o header
        WRITE(fname,*) itime; fname = 'sl'//TRIM(ADJUSTL(fname))
        ! idummy = g(2)%size; g(2)%size = 1
        ! CALL DNS_WRITE_FIELDS(fname, i0, imax,i1,kmax, i2, isize_field, sl, wrk3d)
        ! g(2)%size = idummy
        CALL TLAB_WRITE_ASCII(efile, 'SL_BOUNDARY. To be written in terms of IO_SUBARRAY as in averages.x')
        CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)

! save scalar dissipation as scalar field
        ! WRITE(fname,*) itime; fname = 'chi'//TRIM(ADJUSTL(fname))
        ! CALL IO_WRITE_FIELDS(fname, IO_SCAL, imax,jmax,kmax, i1, field, wrk3d)

! -------------------------------------------------------------------
! Surface PDFs
! -------------------------------------------------------------------
     ELSE IF ( iopt .EQ. 2 ) THEN
        CALL SL_BOUNDARY_VORTICITY_PDF(isl, ith, np, nfield, itxc_size, threshold, buff_nps_u_jmax, &
             u,v,w,s,field, sl, samples, pdf, txc, wrk1d,wrk2d,wrk3d)

! -------------------------------------------------------------------
! Surface JPDFs
! -------------------------------------------------------------------
     ELSE IF ( iopt .GE. 3 ) THEN
        CALL SL_BOUNDARY_VORTICITY_JPDF(iopt, isl, ith, np, nfield, itxc_size, threshold, buff_nps_u_jmax, &
             u,v,w, sl, samples, txc, wrk1d,wrk2d,wrk3d)

     ENDIF

  ENDDO

  CALL TLAB_STOP(0)
END PROGRAM SL_BOUNDARY
