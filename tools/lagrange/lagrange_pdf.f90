#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "LAGRANGE_PDF"

!########################################################################
!# Tool/Library PLOT
!#
!########################################################################
!# HISTORY
!#
!# 2014/08/11 - L. Muessle
!#
!#
!########################################################################
!# DESCRIPTION
!#
!# Post processing particle pdf
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
PROGRAM LAGRANGE_PDF
  
  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE LAGRANGE_GLOBAL
#ifdef USE_MPI
  USE DNS_MPI
#endif
    
  IMPLICIT NONE
#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

! -------------------------------------------------------------------
  TINTEGER  ierr,isize_wrk3d, i
  TREAL, DIMENSION(:),      ALLOCATABLE :: x,y,z, dx,dy,dz
  TREAL, DIMENSION(:),      ALLOCATABLE :: wrk1d,wrk2d, wrk3d
  TREAL, DIMENSION(:,:),      ALLOCATABLE :: txc
  TREAL, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE    :: s  

  TREAL, DIMENSION(:,:),    ALLOCATABLE :: l_q, l_txc, l_hq
  INTEGER(8), DIMENSION(:), ALLOCATABLE :: l_tags
  TLONGINTEGER, DIMENSION(:,:),   ALLOCATABLE         :: particle_bins
  TREAL, DIMENSION(:),   ALLOCATABLE         :: counter_interval
  TREAL y_pdf_max, y_pdf_min

  TINTEGER nitera_first,nitera_last,nitera_save

  CHARACTER*32 inifile
  CHARACTER*64 str, fname
  CHARACTER*128 line
  CHARACTER*32 bakfile


  inifile = 'dns.ini'
  bakfile = TRIM(ADJUSTL(inifile))//'.bak'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  IF ( icalc_particle .EQ. 1 ) THEN
     CALL PARTICLE_READ_GLOBAL('dns.ini')
  ENDIF
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif
!  CALL DNS_READ_LOCAL(inifile) !for nitera stuff

! Get the local information from the dns.ini
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Start',      '0',  nitera_first)
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'End',        '0',  nitera_last )
  CALL SCANINIINT(bakfile, inifile, 'Iteration', 'Restart',    '50', nitera_save )

  number_of_bins = particle_pdf_max/particle_pdf_interval
  inb_particle_txc = 1 

#include "dns_alloc_larrays.h"
  isize_wrk3d = imax*jmax*kmax
  isize_wrk3d = MAX(isize_wrk3d,(imax+1)*jmax*(kmax+1))
  isize_wrk3d = MAX(isize_wrk3d,(jmax*(kmax+1)*inb_lag_total_interp*2))
  isize_wrk3d = MAX(isize_wrk3d,(jmax*(imax+1)*inb_lag_total_interp*2))

  isize_wrk2d = MAX(isize_wrk2d, jmax*inb_lag_total_interp)

  IF (jmax_part .EQ. 1) THEN
     jmax_part   = jmax ! 1 by default
  ENDIF
! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------      
  ALLOCATE(x(imax_total))
  ALLOCATE(y(jmax_total))
  ALLOCATE(z(kmax_total))
  ALLOCATE(dx(imax_total*inb_grid))
  ALLOCATE(dy(jmax_total*inb_grid))
  ALLOCATE(dz(kmax_total*inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d))
  ALLOCATE(wrk3d(isize_wrk3d))

  ALLOCATE(s(imax,jmax,kmax, inb_scal_array))

  IF (ilagrange .EQ. LAG_TYPE_BIL_CLOUD .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_2 &
       .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_3 .OR. ilagrange .EQ. LAG_TYPE_BIL_CLOUD_4) THEN !Allocte memory to read fields
     ALLOCATE(txc(isize_field,3))
     ALLOCATE(l_hq(isize_particle,inb_particle)) !Rubish information. Just to run FIELD_TO_PARTICLE properly
  ENDIF
  ALLOCATE(particle_bins(number_of_bins,3)) 
  ALLOCATE(counter_interval(number_of_bins)) 

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

!#######################################################################
!Setup the file properties
!#######################################################################
  y_pdf_max=y_particle_pdf_pos+0.5*y_particle_pdf_width
  y_pdf_min=y_particle_pdf_pos-0.5*y_particle_pdf_width
  particle_bins=int(0,KIND=8)

  DO i=nitera_first,nitera_last,nitera_save

!#######################################################################
!READ ALL FILES
!#######################################################################
     fname='scal'
     WRITE(str,*) i;  str = TRIM(ADJUSTL(fname))//"."//TRIM(ADJUSTL(str))
     CALL DNS_READ_FIELDS(str, i1, imax,jmax,kmax, inb_scal, i0, isize_wrk3d, s, wrk3d)
!     CALL FI_LIQUIDWATER(ibodyforce, imax,jmax,kmax, body_param, s(:,:,:,1),s(:,:,:,inb_scal_array)) !Update the liquid function
     CALL THERMO_AIRWATER_LINEAR(imax,jmax,kmax, s, s(:,:,:,inb_scal_array), wrk3d)

     WRITE(fname,*) nitera_first; fname = "particle_id."//TRIM(ADJUSTL(fname))
     CALL DNS_READ_PARTICLE_TAGS(fname,l_tags)

     WRITE(fname,*) i; fname = "particle."//TRIM(ADJUSTL(fname))
     CALL DNS_READ_PARTICLE(fname,l_q) ! h_particle only as dummy

! ######################################################################
! Save particle pathlines for particle_pdf
! ######################################################################
     number_of_bins = particle_pdf_max/particle_pdf_interval

     WRITE(fname,*) i; fname = "particle_pdf."//TRIM(ADJUSTL(fname))
     CALL PARTICLE_PDF(fname,s,wrk1d,wrk2d,wrk3d,x,y,z,l_txc,l_tags,l_hq,l_q)

  ENDDO

  CALL DNS_END(0)

  STOP
END PROGRAM
