#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2010/04/01 - J.P. Mellado
!#              Created
!# 2012/12/01 - J.P. Mellado
!#              Not using global variables {imax,jmax,kmax}_total
!#              in routines IO_FIELDS_* any more
!#
!########################################################################
!# DESCRIPTION
!#
!# Extracted from DNS_READ_FIELDS to handle different types of file formats
!#
!########################################################################
!# ARGUMENTS
!#
!# nfield     In      Number of fields in the file
!# iheader    In      Header control flag:
!#                    0 No header
!#                    1 Scalar fields
!#                    2 Flow fields
!#
!########################################################################
SUBROUTINE DNS_WRITE_FIELDS(fname, iheader, nx,ny,nz, nfield, itxc, a, txc)

  USE TLAB_VARS, ONLY : imode_files, imode_eqns
  USE DNS_CONSTANTS, ONLY : efile, lfile
  USE TLAB_VARS, ONLY : itime, rtime
  USE TLAB_VARS, ONLY : visc, froude, rossby, damkohler, prandtl, mach
  USE TLAB_VARS, ONLY : schmidt
  USE THERMO_GLOBAL, ONLY : gama0
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_npro_i, ims_npro_k
#endif

  IMPLICIT NONE

  CHARACTER*(*) fname
  TINTEGER itxc, iheader, nfield, nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,nfield) :: a
  TREAL, DIMENSION(itxc)            :: txc

! -------------------------------------------------------------------
  CHARACTER*32 :: str
  CHARACTER*128 :: line
  TINTEGER nx_total, ny_total, nz_total
  TINTEGER ifield

  TINTEGER isize_max, isize
  PARAMETER(isize_max=20)
  TREAL params(isize_max)

! ###################################################################
#ifdef USE_MPI
  nx_total = nx*ims_npro_i
  ny_total = ny
  nz_total = nz*ims_npro_k
#else
  nx_total = nx
  ny_total = ny
  nz_total = nz
#endif

  line = 'Writing field '//TRIM(ADJUSTL(fname))//' of size'
  WRITE(str,*) nx_total; line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(str))
  WRITE(str,*) ny_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))
  WRITE(str,*) nz_total; line = TRIM(ADJUSTL(line))//'x'//TRIM(ADJUSTL(str))//'...'

  CALL TLAB_WRITE_ASCII(lfile, line)

! ###################################################################
  IF      ( imode_files .EQ. DNS_FILE_RAWARRAY ) THEN
     CALL IO_WRITE_FIELDS_ARRAY(fname, nx,ny,nz, itxc, iheader, nfield, a, txc)

! ###################################################################
  ELSE IF ( imode_files .EQ. DNS_FILE_RAWSPLIT ) THEN
#ifdef USE_MPI
#ifdef USE_MPI_IO
     IF ( itxc .LT. nx*ny*nz ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'DNS_WRITE_FIELDS. Work array size error')
        CALL TLAB_STOP(DNS_ERROR_ALLOC)
     ENDIF
#endif
#endif

! define header info
     isize = 0
     IF      ( iheader .EQ. 1 ) THEN ! Scalar
        isize = isize+1; params(isize) = rtime
        isize = isize+1; params(isize) = visc ! inverse of reynolds
        isize = isize+1+1                     ! prepare space for schmidt and damkohler

     ELSE IF ( iheader .EQ. 2 ) THEN ! Flow
        isize = isize+1; params(isize) = rtime
        isize = isize+1; params(isize) = visc ! inverse of reynolds
        isize = isize+1; params(isize) = froude
        isize = isize+1; params(isize) = rossby
        IF ( imode_eqns .EQ. DNS_EQNS_INTERNAL .OR. imode_eqns .EQ. DNS_EQNS_TOTAL ) THEN
           isize = isize+1; params(isize) = gama0
           isize = isize+1; params(isize) = prandtl
           isize = isize+1; params(isize) = mach
        ENDIF
     ENDIF

     IF ( isize .GT. isize_max ) THEN
        CALL TLAB_WRITE_ASCII(efile, 'DNS_WRITE_FIELDS. Parameters array size error.')
        CALL TLAB_STOP(DNS_ERROR_ALLOC)
     ENDIF

! write data
     DO ifield = 1,nfield
        IF ( iheader .EQ. 1 ) params(isize-1) = schmidt(ifield)   ! Scalar header
        IF ( iheader .EQ. 1 ) params(isize  ) = damkohler(ifield) ! Scalar header
        WRITE(str,'(I2)') ifield
        str=TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(str))
        CALL IO_WRITE_FIELDS_SPLIT(str, iheader, nx,ny,nz,itime, isize, params, a(1,ifield),txc)

     ENDDO

  ELSE IF ( imode_files .EQ. DNS_FILE_NETCDF )  THEN
     ! To be implemented
  ELSE IF ( imode_files .EQ. DNS_NOFILE )       THEN
     ! Do nothing
  ENDIF

  RETURN
END SUBROUTINE DNS_WRITE_FIELDS
