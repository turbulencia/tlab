#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

!########################################################################
!# DESCRIPTION
!#
!# Create random velocity/scalar fields according to a given spectrum and
!# distribution function. 
!#
!# Thermodynamics are initialized to a constant mean value, if needed.
!#
!# For a Gaussian distribution, the covariance matrix can also been 
!# specified.
!#
!########################################################################
PROGRAM INIRAND

  USE DNS_CONSTANTS
  USE DNS_GLOBAL
  USE RAND_LOCAL
#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif
#ifdef USE_FFTW
#include "fftw3.f"
#endif

! -------------------------------------------------------------------
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE, TARGET :: x,y,z
  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE         :: q, s, txc
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE         :: wrk1d,wrk2d,wrk3d

  TINTEGER iq, is, isize_wrk3d, ifourier_type
  TREAL AVG1V3D, dummy

  CHARACTER*32 inifile
  
! ###################################################################
  inifile = 'dns.ini'

  CALL DNS_INITIALIZE

  CALL DNS_READ_GLOBAL(inifile)
  CALL RAND_READ_LOCAL(inifile)
#ifdef USE_MPI
  CALL DNS_MPI_INITIALIZE
#endif

  CALL IO_WRITE_ASCII(lfile,'Initializing random fiels.')

  itime = 0; rtime = C_0_R

! seed for random generator
#ifdef USE_MPI
  seed = seed + ims_pro
#endif
  seed = - ABS(seed)
  
! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
  ALLOCATE(x(g(1)%size,g(1)%inb_grid))
  ALLOCATE(y(g(2)%size,g(2)%inb_grid))
  ALLOCATE(z(g(3)%size,g(3)%inb_grid))

  ALLOCATE(wrk1d(isize_wrk1d*inb_wrk1d))
  ALLOCATE(wrk2d(isize_wrk2d*inb_wrk2d))

  IF ( icalc_flow .EQ. 1 ) ALLOCATE(q(isize_field,inb_flow_array))
  IF ( icalc_scal .EQ. 1 ) ALLOCATE(s(isize_field,inb_scal_array))

  inb_txc = 3
  isize_txc = isize_txc_field*inb_txc
  IF ( inb_txc .GT. 0 ) ALLOCATE(txc(isize_txc_field,inb_txc))

  isize_wrk3d = isize_txc_field
  ALLOCATE(wrk3d(isize_wrk3d))

! -------------------------------------------------------------------
! Read the grid 
! -------------------------------------------------------------------
#include "dns_read_grid.h"

  IF ( ifourier .EQ. 1 ) THEN
     CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  IF ( icalc_flow .EQ. 1 ) THEN ! we need array space
     CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)
  ENDIF

  IF ( jmax_total .EQ. 1 ) THEN; ifourier_type = 2; 
  ELSE;                          ifourier_type = 3; ENDIF

! ###################################################################
  IF ( icalc_flow .EQ. 1 ) THEN
     DO iq = 1,3
        IF ( flag_type .EQ. 1 ) CALL RAND_PDF(imax,jmax,kmax, seed, isymmetric, ipdf, txc(1,2))
        
        IF ( ispectrum .GT. 0 ) THEN
           IF ( flag_type .EQ. 1 ) CALL OPR_FOURIER_F(ifourier_type, imax,jmax,kmax, txc(1,2),txc(1,1), txc(1,3),wrk2d,wrk3d)
           
           CALL RAND_PSD(imax,jmax,kmax, ispectrum, spc_param, flag_type, seed, txc(1,1))
           CALL OPR_FOURIER_B(ifourier_type, imax,jmax,kmax, txc(1,1), txc(1,2), wrk3d)
        ENDIF

! Normalize variance
        dummy = AVG1V3D(imax,jmax,kmax, i1, txc(1,2))
        txc(1:isize_field,2) = txc(1:isize_field,2) - dummy
        dummy = AVG1V3D(imax,jmax,kmax, i2, txc(1,2)); dummy = SQRT(ucov(iq)/dummy)
        q(1:isize_field,iq)  = txc(1:isize_field,2) * dummy

     ENDDO
     
     IF ( ipdf .EQ. 2 ) THEN ! Gaussian PDF
        CALL RAND_COVARIANCE(imax,jmax,kmax, q(:,1),q(:,2),q(:,3), ucov)
     ENDIF
     
     IF ( imode_eqns .EQ. DNS_EQNS_TOTAL .OR. imode_eqns .EQ. DNS_EQNS_INTERNAL ) THEN
        q(:,4) = pbg%mean; q(:,5) = rbg%mean
     ENDIF
     
     CALL DNS_WRITE_FIELDS('flow.rand', i2, imax,jmax,kmax, inb_flow, isize_field, q, txc)
     
  ENDIF
  
  IF ( icalc_scal .EQ. 1 ) THEN
     DO is = 1,inb_scal
        IF ( flag_type .EQ. 1 ) CALL RAND_PDF(imax,jmax,kmax, seed, isymmetric, ipdf, txc(1,2))
        
        IF ( ispectrum .GT. 0 ) THEN
           IF ( flag_type .EQ. 1 ) CALL OPR_FOURIER_F(ifourier_type, imax,jmax,kmax, txc(1,2),txc(1,1), txc(1,3),wrk2d,wrk3d)
           
           CALL RAND_PSD(imax,jmax,kmax, ispectrum, spc_param, flag_type, seed, txc(1,1))
           CALL OPR_FOURIER_B(ifourier_type, imax,jmax,kmax, txc(1,1), txc(1,2), wrk3d)
        ENDIF

! Normalize variance
        dummy = AVG1V3D(imax,jmax,kmax, i1, txc(1,2))
        txc(1:isize_field,2) = txc(1:isize_field,2) - dummy
        dummy = AVG1V3D(imax,jmax,kmax, i2, txc(1,2)); dummy = SQRT(ucov(is)/dummy)
        s(1:isize_field,is)  = txc(1:isize_field,2) * dummy

     ENDDO
     
     CALL DNS_WRITE_FIELDS('scal.rand', i1, imax,jmax,kmax, inb_scal, isize_field, s, txc)
     
  ENDIF
  
  CALL DNS_END(0)

  STOP
END PROGRAM INIRAND
