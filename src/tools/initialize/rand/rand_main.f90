#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

#define C_FILE_LOC "INIRAND"

PROGRAM INIRAND

  USE TLAB_CONSTANTS
  USE TLAB_VARS
  USE TLAB_ARRAYS
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_PROCS
#endif
  USE RAND_LOCAL
#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_pro
#endif

  IMPLICIT NONE

  ! -------------------------------------------------------------------
  TINTEGER iq, is

  ! ###################################################################
  CALL TLAB_START()

  CALL DNS_READ_GLOBAL(ifile)
  CALL RAND_READ_LOCAL(ifile)
#ifdef USE_MPI
  CALL TLAB_MPI_INITIALIZE
#endif

  isize_wrk3d = isize_txc_field

  inb_txc = 3

  CALL TLAB_ALLOCATE(C_FILE_LOC)

  CALL IO_READ_GRID(gfile, g(1)%size,g(2)%size,g(3)%size, g(1)%scale,g(2)%scale,g(3)%scale, x,y,z, area)
  CALL FDM_INITIALIZE(x, g(1), wrk1d)
  CALL FDM_INITIALIZE(y, g(2), wrk1d)
  CALL FDM_INITIALIZE(z, g(3), wrk1d)

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Initializing random fiels.')

#ifdef USE_MPI
  seed = seed + ims_pro         ! seed for random generator
#endif
  seed = - ABS(seed)

  IF ( ifourier .EQ. 1 ) THEN
    CALL OPR_FOURIER_INITIALIZE(txc, wrk1d,wrk2d,wrk3d)
  ENDIF

  CALL OPR_CHECK(imax,jmax,kmax, q, txc, wrk2d,wrk3d)

  itime = 0; rtime = C_0_R

  DO iq = 1,inb_flow
    CALL RAND_FIELD( ucov(iq), q(1,iq), txc(1,1), txc(1,2), txc(1,3), wrk2d,wrk3d )
  ENDDO
  IF ( ipdf .EQ. 2 ) THEN ! Gaussian PDF
    CALL RAND_COVARIANCE(ucov, q(:,1),q(:,2),q(:,3))
  ENDIF
  CALL DNS_WRITE_FIELDS('flow.rand', i2, imax,jmax,kmax, inb_flow, isize_field, q, txc)

  DO is = 1,inb_scal
    CALL RAND_FIELD( ucov(is), s(1,is), txc(1,1), txc(1,2), txc(1,3), wrk2d,wrk3d )
  ENDDO
  CALL DNS_WRITE_FIELDS('scal.rand', i1, imax,jmax,kmax, inb_scal, isize_field, s, txc)

  CALL TLAB_STOP(0)
END PROGRAM INIRAND
