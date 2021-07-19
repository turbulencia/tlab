#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

MODULE TLAB_ARRAYS
  IMPLICIT NONE
  SAVE
  PRIVATE

  TREAL, ALLOCATABLE, PUBLIC :: x(:,:),y(:,:),z(:,:)  ! Grid and associated arrays
  TREAL, ALLOCATABLE, PUBLIC :: q(:,:)                ! Eulerian fields, flow vartiables
  TREAL, ALLOCATABLE, PUBLIC :: s(:,:)                ! Eulerian fields, scalar variables
  TREAL, ALLOCATABLE, PUBLIC :: txc(:,:)              ! Temporary space for Eulerian fields
  TREAL, ALLOCATABLE, PUBLIC :: wrk1d(:,:)            ! Work arrays (scratch space)
  TREAL, ALLOCATABLE, PUBLIC :: wrk2d(:,:)            ! Work arrays (scratch space)
  TREAL, ALLOCATABLE, PUBLIC :: wrk3d(:)              ! Work arrays (scratch space)

  TARGET x, y, z
  TARGET q, s, txc

END MODULE TLAB_ARRAYS

MODULE TLAB_MAIN
  USE DNS_CONSTANTS, ONLY : lfile, efile
  USE TLAB_ARRAYS
  IMPLICIT NONE
  SAVE
  PRIVATE

  CHARACTER*128 str, line

  PUBLIC :: TLAB_ALLOCATE

CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TLAB_ALLOCATE(C_FILE_LOC)
    USE DNS_GLOBAL, ONLY : isize_field, inb_flow_array, inb_scal_array
    USE DNS_GLOBAL, ONLY : isize_txc_field, inb_txc
    USE DNS_GLOBAL, ONLY : isize_wrk1d, isize_wrk2d, isize_wrk3d, inb_wrk1d, inb_wrk2d
    USE DNS_GLOBAL, ONLY : g

    IMPLICIT NONE

    CHARACTER(LEN=*) C_FILE_LOC

    ! -------------------------------------------------------------------
    TINTEGER ierr

    ! ###################################################################
    WRITE(str,*) g(1)%inb_grid; line = 'Allocating array x of size '//TRIM(ADJUSTL(str))//'x'
    WRITE(str,*) g(1)%size; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    CALL TLAB_WRITE_ASCII(lfile,line)
    ALLOCATE(x(g(1)%size,g(1)%inb_grid),stat=ierr)
    IF ( ierr .NE. 0 ) THEN
      CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for x.')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF

    WRITE(str,*) g(2)%inb_grid; line = 'Allocating array y of size '//TRIM(ADJUSTL(str))//'x'
    WRITE(str,*) g(2)%size; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    CALL TLAB_WRITE_ASCII(lfile,line)
    ALLOCATE(y(g(2)%size,g(2)%inb_grid),stat=ierr)
    IF ( ierr .NE. 0 ) THEN
      CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for y.')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF

    WRITE(str,*) g(3)%inb_grid; line = 'Allocating array z of size '//TRIM(ADJUSTL(str))//'x'
    WRITE(str,*) g(3)%size; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    CALL TLAB_WRITE_ASCII(lfile,line)
    ALLOCATE(z(g(3)%size,g(3)%inb_grid),stat=ierr)
    IF ( ierr .NE. 0 ) THEN
      CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for z.')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF

    ! -------------------------------------------------------------------
    IF ( inb_flow_array .GT. 0 ) THEN
      WRITE(str,*) inb_flow_array; line = 'Allocating array flow  of size '//TRIM(ADJUSTL(str))//'x'
      WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
      CALL TLAB_WRITE_ASCII(lfile,line)
      ALLOCATE(q(isize_field,inb_flow_array),stat=ierr)
      IF ( ierr .NE. 0 ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for flow.')
        CALL TLAB_STOP(DNS_ERROR_ALLOC)
      ENDIF
    ENDIF

    IF ( inb_scal_array .GT. 0 ) THEN
      WRITE(str,*) inb_scal_array; line = 'Allocating array scal  of size '//TRIM(ADJUSTL(str))//'x'
      WRITE(str,*) isize_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
      CALL TLAB_WRITE_ASCII(lfile,line)
      ALLOCATE(s(isize_field,inb_scal_array),stat=ierr)
      IF ( ierr .NE. 0 ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for scal.')
        CALL TLAB_STOP(DNS_ERROR_ALLOC)
      ENDIF
    ENDIF

    ! -------------------------------------------------------------------
    IF ( inb_txc .GT. 0 ) THEN
      WRITE(str,*) inb_txc; line = 'Allocating array txc   of size '//TRIM(ADJUSTL(str))//'x'
      WRITE(str,*) isize_txc_field; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
      CALL TLAB_WRITE_ASCII(lfile,line)
      ALLOCATE(txc(isize_txc_field,inb_txc),stat=ierr)
      IF ( ierr .NE. 0 ) THEN
        CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for txc.')
        CALL TLAB_STOP(DNS_ERROR_ALLOC)
      ENDIF
    ENDIF

    WRITE(str,*) isize_wrk3d; line = 'Allocating array wrk3d of size '//TRIM(ADJUSTL(str))
    CALL TLAB_WRITE_ASCII(lfile,line)
    ALLOCATE(wrk3d(isize_wrk3d),stat=ierr)
    IF ( ierr .NE. 0 ) THEN
      CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for wrk3d.')
      CALL TLAB_STOP(DNS_ERROR_ALLOC)
    ENDIF

    ALLOCATE(wrk1d(isize_wrk1d,inb_wrk1d))
    ALLOCATE(wrk2d(isize_wrk2d,inb_wrk2d))

    RETURN
  END SUBROUTINE TLAB_ALLOCATE

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TLAB_START()
#ifdef USE_OPENMP
    USE OMP_LIB
#endif
    USE DNS_GLOBAL, ONLY : istat_min_ver, istat_maj_ver
    USE DNS_GLOBAL, ONLY : dns_omp_numThreads
    USE DNS_GLOBAL, ONLY : imode_verbosity
#ifdef USE_MPI
    USE DNS_MPI
#endif

    IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

    CHARACTER*10 clock(2)

    !########################################################################
    imode_verbosity = 1 ! default value; needed already here

    istat_maj_ver = 1
    istat_min_ver = 0

    ! -------------------------------------------------------------------
    ! Inititalize MPI parallel mode
    ! -------------------------------------------------------------------
#ifdef USE_MPI
#ifdef USE_PSFFT
    CALL MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,ims_nb_thrsupp_provided,ims_err)
    IF ( ims_nb_thrsupp_provided < MPI_THREAD_SERIALIZED ) THEN
      CALL TLAB_WRITE_ASCII(efile,"DNS_START. MPI Library is missing the needed level of thread support for nonblocking communication.")
      CALL TLAB_WRITE_ASCII(efile,"DNS_START. Try MPI_THREAD_FUNNELED after reading the documentation.")
      CALL TLAB_STOP(DNS_ERROR_PSFFT)
    ENDIF
#else
    CALL MPI_INIT(ims_err)
#endif

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ims_npro,ims_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro, ims_err)

    ims_time_min = MPI_WTIME()

    ims_time_trans = C_0_R

#endif

    !########################################################################
    ! First output
    !########################################################################
    CALL DATE_AND_TIME(clock(1),clock(2))
    line='Starting on '//TRIM(ADJUSTL(clock(1)(1:8)))//' at '//TRIM(ADJUSTL(clock(2)))
    CALL TLAB_WRITE_ASCII(lfile,line)

    line='Git-hash '//GITHASH
    CALL TLAB_WRITE_ASCII(lfile,line)

#ifdef USE_MPI
    WRITE(line,*) ims_npro
    line='Number of MPI tasks '//TRIM(ADJUSTL(line))
    CALL TLAB_WRITE_ASCII(lfile,line)

    IF ( ims_npro .EQ. 0 ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'DNS_START. Number of processors is zero.')
      CALL TLAB_STOP(DNS_ERROR_MINPROC)
    ENDIF

#endif

    ! -------------------------------------------------------------------
    ! Inititalize OpenMP mode
    ! -------------------------------------------------------------------
#ifdef USE_OPENMP
    dns_omp_numThreads = omp_get_max_threads()
    WRITE(line,*) dns_omp_numThreads
    line='Number of OMP threads '//TRIM(ADJUSTL(line))
    CALL TLAB_WRITE_ASCII(lfile,line)

#else
    dns_omp_numThreads = 1

#endif

    RETURN
  END SUBROUTINE TLAB_START

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TLAB_STOP(error_code)
    USE DNS_GLOBAL, ONLY : g
    USE DNS_GLOBAL, ONLY : ifourier, fft_plan_fx, fft_plan_bx, fft_plan_fz, fft_plan_bz
#ifdef USE_MPI
    USE DNS_MPI, ONLY : ims_time_min, ims_time_max, ims_time_trans, ims_err
#endif
    IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

    INTEGER, INTENT(IN) :: error_code

    ! ###################################################################
#ifdef USE_FFTW
    IF ( ifourier == 1 ) THEN
      CALL dfftw_destroy_plan(fft_plan_fx)
      CALL dfftw_destroy_plan(fft_plan_bx)
      IF ( g(3)%size > 1 ) THEN
        CALL dfftw_destroy_plan(fft_plan_fz)
        CALL dfftw_destroy_plan(fft_plan_bz)
      END IF
    END IF
#endif

    ! ###################################################################
    IF ( error_code /= 0 ) THEN
      WRITE(line,*) error_code
      line = 'Error code '//TRIM(ADJUSTL(line))//'.'
      CALL TLAB_WRITE_ASCII(efile,line)
    ENDIF

    CALL GETARG(0,line)
    WRITE(line,*) 'Finalizing program '//TRIM(ADJUSTL(line))
    IF ( error_code == 0 ) THEN
      line = TRIM(ADJUSTL(line))//' normally.'
    ELSE
      line = TRIM(ADJUSTL(line))//' abnormally. Check '//TRIM(ADJUSTL(efile))
    ENDIF
    CALL TLAB_WRITE_ASCII(lfile, line)

#ifdef USE_MPI
    ims_time_max = MPI_WTIME()
    WRITE(line,1000) ims_time_max-ims_time_min
    line = 'Time elapse ....................: '//TRIM(ADJUSTL(line))
    CALL TLAB_WRITE_ASCII(lfile, line)

#ifdef PROFILE_ON
    WRITE(line,1000) ims_time_trans
    line = 'Time in array transposition ....: '//TRIM(ADJUST(line))
    CALL TLAB_WRITE_ASCII(lfile, line)
#endif

1000 FORMAT(G_FORMAT_R)

#endif

    CALL TLAB_WRITE_ASCII(lfile, '########################################')

#ifdef USE_MPI
    CALL MPI_FINALIZE(ims_err)
#endif

    RETURN
  END SUBROUTINE TLAB_STOP

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TLAB_WRITE_ASCII(file, line)
    USE DNS_GLOBAL, ONLY : imode_verbosity
#ifdef USE_MPI
    USE DNS_MPI, ONLY : ims_pro
#endif
    IMPLICIT NONE

    CHARACTER*(*), INTENT(IN) :: file, line

    ! -----------------------------------------------------------------------
    CHARACTER*10 clock(2)

    ! #######################################################################
#ifdef USE_MPI
    IF ( ims_pro .EQ. 0 ) THEN
#endif

      IF ( imode_verbosity .GT. 0 ) THEN

        OPEN(UNIT=22, FILE=file, STATUS='unknown',POSITION='APPEND')
        IF      ( imode_verbosity .EQ. 1 ) THEN
          WRITE(22,'(a)') TRIM(ADJUSTL(line))
        ELSE IF ( imode_verbosity .EQ. 2 ) THEN
          CALL DATE_AND_TIME(clock(1),clock(2))
          WRITE(22,'(a)') '['//TRIM(ADJUSTR(clock(2)))//'] '//TRIM(ADJUSTL(line))
        ENDIF
        CLOSE(22)

      ENDIF

#ifndef PARALLEL
      IF ( file .EQ. efile ) THEN
        WRITE(*,*) TRIM(ADJUSTL(line))
      ENDIF
#endif

#ifdef USE_MPI
    ENDIF
#endif

    RETURN
  END SUBROUTINE TLAB_WRITE_ASCII

END MODULE TLAB_MAIN
