#include "types.h"
#include "dns_const.h"
#include "dns_error.h"

MODULE TLAB_PROCS
  USE TLAB_CONSTANTS, ONLY : lfile, efile
#ifdef USE_OPENMP
  USE OMP_LIB
#endif
#ifdef USE_MPI
  USE MPI
  USE TLAB_MPI_VARS
#endif
  IMPLICIT NONE

  SAVE

  PRIVATE

  CHARACTER*128 line

  PUBLIC :: TLAB_START
  PUBLIC :: TLAB_STOP
  PUBLIC :: TLAB_WRITE_ASCII
  PUBLIC :: TLAB_ALLOCATE
  PUBLIC :: TLAB_ALLOCATE_ARRAY1, TLAB_ALLOCATE_ARRAY1_INT, TLAB_ALLOCATE_ARRAY2
#ifdef USE_MPI
  PUBLIC :: TLAB_MPI_PANIC
#endif
CONTAINS

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TLAB_ALLOCATE(C_FILE_LOC)
    USE TLAB_VARS, ONLY : isize_field, inb_flow_array, inb_scal_array
    USE TLAB_VARS, ONLY : isize_txc_field, inb_txc
    USE TLAB_VARS, ONLY : isize_wrk1d, isize_wrk2d, isize_wrk3d, inb_wrk1d, inb_wrk2d
    USE TLAB_VARS, ONLY : g
    USE TLAB_ARRAYS
    IMPLICIT NONE

    CHARACTER(LEN=*) C_FILE_LOC

    ! -------------------------------------------------------------------
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,x,    g(1)%inb_grid, g(1)%size,  'x')
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,y,    g(2)%inb_grid, g(2)%size,  'y')
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,z,    g(3)%inb_grid, g(3)%size,  'z')
    ! -------------------------------------------------------------------
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,q,    inb_flow_array,isize_field,'flow') 
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,s,    inb_scal_array,isize_field,'scal') 
    ! -------------------------------------------------------------------
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,txc,  inb_txc,isize_txc_field,   'txc')
    CALL TLAB_ALLOCATE_ARRAY1(C_FILE_LOC,wrk3d,isize_wrk3d,               'wrk3d')
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,wrk1d,inb_wrk1d,isize_wrk1d,     'wrk1d') 
    CALL TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,wrk2d,inb_wrk2d,isize_wrk2d,     'wrk2d') 

    RETURN
  END SUBROUTINE TLAB_ALLOCATE

  ! ######################################################################
  ! ######################################################################
  SUBROUTINE TLAB_ALLOCATE_ARRAY1(C_FILE_LOC,a,i1,s)

    IMPLICIT NONE

    CHARACTER(LEN=*),                 INTENT(IN   ) :: C_FILE_LOC
    TREAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    TINTEGER,                         INTENT(IN   ) :: i1
    CHARACTER(LEN=*),                 INTENT(IN   ) :: s
    ! --------------------------------------------------------------------
    
    CHARACTER*128                                   :: str, line
    TINTEGER                                        :: ierr
    
    !#####################################################################
    IF ( i1 .LE. 0 ) RETURN 
    
    WRITE(str,*) i1; line = 'Allocating array '//TRIM(ADJUSTL(s))//' of size '//TRIM(ADJUSTL(str))
    CALL TLAB_WRITE_ASCII(lfile,line)
    ALLOCATE(a(i1),stat=ierr)
    IF ( ierr /= 0 ) THEN
       CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//TRIM(ADJUSTL(s)) //'.')
       CALL TLAB_STOP(DNS_ERROR_ALLOC)
    END IF
   
  END SUBROUTINE TLAB_ALLOCATE_ARRAY1

  ! ######################################################################
  ! ######################################################################
  SUBROUTINE TLAB_ALLOCATE_ARRAY1_INT(C_FILE_LOC,a,i1,s)

    IMPLICIT NONE

    CHARACTER(LEN=*),                    INTENT(IN   ) :: C_FILE_LOC
    TINTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    TINTEGER,                            INTENT(IN   ) :: i1
    CHARACTER(LEN=*),                    INTENT(IN   ) :: s
    ! --------------------------------------------------------------------
    
    CHARACTER*128                                      :: str, line
    TINTEGER                                           :: ierr
    
    !#####################################################################
    IF ( i1 .LE. 0 ) RETURN 
    
    WRITE(str,*) i1; line = 'Allocating array '//TRIM(ADJUSTL(s))//' of size '//TRIM(ADJUSTL(str))
    CALL TLAB_WRITE_ASCII(lfile,line)
    ALLOCATE(a(i1),stat=ierr)
    IF ( ierr /= 0 ) THEN
       CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//TRIM(ADJUSTL(s)) //'.')
       CALL TLAB_STOP(DNS_ERROR_ALLOC)
    END IF
   
  END SUBROUTINE TLAB_ALLOCATE_ARRAY1_INT

  ! ######################################################################
  ! ######################################################################
  SUBROUTINE TLAB_ALLOCATE_ARRAY2(C_FILE_LOC,a,i1,i2,s)

    IMPLICIT NONE

    CHARACTER(LEN=*),                   INTENT(IN   ) :: C_FILE_LOC
    TREAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a
    TINTEGER,                           INTENT(IN   ) :: i1, i2
    CHARACTER(LEN=*),                   INTENT(IN   ) :: s
    ! --------------------------------------------------------------------
    
    CHARACTER*128                                     :: str, line
    TINTEGER                                          :: ierr

    !#####################################################################
    IF ( i1 .LE. 0 .OR. i2 .LE. 0 ) RETURN 
    
    WRITE(str,*) i1; line = 'Allocating array '//TRIM(ADJUSTL(s))//' of size '//TRIM(ADJUSTL(str))//'x'
    WRITE(str,*) i2; line = TRIM(ADJUSTL(line))//TRIM(ADJUSTL(str))
    CALL TLAB_WRITE_ASCII(lfile,line)
    ALLOCATE(a(i2,i1),stat=ierr)
    IF ( ierr /= 0 ) THEN
       CALL TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Error while allocating memory space for '//TRIM(ADJUSTL(s))//'.')
       CALL TLAB_STOP(DNS_ERROR_ALLOC)
    END IF
    
  ENd SUBROUTINE TLAB_ALLOCATE_ARRAY2
  
  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TLAB_START()
    USE TLAB_VARS, ONLY : imode_verbosity
    USE TLAB_VARS, ONLY : dns_omp_numThreads

    IMPLICIT NONE

    CHARACTER*10 clock(2)

    !########################################################################
    imode_verbosity = 1 ! default value; needed already here

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
    END IF
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

    line='Git-hash   '//GITHASH
    CALL TLAB_WRITE_ASCII(lfile,line)

    line='Git-branch '//GITBRANCH
    CALL TLAB_WRITE_ASCII(lfile,line)

#ifdef USE_MPI
    WRITE(line,*) ims_npro
    line='Number of MPI tasks '//TRIM(ADJUSTL(line))
    CALL TLAB_WRITE_ASCII(lfile,line)

    IF ( ims_npro .EQ. 0 ) THEN
      CALL TLAB_WRITE_ASCII(efile, 'DNS_START. Number of processors is zero.')
      CALL TLAB_STOP(DNS_ERROR_MINPROC)
    END IF

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
    USE TLAB_VARS, ONLY : g
    USE TLAB_VARS, ONLY : ifourier, fft_plan_fx, fft_plan_bx, fft_plan_fz, fft_plan_bz
    USE TLAB_VARS, ONLY : ivfilter, fft_plan_fy1d, fft_plan_by1d
    IMPLICIT NONE

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
    IF (ivfilter .EQ. 1 ) THEN
      CALL dfftw_destroy_plan(fft_plan_fy1d)
      CALL dfftw_destroy_plan(fft_plan_by1d)
    ENDIF
#endif

    ! ###################################################################
    IF ( error_code /= 0 ) THEN
      WRITE(line,*) error_code
      line = 'Error code '//TRIM(ADJUSTL(line))//'.'
      CALL TLAB_WRITE_ASCII(efile,line)
    END IF

    CALL GETARG(0,line)
    WRITE(line,*) 'Finalizing program '//TRIM(ADJUSTL(line))
    IF ( error_code == 0 ) THEN
      line = TRIM(ADJUSTL(line))//' normally.'
    ELSE
      line = TRIM(ADJUSTL(line))//' abnormally. Check '//TRIM(ADJUSTL(efile))
    END IF
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
    IF ( ims_err .EQ. 0 ) THEN
       CALL MPI_FINALIZE(ims_err)
    ELSE
       CALL MPI_Abort(MPI_COMM_WORLD,ims_err,ims_err)
    ENDIF
#endif
    RETURN
  END SUBROUTINE TLAB_STOP

#ifdef USE_MPI 
  SUBROUTINE TLAB_MPI_PANIC(location,mpi_error_code)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: location
    INTEGER,          INTENT(IN) :: mpi_error_code

    !##############################
    CHARACTER error_string*1024, line*512
    INTEGER error_local, error_len

    CALL MPI_Error_String(mpi_error_code, error_string, error_len, error_local)
    CALL TLAB_WRITE_ASCII(efile,'MPI-ERROR: Source file'//TRIM(ADJUSTL(LOCATION)),.TRUE.)
    CALL TLAB_WRITE_ASCII(efile,error_string,.TRUE.)

    CALL TLAB_STOP(mpi_error_code)
    ! Not supposed to return from this subroutine

  END SUBROUTINE TLAB_MPI_PANIC
#endif

  ! ###################################################################
  ! ###################################################################
  SUBROUTINE TLAB_WRITE_ASCII(file, line, flag_all)
    USE TLAB_VARS, ONLY : imode_verbosity
    IMPLICIT NONE

    CHARACTER*(*),  INTENT(IN)            :: file, line
    LOGICAL,        INTENT(IN), OPTIONAL  :: flag_all

    ! -----------------------------------------------------------------------
    CHARACTER*10 clock(2)

    ! #######################################################################
#ifdef USE_MPI
    IF ( ims_pro == 0 .OR. PRESENT(flag_all) ) THEN
#endif

      IF ( imode_verbosity > 0 ) THEN

        OPEN(UNIT=22, FILE=file, STATUS='unknown',POSITION='APPEND')
        IF      ( imode_verbosity .EQ. 1 ) THEN
          WRITE(22,'(a)') TRIM(ADJUSTL(line))
        ELSE IF ( imode_verbosity .EQ. 2 ) THEN
          CALL DATE_AND_TIME(clock(1),clock(2))
          WRITE(22,'(a)') '['//TRIM(ADJUSTR(clock(2)))//'] '//TRIM(ADJUSTL(line))
        END IF
        CLOSE(22)

      END IF

#ifndef PARALLEL
      IF ( file .EQ. efile ) THEN
        WRITE(*,*) TRIM(ADJUSTL(line))
      END IF
#endif

#ifdef USE_MPI
    END IF
#endif

    RETURN
  END SUBROUTINE TLAB_WRITE_ASCII

END MODULE TLAB_PROCS
