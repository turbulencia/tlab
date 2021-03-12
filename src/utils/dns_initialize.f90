#include "types.h"
#include "dns_error.h"

SUBROUTINE DNS_INITIALIZE
#ifdef USE_OPENMP
  USE OMP_LIB
#endif 
       
  USE DNS_CONSTANTS, ONLY : efile, lfile
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
  CHARACTER*64 line

!########################################################################
  imode_verbosity = 1 ! default value; needed already here

  istat_maj_ver = 1
  istat_min_ver = 0

! -------------------------------------------------------------------
! Inititalize MPI parallel mode
! -------------------------------------------------------------------
#ifdef USE_MPI 
#ifdef USE_PSFFT 
  call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,ims_nb_thrsupp_provided,ims_err)
  IF ( ims_nb_thrsupp_provided < MPI_THREAD_SERIALIZED ) THEN  
     CALL IO_WRITE_ASCII(efile,"DNS_INITIALIZE. MPI Library is missing the needed level of thread support for nonblocking communication.") 
     CALL IO_WRITE_ASCII(efile,"DNS_INITIALIZE. Try MPI_THREAD_FUNNELED after reading the documentation.")
     CALL DNS_STOP(DNS_ERROR_PSFFT) 
  ENDIF
#else 
  call MPI_INIT(ims_err)
#endif 

  call MPI_COMM_SIZE(MPI_COMM_WORLD,ims_npro,ims_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro, ims_err)

  ims_time_min = MPI_WTIME()

  ims_time_trans = C_0_R

#endif

!########################################################################
! First output
!########################################################################
  CALL DATE_AND_TIME(clock(1),clock(2))
  line='Starting on '//TRIM(ADJUSTL(clock(1)(1:8)))//' at '//TRIM(ADJUSTL(clock(2)))
  CALL IO_WRITE_ASCII(lfile,line)

  line='Git-hash '//GITHASH
  CALL IO_WRITE_ASCII(lfile,line)

#ifdef USE_MPI
  WRITE(line,*) ims_npro
  line='Number of MPI tasks '//TRIM(ADJUSTL(line))
  CALL IO_WRITE_ASCII(lfile,line)

  IF ( ims_npro .EQ. 0 ) THEN
     CALL IO_WRITE_ASCII(efile, 'DNS_INITIALIZE. Number of processors is zero.')
     CALL DNS_STOP(DNS_ERROR_MINPROC)
  ENDIF

#endif

! -------------------------------------------------------------------
! Inititalize OpenMP mode
! -------------------------------------------------------------------
#ifdef USE_OPENMP
  dns_omp_numThreads = omp_get_max_threads() 
  WRITE(line,*) dns_omp_numThreads
  line='Number of OMP threads '//TRIM(ADJUSTL(line))
  CALL IO_WRITE_ASCII(lfile,line)

#else
  dns_omp_numThreads = 1

#endif

  RETURN
END SUBROUTINE DNS_INITIALIZE
