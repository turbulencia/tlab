#include "types.h"
#include "dns_error.h"
#include "dns_const.h"
#include "dns_const_mpi.h"

#ifdef USE_PSFFT
#include "nb3dfft_defines.inc" 
#endif

! ######################################################################
! Initialization of PSFFT Library for nonblocking communication 
! ######################################################################
SUBROUTINE DNS_NB3DFFT_INITIALIZE

  USE DNS_CONSTANTS, ONLY : efile
#ifdef USE_PSFFT 
  USE DNS_CONSTANTS, ONLY : lfile
  USE DNS_GLOBAL, ONLY : imax,jmax,kmax
  USE DNS_GLOBAL, ONLY : g
  USE DNS_MPI
  USE NB3DFFT, ONLY : nb3dfft_test_setup, nb3dfft_setup, get_dims 
#endif 

  IMPLICIT NONE

! #######################################################################
#ifdef USE_PSFFT 
  CALL IO_WRITE_ASCII(lfile,'Initialize nonblocking communication.')

  ims_nb_proc_grid = (/ ims_npro_i, ims_npro_k /) 
  CALL NB3DFFT_SETUP(ims_nb_proc_grid, g(1)%size,g(2)%size,g(3)%size, &
       ims_nb_msize)

  CALL GET_DIMS(ims_nb_xsrt,ims_nb_xend,ims_nb_xsiz,1,1) 
  CALL GET_DIMS(ims_nb_ysrt,ims_nb_yend,ims_nb_ysiz,1,2) 
  CALL GET_DIMS(ims_nb_zsrt,ims_nb_zend,ims_nb_zsiz,1,3) 

  IF (       ims_nb_xsrt(1) .EQ. 1 .AND. ims_nb_xend(1) .EQ. g(1)%size &
       .AND. ims_nb_xsiz(2)*ims_nb_xsiz(3) .EQ. ims_size_i(DNS_MPI_I_PARTIAL) ) THEN  
     ! Decomp standing in X okay
  ELSE  
     CALL IO_WRITE_ASCII(efile,'Decomp standing in X-BAD') 
     CALL DNS_STOP(DNS_ERROR_PARPARTITION)
  ENDIF

  IF (      ims_nb_ysrt(1) .EQ. ims_offset_i+1 & 
       .AND.ims_nb_ysrt(2) .EQ. ims_offset_j+1 &
       .AND.ims_nb_ysrt(3) .EQ. ims_offset_k+1 &
       .AND.ims_nb_ysiz(1) .EQ. imax & 
       .AND.ims_nb_ysiz(2) .EQ. jmax &
       .AND.ims_nb_ysiz(3) .EQ. kmax ) THEN  
  ELSE 
     CALL IO_WRITE_ASCII(efile,'Decomp standing in Y--BAD')
     CALL DNS_STOP(DNS_ERROR_PARPARTITION)
  ENDIF

  IF (      ims_nb_zsrt(3) .EQ. 1 .AND. ims_nb_zend(3) .EQ. g(3)%size & 
       .AND.ims_nb_zsiz(1)*ims_nb_zsiz(2) .EQ. ims_size_k(DNS_MPI_K_PARTIAL) ) THEN 
     ! Decomp standing in Z okay 
  ELSE 
     CALL IO_WRITE_ASCII(efile, 'Decomp standing in Z--BAD') 
     CALL DNS_STOP(DNS_ERROR_PARPARTITION)
  ENDIF 

  CALL IO_WRITE_ASCII(lfile,'Checking that NB3DFFT and DNS domain decompositions agree.') 

  CALL nb3dfft_test_setup()

#else
  CALL IO_WRITE_ASCII(efile,'Compiler flag USE_PSFFT needs to be used.') 
  CALL DNS_STOP(DNS_ERROR_PARPARTITION)

#endif

  RETURN
END SUBROUTINE DNS_NB3DFFT_INITIALIZE
