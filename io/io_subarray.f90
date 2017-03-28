#include "types.h" 

#define SIZEOFBYTE 1

SUBROUTINE IO_WRITE_SUBARRAY4(iflag_mode, fname, varname, data, sizes, work)

  USE DNS_TYPES,     ONLY : subarray_dt
  USE DNS_CONSTANTS, ONLY : lfile
#ifdef USE_MPI
  USE DNS_MPI,       ONLY : mpio_aux
#else
  USE DNS_GLOBAL,    ONLY : io_aux
#endif

  IMPLICIT NONE

#ifdef USE_MPI 
#include "mpif.h"
#endif 

  TINTEGER,                                   INTENT(IN)    :: iflag_mode
  CHARACTER*(*),                              INTENT(IN)    :: fname
  TINTEGER,                                   INTENT(IN)    :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
  CHARACTER*32, DIMENSION(sizes(5)),          INTENT(IN)    :: varname
  TREAL,        DIMENSION(sizes(1),sizes(5)), INTENT(IN)    :: data
  REAL(4),      DIMENSION(sizes(1)),          INTENT(INOUT) :: work
  
! -----------------------------------------------------------------------
  TINTEGER iv, isize
  CHARACTER*64 name

#ifdef USE_MPI
  TINTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh, ims_err
#else
  TINTEGER :: ioffset_local
#endif

! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

  isize = ( sizes(3) -sizes(2) ) /sizes(4) +1

#ifdef USE_MPI
  IF ( mpio_aux(iflag_mode)%active ) THEN
#endif

  DO iv = 1,sizes(5)
     name = TRIM(ADJUSTL(fname))
     IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))

     CALL IO_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')

     work(1:isize) = SNGL(data(sizes(2):sizes(3):sizes(4),iv))

#ifdef USE_MPI
     CALL MPI_File_open(mpio_aux(iflag_mode)%communicator, TRIM(ADJUSTL(name)), IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,mpio_fh, ims_err)
     CALL MPI_File_set_view(mpio_fh, mpio_aux(iflag_mode)%offset, MPI_REAL4, mpio_aux(iflag_mode)%subarray, 'native', MPI_INFO_NULL, ims_err) 
     CALL MPI_File_write_all(mpio_fh, work, isize, MPI_REAL4, mpio_status, ims_err) 
     CALL MPI_File_close(mpio_fh, ims_err)  
     
#else
#include "dns_open_file.h"
     ioffset_local = io_aux(iflag_mode)%offset + 1
     WRITE(LOC_UNIT_ID,POS=ioffset_local) work(1:isize)
     CLOSE(LOC_UNIT_ID)
     
#endif
     
  ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_WRITE_SUBARRAY4

!########################################################################
!########################################################################
SUBROUTINE IO_READ_SUBARRAY8(iflag_mode, fname, varname, data, sizes, work)

  USE DNS_TYPES,     ONLY : subarray_dt
  USE DNS_CONSTANTS, ONLY : lfile
#ifdef USE_MPI
  USE DNS_MPI,       ONLY : mpio_aux
#else
  USE DNS_GLOBAL,    ONLY : io_aux
#endif

  IMPLICIT NONE

#ifdef USE_MPI 
#include "mpif.h"
#endif 

  TINTEGER,                                   INTENT(IN)    :: iflag_mode
  CHARACTER*(*),                              INTENT(IN)    :: fname
  TINTEGER,                                   INTENT(IN)    :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
  CHARACTER*32, DIMENSION(sizes(5)),          INTENT(IN)    :: varname
  TREAL,        DIMENSION(sizes(1),sizes(5)), INTENT(OUT)   :: data
  TREAL,        DIMENSION(sizes(1)),          INTENT(INOUT) :: work
  
! -----------------------------------------------------------------------
  TINTEGER iv, isize
  CHARACTER*64 name

#ifdef USE_MPI
  TINTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh, ims_err
#else
  TINTEGER :: ioffset_local
#endif

! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

  isize = ( sizes(3) -sizes(2) ) /sizes(4) +1

#ifdef USE_MPI
  IF ( mpio_aux(iflag_mode)%active ) THEN
#endif

  DO iv = 1,sizes(5)
     name = TRIM(ADJUSTL(fname))
     IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))

     CALL IO_WRITE_ASCII(lfile, 'Reading field '//TRIM(ADJUSTL(name))//'...')
     
#ifdef USE_MPI
     CALL MPI_File_open(mpio_aux(iflag_mode)%communicator, TRIM(ADJUSTL(name)), MPI_MODE_RDONLY,MPI_INFO_NULL,mpio_fh, ims_err)
     CALL MPI_File_set_view(mpio_fh, mpio_aux(iflag_mode)%offset, MPI_REAL8, mpio_aux(iflag_mode)%subarray, 'native', MPI_INFO_NULL, ims_err) 
     CALL MPI_File_read_all(mpio_fh, work, isize, MPI_REAL8, mpio_status, ims_err) 
     CALL MPI_File_close(mpio_fh, ims_err)  
     
#else
#include "dns_open_file.h"
     ioffset_local = io_aux(iflag_mode)%offset + 1
     READ(LOC_UNIT_ID,POS=ioffset_local) work(1:isize)
     CLOSE(LOC_UNIT_ID)
     
#endif
     
     data(sizes(2):sizes(3):sizes(4),iv) = work(1:isize)

  ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_READ_SUBARRAY8

!########################################################################
!########################################################################
SUBROUTINE IO_WRITE_SUBARRAY8(iflag_mode, fname, varname, data, sizes, work)

  USE DNS_TYPES,     ONLY : subarray_dt
  USE DNS_CONSTANTS, ONLY : lfile
#ifdef USE_MPI
  USE DNS_MPI,       ONLY : mpio_aux
#else
  USE DNS_GLOBAL,    ONLY : io_aux
#endif

  IMPLICIT NONE

#ifdef USE_MPI 
#include "mpif.h"
#endif 

  TINTEGER,                                   INTENT(IN)    :: iflag_mode
  CHARACTER*(*),                              INTENT(IN)    :: fname
  TINTEGER,                                   INTENT(IN)    :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
  CHARACTER*32, DIMENSION(sizes(5)),          INTENT(IN)    :: varname
  TREAL,        DIMENSION(sizes(1),sizes(5)), INTENT(IN)    :: data
  TREAL,        DIMENSION(sizes(1)),          INTENT(INOUT) :: work
  
! -----------------------------------------------------------------------
  TINTEGER iv, isize
  CHARACTER*64 name

#ifdef USE_MPI
  TINTEGER :: mpio_status(MPI_STATUS_SIZE), mpio_fh, ims_err
#else
  TINTEGER :: ioffset_local
#endif

! #######################################################################
#define LOC_UNIT_ID 75
#define LOC_STATUS 'unknown'

  isize = ( sizes(3) -sizes(2) ) /sizes(4) +1

#ifdef USE_MPI
  IF ( mpio_aux(iflag_mode)%active ) THEN
#endif

  DO iv = 1,sizes(5)
     name = TRIM(ADJUSTL(fname))
     IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))

     CALL IO_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')

     work(1:isize) = data(sizes(2):sizes(3):sizes(4),iv)

#ifdef USE_MPI
     CALL MPI_File_open(mpio_aux(iflag_mode)%communicator, TRIM(ADJUSTL(name)), IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE),MPI_INFO_NULL,mpio_fh, ims_err)
     CALL MPI_File_set_view(mpio_fh, mpio_aux(iflag_mode)%offset, MPI_REAL8, mpio_aux(iflag_mode)%subarray, 'native', MPI_INFO_NULL, ims_err) 
     CALL MPI_File_write_all(mpio_fh, work, isize, MPI_REAL8, mpio_status, ims_err) 
     CALL MPI_File_close(mpio_fh, ims_err)  
     
#else
#include "dns_open_file.h"
     ioffset_local = io_aux(iflag_mode)%offset + 1
     WRITE(LOC_UNIT_ID,POS=ioffset_local) work(1:isize)
     CLOSE(LOC_UNIT_ID)
     
#endif
     
  ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE IO_WRITE_SUBARRAY8

