#include "types.h"

#define LOC_UNIT_ID i55
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_VISUALS(fname, iformat, nx,ny,nz, nfield, subdomain, field, txc)

  USE TLAB_TYPES,  ONLY : subarray_dt
  USE TLAB_VARS, ONLY : g, isize_txc_field
#ifdef USE_MPI
  USE TLAB_MPI_VARS,    ONLY : ims_pro
  USE TLAB_MPI_PROCS
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif
#include "integers.h"

  TINTEGER iformat, nx,ny,nz, nfield, subdomain(6)
  TREAL, DIMENSION(isize_txc_field,nfield) :: field
  TREAL, DIMENSION(nx*ny*nz,2) :: txc
  CHARACTER*(*) fname

! -------------------------------------------------------------------
  TINTEGER iaux_loc, sizes(5), nx_aux,ny_aux,nz_aux, ifield
  CHARACTER*32 varname(16), name
  TINTEGER iflag_mode

! ###################################################################
  sizes(5) = nfield

  nx_aux = subdomain(2)-subdomain(1)+1
  ny_aux = subdomain(4)-subdomain(3)+1
  nz_aux = subdomain(6)-subdomain(5)+1

  iflag_mode = 0 ! default
  sizes(1) = isize_txc_field     ! array size
  sizes(2) = 1                   ! lower bound
  IF      ( subdomain(2)-subdomain(1)+1 .EQ. g(1)%size .AND. &
            subdomain(6)-subdomain(5)+1 .EQ. 1          ) THEN! xOy plane
     iflag_mode = 1
     sizes(3)   = ny_aux *nx     ! upper bound
     sizes(4)   = 1              ! stride

  ELSE IF ( subdomain(6)-subdomain(5)+1 .EQ. g(3)%size .AND. &
            subdomain(2)-subdomain(1)+1 .EQ. 1          ) THEN! zOy plane
     iflag_mode = 2
     sizes(3)   = ny_aux *nx *nz ! upper bound
     sizes(4)   = nx             ! stride

  ELSE IF ( subdomain(2)-subdomain(1)+1 .EQ. g(1)%size .AND. &
            subdomain(6)-subdomain(5)+1 .EQ. g(3)%size ) THEN
     iflag_mode = 3                                           ! xOy blocks
     sizes(3)   = ny_aux *nx *nz ! upper bound
     sizes(4)   = 1              ! stride

  ENDIF

! ###################################################################
  IF      ( iformat .EQ. 0 ) THEN ! standard scalar format
     iaux_loc = nx*ny*nz
     CALL DNS_WRITE_FIELDS(fname, i1, nx,ny,nz, nfield, iaux_loc, field, txc)

! -------------------------------------------------------------------
  ELSE IF ( iformat .EQ. 1 ) THEN  ! ensight; to be removed
     CALL ENSIGHT_FIELD(fname, i1, nx,ny,nz, nfield, subdomain, field, txc)

! -------------------------------------------------------------------
  ELSE IF ( iformat .EQ. 2 .AND. iflag_mode .GT. 0 ) THEN  ! single precision, using MPI_IO
     IF ( ny_aux .NE. ny ) THEN
        DO ifield = 1,nfield
           CALL REDUCE_BLOCK_INPLACE(nx,ny,nz, i1,subdomain(3),i1, nx,ny_aux,nz, field(1,ifield), txc)
        ENDDO
     ENDIF

     varname = ''
     IF ( nfield .GT. 1 ) THEN
        DO ifield = 1,nfield; WRITE(varname(ifield),*) ifield; varname(ifield) = TRIM(ADJUSTL(varname(ifield))); ENDDO
     ENDIF
     CALL IO_WRITE_SUBARRAY4(iflag_mode, fname, varname, field, sizes, txc)

! -------------------------------------------------------------------
  ELSE                                                     ! single precision, through PE0
     DO ifield = 1,nfield
        IF ( nfield .GT. 1 ) THEN
           WRITE(name,*) ifield; name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(name))
        ELSE
           name = fname
        ENDIF

#ifdef USE_MPI
        IF ( ims_pro .EQ. 0 ) THEN
#endif
#include "dns_open_file.h"

#ifdef USE_MPI
        ENDIF
        CALL DNS_MPI_WRITE_PE0_SINGLE(LOC_UNIT_ID, nx,ny,nz, subdomain, field(1,ifield), txc(1,1), txc(1,2))
        IF ( ims_pro .EQ. 0 ) THEN
#else
           CALL REDUCE_BLOCK_INPLACE(nx,ny,nz, subdomain(1),subdomain(3),subdomain(5), nx_aux,ny_aux,nz_aux, field(1,ifield), txc)
           WRITE(LOC_UNIT_ID) SNGL(field(1:nx_aux*ny_aux*nz_aux,ifield))
#endif
           CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
        ENDIF
#endif
     ENDDO

  ENDIF

  RETURN
END SUBROUTINE IO_WRITE_VISUALS

!########################################################################
! Writing data in Ensight Gold Variable File Format
!########################################################################
SUBROUTINE ENSIGHT_FIELD(name, iheader, nx,ny,nz, nfield, subdomain, field, tmp_mpi)

#ifdef USE_MPI
  USE TLAB_MPI_VARS, ONLY : ims_pro
  USE TLAB_MPI_PROCS
#endif

  implicit NONE

#include "integers.h"

  CHARACTER*(*) name
  TINTEGER, INTENT(IN) :: iheader ! 0 no header; 1 header

  TINTEGER, INTENT(IN) :: nx,ny,nz, nfield, subdomain(6)
  TREAL, DIMENSION(nx,ny,nz,nfield) :: field
  TREAL, DIMENSION(nx*ny*nz,2     ) :: tmp_mpi

! -------------------------------------------------------------------
  CHARACTER*80 line
#ifdef USE_MPI
  TINTEGER ifield
#else
  TINTEGER i,j,k, ifield
#endif

! ###################################################################
! Header in Ensight Gold Variable File Format
! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
#include "dns_open_file.h"

  IF ( iheader .EQ. 1 ) THEN
     line='description line              '
     WRITE(LOC_UNIT_ID) line
     line='part                          '
     WRITE(LOC_UNIT_ID) line
     WRITE(LOC_UNIT_ID) i1
     line='block                         '
     WRITE(LOC_UNIT_ID) line
  ENDIF

#ifdef USE_MPI
  ENDIF
#endif

! ###################################################################
! Body
! ###################################################################
! -------------------------------------------------------------------
! parallel
! -------------------------------------------------------------------
#ifdef USE_MPI
  DO ifield = 1,nfield
     CALL DNS_MPI_WRITE_PE0_SINGLE(LOC_UNIT_ID, nx,ny,nz, subdomain, field, tmp_mpi(1,1), tmp_mpi(1,2))
  END DO

! -------------------------------------------------------------------
! serial
! -------------------------------------------------------------------
#else
  DO ifield = 1,nfield
     DO k = subdomain(5),subdomain(6)
        DO j = subdomain(3),subdomain(4)
           WRITE(LOC_UNIT_ID) (SNGL(field(i,j,k,ifield)),i=subdomain(1),subdomain(2))
        ENDDO
     ENDDO
  END DO

#endif

! ###################################################################
! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif
     CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
  ENDIF
#endif

  RETURN
END SUBROUTINE ENSIGHT_FIELD

!########################################################################
! Writing data in Ensight Gold Geometry File Format
! Note that record-length information has to be avoided
!########################################################################
SUBROUTINE ENSIGHT_GRID(name, nx,ny,nz, subdomain, x,y,z)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER nx,ny,nz, subdomain(6)
  TREAL x(nx), y(ny), z(nz)
  CHARACTER*(*) name

! -------------------------------------------------------------------
  CHARACTER*80 line
  TINTEGER ij

! ###################################################################
#include "dns_open_file.h"

  line='Fortran Binary                '
  WRITE(LOC_UNIT_ID) line

  line='description line 1            '
  WRITE(LOC_UNIT_ID) line
  line='description line 2            '
  WRITE(LOC_UNIT_ID) line
  line='node id off                   '
  WRITE(LOC_UNIT_ID) line
  line='element id off                '
  WRITE(LOC_UNIT_ID) line

  line='part                          '
  WRITE(LOC_UNIT_ID) line
  WRITE(LOC_UNIT_ID) i1
  line='description line              '
  WRITE(LOC_UNIT_ID) line
  line='block rectilinear             '
  WRITE(LOC_UNIT_ID) line
  WRITE(LOC_UNIT_ID) subdomain(2)-subdomain(1)+1,subdomain(4)-subdomain(3)+1,subdomain(6)-subdomain(5)+1
  WRITE(LOC_UNIT_ID) (SNGL(x(ij)),ij=subdomain(1),subdomain(2))
  WRITE(LOC_UNIT_ID) (SNGL(y(ij)),ij=subdomain(3),subdomain(4))
  WRITE(LOC_UNIT_ID) (SNGL(z(ij)),ij=subdomain(5),subdomain(6))

  CLOSE(LOC_UNIT_ID)

  RETURN
END SUBROUTINE ENSIGHT_GRID


! ###################################################################
! ###################################################################
#ifdef USE_MPI

SUBROUTINE VISUALS_MPIO_AUX(opt_format, subdomain)

  USE TLAB_VARS, ONLY : imax,kmax, io_aux
  USE TLAB_MPI_VARS

  IMPLICIT NONE

#include "mpif.h"

  TINTEGER,                 INTENT(IN)  :: opt_format, subdomain(6)

! -----------------------------------------------------------------------
  TINTEGER                :: ndims
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset

! #######################################################################
  io_aux(:)%active = .FALSE.
  io_aux(:)%offset = 0
  IF ( opt_format .EQ. 1 ) io_aux(:)%offset = 244 ! # bytes of ensight header

! ###################################################################
! Saving full vertical xOy planes; using subdomain(5) to define the plane
  IF ( ims_pro_k .EQ. ( (subdomain(5)-1) /kmax) ) io_aux(1)%active = .TRUE.
  io_aux(1)%communicator = ims_comm_x

  ndims = 2
  sizes(1)   = imax *ims_npro_i; sizes(2)   = subdomain(4)-subdomain(3)+1
  locsize(1) = imax;             locsize(2) = subdomain(4)-subdomain(3)+1
  offset(1)  = ims_offset_i;     offset(2)  = 0

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
       MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(1)%subarray, ims_err)
  CALL MPI_Type_commit(io_aux(1)%subarray, ims_err)

! Saving full vertical zOy planes; using subdomain(1) to define the plane
  IF ( ims_pro_i .EQ.  ( (subdomain(1)-1) /imax) ) io_aux(2)%active = .TRUE.
  io_aux(2)%communicator = ims_comm_z

  ndims = 2
                             sizes(1)   = subdomain(4)-subdomain(3)+1; sizes(2)   = kmax *ims_npro_k
                             locsize(1) = subdomain(4)-subdomain(3)+1; locsize(2) = kmax
                             offset(1)  = 0;                           offset(2)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
       MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(2)%subarray, ims_err)
  CALL MPI_Type_commit(io_aux(2)%subarray, ims_err)

! Saving full blocks xOz planes
  io_aux(3)%active = .TRUE.
  io_aux(3)%communicator = MPI_COMM_WORLD

  ndims = 3
  sizes(1)   = imax *ims_npro_i; sizes(2)   = subdomain(4)-subdomain(3)+1; sizes(3)   = kmax *ims_npro_k
  locsize(1) = imax;             locsize(2) = subdomain(4)-subdomain(3)+1; locsize(3) = kmax
  offset(1)  = ims_offset_i;     offset(2)  = 0;                           offset(3)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, &
       MPI_ORDER_FORTRAN, MPI_REAL4, io_aux(3)%subarray, ims_err)
  CALL MPI_Type_commit(io_aux(3)%subarray, ims_err)

  RETURN
END SUBROUTINE VISUALS_MPIO_AUX

#endif
