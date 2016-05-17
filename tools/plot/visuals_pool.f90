#include "types.h"

#define LOC_UNIT_ID i55
#define LOC_STATUS 'unknown'

SUBROUTINE VISUALS_WRITE(fname, itype, iformat, nx,ny,nz, subdomain, field, tmp_mpi)

  USE DNS_TYPES,  ONLY : subarray_structure
  USE DNS_GLOBAL, ONLY : imax_total, kmax_total
#ifdef USE_MPI
  USE DNS_MPI,    ONLY : mpio_aux, ims_pro
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif
#include "integers.h"

  TINTEGER iformat, itype, nx,ny,nz, subdomain(6)
  TREAL, DIMENSION(nx,ny,nz,*) :: field
  TREAL, DIMENSION(nx*ny*nz,2) :: tmp_mpi
  CHARACTER*(*) fname

! -------------------------------------------------------------------
  TINTEGER iaux_loc, nfield, sizes(4), ny_aux, j,k, ifield
  CHARACTER*32 varname(16), name
  TINTEGER iflag_mode

#ifndef USE_MPI
  TYPE(subarray_structure) :: mpio_aux(10)
#endif

! ###################################################################
  IF      ( itype .EQ. 0 ) THEN; nfield = 1;
  ELSE IF ( itype .EQ. 1 ) THEN; nfield = 3;
  ELSE IF ( itype .EQ. 2 ) THEN; nfield = 6; ENDIF
  sizes(5) = nfield

! in-place reduction of data in array field along direction Oy 
! to be done
!  ny_aux = subdomain(4)-subdomain(3)+1
  ny_aux = ny

  iflag_mode = 0 ! default
  sizes(1) = ny_aux *nx *nz  ! array size
  sizes(2) = 1               ! lower bound
  IF      ( subdomain(2)-subdomain(1)+1 .EQ. imax_total .AND. &
            subdomain(6)-subdomain(5)+1 .EQ. 1          ) THEN! xOy plane
     iflag_mode = 1
     sizes(3)   = ny_aux *nx ! upper bound
     sizes(4)   = 1          ! stride
     
  ELSE IF ( subdomain(6)-subdomain(5)+1 .EQ. kmax_total .AND. & 
            subdomain(2)-subdomain(1)+1 .EQ. 1          ) THEN! zOy plane
     iflag_mode = 2
     sizes(3)   = sizes(1)   ! upper bound
     sizes(4)   = nx         ! stride
     
  ELSE IF ( subdomain(2)-subdomain(1)+1 .EQ. imax_total .AND. &
            subdomain(6)-subdomain(5)+1 .EQ. kmax_total ) THEN
     iflag_mode = 3                                           ! xOy blocks
     sizes(3)   = sizes(1)   ! upper bound
     sizes(4)   = 1          ! stride

  ENDIF

! ###################################################################
  IF      ( iformat .EQ. 0 ) THEN ! standard scalar format
     iaux_loc = nx*ny*nz
     CALL DNS_WRITE_FIELDS(fname, i1, nx,ny,nz, nfield, iaux_loc, field, tmp_mpi)

! -------------------------------------------------------------------
  ELSE IF ( iformat .EQ. 1 ) THEN  ! ensight; to be removed
     CALL ENSIGHT_FIELD(fname, i1, nx,ny_aux,nz, nfield, subdomain, field, tmp_mpi)
     
! -------------------------------------------------------------------
  ELSE IF ( iformat .EQ. 2 .AND. iflag_mode .GT. 0 ) THEN  ! single precision, using MPI_IO
     varname = ''
     IF ( nfield .GT. 1 ) THEN
        DO ifield = 1,nfield; WRITE(varname(ifield),*) ifield; varname(ifield) = TRIM(ADJUSTL(varname(ifield))); ENDDO
     ENDIF
     CALL IO_WRITE_SUBARRAY4_NEW(fname, varname, field, sizes, mpio_aux(iflag_mode), tmp_mpi)
        
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
        CALL DNS_MPI_WRITE_PE0_SINGLE(LOC_UNIT_ID, nx,ny,nz, subdomain, field(1,1,1,ifield), tmp_mpi(1,1), tmp_mpi(1,2))
        IF ( ims_pro .EQ. 0 ) THEN
#else
           DO k = subdomain(5),subdomain(6)
              DO j = subdomain(3),subdomain(4)
                 WRITE(LOC_UNIT_ID) (SNGL(field(subdomain(1):subdomain(2),j,k,ifield)))
              ENDDO
           ENDDO
#endif
        CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
        ENDIF
#endif
     ENDDO

  ENDIF
  
  RETURN
END SUBROUTINE VISUALS_WRITE

! ###################################################################
! ###################################################################
#ifdef USE_MPI

SUBROUTINE VISUALS_MPIO_AUX(opt_format, subdomain)

  USE DNS_TYPES,  ONLY : subarray_structure
  USE DNS_GLOBAL, ONLY : imax_total,jmax_total,kmax_total, imax,jmax,kmax
  USE DNS_MPI

  IMPLICIT NONE

#include "mpif.h" 

  TINTEGER,                 INTENT(IN)  :: opt_format, subdomain(6)

! -----------------------------------------------------------------------
  TINTEGER                :: ndims
  TINTEGER, DIMENSION(3)  :: sizes, locsize, offset

! #######################################################################
  mpio_aux(:)%active = .FALSE.
  mpio_aux(:)%offset = 0
  IF ( opt_format .EQ. 1 ) mpio_aux(:)%offset = 244 ! # bytes of ensight header

! ###################################################################
! Saving full vertical xOy planes; using subdomain(5) to define the plane
  IF ( ims_pro_k .EQ. ( subdomain(5) /kmax) ) mpio_aux(1)%active = .TRUE.
  mpio_aux(1)%communicator = ims_comm_x

  ndims = 2
  sizes(1)   = imax_total;   sizes(2)   = subdomain(4)-subdomain(3)+1
  locsize(1) = imax;         locsize(2) = subdomain(4)-subdomain(3)+1
  offset(1)  = ims_offset_i; offset(2)  = 0
  
  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
       MPI_ORDER_FORTRAN, MPI_REAL4, mpio_aux(1)%subarray, ims_err)
  CALL MPI_Type_commit(mpio_aux(1)%subarray, ims_err)

! Saving full vertical zOy planes; using subdomain(1) to define the plane
  IF ( ims_pro_i .EQ.  ( subdomain(1) /imax) ) mpio_aux(2)%active = .TRUE.
  mpio_aux(2)%communicator = ims_comm_z

  ndims = 2
                             sizes(1)   = subdomain(4)-subdomain(3)+1; sizes(2)   = kmax_total 
                             locsize(1) = subdomain(4)-subdomain(3)+1; locsize(2) = kmax 
                             offset(1)  = 0;                           offset(2)  = ims_offset_k

  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
       MPI_ORDER_FORTRAN, MPI_REAL4, mpio_aux(2)%subarray, ims_err)
  CALL MPI_Type_commit(mpio_aux(2)%subarray, ims_err)

! Saving full blocks xOz planes
  mpio_aux(3)%active = .TRUE.
  mpio_aux(3)%communicator = MPI_COMM_WORLD

  ndims = 3
  sizes(1)   = imax_total;   sizes(2)   = subdomain(4)-subdomain(3)+1; sizes(3)   = kmax_total 
  locsize(1) = imax;         locsize(2) = subdomain(4)-subdomain(3)+1; locsize(3) = kmax 
  offset(1)  = ims_offset_i; offset(2)  = 0;                           offset(3)  = ims_offset_k
  
  CALL MPI_Type_create_subarray(ndims, sizes, locsize, offset, & 
       MPI_ORDER_FORTRAN, MPI_REAL4, mpio_aux(3)%subarray, ims_err)
  CALL MPI_Type_commit(mpio_aux(3)%subarray, ims_err)
  
  RETURN
END SUBROUTINE VISUALS_MPIO_AUX

#endif

!########################################################################
! Writing data in Ensight Gold Variable File Format
!########################################################################
SUBROUTINE ENSIGHT_FIELD(name, iheader, nx,ny,nz, nfield, subdomain, field, tmp_mpi)

#ifdef USE_MPI
  USE DNS_MPI, ONLY : ims_pro
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
  TINTEGER i,j,k, ifield

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

