#include "types.h"

SUBROUTINE VISUALS_WRITE(name, itype, iformat, nx,ny,nz, subdomain, field, tmp_mpi)

  IMPLICIT NONE

#include "integers.h"

  TINTEGER iformat, itype, nx,ny,nz, subdomain(6)
  TREAL, DIMENSION(nx,ny,nz,*) :: field
  TREAL, DIMENSION(nx*ny*nz,2) :: tmp_mpi
  CHARACTER*(*) name

! -------------------------------------------------------------------
  TINTEGER iaux_loc, nfield

! ###################################################################
  IF      ( itype .EQ. 0 ) THEN; nfield = 1;
  ELSE IF ( itype .EQ. 1 ) THEN; nfield = 3;
  ELSE IF ( itype .EQ. 2 ) THEN; nfield = 6; ENDIF

  IF      ( iformat .EQ. 0 ) THEN ! standard scalar format
     iaux_loc = nx*ny*nz
     CALL DNS_WRITE_FIELDS(name, i1, nx,ny,nz, nfield, iaux_loc, field, tmp_mpi)

  ELSE IF ( iformat .EQ. 1 ) THEN ! ensight
     CALL ENSIGHT_FIELD(name, i1, nx,ny,nz, nfield, subdomain, field, tmp_mpi)
           
  ELSE IF ( iformat .EQ. 2 ) THEN ! ensight, no header
     CALL ENSIGHT_FIELD(name, i0, nx,ny,nz, nfield, subdomain, field, tmp_mpi)

  ENDIF

  RETURN
END SUBROUTINE VISUALS_WRITE

!########################################################################
! Writing data in Ensight Gold Variable File Format
!########################################################################
#define LOC_UNIT_ID i55
#define LOC_STATUS 'unknown'

SUBROUTINE ENSIGHT_FIELD(name, iheader, nx,ny,nz, nfield, subdomain, field, tmp_mpi)

#ifdef USE_MPI
  USE DNS_MPI
#endif

  IMPLICIT NONE

#include "integers.h"

  CHARACTER*(*) name
  TINTEGER, INTENT(IN) :: iheader ! 0 no header; 1 header

  TINTEGER, INTENT(IN) :: nx,ny,nz, nfield, subdomain(6)
  TREAL, DIMENSION(nx,ny,nz,nfield) :: field
#ifdef USE_MPI
  TREAL, DIMENSION(nx*ny*nz,2     ) :: tmp_mpi
#else
  TREAL, DIMENSION(*)               :: tmp_mpi
#endif

! -------------------------------------------------------------------
  CHARACTER*80 line
  TINTEGER i,j,k, ifield

! ###################################################################        
! Header
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
