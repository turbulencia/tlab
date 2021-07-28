#include "types.h"
#include "dns_error.h"

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/04/26 - J.P. Mellado
!#              PE 0 handling extracted into DNS_MPI_READ_PE0
!#
!########################################################################
!# DESCRIPTION
!#
!# Read scalar restart file. There are USE_MPI and SERIAL modes, and
!# within the USE_MPI mode there are USE_MPI_IO, USE_SPLIT_IO
!# and PE 0 handling.
!#
!# IEEE double precision floating point representation.
!# Unformatted records
!# Byte ordering is big-endian
!# Embedded record information (4-bytes integer)
!#
!########################################################################
!# ARGUMENTS
!#
!# nfield     In      Number of fields in the file
!# iread      In      Control flag, if set to 0 read all fields (array a
!#                    w/ space for nfield), if not read only iread (array
!#                    a only 1 field)
!# iheader    In      Header control flag:
!#                    0 No header
!#                    1 Scalar fields
!#                    2 Flow fields
!#
!########################################################################
#define LOC_UNIT_ID i54
#define LOC_STATUS 'old'

SUBROUTINE IO_READ_FIELDS_ARRAY(name, nx,ny,nz, itxc, iheader, nfield, iread, a, txc)

  USE TLAB_VARS,ONLY : itime, rtime, visc
  USE TLAB_PROCS
#ifdef USE_MPI
  USE TLAB_MPI_VARS,   ONLY : ims_offset_k, ims_npro_i, ims_npro_k, ims_pro, ims_err
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) name
  TINTEGER itxc,iheader,nfield,iread, nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*) :: a
  TREAL, DIMENSION(itxc)       :: txc

! -------------------------------------------------------------------
  TREAL dummy
  TINTEGER nxy, is, iz
  TINTEGER nx_total, ny_total, nz_total

#ifdef USE_MPI

#ifdef USE_MPI_IO
  INTEGER mpio_fh
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp, mpio_locoff, reclen
  INTEGER mpio_locfldsize
  INTEGER status(MPI_STATUS_SIZE)

#else
#ifdef USE_SPLIT_IO
  CHARACTER*32 name_loc
#else
  TREAL, DIMENSION(:), ALLOCATABLE :: tmp
  TINTEGER itmp, iuse_tmp
#endif
#endif

#else
  TINTEGER reclen

#endif

! ###################################################################
#ifdef USE_MPI
  nx_total = nx*ims_npro_i
  ny_total = ny
  nz_total = nz*ims_npro_k
#else
  nx_total = nx
  ny_total = ny
  nz_total = nz
#endif

  nxy = nx*ny

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF      ( iheader .EQ. 1 ) THEN ! Scalar type
     IF ( ims_pro .EQ. 0 ) THEN
        OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted')
        REWIND(LOC_UNIT_ID)
        CALL HEADER_READ_SCAL(LOC_UNIT_ID, i0, itime, rtime,     nx_total,ny_total,nz_total,&
             dummy, nfield)
        CLOSE(LOC_UNIT_ID)
     ENDIF

! Displacement to start of field
     mpio_disp = &
          2*SIZEOFINT + SIZEOFREAL + SIZEOFINT +            &
          2*SIZEOFINT + SIZEOFINT  + SIZEOFINT + SIZEOFINT +&
          2*SIZEOFINT + SIZEOFREAL +                        &
            SIZEOFINT ! add the record length information of first record
     IF ( nfield .GT. 1 ) THEN
        mpio_disp = mpio_disp + 2*SIZEOFINT + SIZEOFINT
     ENDIF

  ELSE IF ( iheader .EQ. 2 ) THEN ! Flow type
     IF ( ims_pro .EQ. 0 ) THEN
        OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted')
        REWIND(LOC_UNIT_ID)
!        CALL HEADER_READ_FLOW(LOC_UNIT_ID, i0, itime, rtime, dt, nx,ny,nz_total,&
        CALL HEADER_READ_FLOW(LOC_UNIT_ID, i0, itime, rtime, dummy, nx_total,ny_total,nz_total,&
             visc, dummy,dummy,dummy,dummy)
        CLOSE(LOC_UNIT_ID)
     ENDIF

! Displacement to start of field
     mpio_disp = &
          2*SIZEOFINT + SIZEOFINT  + SIZEOFREAL + SIZEOFREAL +&
          2*SIZEOFINT + SIZEOFINT  + SIZEOFINT  + SIZEOFINT  +&
          2*SIZEOFINT + SIZEOFREAL + SIZEOFREAL + SIZEOFREAL +&
          2*SIZEOFINT + SIZEOFREAL + SIZEOFREAL +             &
            SIZEOFINT ! add the record length information of first record

  ELSE
     mpio_disp = &
            SIZEOFINT ! add the record length information of first record

  ENDIF
  CALL MPI_BCAST(itime, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(rtime, 1, MPI_REAL8,    0, MPI_COMM_WORLD, ims_err)
  IF ( iheader .EQ. 2 ) THEN ! Flow type
!     CALL MPI_BCAST(dt,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
     CALL MPI_BCAST(visc,  1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)

  mpio_locfldsize = nxy*nz
  mpio_locoff     = nxy
  mpio_locoff     = ims_offset_k*mpio_locoff ! reclen might be of type larger than INT4
  reclen          = nxy
  reclen          = reclen*nz_total
  reclen          = reclen*SIZEOFREAL ! reclen might be of type larger than INT4

  DO is = 1,nfield
     IF ( iread .EQ. 0 ) THEN; iz = is
     ELSE;                     iz = 1; ENDIF

     CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, &
          MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)

     CALL MPI_FILE_READ_AT_ALL(mpio_fh, mpio_locoff, a(1,iz), &
          mpio_locfldsize, MPI_REAL8, status, ims_err)

     mpio_disp = mpio_disp + 2 * SIZEOFINT + reclen

     IF ( is .EQ. iread ) GOTO 11

  ENDDO

11 CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

#else
#ifdef USE_SPLIT_IO
! ###################################################################
! Use SPLITS for restart files
! ###################################################################
  WRITE(name_loc,*) ims_pro
  name_loc=TRIM(ADJUSTL(name))//'.'//TRIM(ADJUSTL(name_loc))

  OPEN(LOC_UNIT_ID,file=name_loc, status=LOC_STATUS,form='unformatted')
  REWIND(LOC_UNIT_ID)

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF      ( iheader .EQ. 1 ) THEN ! Scalar type
     CALL HEADER_READ_SCAL(LOC_UNIT_ID, i0, itime, rtime,     nx_total,ny_total,nz_total,&
          dummy, nfield)
  ELSE IF ( iheader .EQ. 2 ) THEN ! Flow type
!     CALL HEADER_READ_FLOW(LOC_UNIT_ID, i0, itime, rtime, dt, nx,ny,nz_total,&
     CALL HEADER_READ_FLOW(LOC_UNIT_ID, i0, itime, rtime, dummy, nx_total,ny_total,nz_total,&
          visc, dummy,dummy,dummy,dummy)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  DO is = 1,nfield
     IF ( iread .EQ. 0 ) THEN; iz = is
     ELSE;                     iz = 1; ENDIF

     READ(LOC_UNIT_ID) a(:,iz)

     IF ( is .EQ. iread ) GOTO 22

  ENDDO

22 CLOSE(LOC_UNIT_ID)

#else
! ###################################################################
! Let Process 0 handle all IO
! ###################################################################
  IF ( ims_pro .EQ. 0 ) THEN
! Allocate temporary array tmp for USE_MPI PE#0 mode if required
     itmp = kmax*nxy
     IF ( itmp .GT. itxc) THEN; iuse_tmp = 1; ALLOCATE(tmp(itmp))
     ELSE;                      iuse_tmp = 0; ENDIF

! Open file
#include "dns_open_file.h"
     REWIND(LOC_UNIT_ID)

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
     IF      ( iheader .EQ. 1 ) THEN ! Scalar type
        CALL HEADER_READ_SCAL(LOC_UNIT_ID, i1, itime, rtime,     nx_total,ny_total,nz_total,&
             dummy, nfield)
     ELSE IF ( iheader .EQ. 2 ) THEN ! Flow type
!        CALL HEADER_READ_FLOW(LOC_UNIT_ID, i1, itime, rtime, dt, nx,ny,nz_total,&
        CALL HEADER_READ_FLOW(LOC_UNIT_ID, i1, itime, rtime, dummy, nx_total,ny_total,nz_total,&
             visc, dummy,dummy,dummy,dummy)
     ENDIF

  ENDIF

  CALL MPI_BCAST(itime, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  CALL MPI_BCAST(rtime, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
  IF ( iheader .EQ. 2 ) THEN ! Flow type
!     CALL MPI_BCAST(dt,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
     CALL MPI_BCAST(visc,  1, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  DO is = 1,nfield

     IF ( iread .EQ. 0 ) THEN; iz = is
     ELSE;                     iz = 1; ENDIF

     IF ( iuse_tmp .EQ. 1 ) THEN
        CALL DNS_MPI_READ_PE0(LOC_UNIT_ID, i1, nx,ny,nz,nz_total, itmp, a(1,iz), tmp)
     ELSE
        CALL DNS_MPI_READ_PE0(LOC_UNIT_ID, i1, nx,ny,nz,nz_total, itxc, a(1,iz), txc)
     ENDIF

     IF ( is .EQ. iread ) GOTO 33

  ENDDO

33 IF ( ims_pro .EQ. 0 ) THEN
     CLOSE(LOC_UNIT_ID)

     IF ( iuse_tmp .EQ. 1 ) DEALLOCATE(tmp)
  ENDIF

#endif
#endif

#else
! ###################################################################
! Serial case
! ###################################################################
#include "dns_open_file.h"
  REWIND(LOC_UNIT_ID)

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF      ( iheader .EQ. 1 ) THEN ! Scalar type
     CALL HEADER_READ_SCAL(LOC_UNIT_ID, i1, itime, rtime,     nx_total,ny_total,nz_total,&
          dummy, nfield)
  ELSE IF ( iheader .EQ. 2 ) THEN ! Flow type
!     CALL HEADER_READ_FLOW(LOC_UNIT_ID, i1, itime, rtime, dt, nx,ny,nz_total,&
     CALL HEADER_READ_FLOW(LOC_UNIT_ID, i1, itime, rtime, dummy, nx_total,ny_total,nz_total,&
          visc, dummy,dummy,dummy,dummy)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  DO is = 1,nfield
     IF ( iread .EQ. 0 ) THEN; iz = is
     ELSE;                     iz = 1; ENDIF

     READ(LOC_UNIT_ID) reclen
     READ(LOC_UNIT_ID) a(:,iz)
     READ(LOC_UNIT_ID) reclen
     IF ( is .EQ. iread ) GOTO 44
  ENDDO

44 CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_READ_FIELDS_ARRAY

#undef LOC_UNIT_ID
#undef LOC_STATUS

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 1999/01/01 - C. Pantano
!#              Created
!# 2007/04/26 - J.P. Mellado
!#              PE 0 handling extracted into DNS_MPI_READ_PE0
!#
!########################################################################
!# DESCRIPTION
!#
!# Read scalar restart file. There are USE_MPI and SERIAL modes, and
!# within the USE_MPI mode there are USE_MPI_IO, USE_SPLIT_IO
!# and PE 0 handling.
!#
!########################################################################
!# ARGUMENTS
!#
!# nfield     In      Number of fields in the file
!# iheader    In      Header control flag:
!#                    0 No header
!#                    1 Scalar fields
!#                    2 Flow fields
!#
!########################################################################
#define LOC_UNIT_ID i55
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_FIELDS_ARRAY(name, nx,ny,nz, itxc, iheader, nfield, a, txc)

  USE TLAB_VARS, ONLY : itime, rtime, visc, prandtl, schmidt ! header info
  USE TLAB_PROCS
  USE THERMO_GLOBAL, ONLY : gama0
#ifdef USE_MPI
  USE TLAB_MPI_VARS,   ONLY : ims_offset_k, ims_npro_i, ims_npro_k, ims_pro, ims_err
#endif

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) name
  TINTEGER itxc, iheader, nfield, nx,ny,nz
  TREAL, DIMENSION(nx*ny*nz,*) :: a
  TREAL, DIMENSION(itxc)       :: txc

! -------------------------------------------------------------------
  TREAL dummy
  TINTEGER nxy, is
  TINTEGER nx_total, ny_total, nz_total

#ifdef USE_MPI

#ifdef USE_MPI_IO
  INTEGER mpio_fh
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp, mpio_locoff, reclen
  INTEGER mpio_locfldsize
  INTEGER status(MPI_STATUS_SIZE)

#else
#ifdef USE_SPLIT_IO
  CHARACTER*32 name_loc
#else
  TREAL, DIMENSION(:), ALLOCATABLE :: tmp
  TINTEGER itmp, iuse_tmp
#endif
#endif

#endif

! ###################################################################
#ifdef USE_MPI
  nx_total = nx*ims_npro_i
  ny_total = ny
  nz_total = nz*ims_npro_k
#else
  nx_total = nx
  ny_total = ny
  nz_total = nz
#endif

  nxy      = nx*ny
  dummy    = C_0_R

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
! Scalar type
  IF      ( iheader .EQ. 1 ) THEN
     IF ( ims_pro .EQ. 0 ) THEN
        OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted')
        CALL HEADER_WRITE_SCAL(LOC_UNIT_ID, i0, itime, rtime,    nx_total,ny_total,nz_total,&
             schmidt(1), nfield)
        CLOSE(LOC_UNIT_ID)
     ENDIF

! Displacement to start of field
     mpio_disp = &
          2*SIZEOFINT + SIZEOFREAL + SIZEOFINT +            &
          2*SIZEOFINT + SIZEOFINT  + SIZEOFINT + SIZEOFINT +&
          2*SIZEOFINT + SIZEOFREAL
     IF ( nfield .GT. 1 ) THEN
        mpio_disp = mpio_disp + 2*SIZEOFINT + SIZEOFINT
     ENDIF

! Flow type
  ELSE IF ( iheader .EQ. 2 ) THEN
     IF ( ims_pro .EQ. 0 ) THEN
        OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted')
!        CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i0, itime, rtime, dt, nx,ny,nz_total,&
        CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i0, itime, rtime, dummy, nx_total,ny_total,nz_total,&
             visc, gama0, prandtl, dummy, dummy)
        CLOSE(LOC_UNIT_ID)
     ENDIF

! Displacement to start of field
     mpio_disp = &
          2*SIZEOFINT + SIZEOFINT  + SIZEOFREAL + SIZEOFREAL +&
          2*SIZEOFINT + SIZEOFINT  + SIZEOFINT  + SIZEOFINT  +&
          2*SIZEOFINT + SIZEOFREAL + SIZEOFREAL + SIZEOFREAL +&
          2*SIZEOFINT + SIZEOFREAL + SIZEOFREAL

  ELSE
     IF ( ims_pro .EQ. 0 ) THEN
        OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted') ! Create file
        CLOSE(LOC_UNIT_ID)
     ENDIF
     mpio_disp = 0

  ENDIF

! This is an MPI barrier ! do not remove
  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)

  mpio_locfldsize = nxy*nz
  mpio_locoff     = nxy
  mpio_locoff     = ims_offset_k*mpio_locoff ! reclen might be of type larger than INT4
  reclen          = nxy
  reclen          = reclen*nz_total
  reclen          = reclen*SIZEOFREAL ! reclen might be of type larger than INT4

  DO is = 1,nfield
!    Write Reclen
     CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_INTEGER4, &
          MPI_INTEGER4, 'native', MPI_INFO_NULL, ims_err)
     IF ( ims_pro .EQ. 0 ) THEN
        CALL MPI_FILE_WRITE(mpio_fh, reclen, 1, MPI_INTEGER4, status, ims_err)
     ENDIF

!    This is an MPI barrier ! do not remove
     CALL MPI_FILE_SYNC(mpio_fh, ims_err)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)
     CALL MPI_FILE_SYNC(mpio_fh, ims_err)

!    Write Field
     mpio_disp = mpio_disp + SIZEOFINT
     CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, &
          MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
     CALL MPI_FILE_WRITE_AT_ALL(mpio_fh, mpio_locoff, a(1,is),&
          mpio_locfldsize, MPI_REAL8, status, ims_err)

!    Write Reclen
     mpio_disp = mpio_disp + reclen
     CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_INTEGER4, &
          MPI_INTEGER4, 'native', MPI_INFO_NULL, ims_err)
     IF ( ims_pro .EQ. 0 ) THEN
        CALL MPI_FILE_WRITE(mpio_fh, reclen, 1, MPI_INTEGER4, status, ims_err)
     ENDIF

!    This is an MPI barrier ! do not remove
     CALL MPI_FILE_SYNC(mpio_fh, ims_err)
     CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)
     CALL MPI_FILE_SYNC(mpio_fh, ims_err)

     mpio_disp = mpio_disp + SIZEOFINT

  ENDDO

  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)

#else
#ifdef USE_SPLIT_IO
! ###################################################################
! Use SPLITS for restart files
! ###################################################################
  WRITE(name_loc,*) ims_pro
  name_loc=TRIM(ADJUSTL(name))//'.'//TRIM(ADJUSTL(name_loc))

  OPEN(unit=LOC_UNIT_ID,file=name_loc,form='unformatted',status=LOC_STATUS)

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF      ( iheader .EQ. 1 ) THEN ! Scalar type
     CALL HEADER_WRITE_SCAL(LOC_UNIT_ID, i0, itime, rtime,     nx_total,ny_total,nz_total,&
          schmidt(1), nfield)
  ELSE IF ( iheader .EQ. 2 ) THEN ! Flow type
!     CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i0, itime, rtime, dt, nx,ny,nz_total,&
     CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i0, itime, rtime, dummy, nx_total,ny_total,nz_total,&
          visc, gama0, prandtl, dummy, dummy)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  DO is = 1,nfield
     WRITE(LOC_UNIT_ID) a(:,is)
  ENDDO

  CLOSE(LOC_UNIT_ID)

#else
! ###################################################################
! Let Process 0 handle all IO
! ###################################################################
  IF ( ims_pro .EQ. 0 ) THEN
! Allocate temporary array tmp for USE_MPI PE#0 mode if required
     itmp = kmax*nxy
     IF ( itmp .GT. itxc) THEN; iuse_tmp = 1; ALLOCATE(tmp(itmp))
     ELSE;                      iuse_tmp = 0; ENDIF

! Open file
#include "dns_open_file.h"

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
     IF      ( iheader .EQ. 1 ) THEN ! Scalar type
        CALL HEADER_WRITE_SCAL(LOC_UNIT_ID, i1, itime, rtime,     nx_total,ny_total,nz_total,&
             schmidt(1), nfield)
     ELSE IF ( iheader .EQ. 2 ) THEN ! Flow type
!        CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i1, itime, rtime, dt, nx,ny,nz_total,&
        CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i1, itime, rtime, dummy, nx_total,ny_total,nz_total,&
             visc, gama0, prandtl, dummy, dummy)
     ENDIF

  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  DO is = 1, nfield
     IF ( iuse_tmp .EQ. 1 ) THEN
        CALL DNS_MPI_WRITE_PE0(LOC_UNIT_ID, i1, nx,ny,nz,nz_total, itmp, a(1,is), tmp)
     ELSE
        CALL DNS_MPI_WRITE_PE0(LOC_UNIT_ID, i1, nx,ny,nz,nz_total, itxc, a(1,is), txc)
     ENDIF
  ENDDO

  IF ( ims_pro .EQ. 0 ) THEN
     CLOSE(LOC_UNIT_ID)
     IF ( iuse_tmp .EQ. 1 ) THEN
        DEALLOCATE(tmp)
     ENDIF
  ENDIF

#endif
#endif

#else
! ###################################################################
! Serial case
! ###################################################################
#include "dns_open_file.h"

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF      ( iheader .EQ. 1 ) THEN ! Scalar type
     CALL HEADER_WRITE_SCAL(LOC_UNIT_ID, i1, itime, rtime,     nx_total,ny_total,nz_total,&
          schmidt(1), nfield)
  ELSE IF ( iheader .EQ. 2 ) THEN ! Flow type
!     CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i1, itime, rtime, dt, nx,ny,nz_total,&
     CALL HEADER_WRITE_FLOW(LOC_UNIT_ID, i1, itime, rtime, dummy, nx_total,ny_total,nz_total,&
          visc, gama0, prandtl, dummy, dummy)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  DO is = 1,nfield
     WRITE(LOC_UNIT_ID) nx*ny*nz*SIZEOFREAL
     WRITE(LOC_UNIT_ID) a(:,is)
     WRITE(LOC_UNIT_ID) nx*ny*nz*SIZEOFREAL
  ENDDO

  CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_WRITE_FIELDS_ARRAY

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2009/04/06 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Read header of flow file and do some consistency checks
!#
!########################################################################
!# ARGUMENTS
!#
!# itime    Out     Iteration at which it was created
!# rtime    Out     Time at which it was created
!#
!########################################################################
SUBROUTINE HEADER_READ_FLOW(unit, irec, itime, rtime, dt, imax,jmax,kmax, &
     visc, gama, prandtl, dummy1, dummy2)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_PROCS

  IMPLICIT NONE

  TINTEGER unit, itime, imax, jmax, kmax, irec
  TREAL rtime, dt, visc, gama, prandtl, dummy1, dummy2

! -------------------------------------------------------------------
  TINTEGER imaxdum, jmaxdum, kmaxdum
  TREAL time(2), param(6), dummy3
  TINTEGER reclen

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) itime
     READ(unit) rtime
     READ(unit) dt
     READ(unit) reclen
  ELSE
     READ(unit) itime, time
     rtime = time(1)
     dt = time(2)
  ENDIF

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) imaxdum
     READ(unit) jmaxdum
     READ(unit) kmaxdum
     READ(unit) reclen
  ELSE
     READ(unit) imaxdum, jmaxdum, kmaxdum
  ENDIF

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) visc
     READ(unit) gama
     READ(unit) prandtl
     READ(unit) dummy1
     READ(unit) dummy2
     READ(unit) dummy3
     READ(unit) reclen
  ELSE
     READ(unit) param
     visc    = param(1)
     gama    = param(2)
     prandtl = param(3)
     dummy3  = param(4)
     dummy1  = param(5)
     dummy2  = param(6)
  ENDIF

! -------------------------------------------------------------------
! Consistency check
! -------------------------------------------------------------------
  IF ( imaxdum .NE. imax .OR. jmaxdum .NE. jmax .OR. kmaxdum .NE. kmax) THEN
     CLOSE(unit)
     CALL TLAB_WRITE_ASCII(efile, 'HEADER_READ_FLOW: Grid size mismatch')
     CALL TLAB_STOP(DNS_ERROR_DIMGRID)
  ENDIF

  RETURN
END SUBROUTINE HEADER_READ_FLOW

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2009/04/06 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Write header
!#
!# dummy 3 just introduced to keep size of header
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE HEADER_WRITE_FLOW(unit, irec, itime, rtime, dt, imax,jmax,kmax,&
     visc, gama, prandtl, dummy1, dummy2)

  IMPLICIT NONE

  TINTEGER unit, itime, imax, jmax, kmax, irec
  TREAL rtime, dt, visc, gama, prandtl, dummy1, dummy2

! -------------------------------------------------------------------
  TREAL time(2), param(6), dummy3
  TINTEGER reclen

!########################################################################
  time(1) = rtime
  time(2) = dt
  reclen = (SIZEOFINT+2*SIZEOFREAL)
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) itime
     WRITE(unit) rtime
     WRITE(unit) dt
     WRITE(unit) reclen
  ELSE
     WRITE(unit) itime, time
  ENDIF

  reclen = 3*SIZEOFINT
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) imax
     WRITE(unit) jmax
     WRITE(unit) kmax
     WRITE(unit) reclen
  ELSE
     WRITE(unit) imax, jmax, kmax
  ENDIF

  dummy3   = C_0_R
  param(1) = visc
  param(2) = gama
  param(3) = prandtl
  param(4) = dummy3
  param(5) = dummy1
  param(6) = dummy2
  reclen = 6*SIZEOFREAL
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) visc
     WRITE(unit) gama
     WRITE(unit) prandtl
     WRITE(unit) dummy3
     WRITE(unit) dummy1
     WRITE(unit) dummy2
     WRITE(unit) reclen
  ELSE
     WRITE(unit) param
  ENDIF

  RETURN
END SUBROUTINE HEADER_WRITE_FLOW

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2009/04/06 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Read header of scalars file and do some consistency checks
!#
!########################################################################
!# ARGUMENTS
!#
!# itime    Out     Iteration at which it was created
!# rtime    Out     Time at which it was created
!#
!########################################################################
SUBROUTINE HEADER_READ_SCAL(unit, irec, jiter, rtime, imax, jmax, kmax, schmidt, inbsc)

  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_PROCS

  IMPLICIT NONE

  TINTEGER unit, jiter, imax, jmax, kmax, irec
  TINTEGER inbsc
  TREAL rtime, schmidt

! -------------------------------------------------------------------
  TINTEGER imaxdum, jmaxdum, kmaxdum, inbscdum, reclen
  TREAL tmp(1)

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) jiter
     READ(unit) rtime
     READ(unit) reclen
  ELSE
     READ(unit) jiter, tmp
     rtime = tmp(1)
  ENDIF

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) imaxdum
     READ(unit) jmaxdum
     READ(unit) kmaxdum
     READ(unit) reclen
  ELSE
     READ(unit) imaxdum, jmaxdum, kmaxdum
  ENDIF

  IF ( irec .EQ. 1 ) THEN
     READ(unit) reclen
     READ(unit) schmidt
     READ(unit) reclen
  ELSE
     READ(unit) tmp
     schmidt = tmp(1)
  ENDIF

  IF ( inbsc .GT. 1 ) THEN
     IF ( irec .EQ. 1 ) THEN
        READ(unit) reclen
        READ(unit) inbscdum
        READ(unit) reclen
     ELSE
        READ(unit) inbscdum
     ENDIF
  ELSE
     inbscdum = 1
  ENDIF

! -------------------------------------------------------------------
! Consistency check
! -------------------------------------------------------------------
  IF (imaxdum .NE. imax .OR. jmaxdum .NE. jmax .OR. kmaxdum .NE. kmax) THEN
     CLOSE(unit)
     CALL TLAB_WRITE_ASCII(efile, 'HEADER_READ_SCAL. Grid size mismatch')
     CALL TLAB_STOP(DNS_ERROR_DIMGRID)
  ENDIF

  IF ( inbscdum .NE. inbsc ) THEN
     CALL TLAB_WRITE_ASCII(efile, 'HEADER_READ_SCAL. Wrong number of scalars in restart file')
     CALL TLAB_STOP(DNS_ERROR_WRONGSCALNB)
  ENDIF

  RETURN
END SUBROUTINE HEADER_READ_SCAL

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2009/04/06 - J.P. Mellado
!#              Cleaned
!#
!########################################################################
!# DESCRIPTION
!#
!# Write header
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE HEADER_WRITE_SCAL(unit, irec, jiter, rtime, imax,jmax,kmax, schmidt, inbsc)

  IMPLICIT NONE

  TINTEGER unit, jiter, imax, jmax, kmax, irec
  TINTEGER inbsc
  TREAL rtime, schmidt

! -------------------------------------------------------------------
  TREAL tmp(1)
  TINTEGER reclen

!########################################################################
  tmp(1) = rtime
  reclen = SIZEOFINT+SIZEOFREAL
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) jiter
     WRITE(unit) rtime
     WRITE(unit) reclen
  ELSE
     WRITE(unit) jiter, tmp
  ENDIF

  reclen = 3*SIZEOFINT
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) imax
     WRITE(unit) jmax
     WRITE(unit) kmax
     WRITE(unit) reclen
  ELSE
     WRITE(unit) imax, jmax, kmax
  ENDIF

  tmp(1) = schmidt
  reclen = SIZEOFREAL
  IF ( irec .EQ. 1 ) THEN
     WRITE(unit) reclen
     WRITE(unit) schmidt
     WRITE(unit) reclen
  ELSE
     WRITE(unit) tmp
  ENDIF

  IF ( inbsc .GT. 1 ) THEN
     reclen = SIZEOFINT
     IF ( irec .EQ. 1 ) THEN
        WRITE(unit) reclen
        WRITE(unit) inbsc
        WRITE(unit) reclen
     ELSE
        WRITE(unit) inbsc
     ENDIF
  ENDIF

  RETURN
END SUBROUTINE HEADER_WRITE_SCAL
