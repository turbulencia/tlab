!# Compile with
!# mpxlf90_r -d -qsuffix=cpp=f90 -WF,-DUSE_MPI,-DUSE_MPI_IO -q64 -qextname -qsave -qarch=pwr6 -qtune=pwr6 -O3 -qenablevmx -qhot=simd -qessl -o vmpi_io.x vmpi_io.f90
#define SIZEOFREAL 8
#define SIZEOFINT  4

#define TREAL      REAL(8)
#define TINTEGER   INTEGER(4)

!########################################################################
MODULE DNS_MPI
  IMPLICIT NONE
  SAVE

  INTEGER  :: ims_pro,  ims_pro_i,  ims_pro_j,  ims_pro_k
  INTEGER  :: ims_npro, ims_npro_i, ims_npro_j, ims_npro_k
  INTEGER  :: ims_offset_i, ims_offset_j, ims_offset_k
  INTEGER  :: ims_err

END MODULE DNS_MPI

MODULE DNS_GLOBAL
  IMPLICIT NONE
  SAVE

  TINTEGER :: imax_total, jmax_total, kmax_total

END MODULE DNS_GLOBAL

!########################################################################
PROGRAM VMPI_IO

  USE TLabMPI_VARS
  USE TLAB_VARS

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TINTEGER imax, jmax, kmax
  TINTEGER inb_flow, inb_flow_array, inb_scal, inb_scal_array
  TINTEGER ip, it, it_max

  TREAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: q, s
  TREAL, DIMENSION(:),   ALLOCATABLE, SAVE :: wrk3d

  CHARACTER*32 :: name_loc
  TINTEGER ifield

  TINTEGER isize_max, isize
  PARAMETER(isize_max=20)
  TREAL params(isize_max)

  TINTEGER :: i1=1, i2=2

! ###################################################################
  imax_total = 1024
  jmax_total =  512
  kmax_total = 1024

  imax       = 1024; ims_npro_i = imax_total/imax
  jmax       =  512
  kmax       =    8; ims_npro_k = kmax_total/kmax

  inb_flow       = 3
  inb_flow_array = inb_flow
  inb_scal       = 1
  inb_scal_array = inb_scal

  ALLOCATE(q(imax*jmax*kmax, inb_flow_array))
  ALLOCATE(s(imax*jmax*kmax, inb_scal_array))
  ALLOCATE(wrk3d(imax*jmax*kmax))

! ###################################################################
! from DNS_START
  call MPI_INIT(ims_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ims_npro,ims_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro, ims_err)

! ###################################################################
! from TLabMPI_Initialize()
  ims_pro_i = MOD(ims_pro,ims_npro_i) ! Starting at 0
  ims_pro_k =     ims_pro/ims_npro_i  ! Starting at 0

  ims_offset_i = ims_pro_i *imax
  ims_offset_j = 0
  ims_offset_k = ims_pro_k *kmax

!  ims_map_i(1) = ims_pro_k*ims_npro_i
!  DO ip = 2,ims_npro_i
!     ims_map_i(ip) = ims_map_i(ip-1) + 1
!  ENDDO

!  ims_map_k(1) = ims_pro_i
!  DO ip = 2,ims_npro_k
!     ims_map_k(ip) = ims_map_k(ip-1) + ims_npro_i
!  ENDDO

! ###################################################################
  it_max = 20

  DO it = 0,it_max
     IF ( ims_pro .EQ. 0 ) THEN
        OPEN(UNIT=22, FILE='dns.out', STATUS='unknown',POSITION='APPEND')
        WRITE(22,*) 'Iteration ...:', it
        CLOSE(22)
     ENDIF

! velocity
     DO ifield = 1,3
        WRITE(name_loc,'(I2)') ifield
        name_loc='flow.'//TRIM(ADJUSTL(name_loc))
        isize = isize_max
        CALL IO_READ_FIELDS_SPLIT(name_loc, i2, imax,jmax,kmax,it, isize,params, q(1,ifield),wrk3d)
     ENDDO

! scalar
     DO ifield = 1,1
        WRITE(name_loc,'(I2)') ifield
        name_loc='scal.'//TRIM(ADJUSTL(name_loc))
        isize = isize_max
        CALL IO_READ_FIELDS_SPLIT(name_loc, i1, imax,jmax,kmax,it, isize,params, s(1,ifield),wrk3d)
     ENDDO

  ENDDO

  CALL MPI_FINALIZE(ims_err)

  STOP
END PROGRAM VMPI_IO

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2010/04/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Read file. There are USE_MPI (only MPI_IO) and SERIAL modes
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'
#define SIZEOFBYTE 1

SUBROUTINE IO_READ_FIELDS_SPLIT(name, iheader, nx,ny,nz,nt, isize,params, a, wrk)

  USE TLAB_VARS,ONLY : imax_total,jmax_total,kmax_total
#ifdef USE_MPI
  USE TLabMPI_VARS
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) name
  TINTEGER iheader, nx,ny,nz,nt, isize
  TREAL, DIMENSION(isize)            :: params
  TREAL, DIMENSION(nx*ny*nz), TARGET :: a, wrk

! -------------------------------------------------------------------
  TINTEGER header_offset

#ifdef USE_MPI
#ifdef USE_MPI_IO
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE)
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp, mpio_locoff
  TREAL, DIMENSION(:), POINTER :: p_read
  TINTEGER id
#endif

#endif

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( iheader .GT. 0 ) THEN
     IF ( ims_pro .EQ. 0 ) THEN
        OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted',access='stream')
        REWIND(LOC_UNIT_ID)
        CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, imax_total,jmax_total,kmax_total,nt, params)
        CLOSE(LOC_UNIT_ID)
     ENDIF
     CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)

! Displacement to start of field
     mpio_disp = header_offset*SIZEOFBYTE

! Size of array params
     isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL
     CALL MPI_BCAST(params, isize, MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)

  ELSE
     mpio_disp = 0

  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)

  IF ( ims_npro_i .GT. 1 ) THEN; p_read => wrk
  ELSE;                          p_read => a; ENDIF

  mpio_locsize = nx*ny*nz
  mpio_locoff  = mpio_locsize         ! mpio_locoff might be of type larger than INT4
  mpio_locoff  = ims_pro*mpio_locoff  ! mpio_locoff might be of type larger than INT4
  CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_FILE_READ_AT_ALL(mpio_fh, mpio_locoff, p_read, mpio_locsize, MPI_REAL8, status, ims_err)

!  IF ( ims_npro_i .GT. 1 ) THEN
!     id  = TLabMPI_I_PARTIAL
!     CALL TLabMPI_TRPB_I(p_read, a, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
!  ENDIF

  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
  NULLIFY(p_read)

#endif

#else
! ###################################################################
! Serial case
! ###################################################################
  OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted',access='stream')
  REWIND(LOC_UNIT_ID)

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( iheader .GT. 0 ) THEN
     CALL IO_READ_HEADER(LOC_UNIT_ID, header_offset, imax_total,jmax_total,kmax_total,nt, params)
! Size of array params
     isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  READ(LOC_UNIT_ID) a

  CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_READ_FIELDS_SPLIT

#undef LOC_UNIT_ID
#undef LOC_STATUS

!########################################################################
!# Tool/Library DNS
!#
!########################################################################
!# HISTORY
!#
!# 2010/04/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

SUBROUTINE IO_WRITE_FIELDS_SPLIT(name, iheader, nx,ny,nz,nt, isize,params, a, wrk)

  USE TLAB_VARS, ONLY : imax_total,jmax_total,kmax_total
#ifdef USE_MPI
  USE TLabMPI_VARS
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) name
  TINTEGER iheader, isize, nx,ny,nz,nt
  TREAL, DIMENSION(isize)            :: params
  TREAL, DIMENSION(nx*ny*nz), TARGET :: a, wrk

! -------------------------------------------------------------------
#ifdef USE_MPI
#ifdef USE_MPI_IO
  INTEGER mpio_fh, mpio_locsize, status(MPI_STATUS_SIZE), header_offset
  INTEGER(KIND=MPI_OFFSET_KIND) mpio_disp, mpio_locoff
  TREAL, DIMENSION(:), POINTER :: p_write
  TINTEGER id
#endif

#endif

#ifdef USE_MPI
#ifdef USE_MPI_IO
! ###################################################################
! Use MPI_IO for restart files
! ###################################################################
! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( ims_pro .EQ. 0 ) THEN
     OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted',access='stream')
     IF ( iheader .GT. 0 ) THEN
        CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, imax_total,jmax_total,kmax_total,nt, params)

! Displacement to start of field
        header_offset = 5*SIZEOFINT + isize*SIZEOFREAL

     ELSE
        header_offset = 0

     ENDIF
     CLOSE(LOC_UNIT_ID)
  ENDIF

  CALL MPI_BCAST(header_offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
  mpio_disp = header_offset*SIZEOFBYTE

! This is an MPI barrier ! do not remove
  CALL MPI_BARRIER(MPI_COMM_WORLD, ims_err)

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  IF ( ims_npro_i .GT. 1 ) THEN
!     id  = TLabMPI_I_PARTIAL
!     CALL TLabMPI_TRPF_I(a, wrk, ims_ds_k(1,id), ims_dr_k(1,id), ims_ts_k(1,id), ims_tr_k(1,id))
     p_write => wrk
  ELSE
     p_write => a
  ENDIF

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, name, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)

  mpio_locsize = nx*ny*nz
  mpio_locoff  = mpio_locsize         ! reclen might be of type larger than INT4
  mpio_locoff  = ims_pro*mpio_locoff  ! reclen might be of type larger than INT4

  CALL MPI_FILE_SET_VIEW(mpio_fh, mpio_disp, MPI_REAL8, MPI_REAL8, 'native', MPI_INFO_NULL, ims_err)
  CALL MPI_FILE_WRITE_AT_ALL(mpio_fh, mpio_locoff, p_write, mpio_locsize, MPI_REAL8, status, ims_err)
  CALL MPI_FILE_CLOSE(mpio_fh, ims_err)
  NULLIFY(p_write)

#endif

#else
! ###################################################################
! Serial case
! ###################################################################
  OPEN(LOC_UNIT_ID,file=name,status=LOC_STATUS,form='unformatted',access='stream')

! -------------------------------------------------------------------
! header
! -------------------------------------------------------------------
  IF ( iheader .GT. 0 ) THEN
     CALL IO_WRITE_HEADER(LOC_UNIT_ID, isize, imax_total,jmax_total,kmax_total,nt, params)
  ENDIF

! -------------------------------------------------------------------
! fields
! -------------------------------------------------------------------
  WRITE(LOC_UNIT_ID) a

  CLOSE(LOC_UNIT_ID)

#endif

  RETURN
END SUBROUTINE IO_WRITE_FIELDS_SPLIT

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2010/04/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE IO_READ_HEADER(unit, offset, nx,ny,nz,nt, params)

  IMPLICIT NONE

  TINTEGER unit, offset, nx,ny,nz,nt
  TREAL, DIMENSION(*) :: params

! -------------------------------------------------------------------
  TINTEGER isize, nx_loc, ny_loc, nz_loc, nt_loc

!########################################################################
  READ(unit) offset, nx_loc, ny_loc, nz_loc, nt_loc

  isize = offset - 5*SIZEOFINT
!  IF ( isize .GT. 0 .AND. MOD(isize,SIZEOFREAL) .EQ. 0 ) THEN
     isize = isize/SIZEOFREAL
     READ(unit) params(1:isize)

!  ELSE
!     CALL TLab_Write_ASCII(efile,'IO_READ_HEADER. Header format incorrect.')
!     CALL TLab_Stop(DNS_ERROR_RECLEN)

!  ENDIF

! Check
!  IF ( nx .NE. nx_loc .OR. ny .NE. ny_loc .OR. nz .NE. nz_loc ) THEN
!     CLOSE(unit)
!     CALL TLab_Write_ASCII(efile, 'IO_READ_HEADER: Grid size mismatch')
!     CALL TLab_Stop(DNS_ERROR_DIMGRID)
!  ENDIF

!  IF ( nt .NE. nt_loc ) THEN
!     CALL TLab_Write_ASCII(wfile, 'IO_READ_HEADER: ItNumber size mismatch')
!     nt = nt_loc
!  ENDIF

  RETURN
END SUBROUTINE IO_READ_HEADER

!########################################################################
!# Tool/Library IO
!#
!########################################################################
!# HISTORY
!#
!# 2010/04/04 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS
!#
!########################################################################
SUBROUTINE IO_WRITE_HEADER(unit, isize, nx,ny,nz,nt, params)

  IMPLICIT NONE

  TINTEGER unit, isize, nx,ny,nz,nt
  TREAL, DIMENSION(isize) :: params

! -------------------------------------------------------------------
  TINTEGER offset

!########################################################################
  offset = 5*SIZEOFINT + isize*SIZEOFREAL

  WRITE(unit) offset, nx, ny, nz, nt
  WRITE(unit) params(1:isize)

  RETURN
END SUBROUTINE IO_WRITE_HEADER
