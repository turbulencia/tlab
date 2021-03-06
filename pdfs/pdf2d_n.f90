#include "types.h"
#include "dns_error.h"

#define NVARS_LOC 16

!########################################################################
!#
!# Array gate contains the global intermittency conditioning field.
!# A last j-plane is added containing the PDF constructed using all the volume.
!#
!# igate     In    Gate level. If 0, no intermittency considered
!# nv        In    Number of variables
!# ibc       In    BCs: 0 homogeneous interval
!#                      1 local interval
!#                      2 local interval, analysis and drop no point
!#                      3 local interval, analysis and drop left point
!#                      4 local interval, analysis and drop right point
!#                      5 local interval, analysis and drop both points
!#
!########################################################################
SUBROUTINE PDF2D_N(fname, varname, igate, nx,ny,nz, &
  nv, nbins, ibc, amin,amax, y, gate, data, pdf, wrk1d)

  USE DNS_TYPES,      ONLY : pointers_dt
  USE DNS_CONSTANTS,  ONLY : efile, lfile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname, varname(nv)
  TINTEGER,           INTENT(IN)    :: nx,ny,nz, nv, nbins, ibc(nv)
  TREAL,              INTENT(IN)    :: amin(nv), amax(nv)
  TREAL,              INTENT(IN)    :: y(ny)
  TREAL,              INTENT(OUT)   :: pdf(nbins+2,ny+1,nv)
  TREAL,              INTENT(INOUT) :: wrk1d(nbins)
  INTEGER(1),         INTENT(IN)    :: gate(*), igate
  TYPE(pointers_dt),  INTENT(IN)    :: data(nv)

  ! -------------------------------------------------------------------
  TINTEGER j, ip, nplim, iv, ibc_loc
  TREAL plim

  CHARACTER*512 line1
  CHARACTER*64 name

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  IF ( NVARS_LOC .LT. nv ) THEN
    CALL IO_WRITE_ASCII(efile, 'PDF2D_N. Aux array size too small')
    CALL DNS_STOP(DNS_ERROR_WRKSIZE)
  ENDIF

  ! threshold in the PDF analysis
  ! To be reviewed for small sample sizes because this value might never achieved
  plim = C_1EM4_R

  DO iv = 1,nv
    ! ###################################################################
    ! PDF calculation of 1 variable along planes
    ! ###################################################################
    DO j = 1,ny
      IF ( igate .EQ. 0 ) THEN
        CALL PDF1V2D(  ibc(iv), nx,ny,nz, j,        amin(iv), amax(iv),       data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
      ELSE
        CALL PDF1V2D1G(ibc(iv), nx,ny,nz, j, igate, amin(iv), amax(iv), gate, data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
      ENDIF

      IF ( ibc(iv) .GT. 1 ) THEN ! threshold for analysis set s.t. single points are removed
        ibc_loc = ibc(iv)-2
        CALL PDF_ANALIZE(nbins, ibc_loc, pdf(1,j,iv), plim, amin(iv), amax(iv), nplim)
        IF ( igate .EQ. 0 ) THEN
          CALL PDF1V2D(  i0, nx,ny,nz, j,        amin(iv), amax(iv),       data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
        ELSE
          CALL PDF1V2D1G(i0, nx,ny,nz, j, igate, amin(iv), amax(iv), gate, data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
        ENDIF
      ENDIF

    ENDDO

    ! ###################################################################
    ! PDF calculation of 1 variable in 3D space
    ! ###################################################################
    j = ny +1

    IF ( igate .EQ. 0 ) THEN
      CALL PDF1V3D(  ibc(iv), nx,ny,nz,        amin(iv), amax(iv),       data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
    ELSE
      CALL PDF1V3D1G(ibc(iv), nx,ny,nz, igate, amin(iv), amax(iv), gate, data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
    ENDIF

    IF ( ibc(iv) .GT. 1 ) THEN ! threshold for analysis set s.t. single points are removed
      ibc_loc = ibc(iv)-2
      CALL PDF_ANALIZE(nbins, ibc_loc, pdf(1,j,iv), plim, amin(iv), amax(iv), nplim)
      IF ( igate .EQ. 0 ) THEN
        CALL PDF1V3D(  i0, nx,ny,nz,        amin(iv), amax(iv),       data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
      ELSE
        CALL PDF1V3D1G(i0, nx,ny,nz, igate, amin(iv), amax(iv), gate, data(iv)%field, nbins, pdf(1,j,iv), wrk1d)
      ENDIF
    ENDIF

  ENDDO

  ! ###################################################################
  ! Save to disk
  ! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
    DO iv = 1,nv
      name = TRIM(ADJUSTL(fname))
      IF ( varname(iv) .NE. '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(varname(iv)))

      CALL IO_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')
#include "dns_open_file.h"
      WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(pdf(:,:,iv))
      CLOSE(LOC_UNIT_ID)

    ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

END SUBROUTINE PDF2D_N

!########################################################################
!########################################################################
SUBROUTINE JPDF2D( fname, nx,ny,nz, nbins, y, data, pdf, wrk2d )

  USE DNS_TYPES,      ONLY : pointers_dt
  USE DNS_CONSTANTS,  ONLY : efile, lfile

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname
  TINTEGER,           INTENT(IN)    :: nx,ny,nz, nbins(2)
  TREAL,              INTENT(IN)    :: y(ny)
  TREAL,              INTENT(OUT)   :: pdf(nbins(1)*nbins(2)+4,ny+1)
  TREAL,              INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))
  TYPE(pointers_dt),  INTENT(IN)    :: data(2)

  ! -------------------------------------------------------------------
  TINTEGER j
  CHARACTER*64 name

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  DO j = 1,ny
    CALL PDF2V2D( nx,ny,nz, j, data(1)%field,data(2)%field, nbins, pdf(1,j), wrk2d )
  ENDDO

  ! PDF calculation of 1 variable in 3D space
  j = ny +1
  CALL PDF2V3D( nx,ny,nz, data(1)%field,data(2)%field, nbins, pdf(1,j), wrk2d )

  ! ###################################################################
  ! Save to disk
  ! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
    name = TRIM(ADJUSTL(fname))

    CALL IO_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')
#include "dns_open_file.h"
    WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(pdf(:,:))
    CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

END SUBROUTINE JPDF2D
