#include "types.h"

!########################################################################
!#
!# Calcualte {pdf(u_i), i=1,...,nv] in ny planes. (Histograms, not normalized.)
!# A last j-plane is added with the PDF in all the volume.
!#
!# ibc       In    BCs: 0 homogeneous interval
!#                      1 local interval
!#                      2 local interval, analysis and drop no point
!#                      3 local interval, analysis and drop left point
!#                      4 local interval, analysis and drop right point
!#                      5 local interval, analysis and drop both points
!#
!########################################################################
SUBROUTINE PDF1V_N( fname, time, nx,ny,nz, nv, nbins, ibc, umin,umax,u, igate,gate, y, pdf, wrk1d )

  USE TLAB_TYPES,      ONLY : pointers_dt
  USE DNS_CONSTANTS,  ONLY : lfile
  USE TLAB_PROCS
  USE PDFS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*),     INTENT(IN   ) :: fname
  TREAL,             INTENT(IN   ) :: time
  TINTEGER,          INTENT(IN   ) :: nx,ny,nz, nv, nbins, ibc(nv)
  TREAL,             INTENT(IN   ) :: umin(nv),umax(nv)            ! Random variables
  TYPE(pointers_dt), INTENT(IN   ) :: u(nv)
  INTEGER(1),        INTENT(IN   ) :: gate(*), igate               ! discrete conditioning criteria
  TREAL,             INTENT(IN   ) :: y(ny)                        ! heights of each plane
  TREAL,             INTENT(  OUT) :: pdf(nbins+2,ny+1,nv)         ! last 2 bins contain the interval bounds
  TREAL,             INTENT(INOUT) :: wrk1d(nbins)

  ! -------------------------------------------------------------------
  TINTEGER iv, j, nplim, ibc_loc
  TREAL plim, umin_loc, umax_loc

  CHARACTER*64 name

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  plim = C_1EM4_R                 ! relative threshold in PDF analysis; adapt to sample sizeo

  DO iv = 1,nv

    DO j = 1,ny                   ! calculation in planes
      IF ( igate == 0 ) THEN
        CALL PDF1V2D(  ibc(iv), nx,ny,nz, j,             umin(iv),umax(iv),u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
      ELSE
        CALL PDF1V2D1G(ibc(iv), nx,ny,nz, j, igate,gate, umin(iv),umax(iv),u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
      END IF

      IF ( ibc(iv) > 1 ) THEN
        ibc_loc = ibc(iv)-2; umin_loc = umin(iv); umax_loc = umax(iv)
        CALL PDF_ANALIZE(ibc_loc, nbins,pdf(1,j,iv), umin_loc,umax_loc, plim,nplim)
        IF ( igate == 0 ) THEN
          CALL PDF1V2D(  i0, nx,ny,nz, j,             umin_loc,umax_loc,u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
        ELSE
          CALL PDF1V2D1G(i0, nx,ny,nz, j, igate,gate, umin_loc,umax_loc,u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
        END IF
      END IF

    END DO

    IF ( ny > 1 ) THEN            ! calculation in whole volume, saved as plane j=ny+1
      IF ( igate == 0 ) THEN
        CALL PDF1V2D(  ibc(iv), nx*ny,1,nz, 1,             umin(iv),umax(iv),u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
      ELSE
        CALL PDF1V2D1G(ibc(iv), nx*ny,1,nz, 1, igate,gate, umin(iv),umax(iv),u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
      END IF

      IF ( ibc(iv) > 1 ) THEN
        ibc_loc = ibc(iv)-2; umin_loc = umin(iv); umax_loc = umax(iv)
        CALL PDF_ANALIZE(ibc_loc, nbins,pdf(1,j,iv), umin_loc,umax_loc, plim,nplim)
        IF ( igate == 0 ) THEN
          CALL PDF1V2D(  i0, nx*ny,1,nz, 1,             umin_loc,umax_loc,u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
        ELSE
          CALL PDF1V2D1G(i0, nx*ny,1,nz, 1, igate,gate, umin_loc,umax_loc,u(iv)%field, nbins,pdf(1,j,iv), wrk1d)
        END IF
      END IF

    END IF

  END DO

  ! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro == 0 ) THEN
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
    DO iv = 1,nv
      name = TRIM(ADJUSTL(fname))
      IF ( u(iv)%tag /= '' ) name = TRIM(ADJUSTL(fname))//'.'//TRIM(ADJUSTL(u(iv)%tag))
      CALL TLAB_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')
#include "dns_open_file.h"
      IF ( ny > 1 ) THEN
        WRITE(LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:,:,iv))
      ELSE
        WRITE(LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:,1,iv))
      END IF
      CLOSE(LOC_UNIT_ID)
    END DO

#ifdef USE_MPI
  END IF
#endif

  RETURN

END SUBROUTINE PDF1V_N

!########################################################################
!########################################################################
SUBROUTINE PDF2V( fname, time, nx,ny,nz, nbins, u,v, y, pdf, wrk2d )

  USE DNS_CONSTANTS,  ONLY : lfile
  USE TLAB_PROCS
  USE PDFS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*), INTENT(IN   ) :: fname
  TREAL,         INTENT(IN   ) :: time
  TINTEGER,      INTENT(IN   ) :: nx,ny,nz, nbins(2)
  TREAL,         INTENT(IN   ) :: u(nx*ny*nz), v(nx*ny*nz)
  TREAL,         INTENT(IN   ) :: y(ny)
  TREAL,         INTENT(  OUT) :: pdf(nbins(1)*nbins(2) +2 +2*nbins(1),ny+1)
  TREAL,         INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2))

  ! -------------------------------------------------------------------
  TINTEGER j
  CHARACTER*64 name

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  DO j = 1,ny               ! calculation in planes
    CALL PDF2V2D( nx,ny,nz, j, u,v, nbins,pdf(1,j), wrk2d )
  END DO

  IF ( ny > 1 ) THEN        ! calculation in whole volume, saved as plane ny+1
    CALL PDF2V2D( nx*ny,1,nz, 1, u,v, nbins,pdf(1,j), wrk2d )
  END IF

  ! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro == 0 ) THEN
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
    name = TRIM(ADJUSTL(fname))
    CALL TLAB_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')
#include "dns_open_file.h"
    IF ( ny > 1 ) THEN
      WRITE(LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:,:))
    ELSE
      WRITE(LOC_UNIT_ID) SNGL(time), ny, nbins, SNGL(y(:)), SNGL(pdf(:,1))
    END IF
    CLOSE(LOC_UNIT_ID)
#ifdef USE_MPI
  END IF
#endif

  RETURN

END SUBROUTINE PDF2V
