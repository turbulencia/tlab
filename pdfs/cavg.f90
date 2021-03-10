#include "types.h"

!########################################################################
!#
!# Calcualte {<a|u_i>, i=1,...,nv] in ny planes.
!# A last j-plane is added with the average in all the volume.
!# Derived from PDF1V_N
!#
!########################################################################
SUBROUTINE CAVG1V_N( fname, varname, nx,ny,nz, nv, nbins, ibc, umin,umax,u, igate,gate, a, y, avg, wrk1d )

  USE DNS_TYPES,      ONLY : pointers_dt
  USE DNS_CONSTANTS,  ONLY : lfile
  USE PDFS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname, varname(nv)
  TINTEGER,           INTENT(IN)    :: nx,ny,nz, nv, nbins, ibc(nv) ! ibc=0 for external interval, 1 for local
  TREAL,              INTENT(IN)    :: umin(nv),umax(nv)            ! Random variables
  TYPE(pointers_dt),  INTENT(IN)    :: u(nv)
  INTEGER(1),         INTENT(IN)    :: gate(*), igate               ! discrete conditioning criteria
  TREAL,              INTENT(IN)    :: a(nx*ny*nz)
  TREAL,              INTENT(IN)    :: y(ny)                        ! heights of each plane
  TREAL,              INTENT(OUT)   :: avg(nbins+2,ny+1,nv)         ! last 2 bins contain the interval bounds
  TREAL,              INTENT(INOUT) :: wrk1d(nbins,2)

  ! -------------------------------------------------------------------
  TINTEGER j, iv
  CHARACTER*64 name

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  DO iv = 1,nv

    DO j = 1,ny             ! calculation in planes
      IF ( igate .EQ. 0 ) THEN
        ! CALL CAVG1V2D(  ibc(iv), nx,ny,nz, j,             umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
        CALL PDF1V2D(  ibc(iv), nx,ny,nz, j,             umin(iv),umax(iv),u(iv)%field, nbins,avg(1,j,iv), wrk1d, a,wrk1d(1,2))
      ELSE
        ! CALL CAVG1V2D1G(ibc(iv), nx,ny,nz, j, igate,gate, umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
        CALL PDF1V2D1G(ibc(iv), nx,ny,nz, j, igate,gate, umin(iv),umax(iv),u(iv)%field, nbins,avg(1,j,iv), wrk1d, a,wrk1d(1,2))
      ENDIF
      avg(1:nbins,j,iv) = wrk1d(1:nbins,2)
    ENDDO

    IF ( ny .GT. 1 ) THEN   ! calculation in whole volume, saved as plane ny+1
      j = ny +1

      IF ( igate .EQ. 0 ) THEN
        ! CALL CAVG1V2D(  ibc(iv), nx*ny,1,nz, 1,             umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
        CALL PDF1V2D(  ibc(iv), nx*ny,1,nz, 1,             umin(iv),umax(iv),u(iv)%field, nbins,avg(1,j,iv), wrk1d, a,wrk1d(1,2))
      ELSE
        ! CALL CAVG1V2D1G(ibc(iv), nx*ny,1,nz, 1, igate,gate, umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
        CALL PDF1V2D1G(ibc(iv), nx*ny,1,nz, 1, igate,gate, umin(iv),umax(iv),u(iv)%field, nbins,avg(1,j,iv), wrk1d, a,wrk1d(1,2))
      ENDIF
      avg(1:nbins,j,iv) = wrk1d(1:nbins,2)

    ENDIF

  ENDDO

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
      IF ( ny .GT. 1 ) THEN
        WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(avg(:,:,iv))
      ELSE
        WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(avg(:,1,iv))
      ENDIF
      CLOSE(LOC_UNIT_ID)
    ENDDO

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

END SUBROUTINE CAVG1V_N

!########################################################################
!########################################################################
SUBROUTINE CAVG2V( fname, nx,ny,nz, nbins, u,v, a, y, avg, wrk2d )

  USE DNS_CONSTANTS,  ONLY : lfile
  USE PDFS

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname
  TINTEGER,           INTENT(IN)    :: nx,ny,nz, nbins(2)
  TREAL,              INTENT(IN)    :: u(nx*ny*nz), v(nx*ny*nz), a(nx*ny*nz)
  TREAL,              INTENT(IN)    :: y(ny)
  TREAL,              INTENT(OUT)   :: avg(nbins(1)*nbins(2) +2 +2*nbins(1),ny+1)
  TREAL,              INTENT(INOUT) :: wrk2d(nbins(1)*nbins(2),2)

  ! -------------------------------------------------------------------
  TINTEGER j
  CHARACTER*64 name

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Calculating '//TRIM(ADJUSTL(fname))//'...')

  DO j = 1,ny             ! calculation in planes
    ! CALL CAVG2V2D( nx,ny,nz, j, u,v, a, nbins,wrk2d,avg(1,j), wrk2d(1,2) )
    CALL PDF2V2D( nx,ny,nz, j, u,v, nbins,avg(1,j), wrk2d, a,wrk2d(1,2) )
    avg(1:nbins(1)*nbins(2),j) = wrk2d(1:nbins(1)*nbins(2),2)
  ENDDO

  IF ( ny .GT. 1 ) THEN   ! calculation in whole volume, saved as plane ny+1
    j = ny +1
    ! CALL CAVG2V2D( nx*ny,1,nz, 1, u,v, a, nbins,wrk2d,avg(1,j), wrk2d(1,2) )
    CALL PDF2V2D( nx*ny,1,nz, 1, u,v, nbins,avg(1,j), wrk2d, a,wrk2d(1,2) )
    avg(1:nbins(1)*nbins(2),j) = wrk2d(1:nbins(1)*nbins(2),2)
  ENDIF

  ! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
    name = TRIM(ADJUSTL(fname))
    CALL IO_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')
#include "dns_open_file.h"
    IF ( ny .GT. 1 ) THEN
      WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(avg(:,:))
    ELSE
      WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(avg(:,1))
    ENDIF
    CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

END SUBROUTINE CAVG2V
