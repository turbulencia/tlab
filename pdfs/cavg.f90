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

  IMPLICIT NONE

#include "integers.h"
#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER*(*) fname, varname(nv)
  TINTEGER,           INTENT(IN)    :: nx,ny,nz, nv, nbins, ibc ! ibc=0 for external interval, 1 for local
  TREAL,              INTENT(IN)    :: umin(nv),umax(nv)        ! Random variables
  TYPE(pointers_dt),  INTENT(IN)    :: u(nv)
  INTEGER(1),         INTENT(IN)    :: gate(*), igate           ! discrete conditioning criteria
  TREAL,              INTENT(IN)    :: a(nx*ny*nz)
  TREAL,              INTENT(IN)    :: y(ny)                    ! heights of each plane
  TREAL,              INTENT(OUT)   :: avg(nbins+2,ny+1,nv)     ! last 2 bins contain the interval bounds
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

    DO j = 1,ny ! calculation in planes
      IF ( igate .EQ. 0 ) THEN
        CALL CAVG1V2D(  ibc, nx,ny,nz, j,             umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
      ELSE
        CALL CAVG1V2D1G(ibc, nx,ny,nz, j, igate,gate, umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
      ENDIF
    ENDDO

    j = ny +1     ! calculation in whole volume, saved as plane ny+1
    IF ( igate .EQ. 0 ) THEN
      CALL CAVG1V3D(  ibc, nx,ny,nz,             umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
    ELSE
      CALL CAVG1V3D1G(ibc, nx,ny,nz, igate,gate, umin(iv),umax(iv),u(iv)%field, a, nbins,wrk1d,avg(1,j,iv), wrk1d(1,2))
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
      WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(avg(:,:,iv))
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

  DO j = 1,ny   ! avg calculation along planes
    CALL CAVG2V2D( nx,ny,nz, j, u,v, a, nbins,wrk2d,avg(1,j), wrk2d(1,2) )
  ENDDO

  j = ny +1     ! avg calculation in 3D space
  CALL CAVG2V3D( nx,ny,nz, u,v, a, nbins,wrk2d,avg(1,j), wrk2d(1,2) )

  ! ###################################################################
#ifdef USE_MPI
  IF ( ims_pro .EQ. 0 ) THEN
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
    name = TRIM(ADJUSTL(fname))
    CALL IO_WRITE_ASCII(lfile, 'Writing field '//TRIM(ADJUSTL(name))//'...')
#include "dns_open_file.h"
    WRITE(LOC_UNIT_ID) ny, nbins, SNGL(y(:)), SNGL(avg(:,:))
    CLOSE(LOC_UNIT_ID)

#ifdef USE_MPI
  ENDIF
#endif

  RETURN

END SUBROUTINE CAVG2V
