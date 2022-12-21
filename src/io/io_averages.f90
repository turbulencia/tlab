#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!#
!# Write N profiles in netCDF or ASCII format
!#
!########################################################################
SUBROUTINE IO_WRITE_AVERAGES( fname, itime,rtime, ny,nv,ng, y, varnames, groupnames, avg )

  USE TLAB_CONSTANTS, ONLY : lfile,efile
  USE TLAB_PROCS
#ifdef USE_MPI
  USE MPI
#endif
#ifdef USE_NETCDF
  USE NETCDF
#endif
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN   ) :: fname, varnames(ng), groupnames(ng)
  TINTEGER,         INTENT(IN   ) :: itime
  TREAL,            INTENT(IN   ) :: rtime
  TINTEGER,         INTENT(IN   ) :: ny,nv,ng
  TREAL,            INTENT(IN   ) :: y(ny)
  TREAL,            INTENT(IN   ) :: avg(ny,nv)

  ! -------------------------------------------------------------------
  INTEGER iv,ig
#ifdef USE_NETCDF
  INTEGER nvg
  INTEGER fid, dtid,dyid, tid,yid,itid, vid(nv)
  CHARACTER(LEN=32) vname(nv), gname(nv)
#else
  INTEGER j
  CHARACTER*1300 line
#endif

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL TLAB_WRITE_ASCII(lfile,'Writing '//TRIM(ADJUSTL(fname))//'...')

#ifdef USE_MPI
  IF ( ims_pro == 0 ) THEN
#endif

    ! -----------------------------------------------------------------------
    ! Using NetCDF format
#ifdef USE_NETCDF
    iv = 1
    DO ig = 1,ng
      nvg = nv -iv +1
      CALL LIST_STRING( varnames(ig), nvg, vname(iv:nv))
      gname(iv:iv+nvg-1) = groupnames(ig)
      iv = iv +nvg
      IF ( iv > nv ) EXIT
    ENDDO

    CALL NC_CHECK( NF90_CREATE( TRIM(ADJUSTL(fname))//'.nc', NF90_NETCDF4, fid ) )

    CALL NC_CHECK( NF90_DEF_DIM( fid, "t", NF90_UNLIMITED, dtid) )
    CALL NC_CHECK( NF90_DEF_DIM( fid, "y", ny,             dyid) )

    CALL NC_CHECK( Nf90_DEF_VAR( fid, "t",  NF90_FLOAT, (/ dtid /), tid) )
    CALL NC_CHECK( Nf90_DEF_VAR( fid, "y",  NF90_FLOAT, (/ dyid /), yid) )
    CALL NC_CHECK( Nf90_DEF_VAR( fid, "it", NF90_INT,   (/ dtid /),itid) )
    DO iv = 1,nv
      CALL NC_CHECK( Nf90_DEF_VAR( fid, TRIM(ADJUSTL(vname(iv))), NF90_FLOAT, (/ dyid, dtid /), vid(iv)) )
      CALL NC_CHECK( NF90_PUT_ATT( fid, vid(iv), 'group', gname(iv) ) )
    END DO

    CALL NC_CHECK( NF90_ENDDEF( fid ) )

    CALL NC_CHECK( NF90_PUT_VAR( fid, tid, SNGL(rtime) ) )
    CALL NC_CHECK( NF90_PUT_VAR( fid,itid, itime ) )
    CALL NC_CHECK( NF90_PUT_VAR( fid, yid, SNGL(y)     ) )
    DO iv = 1,nv
      CALL NC_CHECK( NF90_PUT_VAR( fid, vid(iv), SNGL(avg(1:ny,iv)) ) )
    END DO

    CALL NC_CHECK( NF90_CLOSE( fid ) )

    ! -----------------------------------------------------------------------
    ! Using ASCII format
#else
#define L_AVGMAX 250

    IF ( L_AVGMAX < nv ) THEN
      CALL TLAB_WRITE_ASCII(efile,'IO_WRITE_AVERAGES. Not enough space in format definition.')
      CALL TLAB_STOP(DNS_ERROR_AVGTMP)
    END IF

#define LOC_UNIT_ID 23
#define LOC_STATUS 'unknown'
#ifdef USE_RECLEN
    OPEN(UNIT=LOC_UNIT_ID, RECL=1050, FILE=fname, STATUS=LOC_STATUS) ! this is probably outdated
#else
    OPEN(UNIT=LOC_UNIT_ID, FILE=fname, STATUS=LOC_STATUS)
#endif

    WRITE(LOC_UNIT_ID, '(A8,E14.7E3)') 'RTIME = ', rtime
    line = ''
    DO ig = 1,ng
      line = TRIM(ADJUSTL(line))//' '//TRIM(ADJUSTL(varnames(ig)))
      WRITE(LOC_UNIT_ID,1010) 'GROUP = '//TRIM(ADJUSTL(groupnames(ig)))//' '//TRIM(ADJUSTL(varnames(ig)))
    ENDDO
    WRITE(LOC_UNIT_ID,1010) 'I J Y '//TRIM(ADJUSTL(line)) ! Independent, dependent variables

    DO j = 1,ny
      WRITE(LOC_UNIT_ID,1020) 1, j, y(j), (avg(j,iv),iv=1,nv)
    END DO

    CLOSE(LOC_UNIT_ID)

1010 FORMAT(A)
1020 FORMAT(I5,(1X,I5),L_AVGMAX(1X,G_FORMAT_R))

#endif

#ifdef USE_MPI
  END IF
#endif

  RETURN
END SUBROUTINE IO_WRITE_AVERAGES

! ###################################################################
#ifdef USE_NETCDF
SUBROUTINE NC_CHECK( status )
  USE TLAB_CONSTANTS, ONLY : efile
  USE TLAB_PROCS
  USE NETCDF

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: status

  IF ( status /= nf90_noerr ) THEN
    CALL TLAB_WRITE_ASCII(efile,'NETCDF error signal '//TRIM(ADJUSTL(NF90_STRERROR(status))))
    CALL TLAB_STOP(DNS_ERROR_UNDEVELOP)
  END IF

  RETURN
END SUBROUTINE NC_CHECK
#endif
