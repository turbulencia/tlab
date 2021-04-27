#include "types.h"
#include "dns_error.h"

!########################################################################
!#
!# Write N profiles in netCDF format
!#
!########################################################################
SUBROUTINE IO_WRITE_AVERAGES( fname, itime,rtime, nv,ny, y, varname, vargroup, var )

  USE DNS_CONSTANTS, ONLY : lfile
  USE NETCDF

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  CHARACTER(LEN=*), INTENT(IN   ) :: fname, varname(nv), vargroup(nv) 
  TINTEGER,         INTENT(IN   ) :: itime
  TREAL,            INTENT(IN   ) :: rtime
  TINTEGER,         INTENT(IN   ) :: ny,nv
  TREAL,            INTENT(IN   ) :: y(ny)
  TREAL,            INTENT(IN   ) :: var(ny,nv)

  ! -------------------------------------------------------------------
  INTEGER iv
  INTEGER fid, dtid,dyid, tid,yid,itid, vid(nv)

#ifdef USE_MPI
  INTEGER ims_pro, ims_err
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,ims_pro,ims_err)
#endif

  ! ###################################################################
  CALL IO_WRITE_ASCII(lfile,'Writing '//TRIM(ADJUSTL(fname))//'...')

#ifdef USE_MPI
  IF ( ims_pro == 0 ) THEN
#endif
    CALL NC_CHECK( NF90_CREATE( TRIM(ADJUSTL(fname))//'.nc', NF90_NETCDF4, fid ) )

    CALL NC_CHECK( NF90_DEF_DIM( fid, "t", NF90_UNLIMITED, dtid) )
    CALL NC_CHECK( NF90_DEF_DIM( fid, "y", ny,             dyid) )

    CALL NC_CHECK( Nf90_DEF_VAR( fid, "t",  NF90_FLOAT, (/ dtid /), tid) )
    CALL NC_CHECK( Nf90_DEF_VAR( fid, "y",  NF90_FLOAT, (/ dyid /), yid) )
    CALL NC_CHECK( Nf90_DEF_VAR( fid, "it", NF90_INT,   (/ dtid /),itid) )
    DO iv = 1,nv
      CALL NC_CHECK( Nf90_DEF_VAR( fid, TRIM(ADJUSTL(varname(iv))), NF90_FLOAT, (/ dyid, dtid /), vid(iv)) )
      CALL NC_CHECK( NF90_PUT_ATT( fid, vid(iv), 'group', vargroup(iv) ) )
    END DO

    CALL NC_CHECK( NF90_ENDDEF( fid ) )

    CALL NC_CHECK( NF90_PUT_VAR( fid, tid, SNGL(rtime) ) )
    CALL NC_CHECK( NF90_PUT_VAR( fid,itid, itime ) )
    CALL NC_CHECK( NF90_PUT_VAR( fid, yid, SNGL(y)     ) )
    DO iv = 1,nv
      CALL NC_CHECK( NF90_PUT_VAR( fid, vid(iv), SNGL(var(1:ny,iv)) ) )
    END DO

    CALL NC_CHECK( NF90_CLOSE( fid ) )

#ifdef USE_MPI
  END IF
#endif

  RETURN
END SUBROUTINE IO_WRITE_AVERAGES

! ###################################################################
SUBROUTINE NC_CHECK( status )

  USE DNS_CONSTANTS, ONLY : efile
  USE NETCDF

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: status

  IF ( status /= nf90_noerr ) THEN
    CALL IO_WRITE_ASCII(efile,'NETCDF error signal '//TRIM(ADJUSTL(NF90_STRERROR(status))))
    CALL DNS_STOP(DNS_ERROR_UNDEVELOP)
  END IF

  RETURN
END SUBROUTINE NC_CHECK
