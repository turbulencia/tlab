!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2007/01/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE STATS_SPATIAL_LAYER(x,y,z,dx,dy,dz, vaux, txc, wrk1d,wrk2d,wrk3d)

#include "types.h"
#include "dns_const.h"
#include "avgij_map.h"

  USE DNS_GLOBAL
  USE DNS_LOCAL
#ifdef USE_MPI
  USE DNS_MPI
#endif
#ifdef LES
  USE LES_GLOBAL, ONLY : iles
#endif

  IMPLICIT NONE

#ifdef USE_MPI
#include "mpif.h"
#endif

  TREAL, DIMENSION(*) :: x, y, z, dx, dy, dz
  TREAL, DIMENSION(*) :: txc, vaux, wrk1d, wrk2d, wrk3d

! -----------------------------------------------------------------------
  TINTEGER is

! #######################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING STATS_SPATIAL_LAYER' )
#endif

! #######################################################################
! Averages
! #######################################################################
  IF ( fstavg .EQ. 1 .AND. frunstat .EQ. 1 ) THEN
#ifdef USE_MPI
     IF ( ims_pro .EQ. 0 ) THEN
#endif
        CALL AVG_FLOW_SPATIAL_LAYER(isize_txc, buff_u_jmin,buff_u_jmax, &
             x,y,dy, vaux(vindex(VA_MEAN_WRK)), txc, wrk1d,wrk2d)
        
        IF ( icalc_scal .EQ. 1 ) THEN
           DO is = 1,inb_scal
              CALL AVG_SCAL_SPATIAL_LAYER(is, isize_txc, buff_u_jmin,buff_u_jmax, &
                   x,y,dy, vaux(vindex(VA_MEAN_WRK)), &
                   vaux(vindex(VA_MEAN_WRK)+MA_MOMENTUM_SIZE+MS_SCALAR_SIZE*(is-1)), txc, wrk1d)
           ENDDO
        ENDIF

#ifdef LES               
        IF ( iles .EQ. 1 ) THEN
           CALL LES_AVG_SPATIAL_LAYER(isize_txc, x,y, vaux(vindex(VA_MEAN_WRK)), txc, wrk1d,wrk2d)
        ENDIF
#endif

#ifdef USE_MPI
     ENDIF
#endif
  ENDIF

#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'LEAVING STATS_SPATIAL_LAYER' )
#endif

  RETURN
END SUBROUTINE STATS_SPATIAL_LAYER
