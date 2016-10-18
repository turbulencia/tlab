#include "types.h"
#include "dns_const.h"
#include "avgij_map.h"

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
SUBROUTINE STATS_SPATIAL_LAYER(vaux, txc, wrk1d,wrk2d,wrk3d)

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

  TREAL, DIMENSION(*) :: txc, vaux, wrk1d, wrk2d, wrk3d

! -----------------------------------------------------------------------
  TINTEGER is

! Pointers to existing allocated space
  TREAL, DIMENSION(:), POINTER :: x,y,z, dx,dy,dz

! #######################################################################
#ifdef TRACE_ON
  CALL IO_WRITE_ASCII(tfile, 'ENTERING STATS_SPATIAL_LAYER' )
#endif

! Define pointers
  x => g(1)%nodes; dx => g(1)%jac(:,1)
  y => g(2)%nodes; dy => g(2)%jac(:,1)
  z => g(3)%nodes; dz => g(3)%jac(:,1)

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
