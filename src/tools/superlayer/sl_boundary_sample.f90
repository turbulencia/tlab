!########################################################################
!# Tool/Library SUPERLAYER
!#
!########################################################################
!# HISTORY
!#
!# 2007/09/14 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Sampling the fields given in array b along the surface given by sl
!# Data stored in array c
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
SUBROUTINE SL_BOUNDARY_SAMPLE(imax,jmax,kmax, nfield_loc, nfield, y, sl, b, c)
  
  IMPLICIT NONE

#include "types.h"

  TINTEGER imax, jmax, kmax, nfield_loc, nfield
  TREAL y(jmax)
  TREAL sl(imax,kmax)
  TREAL b(imax,jmax,kmax,nfield_loc)
  TREAL c(nfield,imax,kmax)

! -------------------------------------------------------------------
  TREAL dy_u, dy_loc
  TINTEGER i, k, ifield, jm

! ###################################################################
  dy_u = y(2)-y(1)

! ###################################################################
! Loop on the points
! ###################################################################
  DO k = 1,kmax
     DO i = 1,imax
        jm = INT((sl(i,k)-y(1))/dy_u) + 1
        dy_loc = sl(i,k)-y(jm)
        DO ifield = 1,nfield_loc
           c(ifield,i,k) = b(i,jm,k,ifield) + &
                (b(i,jm+1,k,ifield)-b(i,jm,k,ifield))/(y(jm+1)-y(jm))*dy_loc
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE SL_BOUNDARY_SAMPLE
