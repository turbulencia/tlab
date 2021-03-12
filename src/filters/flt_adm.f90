#include "types.h"

!########################################################################
!# HISTORY
!#
!# 2008/11/25 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Explicit filter after 3 levels of deconvolution
!#
!# uf = [ I + (I-G) + (I-G)*(I-G) ]*G*u = (3G - 3G^2 + G^3)*u 
!#    = G*( 3u - 3G*u + G*G*u )
!#
!########################################################################
SUBROUTINE FLT_ADM(imax, jkmax, periodic, a, u, uf, tmp)

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER imax, jkmax
  TREAL, DIMENSION(jkmax,imax) :: u, uf, tmp
  TREAL, DIMENSION(imax,5)     :: a

! -------------------------------------------------------------------
  TINTEGER i, jk

! #######################################################################
  CALL FLT_E4(imax, jkmax, periodic, u,  uf,  a)
  CALL FLT_E4(imax, jkmax, periodic, uf, tmp, a)

  DO i = 1,imax
     DO jk = 1,jkmax
        tmp(jk,i) = tmp(jk,i) + C_3_R*( u(jk,i) - uf(jk,i) )
     ENDDO
  ENDDO

  CALL FLT_E4(imax, jkmax, periodic, tmp, uf, a)
  
  RETURN
END SUBROUTINE FLT_ADM
