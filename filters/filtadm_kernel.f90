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
SUBROUTINE FILTADM_KERNEL(imax, jkmax, periodic, u, uf, tmp, a)

  IMPLICIT NONE

  LOGICAL periodic
  TINTEGER imax, jkmax
  TREAL, DIMENSION(jkmax,imax) :: u, uf, tmp
  TREAL, DIMENSION(imax,5)     :: a

! -------------------------------------------------------------------
  TINTEGER i, jk

! #######################################################################
  CALL FILT4E_KERNEL(imax, jkmax, periodic, u,  uf,  a)
  CALL FILT4E_KERNEL(imax, jkmax, periodic, uf, tmp, a)

  DO i = 1,imax
     DO jk = 1,jkmax
        tmp(jk,i) = tmp(jk,i) + C_3_R*( u(jk,i) - uf(jk,i) )
     ENDDO
  ENDDO

  CALL FILT4E_KERNEL(imax, jkmax, periodic, tmp, uf, a)
  
  RETURN
END SUBROUTINE FILTADM_KERNEL
