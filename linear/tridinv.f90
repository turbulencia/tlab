!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2008/08/13 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# Calculate the inverse (lower-triangular matrix) of the unitary 
!# bidiagonal matrix with subdiagonal given by a
!#
!########################################################################
!# ARGUMENTS 
!#
!# a    In    Subdiagonal. First element corresponds to second raw !
!#
!########################################################################
SUBROUTINE TRIDINV(size,a,inv)

  IMPLICIT NONE

#include "types.h"

  TINTEGER size
  TREAL a(size), inv(size,size)

! -----------------------------------------------------------------------
  TINTEGER icol, iraw

! #######################################################################
  DO icol = 1,size-1
     inv(icol+1,icol) =-a(icol)
     DO iraw = icol+2,size
        inv(iraw,icol) =-inv(iraw-1,icol)*a(iraw-1)
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE TRIDINV
