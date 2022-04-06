#include "types.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   IO of geometry field eps as int1 with header.
!#  
!#    
!########################################################################
!# ARGUMENTS 
!#
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!########################################################################

subroutine IBM_IO_READ_GEOMETRY(wrk3d)
  
  use TLAB_VARS,                   only : imax,jmax,kmax, isize_field
  use DNS_IBM,                     only : eps
  use TLAB_CONSTANTS,              only : dp
  use IO_FIELDS
  use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

  implicit none

#include "integers.h"
  
  TREAL, dimension(isize_field), target, intent(inout) ::  wrk3d

  integer(1), pointer                                  :: int_wrk(:) => null()
  TINTEGER,   parameter                                :: params_size = 1
  TREAL                                                :: params(params_size)
  character(len=32)                                    :: fname
  TINTEGER                                             :: isize

  ! ================================================================== !
  
  ! pass memory address from double precision array to int1 array
  call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax*jmax*kmax])
  int_wrk(:) = int(wrk3d(:), 1)
  
  ! file name
  fname = 'eps0.1'
  isize = i0 ! header without params

  ! read eps field as int(1)
  call IO_READ_FIELD_INT1(fname, i1, imax,jmax,kmax, i0, isize, params, int_wrk)

  ! type casting
  eps(:) = real(int_wrk(:), dp)

  return
end subroutine IBM_IO_READ_GEOMETRY

!#####################read###################################################

subroutine IBM_IO_WRITE_GEOMETRY(wrk3d)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : isize_field, imax,jmax,kmax
  use IO_FIELDS

  implicit none

#include "integers.h"
  
  integer(1), dimension(isize_field), intent(inout) :: wrk3d

  TINTEGER, parameter                               :: param_size = 1
  TINTEGER                                          :: isize
  TREAL                                             :: params(param_size)
  character(len=32)                                 :: fname
  ! ================================================================== !

  ! dp to int1
  wrk3d(:)  = int(eps(:), 1) 
 
  ! header (offset, nx, ny, nz, nt == 20 bits)
  isize = i0 ! header without params
  
  ! name
  fname = 'eps0.1'
  
  ! write eps field as int(1) 
  call IO_WRITE_FIELD_INT1(fname, i1, imax,jmax,kmax, i0, isize, params, wrk3d)

  return
end subroutine IBM_IO_WRITE_GEOMETRY

!########################################################################