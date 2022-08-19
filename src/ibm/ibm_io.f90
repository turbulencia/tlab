#include "types.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#   IO of geometry field eps.
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

subroutine IBM_IO_READ_INT_GEOMETRY(wrk3d)
  
    use TLAB_TYPES,              only : dp
    use TLAB_VARS,                   only : imax,jmax,kmax, isize_field
  use IBM_VARS,                    only : eps, eps_name
  use IO_FIELDS
  use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

  implicit none

#include "integers.h"
  
  TREAL, dimension(isize_field), target, intent(inout) ::  wrk3d

  integer(1), pointer                                  :: int_wrk(:) => null()
  TINTEGER,   parameter                                :: params_size = 1
  TREAL                                                :: params(params_size)
  TINTEGER                                             :: isize
  ! ================================================================== !
  
  ! pass memory address from double precision array to int1 array
  call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax*jmax*kmax])
  int_wrk(:) = int(wrk3d(:), 1)
  
  ! header without params
  isize = i0 

  ! read eps field as int(1)
  call IO_READ_FIELD_INT1(eps_name, i1, imax,jmax,kmax, i0, isize, params, int_wrk)

  ! type casting
  eps(:) = real(int_wrk(:), dp)

  return
end subroutine IBM_IO_READ_INT_GEOMETRY

!########################################################################

subroutine IBM_IO_WRITE_INT_GEOMETRY(wrk3d)
  
  use IBM_VARS,  only : eps, eps_name
  use TLAB_VARS, only : isize_field, imax,jmax,kmax
  use IO_FIELDS

  implicit none

#include "integers.h"
  
  integer(1), dimension(isize_field), intent(inout) :: wrk3d

  TINTEGER, parameter                               :: param_size = 1
  TINTEGER                                          :: isize
  TREAL                                             :: params(param_size)
  ! ================================================================== !

  ! dp to int1
  wrk3d(:)  = int(eps(:), 1) 
 
  ! header (offset, nx, ny, nz, nt == 20 byte)
  isize = i0 ! header without params
  
  ! write eps field as int(1) 
  call IO_WRITE_FIELD_INT1(eps_name, i1, imax,jmax,kmax, i0, isize, params, wrk3d)

  return
end subroutine IBM_IO_WRITE_INT_GEOMETRY

!########################################################################

subroutine IBM_IO_WRITE_BIT_GEOMETRY(wrk3d)
  
  use IBM_VARS,  only : eps, eps_name
  use TLAB_VARS, only : isize_field, imax,jmax,kmax
  use IO_FIELDS

  implicit none

#include "integers.h"
  
  integer(1), dimension(isize_field), target, intent(inout) :: wrk3d

  integer(1), dimension(:),           pointer               :: eps_bit 

  TINTEGER, parameter                                       :: param_size = 1
  TINTEGER                                                  :: isize, bsize_field, imax_bit
  TREAL                                                     :: params(param_size)
  ! ================================================================== !
  
  ! size of bit-array
  bsize_field = isize_field / i8 
  imax_bit    = imax / i8 ! already checked in IBM_READ_CONSISTENCY_CHECK if possible

  ! assign to scratch
  eps_bit => wrk3d(1:bsize_field)
  
  ! dp real to bitwise int1
  call IBM_IO_R2B(isize_field, bsize_field, eps, eps_bit)
 
  ! header (offset, nx, ny, nz, nt == 20 byte)
  isize = i0 ! header without params

  ! write bitwise eps_bit field as int(1) 
  call IO_WRITE_FIELD_INT1(eps_name, i1, imax_bit,jmax,kmax, i0, isize, params, eps_bit)

  return
end subroutine IBM_IO_WRITE_BIT_GEOMETRY

!########################################################################

subroutine IBM_IO_READ_BIT_GEOMETRY(wrk3d)
  
    use TLAB_TYPES,              only : dp
    use IBM_VARS,                    only : eps, eps_name
  use TLAB_VARS,                   only : isize_field, imax,jmax,kmax
  use IO_FIELDS
  use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

  implicit none

#include "integers.h"
  
  TREAL, dimension(isize_field), target, intent(inout) ::  wrk3d

  integer(1), pointer                                  :: int_wrk(:) => null()
  TINTEGER,   parameter                                :: params_size = 1
  TINTEGER                                             :: isize, bsize_field, imax_bit
  TREAL                                                :: params(params_size)
  ! ================================================================== !
  
  ! size of bit-array
  bsize_field = isize_field / i8 
  imax_bit    = imax / i8 ! already checked in IBM_READ_CONSISTENCY_CHECK if possible

  ! pass memory address from double precision array to int1 array
  call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax_bit*jmax*kmax])
  int_wrk(:) = int(wrk3d(:), 1)
  
  ! header without params
  isize = i0 

  ! read eps field as int(1)
  call IO_READ_FIELD_INT1(eps_name, i1, imax_bit,jmax,kmax, i0, isize, params, int_wrk)

  ! dp real to bitwise int1
  call IBM_IO_B2R(bsize_field, isize_field, int_wrk, eps)

  ! type casting
  ! eps(:) = real(int_wrk(:), dp)

  return
end subroutine IBM_IO_READ_BIT_GEOMETRY

!########################################################################

subroutine IBM_IO_R2B(rsize,bsize,r,b)
  
  implicit none

#include "integers.h"
  
  TINTEGER,                     intent(in)  :: rsize
  TINTEGER,                     intent(in)  :: bsize
  TREAL,      dimension(rsize), intent(in)  :: r 
  integer(1), dimension(bsize), intent(out) :: b

  TINTEGER                                  :: i, ip, ib 
  ! ================================================================== !
  
  ! initialize bit-array
  b(1:bsize) = int(0,1)
  
  ! convert
  do i = 1, bsize
    ip = (i - i1) * i8 
    do ib = 1, 8
      if ( r(ip+ib) > C_0_R ) then 
        b(i) = ibset(b(i), ib-i1)
      end if
    end do
  end do

  return 
end subroutine IBM_IO_R2B

!########################################################################

subroutine IBM_IO_B2R(bsize,rsize,b,r)
  
  implicit none

#include "integers.h"
  
  TINTEGER,                     intent(in)  :: bsize
  TINTEGER,                     intent(in)  :: rsize
  integer(1), dimension(bsize), intent(in)  :: b
  TREAL,      dimension(rsize), intent(out) :: r 

  TINTEGER                                  :: i, ip, ib 
  ! ================================================================== !
  
  ! initialize real-array
  r(:) = C_0_R

  ! convert
  do i = 1, bsize
    ip = (i - i1) * i8 
    do ib = 1, 8
      if ( btest(b(i), ib-i1) ) then
        r(ip+ib) = C_1_R
      end if
    end do
  end do

  return 
end subroutine IBM_IO_B2R

!########################################################################