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
  
  use TLAB_CONSTANTS,              only : wp, wi
  use TLAB_VARS,                   only : imax,jmax,kmax, isize_field
  use IBM_VARS,                    only : eps, eps_name
  use IO_FIELDS
  use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

  implicit none

  real(wp), dimension(isize_field), target, intent(inout) ::  wrk3d

  integer(1),                       pointer               :: int_wrk(:) => null()
  integer(wi),                      parameter             :: params_size = 1
  real(wp)                                                :: params(params_size)
  integer(wi)                                             :: isize
  ! ================================================================== !
  
  ! pass memory address from double precision array to int1 array
  call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax*jmax*kmax])
  int_wrk(:) = int(wrk3d(:), 1)
  
  ! header without params
  isize = 0 

  ! read eps field as int(1)
  call IO_READ_FIELD_INT1(eps_name, 1, imax,jmax,kmax, 0, isize, params, int_wrk)

  ! type casting
  eps(:) = real(int_wrk(:), wp)

  return
end subroutine IBM_IO_READ_INT_GEOMETRY

!########################################################################

subroutine IBM_IO_WRITE_INT_GEOMETRY(wrk3d)
  
  use IBM_VARS,       only : eps, eps_name
  use TLAB_VARS,      only : isize_field, imax,jmax,kmax
  use TLAB_CONSTANTS, only : wp, wi
  use IO_FIELDS

  implicit none

  integer(1),  dimension(isize_field), intent(inout) :: wrk3d

  integer(wi), parameter                             :: param_size = 1
  integer(wi)                                        :: isize
  real(wp)                                           :: params(param_size)
  ! ================================================================== !

  ! wp to int1
  wrk3d(:) = int(eps(:), 1) 
 
  ! header (offset, nx, ny, nz, nt == 20 byte)
  isize = 0 ! header without params
  
  ! write eps field as int(1) 
  call IO_WRITE_FIELD_INT1(eps_name, 1, imax,jmax,kmax, 0, isize, params, wrk3d)

  return
end subroutine IBM_IO_WRITE_INT_GEOMETRY

!########################################################################

subroutine IBM_IO_WRITE_BIT_GEOMETRY(wrk3d)
  
  use IBM_VARS,       only : eps, eps_name
  use TLAB_VARS,      only : isize_field, imax,jmax,kmax
  use TLAB_CONSTANTS, only : wp, wi
  use IO_FIELDS

  implicit none

  integer(1), dimension(isize_field), target, intent(inout) :: wrk3d

  integer(1), dimension(:),           pointer               :: eps_bit 

  integer(wi), parameter                                    :: param_size = 1
  integer(wi)                                               :: isize, bsize_field, imax_bit
  real(wp)                                                  :: params(param_size)
  ! ================================================================== !
  
  ! size of bit-array
  bsize_field = isize_field / 8 
  imax_bit    = imax / 8 ! already checked in IBM_READ_CONSISTENCY_CHECK if possible

  ! assign to scratch
  eps_bit => wrk3d(1:bsize_field)
  
  ! wp real to bitwise int1
  call IBM_IO_R2B(isize_field, bsize_field, eps, eps_bit)
 
  ! header (offset, nx, ny, nz, nt == 20 byte)
  isize = 0 ! header without params

  ! write bitwise eps_bit field as int(1) 
  call IO_WRITE_FIELD_INT1(eps_name, 1, imax_bit,jmax,kmax, 0, isize, params, eps_bit)

  return
end subroutine IBM_IO_WRITE_BIT_GEOMETRY

!########################################################################

subroutine IBM_IO_READ_BIT_GEOMETRY(wrk3d)
  
  use IBM_VARS,                    only : eps, eps_name
  use TLAB_VARS,                   only : isize_field, imax,jmax,kmax
  use TLAB_CONSTANTS,              only : wp, wi
  use IO_FIELDS
  use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

  implicit none

  real(wp), dimension(isize_field), target, intent(inout) ::  wrk3d

  integer(1),    pointer                                  :: int_wrk(:) => null()
  integer(wi),   parameter                                :: params_size = 1
  integer(wi)                                             :: isize, bsize_field, imax_bit
  real(wp)                                                :: params(params_size)
  ! ================================================================== !
  
  ! size of bit-array
  bsize_field = isize_field / 8 
  imax_bit    = imax / 8 ! already checked in IBM_READ_CONSISTENCY_CHECK if possible

  ! pass memory address from double precision array to int1 array
  call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax_bit*jmax*kmax])
  int_wrk(:) = int(wrk3d(:), 1)
  
  ! header without params
  isize = 0 

  ! read eps field as int(1)
  call IO_READ_FIELD_INT1(eps_name, 1, imax_bit,jmax,kmax, 0, isize, params, int_wrk)

  ! wp real to bitwise int1
  call IBM_IO_B2R(bsize_field, isize_field, int_wrk, eps)

  ! type casting
  ! eps(:) = real(int_wrk(:), dp)

  return
end subroutine IBM_IO_READ_BIT_GEOMETRY

!########################################################################

subroutine IBM_IO_R2B(rsize,bsize,r,b)

  use TLAB_CONSTANTS, only : wp, wi

  implicit none

  integer(wi),                   intent(in)  :: rsize
  integer(wi),                   intent(in)  :: bsize
  real(wp),    dimension(rsize), intent(in)  :: r 
  integer(1),  dimension(bsize), intent(out) :: b

  integer(wi)                                :: i, ip, ib 
  ! ================================================================== !
  
  ! initialize bit-array
  b(1:bsize) = int(0,1)
  
  ! convert
  do i = 1, bsize
    ip = (i - 1) * 8 
    do ib = 1, 8
      if ( r(ip+ib) > 0.0_wp ) then 
        b(i) = ibset(b(i), ib-1)
      end if
    end do
  end do

  return 
end subroutine IBM_IO_R2B

!########################################################################

subroutine IBM_IO_B2R(bsize,rsize,b,r)

  use TLAB_CONSTANTS, only : wp, wi
  
  implicit none
  
  integer(wi),                   intent(in)  :: bsize
  integer(wi),                   intent(in)  :: rsize
  integer(1),  dimension(bsize), intent(in)  :: b
  real(wp),    dimension(rsize), intent(out) :: r 

  integer(wi)                                :: i, ip, ib 
  ! ================================================================== !
  
  ! initialize real-array
  r(:) = 0.0_wp

  ! convert
  do i = 1, bsize
    ip = (i - 1) * 8 
    do ib = 1, 8
      if ( btest(b(i), ib-1) ) then
        r(ip+ib) = 1.0_wp
      end if
    end do
  end do

  return 
end subroutine IBM_IO_B2R

!########################################################################