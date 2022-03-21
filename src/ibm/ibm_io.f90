#include "types.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/XX/XX - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#  
!#  
!#    
!#    
!#
!# 
!########################################################################
!# ARGUMENTS 
!#
!# 
!#                            
!#                           
!#
!########################################################################
!# REQUIREMENTS
!#
!# 
!#                            
!#                           
!#
!########################################################################
subroutine IBM_IO_READ(wrk3d)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : imax,jmax,kmax
  use IO_FIELDS

  implicit none

#include "integers.h"
  
  TREAL, dimension(*), intent(inout) ::  wrk3d

  character(len=32)                  :: fname

  ! ================================================================== !
  ! read eps field (filename 'eps0.1')
  write(fname,*) i0; 
  fname = trim(adjustl('eps'))//trim(adjustl(fname))
  
  call IO_READ_FIELDS(fname, IO_FLOW, imax,jmax,kmax, i1, i1, eps, wrk3d)
  
  return
end subroutine IBM_IO_READ
!########################################################################
! subroutine IBM_IO_READ_INT1(a, wrk3d)
  
!   use TLAB_VARS,                   only : imax,jmax,kmax
!   use TLAB_CONSTANTS,              only : dp
!   use IO_FIELDS
!   use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

!   implicit none

! #include "integers.h"
  
!   TREAL, dimension(imax*jmax*kmax), target, intent(  out) ::  a
!   TREAL, dimension(imax*jmax*kmax), target, intent(inout) ::  wrk3d

!   integer(1), pointer                                     :: int_wrk(:) => null()
!   real(dp),   pointer                                     :: dp_a(:)    => null()
!   TINTEGER,   parameter                                   :: isize=1
!   TREAL                                                   :: params(isize)
!   ! character(len=32)                                    :: fname

!   ! ================================================================== !
!   ! write eps field (filename 'eps0.1') as int(1) without header
  
!   ! write(fname,*) i0; 
!   ! fname = trim(adjustl('eps'))//trim(adjustl(fname))

!   ! Pass memory address from double precision array to int1 array
!   call c_f_pointer(c_loc(wrk3d), int_wrk, shape=[imax*jmax*kmax])
!   int_wrk(:) = int(wrk3d(:), 1)
  
!   call IO_WRITE_FIELD_INT1('eps0.1', i0, imax,jmax,kmax, i0, isize, params, int_wrk)

!   ! Pass memory address from int1 array to double precision
!   call c_f_pointer(c_loc(wrk3d), dp_a, shape=[imax*jmax*kmax])
!   a(:) = real(dp_a(:), dp)

!   return
! end subroutine IBM_IO_READ_INT1
!########################################################################
subroutine IBM_IO_WRITE(wrk3d)
  
  use DNS_IBM,   only : eps
  use TLAB_VARS, only : imax,jmax,kmax
  use IO_FIELDS

  implicit none

#include "integers.h"
  
  TREAL, dimension(*), intent(inout) ::  wrk3d

  character(len=32)                  :: fname

  ! ================================================================== !
  ! read eps field (filename 'eps0.1')
  write(fname,*) i0; 
  fname = trim(adjustl('eps'))//trim(adjustl(fname))
  
  call IO_WRITE_FIELDS(fname, IO_FLOW, imax,jmax,kmax, i1, eps, wrk3d)

  return
end subroutine IBM_IO_WRITE
!########################################################################
subroutine IBM_IO_WRITE_INT1(a)
  
  use TLAB_VARS,                   only : imax,jmax,kmax
  use IO_FIELDS
  use, intrinsic :: ISO_C_binding, only : c_f_pointer, c_loc

  implicit none

#include "integers.h"
  
  TREAL, dimension(imax*jmax*kmax), target, intent(in) ::  a

  integer(1), pointer                                  :: int_a(:) => null()
  TINTEGER, parameter                                  :: isize=1
  TREAL                                                :: params(isize)
  ! character(len=32)                                    :: fname

  ! ================================================================== !
  ! write eps field (filename 'eps0.1') as int(1) without header

  ! write(fname,*) i0; 
  ! fname = trim(adjustl('eps'))//trim(adjustl(fname))

  ! Pass memory address from double precision array to int1 array
  call c_f_pointer(c_loc(a), int_a, shape=[imax*jmax*kmax])
  int_a(:) = int(a(:), 1)

  call IO_WRITE_FIELD_INT1('eps0.1', i0, imax,jmax,kmax, i0, isize, params, int_a)

  return
end subroutine IBM_IO_WRITE_INT1