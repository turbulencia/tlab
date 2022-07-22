#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/03/25 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#  Cubic spline function with different BCs.
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

subroutine CUBIC_SPLINE(bc, bcval, norg, nint, xorg, yorg, xint, yint, wrk)

  implicit none
  
#include "integers.h"

  TINTEGER, dimension(2     ),  intent(in)   :: bc         ! boundary condition for both endpoints
                                                           ! bc=0 periodic
                                                           ! bc=1 clamped
                                                           ! bc=2 fixed1
                                                           ! bc=3 natural
                                                           ! bc=4 fixed2
  TREAL,    dimension(2     ),  intent(in)   :: bcval      ! boundary values of 1st or 2nd deriv. at endpoints 
  TINTEGER,                     intent(in)   :: norg, nint ! data sizes
  TREAL,    dimension(norg  ),  intent(in)   :: xorg, yorg ! original data
  TREAL,    dimension(nint  ),  intent(in)   :: xint       ! interpolated
  TREAL,    dimension(nint  ),  intent(out)  :: yint       ! interpolated
  TREAL,    dimension(norg,11), intent(inout):: wrk        ! scratch

  target                                     :: wrk

  ! -------------------------------------------------------------------
  TREAL,    dimension(:),       pointer      :: rhs                ! rhs and solution
  TREAL,    dimension(:),       pointer      :: dx                 ! step size
  TREAL,    dimension(:),       pointer      :: a, b, c, d         ! spline coeff
  TREAL,    dimension(:),       pointer      :: aa, bb, cc, dd, ee ! diagonals 
  TINTEGER                                   :: i 
  TREAL                                      :: wrk_dummy 

  ! ###################################################################
  ! assigning pointers to scratch
  wrk(:,:) = C_0_R        
  dx  => wrk(1:norg-1,1 ); rhs => wrk(1:norg  ,2 )
  a   => wrk(1:norg-1,3 ); b   => wrk(1:norg-1,4 ); c  => wrk(1:norg-1,5); d => wrk(1:norg-1,6)
  aa  => wrk(1:norg  ,7 ); bb  => wrk(1:norg  ,8 ); cc => wrk(1:norg  ,9)
  dd  => wrk(1:norg-1,10); ee  => wrk(1:norg-1,11)

  ! step size (no need to be uniform, even for periodic BCs)
  do i = 1,norg-i1
    dx(i) = xorg(i+1) - xorg(i)
  end do

  ! check input data
#ifdef _DEBUG
  call CUBIC_SPLINE_CHECK_INPUT(bc, bcval, norg, nint, xorg, xint, dx)
#endif

  ! build diagonals and rhs
  call CUBIC_SPLINE_LHS(bc,norg,dx,aa,bb,cc)
  call CUBIC_SPLINE_RHS(bc,bcval,norg,dx,yorg,rhs)

  ! solve tridiagonal matrix with thomas algorithm
  if ( bc(1) > 0 .and. bc(2) > 0 ) then
    call TRIDFS(norg,   aa,bb,cc    )
    call TRIDSS(norg,i1,aa,bb,cc,rhs)
  else ! if periodic BCs
    call TRIDPFS(norg-1,   aa,bb,cc,dd,ee              )
    call TRIDPSS(norg-1,i1,aa,bb,cc,dd,ee,rhs,wrk_dummy)
    rhs(norg) = rhs(1)
  end if

  ! compute spline coefficients
  call CUBIC_SPLINE_COEFF(norg,dx,yorg,rhs,a,b,c,d)

  ! compute yint of cubic spline interpolation
  call CUBIC_SPLINE_FUNC(norg,nint,a,b,c,d,xorg,xint,yint)

  ! disassociate pointers
  nullify(dx,a,b,c,d,aa,bb,cc,rhs)
  
  return
end subroutine CUBIC_SPLINE

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_CHECK_INPUT(bc, bcval,norg, nint, xorg, xint, dx)

  use TLAB_CONSTANTS, only: efile
  use TLAB_PROCS
  
  implicit none
  
  TINTEGER, dimension(2     ), intent(in):: bc
  TREAL,    dimension(2     ), intent(in):: bcval
  TINTEGER,                    intent(in):: norg, nint
  TREAL,    dimension(norg  ), intent(in):: xorg
  TREAL,    dimension(nint  ), intent(in):: xint
  TREAL,    dimension(norg-1), intent(in):: dx
  
  ! -------------------------------------------------------------------
  TINTEGER                               :: i
  
  ! ###################################################################
  ! check data
  if ( norg < 3 ) then 
    call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. At least three data points needed for cubic spline interpolation.')
    call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
  end if 
  if ( (xint(1) < xorg(1)) .or. (xint(nint) > xorg(norg)) ) then
    call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. No extrapolation, just interpolation - check borders of xint.')
    call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
  end if
  do i = 1, norg-1
    if ( dx(i) < 0 ) then
      call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. x must be strictly increasing.')
      call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
    end if
  end do

  ! check boundary conditions
  if (( bc(1) < 0 .or. bc(1) > 5 ) .or. ( bc(2) < 0 .or. bc(2) > 5 )) then
    call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. Wrong choice of boundary conditions.')
    call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
  end if 
  if (( bc(1) == CS_BCS_PERIODIC .and. bc(2) /= CS_BCS_PERIODIC ) .or. &
      ( bc(2) == CS_BCS_PERIODIC .and. bc(1) /= CS_BCS_PERIODIC )) then
    call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. Periodic boundary conditions only possible for both endpoints.')
    call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
  end if 
  do i = 1,2
    if (( bc(i) == CS_BCS_PERIODIC .and. bcval(i) /= 0 ) .or. &
        ( bc(i) == CS_BCS_CLAMPED  .and. bcval(i) /= 0 ) .or. &
        ( bc(i) == CS_BCS_NATURAL  .and. bcval(i) /= 0 )) then
      call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. Wrong choice of combination of boundary condition and derivative values.')
      call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
    endif
  end do

  return
end subroutine CUBIC_SPLINE_CHECK_INPUT

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_LHS(bc, n, dx, a, b, c)
  
  implicit none

  TINTEGER, dimension(2  ), intent(in) :: bc    
  TINTEGER,                 intent(in) :: n
  TREAL,    dimension(n-1), intent(in) :: dx
  TREAL,    dimension(n  ), intent(out):: a,b,c
  
  ! -------------------------------------------------------------------
  TINTEGER                            :: i
  
  ! ###################################################################
  ! left boundary point
  select case( bc(1) )
  case( CS_BCS_CLAMPED, CS_BCS_FIXED_1 )
    a(1) = C_0_R
    c(1) = C_1_R
  case( CS_BCS_NATURAL, CS_BCS_FIXED_2 )
    a(1) = C_0_R
    c(1) = C_0_R
  case( CS_BCS_PERIODIC )
    a(1) = dx(n-1) / (dx(n-1) + dx(1))
    c(1) = dx(  1) / (dx(n-1) + dx(1))
  end select
  
  ! inner points
  b(1) = C_2_R
  do i = 2,n-1
    a(i) = dx(i-1) / (dx(i-1) + dx(i))
    b(i) = C_2_R
    c(i) = dx(i  ) / (dx(i-1) + dx(i))
  end do
  b(n) = C_2_R
  
  ! right boundary point
  select case( bc(2) )
  case( CS_BCS_CLAMPED, CS_BCS_FIXED_1 )
    a(n) = C_1_R
    c(n) = C_0_R
  case( CS_BCS_NATURAL, CS_BCS_FIXED_2 )
    a(n) = C_0_R
    c(n) = C_0_R
  case( CS_BCS_PERIODIC)
    a(n-1) = dx(n-2) / (dx(n-1) + dx(n-2))
    c(n-1) = dx(n-1) / (dx(n-1) + dx(n-2))
  end select
  
  return
end subroutine CUBIC_SPLINE_LHS

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_RHS(bc, bcval, n, dx, y, rhs)
  
  implicit none
  
  TINTEGER, dimension(2  ), intent(in) :: bc    
  TREAL,    dimension(2  ), intent(in) :: bcval  
  TINTEGER,                 intent(in) :: n
  TREAL,   dimension(n-1),  intent(in) :: dx
  TREAL,   dimension(n  ),  intent(in) :: y
  TREAL,   dimension(n  ),  intent(out):: rhs
  
  ! -------------------------------------------------------------------
  TINTEGER                            :: i
  
  ! ###################################################################
  ! left boundary point
  select case( bc(1) )
  case( CS_BCS_CLAMPED )
    rhs(1) = (C_6_R / dx(  1)) *               (y(  2) - y(1)) / dx(  1)
  case( CS_BCS_FIXED_1 )
    rhs(1) = (C_6_R / dx(  1)) * (-bcval(1) + ((y(  2) - y(1)) / dx(  1)))
  case( CS_BCS_NATURAL )
    rhs(1) = C_0_R
  case( CS_BCS_FIXED_2 )
    rhs(1) = C_2_R * bcval(1)
  case( CS_BCS_PERIODIC)
    rhs(1) = (y(2) - y(1)) / dx(1) - (y(1) - y(n-1)) / dx(n-1)
    rhs(1) = C_6_R * rhs(1) / (dx(1) + dx(n-1))
  end select

  ! inner points
  do i = 2,n-1
    rhs(i) = (y(i+1) - y(i)) / dx(i) - (y(i) - y(i-1)) / dx(i-1)
    rhs(i) = C_6_R * rhs(i) / (dx(i) + dx(i-1))
  end do
  
  ! right boundary point
  select case( bc(2) )
  case( CS_BCS_CLAMPED )
    rhs(n) = (C_6_R / dx(n-1)) *               (y(n-1) - y(n)) / dx(n-1)
  case( CS_BCS_FIXED_1 )
    rhs(n) = (C_6_R / dx(n-1)) * ( bcval(2) + ((y(n-1) - y(n)) / dx(n-1)))
  case( CS_BCS_NATURAL )
    rhs(n) = C_0_R
  case( CS_BCS_FIXED_2 )
    rhs(n) = C_2_R * bcval(2)
  case( CS_BCS_PERIODIC )
    rhs(n-1) = (y(n) - y(n-1)) / dx(n-1) - (y(n-1) - y(n-2)) / dx(n-2)
    rhs(n-1) = C_6_R * rhs(n-1) / (dx(n-1) + dx(n-2))
  end select

  return
end subroutine CUBIC_SPLINE_RHS

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_COEFF(n, dx, y, rhs, a, b, c, d)
  
  implicit none
  
  TINTEGER,                intent(in) :: n
  TREAL,   dimension(n-1), intent(in) :: dx
  TREAL,   dimension(n  ), intent(in) :: y
  TREAL,   dimension(n  ), intent(in) :: rhs
  TREAL,   dimension(n-1), intent(out):: a,b,c,d
  
  ! -------------------------------------------------------------------
  TINTEGER                            :: i
  
  ! ###################################################################
  do i = 1,n-1
    a(i) = (rhs(i+1) - rhs(i)) / (C_6_R * dx(i))
    b(i) = rhs(i) / C_2_R
    c(i) = (y(i+1) - y(i)) / dx(i) - (rhs(i+1) + C_2_R*rhs(i)) * (dx(i) / C_6_R)
    d(i) = y(i)
  end do

  return
end subroutine CUBIC_SPLINE_COEFF

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_FUNC(norg, nint, a, b, c, d, xorg, xint, yint)
  
  implicit none
  
  TINTEGER,                   intent(in) :: norg,nint
  TREAL,   dimension(norg-1), intent(in) :: a,b,c,d
  TREAL,   dimension(norg  ), intent(in) :: xorg
  TREAL,   dimension(nint  ), intent(in) :: xint
  TREAL,   dimension(nint  ), intent(out):: yint
  
  ! -------------------------------------------------------------------
  TINTEGER                               :: i, idx
  TREAL                                  :: z
  
  ! ###################################################################

  do i = 1,nint
    call CUBIC_SPLINE_BISECT(norg,xorg,xint(i),idx)
    z       = xint(i) - xorg(idx)
    yint(i) = a(idx)*z**3 + b(idx)*z**2 + c(idx)*z + d(idx) 
  end do 

  return
end subroutine CUBIC_SPLINE_FUNC

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_BISECT(n, a, x, idx)
  
  implicit none

#include "integers.h"
  
  TINTEGER,               intent(in) :: n   ! size of array a
  TREAL,    dimension(n), intent(in) :: a   ! ordered array
  TREAL,                  intent(in) :: x   ! element to insert
  TINTEGER,               intent(out):: idx ! index for insertion
  
  ! -------------------------------------------------------------------
  TINTEGER                           :: mid, lo, hi
  
  ! ###################################################################

  lo = i1; hi = n

  do while ( lo < hi )
    mid = (lo + hi) / i2
    if (x < a(mid)) then
      hi = mid 
    else
      lo = mid + i1
    end if
  end do      

  idx = lo - i1

  return
end subroutine CUBIC_SPLINE_BISECT