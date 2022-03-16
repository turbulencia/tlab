#include "types.h"
#include "dns_error.h"

!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/03/15 - J. Kostelecky
!#              Created
!#
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#  Cubic spline function with natural BCs (second derivative at 
!#  curve ends are zero).
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

subroutine CUBIC_SPLINE(norg, nint, xorg, yorg, xint, yint, wrk)

  implicit none
  
#include "integers.h"

  TINTEGER,                    intent(in)   :: norg, nint ! data sizes
  TREAL,    dimension(norg  ), intent(in)   :: xorg, yorg ! original data
  TREAL,    dimension(nint  ), intent(in)   :: xint       ! interpolated
  TREAL,    dimension(nint  ), intent(out)  :: yint       ! interpolated
  TREAL,    dimension(norg,9), intent(inout):: wrk        ! scratch
  target                                    :: wrk

  ! ! -------------------------------------------------------------------
  TREAL,    dimension(:),      pointer      :: rhs        ! rhs and solution
  TREAL,    dimension(:),      pointer      :: dx         ! step size
  TREAL,    dimension(:),      pointer      :: a, b, c, d ! spline coeff
  TREAL,    dimension(:),      pointer      :: aa, bb, cc ! diagonals 
  TINTEGER                                  :: i, ndx 

  ! ###################################################################
  ! assigning pointers to scratch
  wrk(:,:) = C_0_R
  dx => wrk(1:norg-1,1); rhs => wrk(1:norg  ,2)
  a  => wrk(1:norg-1,3); b   => wrk(1:norg-1,4); c  => wrk(1:norg-1,5); d => wrk(1:norg-1,6)
  aa => wrk(1:norg-2,7); bb  => wrk(1:norg-2,8); cc => wrk(1:norg-2,9)
  
  ! step size
  ndx = norg - i1 ! size of dx
  do i = 1,ndx
    dx(i) = xorg(i+1) - xorg(i)
  end do

  ! check input data
#ifdef _DEBUG
  call CUBIC_SPLINE_CHECK_INPUT(norg, nint, xorg, yorg, xint, yint, dx)
#endif

  ! lhs, rhs (no Bcs, added later)
  call CUBIC_SPLINE_LHS(ndx,dx,aa,bb,cc)
  call CUBIC_SPLINE_RHS(ndx,dx,yorg,rhs(2:ndx))

  ! solve tridiagonal matrix with thomas algorithm
  call TRIDFS(ndx-i1,   aa,bb,cc           )
  call TRIDSS(ndx-i1,i1,aa,bb,cc,rhs(2:ndx))

  ! add BCs for natural splines
  rhs(1   ) = C_0_R
  rhs(norg) = C_0_R
  
  ! compute spline coefficients
  call CUBIC_SPLINE_COEFF(ndx,dx,yorg,rhs,a,b,c,d)

  ! compute yint of cubic spline interpolation
  call CUBIC_SPLINE_FUNC(norg,nint,a,b,c,d,xorg,xint,yint)

  nullify(dx,a,b,c,d,aa,bb,cc,rhs)
  
  return
end subroutine CUBIC_SPLINE

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_CHECK_INPUT(norg, nint, xorg, yorg, xint, yint, dx)

  use TLAB_CONSTANTS, only: efile
  use TLAB_PROCS
  
  implicit none
  
  TINTEGER,                    intent(in):: norg, nint
  TREAL,    dimension(norg  ), intent(in):: xorg, yorg
  TREAL,    dimension(nint  ), intent(in):: xint
  TREAL,    dimension(nint  ), intent(in):: yint
  TREAL,    dimension(norg-1), intent(in):: dx
  
  ! -------------------------------------------------------------------
  TINTEGER                               :: i
  
  ! ###################################################################

  if ( norg < 3 ) then 
    call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. At least three data points needed for cubic spline interpolation.')
    call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
  end if 

  if ( (xint(1) < xorg(1)) .or. (xint(nint) > xorg(norg)) ) then
    call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. No extrapolation, just interpolation - check boarders of xint.')
    call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
  end if

  do i = 1, norg-1
    if ( dx(i) < 0 ) then
      call TLAB_WRITE_ASCII(efile, 'CUBIC SPLINE. x must be strictly increasing.')
      call TLAB_STOP(DNS_ERROR_CUBIC_SPLINE)
    end if
  end do

  return
end subroutine CUBIC_SPLINE_CHECK_INPUT

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_LHS(n, dx, a, b, c)
  
  implicit none
  
  TINTEGER,                intent(in) :: n
  TREAL,   dimension(n  ), intent(in) :: dx
  TREAL,   dimension(n-1), intent(out):: a,b,c
  
  ! -------------------------------------------------------------------
  TINTEGER                            :: i
  
  ! ###################################################################
  a(1) = C_0_R
  b(1) = C_2_R
  
  do i = 2,n-1
    a(i  ) = dx(i) / (dx(i  ) + dx(i+1))
    b(i  ) = C_2_R
    c(i-1) = dx(i) / (dx(i-1) + dx(i  ))
  end do
  
  c(n-1) = C_0_R

  return
end subroutine CUBIC_SPLINE_LHS

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_RHS(n, dx, y, rhs)
  
  implicit none
  
  TINTEGER,                intent(in) :: n
  TREAL,   dimension(n  ), intent(in) :: dx
  TREAL,   dimension(n+1), intent(in) :: y
  TREAL,   dimension(n-1), intent(out):: rhs
  
  ! -------------------------------------------------------------------
  TINTEGER                            :: i
  
  ! ###################################################################
  do i = 2,n
    rhs(i-1) = (y(i+1) - y(i)) / dx(i) - (y(i) - y(i-1)) / dx(i-1)
    rhs(i-1) = C_6_R * rhs(i-1) / (dx(i) + dx(i-1))
  end do

  return
end subroutine CUBIC_SPLINE_RHS

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_COEFF(n, dx, y, rhs, a, b, c, d)
  
  implicit none
  
  TINTEGER,                intent(in) :: n
  TREAL,   dimension(n  ), intent(in) :: dx
  TREAL,   dimension(n+1), intent(in) :: y
  TREAL,   dimension(n-1), intent(in) :: rhs
  TREAL,   dimension(n  ), intent(out):: a,b,c,d
  
  ! -------------------------------------------------------------------
  TINTEGER                            :: i
  
  ! ###################################################################
  do i = 1,n
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
  TINTEGER                           :: i, mid, lo, hi
  
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