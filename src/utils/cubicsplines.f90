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

  use TLab_Constants, only : wi, wp

  implicit none
  
  integer(wi), dimension(2     ),  intent(in)   :: bc         ! boundary condition for both endpoints
                                                                ! bc=0 periodic
                                                                ! bc=1 clamped
                                                                ! bc=2 fixed1
                                                                ! bc=3 natural
                                                                ! bc=4 fixed2
  real(wp),    dimension(2     ),  intent(in)   :: bcval      ! boundary values of 1st or 2nd deriv. at endpoints 
  integer(wi),                     intent(in)   :: norg, nint ! data sizes
  real(wp),    dimension(norg  ),  intent(in)   :: xorg, yorg ! original data
  real(wp),    dimension(nint  ),  intent(in)   :: xint       ! interpolated
  real(wp),    dimension(nint  ),  intent(out)  :: yint       ! interpolated
  real(wp),    dimension(norg,11), intent(inout):: wrk        ! scratch

  target                                        :: wrk

  ! -------------------------------------------------------------------
  real(wp),    dimension(:),       pointer      :: rhs                ! rhs and solution
  real(wp),    dimension(:),       pointer      :: dx                 ! step size
  real(wp),    dimension(:),       pointer      :: a, b, c, d         ! spline coeff
  real(wp),    dimension(:),       pointer      :: aa, bb, cc, dd, ee ! diagonals 
  integer(wi)                                   :: i 
  real(wp)                                      :: wrk_dummy 

  ! ###################################################################
  ! assigning pointers to scratch
  wrk(:,:) = 0.0_wp        
  dx  => wrk(1:norg-1,1 ); rhs => wrk(1:norg  ,2 )
  a   => wrk(1:norg-1,3 ); b   => wrk(1:norg-1,4 ); c  => wrk(1:norg-1,5); d => wrk(1:norg-1,6)
  aa  => wrk(1:norg  ,7 ); bb  => wrk(1:norg  ,8 ); cc => wrk(1:norg  ,9)
  dd  => wrk(1:norg-1,10); ee  => wrk(1:norg-1,11)

  ! step size (no need to be uniform, even for periodic BCs)
  do i = 1,norg-1
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
    call TRIDSS(norg,1,aa,bb,cc,rhs)
  else ! if periodic BCs
    call TRIDPFS(norg-1,   aa,bb,cc,dd,ee              )
    call TRIDPSS(norg-1,1,aa,bb,cc,dd,ee,rhs,wrk_dummy)
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

  use TLab_Constants, only: efile, wi, wp
  use TLab_WorkFlow
  
  implicit none
  
  integer(wi), dimension(2     ), intent(in):: bc
  real(wp),    dimension(2     ), intent(in):: bcval
  integer(wi),                    intent(in):: norg, nint
  real(wp),    dimension(norg  ), intent(in):: xorg
  real(wp),    dimension(nint  ), intent(in):: xint
  real(wp),    dimension(norg-1), intent(in):: dx
  
  ! -------------------------------------------------------------------
  integer(wi)                               :: i
  
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

  use TLab_Constants, only : wi, wp
  
  implicit none

  integer(wi), dimension(2  ), intent(in) :: bc    
  integer(wi),                 intent(in) :: n
  real(wp),    dimension(n-1), intent(in) :: dx
  real(wp),    dimension(n  ), intent(out):: a,b,c
  
  ! -------------------------------------------------------------------
  integer(wi)                             :: i
  
  ! ###################################################################
  ! left boundary point
  select case( bc(1) )
  case( CS_BCS_CLAMPED, CS_BCS_FIXED_1 )
    a(1) = 0.0_wp
    c(1) = 1.0_wp
  case( CS_BCS_NATURAL, CS_BCS_FIXED_2 )
    a(1) = 0.0_wp
    c(1) = 0.0_wp
  case( CS_BCS_PERIODIC )
    a(1) = dx(n-1) / (dx(n-1) + dx(1))
    c(1) = dx(  1) / (dx(n-1) + dx(1))
  end select
  
  ! inner points
  b(1) = 2.0_wp
  do i = 2,n-1
    a(i) = dx(i-1) / (dx(i-1) + dx(i))
    b(i) = 2.0_wp
    c(i) = dx(i  ) / (dx(i-1) + dx(i))
  end do
  b(n) = 2.0_wp
  
  ! right boundary point
  select case( bc(2) )
  case( CS_BCS_CLAMPED, CS_BCS_FIXED_1 )
    a(n) = 1.0_wp
    c(n) = 0.0_wp
  case( CS_BCS_NATURAL, CS_BCS_FIXED_2 )
    a(n) = 0.0_wp
    c(n) = 0.0_wp
  case( CS_BCS_PERIODIC)
    a(n-1) = dx(n-2) / (dx(n-1) + dx(n-2))
    c(n-1) = dx(n-1) / (dx(n-1) + dx(n-2))
  end select
  
  return
end subroutine CUBIC_SPLINE_LHS

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_RHS(bc, bcval, n, dx, y, rhs)

  use TLab_Constants, only : wi, wp
  
  implicit none
  
  integer(wi), dimension(2  ), intent(in) :: bc    
  real(wp),    dimension(2  ), intent(in) :: bcval  
  integer(wi),                 intent(in) :: n
  real(wp),   dimension(n-1),  intent(in) :: dx
  real(wp),   dimension(n  ),  intent(in) :: y
  real(wp),   dimension(n  ),  intent(out):: rhs
  
  ! -------------------------------------------------------------------
  integer(wi)                             :: i
  
  ! ###################################################################
  ! left boundary point
  select case( bc(1) )
  case( CS_BCS_CLAMPED )
    rhs(1) = (6.0_wp / dx(  1)) *               (y(  2) - y(1)) / dx(  1)
  case( CS_BCS_FIXED_1 )
    rhs(1) = (6.0_wp / dx(  1)) * (-bcval(1) + ((y(  2) - y(1)) / dx(  1)))
  case( CS_BCS_NATURAL )
    rhs(1) = 0.0_wp
  case( CS_BCS_FIXED_2 )
    rhs(1) = 2.0_wp * bcval(1)
  case( CS_BCS_PERIODIC)
    rhs(1) = (y(2) - y(1)) / dx(1) - (y(1) - y(n-1)) / dx(n-1)
    rhs(1) = 6.0_wp * rhs(1) / (dx(1) + dx(n-1))
  end select

  ! inner points
  do i = 2,n-1
    rhs(i) = (y(i+1) - y(i)) / dx(i) - (y(i) - y(i-1)) / dx(i-1)
    rhs(i) = 6.0_wp * rhs(i) / (dx(i) + dx(i-1))
  end do
  
  ! right boundary point
  select case( bc(2) )
  case( CS_BCS_CLAMPED )
    rhs(n) = (6.0_wp / dx(n-1)) *               (y(n-1) - y(n)) / dx(n-1)
  case( CS_BCS_FIXED_1 )
    rhs(n) = (6.0_wp / dx(n-1)) * ( bcval(2) + ((y(n-1) - y(n)) / dx(n-1)))
  case( CS_BCS_NATURAL )
    rhs(n) = 0.0_wp
  case( CS_BCS_FIXED_2 )
    rhs(n) = 2.0_wp * bcval(2)
  case( CS_BCS_PERIODIC )
    rhs(n-1) = (y(n) - y(n-1)) / dx(n-1) - (y(n-1) - y(n-2)) / dx(n-2)
    rhs(n-1) = 6.0_wp * rhs(n-1) / (dx(n-1) + dx(n-2))
  end select

  return
end subroutine CUBIC_SPLINE_RHS

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_COEFF(n, dx, y, rhs, a, b, c, d)

  use TLab_Constants, only : wi, wp
  
  implicit none
  
  integer(wi),                intent(in) :: n
  real(wp),   dimension(n-1), intent(in) :: dx
  real(wp),   dimension(n  ), intent(in) :: y
  real(wp),   dimension(n  ), intent(in) :: rhs
  real(wp),   dimension(n-1), intent(out):: a,b,c,d
  
  ! -------------------------------------------------------------------
  integer(wi)                            :: i
  
  ! ###################################################################
  do i = 1,n-1
    a(i) = (rhs(i+1) - rhs(i)) / (6.0_wp * dx(i))
    b(i) = rhs(i) / 2.0_wp
    c(i) = (y(i+1) - y(i)) / dx(i) - (rhs(i+1) + 2.0_wp*rhs(i)) * (dx(i) / 6.0_wp)
    d(i) = y(i)
  end do

  return
end subroutine CUBIC_SPLINE_COEFF

!########################################################################
!########################################################################

subroutine CUBIC_SPLINE_FUNC(norg, nint, a, b, c, d, xorg, xint, yint)

  use TLab_Constants, only : wi, wp
  
  implicit none
  
  integer(wi),                   intent(in) :: norg,nint
  real(wp),   dimension(norg-1), intent(in) :: a,b,c,d
  real(wp),   dimension(norg  ), intent(in) :: xorg
  real(wp),   dimension(nint  ), intent(in) :: xint
  real(wp),   dimension(nint  ), intent(out):: yint
  
  ! -------------------------------------------------------------------
  integer(wi)                               :: i, idx
  real(wp)                                  :: z
  
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

  use TLab_Constants, only : wi, wp
  
  implicit none

  integer(wi),               intent(in) :: n   ! size of array a
  real(wp),    dimension(n), intent(in) :: a   ! ordered array
  real(wp),                  intent(in) :: x   ! element to insert
  integer(wi),               intent(out):: idx ! index for insertion
  
  ! -------------------------------------------------------------------
  integer(wi)                           :: mid, lo, hi
  
  ! ###################################################################

  lo = 1; hi = n

  do while ( lo < hi )
    mid = (lo + hi) / 2
    if (x < a(mid)) then
      hi = mid 
    else
      lo = mid + 1
    end if
  end do      

  idx = lo - 1

  return
end subroutine CUBIC_SPLINE_BISECT