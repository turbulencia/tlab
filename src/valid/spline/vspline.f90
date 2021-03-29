program vspline

  IMPLICIT NONE
  
#include "types.h"
#include "integers.h"

! ###################################################################
  TINTEGER, parameter          :: imax    = 10                      ! number of data points m > k
  TINTEGER, parameter          :: k       = 3                       ! degree of the spline k=[1,5]! cubic splines (k=3)
  TINTEGER, parameter          :: nest    = imax+2*k                ! number of knots of spline >=2*k+2, always large enough is nest=m+k+1,
  TINTEGER, parameter          :: lwrk    = imax*(k+1)+nest*(8+5*k) ! smallest dimension of working arrays >=imax*(k+1)+nest*(7+3*k) 
  TINTEGER, parameter          :: imax_sp = 28                      ! number of spline data points, should be (imax+x*(imax-1)) with x in N
  TINTEGER                     :: iopt, n, ier 
  TINTEGER                     :: i
  TINTEGER                     :: c1, c2, c3, c11, c22, c33         ! for time messurement
  TREAL                        :: start, finish                     ! for time messurement
  TREAL                        :: xb, xe, s, fp
! data arrays
  TREAL,    dimension(imax)    :: x, y, w 
  TREAL,    dimension(nest)    :: t, c
! data arrays for spline
  TREAL,    dimension(imax_sp) :: x_sp, y_sp
! working arrays
  TREAL,    dimension(lwrk)    :: wrk  
  TINTEGER, dimension(lwrk)    :: iwrk

  character(len=64)            :: line
! ###################################################################
! ###################################################################
! create data points and set spline parameters

! x-position of data x(i) and weights of points w(i)
  do i = 1,imax
    x(i) = i            ! x-values - must be strictly increasing
    w(i) = 1.0          ! weights of data points - weights are all equal here
  end do 
  call RANDOM_NUMBER(y) ! random y-values in the range of [0,1]

! write out random data points
  write(*,*) '================================================================================'
  do i = 1,imax
    write(*,*) 'random data point (x,y)   : ', x(i),y(i)
  end do
  write(*,*) '================================================================================'

! set interval for spline approximation
  xb = x(1);  xe = x(imax)

! create equally spaced array to evaluate spline with boundaries [xb,xe]
  do i=1, imax_sp
    x_sp(i) = xb + (xe - xb) * (i - 1) / (imax_sp - 1)
  end do

! parameters
  iopt = 0    ! (iopt=0 or 1) smoothing spline, weighted least-squares spline (iopt=-1)
  s    = 0    ! control the tradeoff between closeness of fit and smoothness
!s    : (in case iopt>=0) s must specify the smoothing factor.
!t    : array,length n, which contains the position of the knots.
!n    : integer, giving the total number of knots of s(x).
!c    : array,length n, which contains the b-spline coefficients.
!k    : integer, giving the degree of s(x).
!x    : array,length m, which contains the points where s(x) must
!fp   : contains the weighted sum of squared residuals of the spline approximation
!ier  : ier contains a non-positive value on exit [-2,-1,0], if error ier=[1,2,3,10]

! ###################################################################
! compute spline function with time function
  call system_clock(c1,c2,c3)
  call cpu_time(start)
  write(*,*) 'System clock start        : ', c1, c2, c3  
  write(*,*) 'CPU     time start        : ', start 
  write(*,*) '================================================================================'

! evaluation of spline function [curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)]
  call curfit(iopt,imax,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  if ( ier /= 0 .AND. ier /= -1 ) then
    write(line, *) 'Spline. Curfit error code = ', ier! call splev(t,n,c,k,x,y,m,ier)

    stop
  end if

! evaluation of the spline [call splev(t,n,c,k,x,y,m,ier)] 
! function to evaluate a B-spline or its derivatives
    !### input parameters:
    !t    : array,length n, which contains the position of the knots.
    !n    : integer, giving the total number of knots of s(x).
    !c    : array,length n, which contains the b-spline coefficients.
    !k    : integer, giving the degree of s(x).
    !x    : array,length m, which contains the points where s(x) must be evaluated.
    !m    : integer, giving the number of points where s(x) must be evaluated.
    !### output parameter:
    !y    : array,length m, giving the value of s(x) at the different points.
    !ier  : error flag: ier = 0 : normal return, ier =10 : invalid input data (see restrictions)
  call splev(t,n,c,k,x_sp,y_sp,imax_sp,ier) 
  if ( ier /= 0 ) then
    write(line, *) 'Spline. Splev error code = ', ier
    stop
  end if

  call system_clock(c11,c22,c33)
  call cpu_time(finish)
  write(*,*) 'System clock finish       : ', c1, c2, c3  
  write(*,*) 'CPU     time finfish      : ', start
  write(*,*) '================================================================================'
  write(*,*) 'time consumption sys_clock: ', c11-c1  
  write(*,*) 'time consumption cpu_clock: ', finish-start 
  write(*,*) '================================================================================'

! write out interpolated spline points 
  do i = 1,imax_sp
    write(*,*) 'spline point (x,y)        : ', x_sp(i),y_sp(i)
  end do
  write(*,*) '================================================================================'

! write out spline results to xxx.txt
  write(*,*) 'smoothing spline degree   : k   =', k
  write(*,*) 'smoothing factor          : s   =', s
  write(*,*) 'sum squared residuals     : fp  =', fp 
  write(*,*) 'total number of knots     : n   =', n 
  write(*,*) '================================================================================'
  do i=1,n-k-1
    write(*,*) 'b-spline coefficient      : ', c(i)
  end do
  write(*,*) '================================================================================'
  do i=1,nest
    write(*,*) 'b-spline knot position    : ', t(i)
  end do
  write(*,*) '================================================================================'

! ###################################################################
! write out data to files

! data points to xxx.txt
  open(11, file = 'data_points.txt')  
   do i=1,imax
      write(11,*) x(i), y(i)
   end do  
  close(11)

! data points interpolated to xxx.txt
  open(13, file = 'data_spline.txt')  
   do i=1,imax_sp
      write(13,*) x_sp(i), y_sp(i)
   end do  
  close(13)

! ###################################################################
! ###################################################################
  stop
end program vspline