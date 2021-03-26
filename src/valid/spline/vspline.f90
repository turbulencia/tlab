program vspline

  IMPLICIT NONE
  
#include "types.h"
#include "integers.h"

! ###################################################################
  TINTEGER, parameter          :: imax    = 10                      ! number of data points
  TINTEGER, parameter          :: k       = 5                       ! degree of the spline k=[1,5]
  TINTEGER, parameter          :: nest    = imax+2*k                ! number of knots of spline >=2*k+2
  TINTEGER, parameter          :: lwrk    = imax*(k+1)+nest*(8+5*k) ! smallest dimension of working arrays >=imax*(k+1)+nest*(7+3*k) 
  TINTEGER, parameter          :: imax_sp = 39                      ! number of spline data points
  TINTEGER                     :: iopt, n, ier 
  TINTEGER                     :: i
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
    x(i) = i            ! x-values
    w(i) = 1.0          ! weights of data points
  end do 
  call RANDOM_NUMBER(y) ! random y-values in the range of [0,1]

! write out random data points
  write(*,*) '================================================================================'
  do i = 1,imax
    write(*,*) 'random data point (x,y)   : ', x(i),y(i)
  end do
  write(*,*) '================================================================================'

! set boundaries
  xb = x(1);  xe = x(imax)

! parameters
  iopt = 0    ! smoothing spline
  s    = 0    ! control the tradeoff between closeness of fit and smoothness
!  n   = no need to specity for iopt=0
!  t   = contains the knots of the spline
!  c   = contains the coefficients in the b-spline
!  fp  = contains the weighted sum of squared residuals of the spline approximation
!  ier = ier contains a non-positive value on exit [-2,-1,0], if error ier=[1,2,3,10]

! ###################################################################
! compute spline function

! evaluation of spline function
  call curfit(iopt,imax,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
  if ( ier /= 0 .AND. ier /= -1 ) then
    write(line, *) 'Spline. Curfit error code = ', ier
    stop
  end if

! create equally spaced array to evaluate spline with boundaries [xb,xe]
  do i=1, imax_sp
    x_sp(i) = xb + (xe - xb) * (i - 1) / (imax_sp - 1)
  end do

! evaluation of the spline
  call splev(t,n,c,k,x_sp,y_sp,imax_sp,ier)
  if ( ier /= 0 ) then
    write(line, *) 'Spline. Splev error code = ', ier
    stop
  end if

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