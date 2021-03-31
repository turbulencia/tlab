program vspline

  implicit none

#include "types.h"
#include "integers.h"

! ###################################################################
! define spline parameters here  
  TINTEGER, parameter          :: imax    = 10                      ! number of data points imax >= 10
  TINTEGER, parameter          :: k       = 3                       ! degree of the spline k=[1,5]
  TINTEGER, parameter          :: mesh    = 2                       ! mesh refinement factor (mesh=1 for x=x_sp)
! do not change the following
  TINTEGER, parameter          :: nest    = imax+2*k                ! number of knots of spline >=2*k+2, always large enough is nest=m+k+1,
  TINTEGER, parameter          :: lwrk    = imax*(k+1)+nest*(8+5*k) ! smallest dimension of working arrays >=imax*(k+1)+nest*(7+3*k) 
  TINTEGER, parameter          :: imax_sp = (imax+(mesh-1)*(imax-1))! number of spline data points
  TINTEGER                     :: iopt, n, ier
  TINTEGER                     :: i
  TREAL                        :: xb, xe, s, fp
! validation routine
  TREAL                        :: res
! time messurement
  TINTEGER                     :: c1, c2, c3, c11, c22, c33         
! data arrays
  TREAL,    dimension(imax)    :: x, y, w 
  TREAL,    dimension(nest)    :: t, c
! data arrays for spline
  TREAL,    dimension(imax_sp) :: x_sp, y_sp
! working arrays
  TREAL,    dimension(lwrk)    :: wrk  
  TINTEGER, dimension(lwrk)    :: iwrk
! ###################################################################
! create data points: 
! x-position of data x(i) and weights of points w(i)
  do i = 1,imax
    x(i) = i            ! x-values - must be strictly increasing
    w(i) = 1.0          ! weights of data points - weights are all equal here
  end do
  
! y-postion of data y(i) 
  call RANDOM_NUMBER(y) ! random y-values in the range of [0,1]

! set interval for spline approximation
  xb = x(1);  xe = x(imax)

! create equally spaced array to evaluate spline with boundaries [xb,xe]
  do i=1, imax_sp
    x_sp(i) = xb + (xe - xb) * (i - 1) / (imax_sp - 1)
  end do

  write(*,*) '================================================================================'
  write(*,*) '================================================================================'
  write(*,*) 'Validation routine of Spline library "fitpack"'
  write(*,*) '================================================================================'
  write(*,*) '================================================================================'

! ###################################################################
! compute spline function multiple times and messure time
  call system_clock(c1,c2,c3)
  i=1
  do while (i<10000)
    i = i+1
  ! evaluation of spline function [curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)]
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
    call curfit(iopt,imax,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
    if ( ier /= 0 .AND. ier /= -1 ) then
      write(*,*) 'Spline. Curfit error code = ', ier! call splev(t,n,c,k,x,y,m,ier)
      stop
    end if
  ! evaluation of the spline [call splev(t,n,c,k,x,y,m,ier)] function to evaluate a B-spline or its derivatives
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
      write(*,*) 'Spline. Splev error code = ', ier
      stop
    end if
  end do
! 
  call system_clock(c11,c22,c33)
   write(*,30) i, imax
  30 format(1X,'Loop over spline generation ', I5, '-times',' for each', I3,' points')
  write(*,'(1X,A,I4)') 'Time consumption (ms)     : ', c11-c1  
  write(*,*) '================================================================================'

! ###################################################################
! first validation routine: check if spline points are on original data points
! compute residual and check
res = 0.0
do i = 1, imax
  res = res + abs(y(i) - y_sp(1 + (i-1)*mesh))
end do

write(*,*) 'Validation: compute total sum of residua sum(delta_y) of spline and data points!'

if (res <= 1e-10) then
  write(*,*) 'Validation of spline routine SUCCESSFUL (res<1e-10)!'
  write(*,'(1X,A,ES9.2)') 'Total residua             :', res
  write(*,*) '================================================================================'
else
  write(*,*) 'Validation of spline routine FAILED (res>1e-10)!'
  write(*,'(1X,A,ES9.2)') 'Total residua             :', res
  write(*,*) '================================================================================' 
end if 

! write out spline results to xxx.txt
  write(*,'(1X,A,I4)') 'number of data points     : imax    = ', imax
  write(*,'(1X,A,I4)') 'number of spline points   : imax_sp = ', imax_sp
  write(*,'(1X,A,I4)') 'smoothing spline degree   : k       = ', k
  write(*,'(1X,A,F4.2)') 'smoothing factor          : s       = ', s
  write(*,'(1X,A,F4.2)') 'sum squared residuals     : fp      = ', fp 
  write(*,'(1X,A,I4)') 'total number of knots     : n       = ', n 
  write(*,*) '================================================================================'
! ###################################################################
! write out data to files
  open(11, file = 'data_points.txt')  
   do i=1,imax
      write(11,'(F18.14,F18.14)') x(i), y(i)
   end do  
  close(11)

  open(13, file = 'data_spline.txt')  
   do i=1,imax_sp
      write(13,'(F18.14,F18.14)') x_sp(i), y_sp(i)
   end do  
  close(13)
! ###################################################################
! ! write out random data points
!   write(*,*) '================================================================================'
!   do i = 1,imax
!     write(*,*) 'random data point (x,y)   : ', x(i),y(i)
!   end do
!   write(*,*) '================================================================================'

! ! write out interpolated spline points 
!   do i = 1,imax_sp
!     write(*,*) 'spline point (x,y)        : ', x_sp(i),y_sp(i)
!   end do
!   write(*,*) '================================================================================'

!   do i=1,n-k-1
!     write(*,*) 'b-spline coefficient      : ', c(i)
!   end do
!   write(*,*) '================================================================================'
!   do i=1,nest
!     write(*,*) 'b-spline knot position    : ', t(i)
!   end do
!   write(*,*) '================================================================================'
! ###################################################################
  stop
end program vspline