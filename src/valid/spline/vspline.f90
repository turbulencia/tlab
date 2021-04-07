program vspline

  implicit none

#include "types.h"
#include "integers.h"

! ###################################################################
! define spline parameters here  
  TINTEGER, parameter          :: imax    = 10                      ! number of data points imax >= 10
  TINTEGER, parameter          :: k       = 3                       ! degree of the spline k=[1,5]
  TINTEGER, parameter          :: mesh    = 10                      ! mesh refinement factor (mesh=1 for x=x_new)
! do not change the following
  TINTEGER, parameter          :: nest    = imax+2*k                ! number of knots of spline >=2*k+2, always large enough is nest=m+k+1,
  TINTEGER, parameter          :: lwrk    = imax*(k+1)+nest*(8+5*k) ! smallest dimension of working arrays >=imax*(k+1)+nest*(7+3*k) 
  TINTEGER, parameter          :: imax_new= (imax+(mesh-1)*(imax-1))! number of spline data points
  TINTEGER                     :: iopt, n, ier
  TINTEGER                     :: i
  TREAL                        :: xb, xe, s, fp
! validation routine
  TREAL                        :: res_2, res_inf
! time messurement
  TINTEGER                     :: c1, c2, c3, c11, c22, c33         
! data arrays 
! (better one working array for w,t,c,wrk,iwrk, see opr_interpolate_pool.f90)
  TREAL,    dimension(imax)    :: x, y, w
  TREAL,    dimension(nest)    :: t, c
! data arrays for spline
  TREAL,    dimension(imax_new):: x_new, y_sp, y_new, delta 
! working arrays
  TREAL,    dimension(lwrk)    :: wrk  
  TINTEGER, dimension(lwrk)    :: iwrk
! ###################################################################
! create data points with sin(x) as exact function 
  do i = 1,imax
    x(i) = (i-1) * (2*C_PI_R / (imax-1)) ! x-values - must be strictly increasing
    w(i) = 1.0                           ! weights of data points, here: all equal 
    y(i) = sin(x(i))                     ! y-values
  end do

! set interval for spline approximation
  xb = x(1);  xe = x(imax)

! create equally spaced array to evaluate spline with boundaries [xb,xe]
  do i = 1,imax_new
    x_new(i) = xb + (xe - xb) * (i - 1) / (imax_new - 1)
    y_new(i) = sin(x_new(i))
  end do

! ! in case of radom data points uncomment the following (no exact solution, validation not possible!!!)
!   call RANDOM_NUMBER(y)
!   y_new = C_0_R
!   do i = 1,imax
!     x(i) = i
!     y_new(1+(i-1)*mesh) = y(i)            
!   end do
!   xb = x(1);  xe = x(imax)
!   do i = 1,imax_new
!     x_new(i) = xb + (xe - xb) * (i - 1) / (imax_new - 1)
!   end do

  write(*,*) '================================================================================'
  write(*,*) '================================================================================'
  write(*,*) 'Validation routines of spline library "fitpack"'
  write(*,*) 'Here: Interpolating spline (data points == spline points) of order k = [1,..,5]'
  write(*,*) '================================================================================'
  write(*,*) '================================================================================'
! ###################################################################
! spline function parameter
  iopt = 0    ! (iopt=0 or 1) smoothing spline, weighted least-squares spline (iopt=-1)
  s    = 0.0  ! control the tradeoff between closeness of fit and smoothness

! compute spline function multiple times and messure time
  call system_clock(c1,c2,c3)
  i=1
  do while (i<=1000)
    i = i+1
  ! evaluation of spline function [curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)]
      !s    : (in case iopt>=0) s must specify the smoothing factor. 
      !       if s=0 then interpolation spline, data points coincident with spline points
      !t    : array,length n, which contains the position of the knots.
      !n    : integer, giving the total number of knots of s(x).
      !c    : array,length n, which contains the b-spline coefficients.
      !k    : integer, giving the degree of s(x).
      !x    : array,length m, which contains the points where s(x) must
      !fp   : contains the weighted sum of squared residuals of the spline approximation
      !ier  : ier contains a non-positive value on exit [-2,-1,0], if error ier=[1,2,3,10]
    call curfit(iopt,imax,x,y,w,xb,xe,k,s,nest,n,t,c,fp,wrk,lwrk,iwrk,ier)
    if ( ier /= 0 .AND. ier /= -1 ) then
      write(*,*) 'Spline. Curfit error code = ', ier
      stop
    end if
  ! evaluation of the spline [call splev(t,n,c,k,x,y,m,ier)] function to evaluate a B-spline or its derivatives
      !###### input parameters:
      !t    : array,length n, which contains the position of the knots.
      !n    : integer, giving the total number of knots of s(x).
      !c    : array,length n, which contains the b-spline coefficients.
      !k    : integer, giving the degree of s(x).
      !x    : array,length m, which contains the points where s(x) must be evaluated.
      !m    : integer, giving the number of points where s(x) must be evaluated.
      !###### output parameter:
      !y    : array,length m, giving the value of s(x) at the different points.
      !ier  : error flag: ier = 0 : normal return, ier =10 : invalid input data (see restrictions)
    call splev(t,n,c,k,x_new,y_sp,imax_new,ier) 
    if ( ier /= 0 ) then
      write(*,*) 'Spline. Splev error code = ', ier
      stop
    end if
  end do

  call system_clock(c11,c22,c33)
  write(*,*) '0. Validation routine: Messure time consumption'
  write(*,30) i-1, imax
  30 format(1X,'Loop over spline generation ', I5, '-times',' for ', I3,' points each')
  write(*,'(1X,A,I4)') 'Time consumption (ms)     : ', c11-c1  
  write(*,*) '================================================================================'
! ###################################################################
! 1. Validation
  write(*,*) '1. Validation routine: Check if data points coincident with spline points'
  do i = 1,imax
    if ((abs(y(i) - y_sp(1 + (i-1)*mesh))) <= 1e-10) then
      write(*,40) i, x(i), y(i), (1 + (i-1)*mesh), abs(y(i) - y_sp(1 + (i-1)*mesh))
      40 format(1X,'Data point ', I3, ' with ', '(', F4.2, ',', F4.2,')', ' is on spline point ', I3,', residua = ' ,ES9.2)
    else
      write(*,50) i, x(i), y(i), (1 + (i-1)*mesh), x_new(1 + (i-1)*mesh), y_sp(1 + (i-1)*mesh), abs(y(i) - y_sp(1 + (i-1)*mesh))
      50 format(1X,'Data point ', I3, ' with ', '(', F4.2, ',', F4.2,')', ' is not on spline point ', I3, ' with ', '(', F4.2, ',', F4.2,')',', residua = ' ,ES9.2)
    end if 
  end do

! 2. Validation
  write(*,*) '================================================================================' 
  write(*,*) '2. Validation routine: Error norms of exact solution and spline points'

  do i = 1 , imax_new
    delta(i) = y_new(i) - y_sp(i)
  end do  
  res_2   = norm2(delta)
  res_inf = maxval(abs(delta))

  write(*,'(1X,A,ES9.2)') 'L_2   error norm          :', res_2
  write(*,'(1X,A,ES9.2)') 'L_inf error norm          :', res_inf
  write(*,*) '================================================================================'

! write out spline results to xxx.txt
  write(*,*) '3. Validation routine: Spline interpolation results and parameters'
  write(*,'(1X,A,I4)') 'number of data points     : imax     = ', imax
  write(*,'(1X,A,I4)') 'mesh refinement factor    : mesh     = ', mesh 
  write(*,'(1X,A,I4)') 'number of spline points   : imax_new = ', imax_new
  write(*,'(1X,A,I4)') 'smoothing spline degree   : k        = ', k
  write(*,'(1X,A,F4.2)') 'smoothing factor          : s        = ', s
  write(*,'(1X,A,F4.2)') 'sum squared residuals     : fp       = ', fp 
  write(*,'(1X,A,I4)') 'total number of knots     : n        = ', n 
  write(*,*) '================================================================================'
! ###################################################################
! write out data to files
  open(11, file = 'data_points.txt')  
   do i=1,imax
      write(11,'(F18.14,F18.14)') x(i), y(i)
   end do  
  close(11)

  open(12, file = 'data_exact.txt')
  open(13, file = 'data_spline.txt')  
   do i=1,imax_new
      write(12,'(F18.14,F18.14)') x_new(i), y_new(i)
      write(13,'(F18.14,F18.14)') x_new(i), y_sp(i)
   end do  
  close(12)
  close(13)
! ###################################################################
! ! write out random data points
!   write(*,*) '================================================================================'
!   do i = 1,imax
!     write(*,*) 'random data point (x,y)   : ', x(i),y(i)
!   end do
!   write(*,*) '================================================================================'

! ! write out interpolated spline points 
!   do i = 1,imax_new
!     write(*,*) 'spline point (x,y)        : ', x_new(i),y_sp(i)
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