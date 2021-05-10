program vspline

  implicit none

#include "types.h"
#include "integers.h"

! ###################################################################
! define spline parameters here  
  TINTEGER, parameter          :: imax    = 10                      ! number of data points imax >= 10
  TINTEGER, parameter          :: k       = 3                       ! degree of the spline k=[1,5]
  TINTEGER, parameter          :: mesh    = 15                      ! mesh refinement factor (mesh=1 for x=x_new)
! do not change the following
  TINTEGER, parameter          :: nest    = imax+2*k                ! number of knots of spline >=2*k+2, always large enough is nest=m+k+1,
  TINTEGER, parameter          :: lwrk    = imax*(k+1)+nest*(8+5*k) ! smallest dimension of working arrays >=imax*(k+1)+nest*(7+3*k) 
  TINTEGER, parameter          :: imax_new= (imax+(mesh-1)*(imax-1))! number of spline data points
  TINTEGER                     :: iopt, n, ier
  TINTEGER                     :: i, j
  TREAL                        :: xb, xe, s, fp
! validation routine
  TREAL                        :: res_2, res_inf
! time messurement
  TINTEGER                     :: iter, c1, c2, c3, c11, c22, c33   
  TINTEGER, parameter          :: start=0, end=1000, step=250
  TREAL,    dimension((end-start)/step+1) :: delta_t
! data arrays                                                       ! better one working array for w,t,c,wrk,iwrk, see opr_interpolate_pool.f90
  TREAL,    dimension(imax)    :: x, y, w
  TREAL,    dimension(nest)    :: t, c
! data arrays for spline
  TREAL,    dimension(imax_new):: x_new, y_sp, y_new, delta 
! working arrays
  TREAL,    dimension(lwrk)    :: wrk  
  TINTEGER, dimension(lwrk)    :: iwrk
! ###################################################################
! Testrun - scaling -- define spline parameters here  
  TINTEGER, parameter          :: imax2_new = 2000, jmax2 = 2000
  TINTEGER, parameter          :: imax2     = imax2_new - imax2_new / 20 * 6
  TINTEGER, parameter          :: k2        = 3                           
! do not change the following
  TINTEGER, parameter          :: nest2    = imax2_new+2*k2                 
  TINTEGER, parameter          :: lwrk2    = imax2_new*(k2+1)+nest2*(8+5*k2) 
  TREAL                        :: xb2, xe2
! data arrays
  TREAL, dimension(imax2)      :: x2, w2
  TREAL, dimension(imax2,jmax2):: y2c                               ! for column access
  TREAL, dimension(jmax2,imax2):: y2r                               ! for row    access
  TREAL, dimension(nest2)      :: t2, c_2
! data arrays for spline
  TREAL, dimension(imax2_new)       :: x2_new, epsi = 1             ! epsi: 0 = solid, 1 = fluid 
  TREAL, dimension(imax2_new,jmax2) :: y2c_new    
  TREAL, dimension(jmax2,imax2_new) :: y2r_new    
! working arrays
  TREAL,    dimension(lwrk2)        :: wrk2  
  TINTEGER, dimension(lwrk2)        :: iwrk2
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

! scaling of spline functions
  write(*,*) '0. Validation routine: Messure time consumption'
  j = 1
  do iter = start, end, step
    call system_clock(c1,c2,c3)
    i=1
    do while (i<=iter)
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
    ! time consumption + save
    call system_clock(c11,c22,c33)
    write(*,30) i-1, c11 - c1 
    30 format(1X, I5,' iterations with delta_t (ms) : ',I4) 
    delta_t(j) = c11 - c1
    j = j+1
  end do
  write(*,*) '================================================================================'
 
  ! write out scaling 
  open(10, file = 'data_scaling.txt')  
   do i=1,(end-start)/step+1
      write(10,'(I10,I10)') start+(i-1)*step, int(delta_t(i))
   end do  
  close(10)
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

! write out spline results
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
! ## DEBUGGING
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
! ## acces larger data - messure scaling for column/row access ######
! ###################################################################
! create mask / epsilon field (Interface points on grid, but here with value 1 (7,14, 27,34 ...))
  do i = 10,imax2_new,20
    do j = 0,5              ! gap size, 6 solid points, 12 fluid points, 2 interface points
      epsi(i+j-2) = 0       
    end do
  end do

! create grid
  x2_new = (/ (i, i = 1, imax2_new, 1) /)
  x2     = pack(x2_new , epsi>0)
  w2     = (/ (1, i = 1, imax2, 1) /) ! weights of data points, all equal here

! set interval for spline approximation
  xb2 = x2(1);  xe2 = x2(imax2)

! create random data field, insert zeros for interface (for row and column access)
  call RANDOM_NUMBER(y2c)  
  do i = 1,imax2-1
    if(( x2(i+1) - x2(i) ) > 1 ) then
      do j = 1,jmax2
        y2c(i,j)   = 0
        y2c(i+1,j) = 0
      end do
    end if
  end do
  y2r = transpose(y2c)

! evaluate global spline function 
  iopt = 0    
  s    = 0.0  

! ! for 1 - line 
  write(*,*) '================================================================================'
  write(*,*) '4. Time Messurement - global Spline generation - large field'

! column access (should be fastest)
  call curfit(iopt,imax2,x2,y2c(:,1),w2,xb2,xe2,k2,s,nest2,n,t2,c_2,fp,wrk2,lwrk2,iwrk2,ier)
  call splev(t2,n,c_2,k2,x2_new,y2c_new(:,1),imax2_new,ier)
  ! start time messurement here
  call system_clock(c1,c2,c3)
  do j = 2,jmax2              ! skip first one
    call curfit(iopt,imax2,x2,y2c(:,j),w2,xb2,xe2,k2,s,nest2,n,t2,c_2,fp,wrk2,lwrk2,iwrk2,ier)
    if ( ier /= 0 .AND. ier /= -1 ) then
      write(*,*) 'Spline. Curfit error code = ', ier
      stop
    end if
    call splev(t2,n,c_2,k2,x2_new,y2c_new(:,j),imax2_new,ier) 
    if ( ier /= 0 ) then
      write(*,*) 'Spline. Splev error code = ', ier
      stop
    end if
  end do
  call system_clock(c11,c22,c33)
  write(*,70) c11 - c1 
  70 format(1X,'Delta_t for field - column access (ms) : ', I4) 

! row access
  call curfit(iopt,imax2,x2,y2r(1,:),w2,xb2,xe2,k2,s,nest2,n,t2,c_2,fp,wrk2,lwrk2,iwrk2,ier)
  call splev(t2,n,c_2,k2,x2_new,y2r_new(1,:),imax2_new,ier) 
  ! start time messurement here
  call system_clock(c1,c2,c3)
  do j = 2,jmax2
    call curfit(iopt,imax2,x2,y2r(j,:),w2,xb2,xe2,k2,s,nest2,n,t2,c_2,fp,wrk2,lwrk2,iwrk2,ier)
    if ( ier /= 0 .AND. ier /= -1 ) then
      write(*,*) 'Spline. Curfit error code = ', ier
      stop
    end if
    call splev(t2,n,c_2,k2,x2_new,y2r_new(j,:),imax2_new,ier) 
    if ( ier /= 0 ) then
      write(*,*) 'Spline. Splev error code = ', ier
      stop
    end if
  end do
  call system_clock(c11,c22,c33)
  write(*,80) c11 - c1 
  80 format(1X,'Delta_t for field - row    access (ms) : ', I4) 

! write out spline results
  write(*,*) '================================================================================'
  write(*,*) '5. Validation routine: Spline interpolation results and parameters'
  write(*,'(1X,A,I4)') 'number of data points     : imax     = ', imax2
  write(*,'(1X,A,I4)') 'number of lines           : jmax     = ', jmax2
  write(*,'(1X,A,I4)') 'number of spline points   : imax_new = ', imax2_new
  write(*,'(1X,A,I4)') 'smoothing spline degree   : k        = ', k2
  write(*,'(1X,A,F4.2)') 'smoothing factor          : s        = ', s
  write(*,'(1X,A,F4.2)') 'sum squared residuals     : fp       = ', fp 
  write(*,'(1X,A,I4)') 'total number of knots     : n        = ', n 
  write(*,*) '================================================================================'
! ###################################################################
! write out data to files, one line of field y2r, y2c
  open(14, file = 'data_points_large.txt')  
    do i=1,imax2
      write(14,'(F20.14,F20.14)') x2(i), y2c(i,1)
    end do  
  close(14)

  open(15, file = 'data_spline_large_row.txt')
  open(16, file = 'data_spline_large_column.txt')  
    do i=1,imax2_new
      write(15,'(F20.14,F20.14)') x2_new(i), y2c_new(i,1)
      write(16,'(F20.14,F20.14)') x2_new(i), y2r_new(1,i)
    end do  
  close(15)
  close(16)  
! ###################################################################
! ## DEBUGGING
! ###################################################################
  ! write(*,*) '================================================================================'
  ! do i = 1,imax2_new
  !   write(*,*) 'x2new  data point (x)   : ', i, x2_new(i)
  ! end do
  ! write(*,*) '================================================================================'
  ! do i = 1,imax2_new
  !   write(*,*) 'epsilon field     (x)   : ', i, epsi(i)
  ! end do
  ! write(*,*) '================================================================================'
  ! do i = 1,imax2
  !   write(*,*) 'x2gap data point  (x)   : ', i, x2(i)
  ! end do
  ! write(*,*) '================================================================================'
  ! do i = 1,imax2
  !   write(*,*) 'w2                (x)   : ', i, w2(i)
  ! end do
  ! write(*,*) '================================================================================'
  ! do i = 1,imax2
  !   write(*,*) 'y2                (x)   : ', i, y2c(i,1) , y2r(1,i)
  ! end do
  ! write(*,*) '================================================================================'
  ! do i = 1,imax2
  !   write(*,*) 'y2_new            (x)   : ', i, y2c_new(i,1) , y2r_new(1,i)
  ! end do
  ! write(*,*) '================================================================================'
  ! ! write(*,*) 'y2                (x)   : ', y2(:,1)
  ! write(*,*) '================================================================================'
  stop
end program vspline