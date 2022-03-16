#include "types.h"

!########################################################################
!# Valid
!#
!########################################################################
!# HISTORY
!#
!# 2022/03/16 - J. Kostelecky
!#              Created           
!#
!########################################################################
!# DESCRIPTION
!#
!# Validate cubic splines.
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################

program CSPLINE

  implicit none
 
#include "integers.h"
  
  ! define spline parameters here  
  TINTEGER, parameter          :: imax     = 100                       ! number of data points imax >= 10
  TINTEGER, parameter          :: mesh     = 10                       ! mesh refinement factor (mesh=1 for x_new=x)
  TINTEGER, parameter          :: imax_new = (imax+(mesh-1)*(imax-1)) ! number of spline data points
  TINTEGER                     :: i, j, n
  TREAL                        :: xb, xe
  ! validation routine
  TREAL                        :: res_2, res_inf
  ! time messurement
  TINTEGER                     :: iter, c1, c2, c3, c11, c22, c33   
  TINTEGER, parameter          :: start=0, end=50000, step=5000
  TREAL,    dimension((end-start)/step+1) :: delta_t
  ! data arrays                                                 
  TREAL,    dimension(imax)    :: x, y
  ! data arrays for spline
  TREAL,    dimension(imax_new):: x_new, y_sp, y_new, delta
  ! working arrays
  TREAL,    dimension(9*imax)    :: wrk  
! ###################################################################
! Initialize
  do i = 1,imax
    x(i) = (i-C_1_R) * (C_2_R*C_PI_R / (imax-C_1_R)) ! x-values - must be strictly increasing
    y(i) = sin(x(i))                                 ! y-values
  end do

! set interval for spline approximation
  xb = x(1);  xe = x(imax)

! create equally spaced array to evaluate spline with boundaries [xb,xe]
  do i = 1,imax_new
    x_new(i) = xb + (xe - xb) * (i - C_1_R) / (imax_new - C_1_R)
    y_new(i) = sin(x_new(i))
  end do

! cubic spline function
  call CUBIC_SPLINE(imax, imax_new, x, y, x_new, y_sp, wrk)

  write(*,*) '================================================================================'
  write(*,*) '================  Validation routines for cubic splines ========================'
  write(*,*) '================================================================================'
! ###################################################################
! 0. Validation - scaling of spline functions
  write(*,*) '0. Validation routine: Messure time consumption'
  j = 1
  do iter = start, end, step
    call system_clock(c1,c2,c3)
    i = 1
    do while (i<=iter)
      i = i + 1
        call CUBIC_SPLINE(imax, imax_new, x, y, x_new, y_sp, wrk)
    end do
    ! time consumption + save
    call system_clock(c11,c22,c33)
    write(*,30) i-1, c11 - c1 
    30 format(1X, I5,' iterations with delta_t (ms) : ',I4) 
    delta_t(j) = c11 - c1
    j = j + 1
  end do

! 1. Validation
  write(*,*) '================================================================================' 
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

! IO - error and function values
  open(20, file='cspline.dat')
  do i = 1,imax_new
    write(20,1000) x_new(i), y_new(i), y_sp(i), y_new(i) - y_sp(i) 
  end do
  close(20)
  1000 format(6(1x,e16.10))

  stop

end program CSPLINE