#include "types.h"

!########################################################################
!# The stretching factor is given by a hyperbolic tangent function, which means (see manual)
!#
!# dx/ds = [ 1 + EXP[(s-s_1)/delta_1] ] **(delta_1*f1/h_0)
!#        *[ 1 + EXP[(s-s_2)/delta_2] ] **(delta_2*f2/h_0)
!#        *...
!#
!# Using compact schemes to integrate this equation
!########################################################################
SUBROUTINE BLD_EXP(idir, x, imax, scalex, w)
  USE GRID_LOCAL
  IMPLICIT NONE

  TINTEGER, INTENT(IN   ) :: idir, imax
  TREAL,    INTENT(INOUT) :: x(imax), scalex, w(imax,7)

  ! -----------------------------------------------------------------------
  TREAL st(3), df(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
  TINTEGER i, iloc, iseg, im, ibc, i1
  TREAL ds

  ! #######################################################################
  iseg = 1                      ! Assumes just 1 segment

  iloc = 1
  IF ( g_build(idir)%mirrored ) iloc = imax/2 ! mirrored case; first point in array is imax/2

  ! create uniform reference grid s
  ds = scalex /M_REAL(imax-iloc)
  DO i=iloc,imax
    x(i) = M_REAL(i-iloc) *ds
  ENDDO

  ! -----------------------------------------------------------------------
  ! define local variables for readability below; strides of 3
  st    = g_build(idir)%vals(1::3,iseg)     ! transition point in uniform grid
  df    = g_build(idir)%vals(2::3,iseg) /ds ! normalized maximum stretching factor
  delta = g_build(idir)%vals(3::3,iseg)

  ! create mapping from grid s to grid x
#define jac(j)  w(j,6)
#define rhs(j)  w(j,7)
  rhs(iloc:imax) = C_1_R
  DO im = 1,3                           ! 3 modes at most
    IF ( ABS(delta(im)) > C_0_R ) THEN
      rhs(iloc:imax) = rhs(iloc:imax) *( EXP((x(iloc:imax)-st(im))/delta(im)) +C_1_R )**(df(im)*delta(im))
    END IF
  END DO
  jac(iloc:imax) = ds
  ibc = 1           ! Boundary condition at the first node
  i1 = 1            ! One equation to be solved
  CALL INT_C1N6_LHS(imax-iloc+1,     ibc,             w(iloc,1),w(iloc,2),w(iloc,3),w(iloc,4),w(iloc,5))
  CALL INT_C1N6_RHS(imax-iloc+1, i1, ibc,  jac(iloc), rhs(iloc), x(iloc))
  CALL PENTADFS(imax-iloc,     w(iloc+1,1),w(iloc+1,2),w(iloc+1,3),w(iloc+1,4),w(iloc+1,5))
  CALL PENTADSS(imax-iloc, i1, w(iloc+1,1),w(iloc+1,2),w(iloc+1,3),w(iloc+1,4),w(iloc+1,5), x(iloc+1))
  x(iloc) = C_0_R   ! Boundary condition
#undef jac
#undef rhs

  ! correct value of scale
  scalex = x(imax)

  RETURN
END SUBROUTINE BLD_EXP
