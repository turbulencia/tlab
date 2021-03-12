!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2000/08/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# 3rd-order B-Spline interpolation for periodic sequence f of imax points,
!# x, between 0 and 1, gives the relative position of the desired point 
!# in the interval (i,i+1)
!#
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
#include "types.h"
  
FUNCTION BSPLINES3P(f, imax, i, x)
  
  IMPLICIT NONE
  
  TINTEGER imax, i
  TREAL f(imax), x
  TREAL BSPLINES3P
  
! -----------------------------------------------------------------------
  TREAL b1, b2, b3, b4
  TINTEGER i1, i3, i4, iw
  
! #######################################################################
  b1 = (C_1_R-x)*(C_1_R-x)*(C_1_R-x)
  b2 = C_4_R           - C_6_R*x*x + C_3_R*x*x*x
  b3 = C_1_R + C_3_R*x + C_3_R*x*x - C_3_R*x*x*x
  b4 = x*x*x
  
  iw = imax+1-i
  i1 = imax-MOD(iw,imax)
  i3 = MOD(i,imax)+1
  iw = i+1
  i4 = MOD(iw,imax)+1
  
  BSPLINES3P = ( b1*f(i1) + b2*f(i) + b3*f(i3) + b4*f(i4) ) /C_6_R
  
END FUNCTION BSPLINES3P
      
!########################################################################
!# Non-periodic BCs
!# BCs set by ghost cell with function value equal the first/last point
!########################################################################
FUNCTION BSPLINES3(f, imax, i, x)
  
  IMPLICIT NONE
  
  TINTEGER imax, i
  TREAL f(imax), x
  TREAL BSPLINES3
  
! -----------------------------------------------------------------------
  TREAL b1, b2, b3, b4
  TINTEGER i1, i3, i4
  
! #######################################################################
  b1 = (C_1_R-x)*(C_1_R-x)*(C_1_R-x)
  b2 = C_4_R           - C_6_R*x*x + C_3_R*x*x*x
  b3 = C_1_R + C_3_R*x + C_3_R*x*x - C_3_R*x*x*x
  b4 = x*x*x
  
  i1 = MAX(i-1,1)
  i3 = i + 1 
  i4 = MIN(i3+1,imax)

  BSPLINES3 = ( b1*f(i1) + b2*f(i) + b3*f(i3) + b4*f(i4) ) /C_6_R
  
END FUNCTION BSPLINES3
!########################################################################
!# Tool/Library
!#
!########################################################################
!# HISTORY
!#
!# 2000/02/01 - J.P. Mellado
!#              Created
!#
!########################################################################
!# DESCRIPTION
!#
!# 3rd-order B-Spline interpolation for periodic sequence f of imax points
!# in non-uniform case.
!# Cox-de Boor algorithm (Pozrikidis, p409)
!########################################################################
!# ARGUMENTS 
!#
!########################################################################
FUNCTION BSPLINES3_NU(f, t, imax, i, tint)

  IMPLICIT NONE

  TINTEGER imax, i
  TREAL, DIMENSION(imax) :: f, t
  TREAL tint
  TREAL BSPLINES3_NU

! -----------------------------------------------------------------------
  TREAL b1, b2, b3, b4, aux1, aux2, aux3
  TREAL t1,t2,t3,t4,t5,dt,dt1,dt2,dt3
  TINTEGER i1, i3, i4, iw

! #######################################################################
! Coefficients

  dt  =(t(imax)-t(1))/M_REAL(imax-1)
  dt1 = t(2)-t(1)
  dt2 = t(3)-t(1)
  dt3 = t(imax)-t(imax-1)

  IF      ( i .EQ. 1 ) THEN
     t1 = t(1) - dt - dt3
     t2 = t(1) - dt
     t3 = t(2)
     t4 = t(3)
     t5 = t(4)
  ELSE IF ( i .EQ. 2 ) THEN
     t1 = t(1) - dt 
     t2 = t(1)
     t3 = t(3)
     t4 = t(4)
     t5 = t(5)
  ELSE IF ( i .EQ. imax-2 ) THEN
     t1 = t(imax-4)
     t2 = t(imax-3)
     t3 = t(imax-1)
     t4 = t(imax)
     t5 = t(imax) + dt
  ELSE IF ( i .EQ. imax-1 ) THEN
     t1 = t(imax-3)
     t2 = t(imax-2)
     t3 = t(imax)
     t4 = t(imax) + dt
     t5 = t(imax) + dt + dt1
  ELSE IF ( i .EQ. imax ) THEN
     t1 = t(imax-2)
     t2 = t(imax-1)
     t3 = t(imax) + dt
     t4 = t(imax) + dt + dt1 
     t5 = t(imax) + dt + dt2
  ELSE 
     t1 = t(i-2)
     t2 = t(i-1)
     t3 = t(i+1)
     t4 = t(i+2)
     t5 = t(i+3)
  ENDIF

  aux1 = ( (tint-t2)*(t3  -tint) / (t3-t2) + (t4-tint)*(tint-t(i)) / (t4-t(i)) ) / &
         ((t4-t2)*(t3-t(i)))
  aux2 = ((t3-t1)*(t3-t2)*(t3-t(i)))
  aux3 = ((t5-t(i))*(t4-t(i))*(t3-t(i)))

  b1 =  (t3-tint )**C_3_R           /aux2
  b2 =  (t3-tint )**C_2_R*(tint-t1) /aux2 + aux1*(t4-tint)
  b3 = (tint-t(i))**C_2_R*(t5-tint) /aux3 + aux1*(tint-t2)
  b4 = (tint-t(i))**C_3_R           /aux3

! Interpolation

  iw = imax+1-i
  i1 = imax-MOD(iw,imax)
  i3 = MOD(i,imax)+1
  iw = i+1
  i4 = MOD(iw,imax)+1

  BSPLINES3_NU = b1*f(i1) + b2*f(i) + b3*f(i3) + b4*f(i4) 

END FUNCTION BSPLINES3_NU
