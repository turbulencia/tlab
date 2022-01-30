#include "types.h"

MODULE GRID_LOCAL
  IMPLICIT NONE
  SAVE

  TINTEGER, PARAMETER :: GTYPE_UNIFORM = 0
  TINTEGER, PARAMETER :: GTYPE_TANH    = 5
  TINTEGER, PARAMETER :: GTYPE_EXP     = 6

  TINTEGER, PARAMETER :: MAX_PARAMES = 9
  TINTEGER, PARAMETER :: MAX_OPTIONS = 5
  TINTEGER, PARAMETER :: MAX_SEGMENT =10

  TYPE grid_build_dt                          ! Information to construct grid in one direction
    SEQUENCE
    TINTEGER nseg                             ! number of segments in this direction
    LOGICAL mirrored                          ! It true, mirror the grid
    TREAL fixed_scale                         ! If positive, rescale grid to this value
    TINTEGER SIZE(MAX_SEGMENT)                ! number of points in each segment
    TREAL    END(MAX_SEGMENT)                 ! physical end of each segment
    TINTEGER opts(MAX_OPTIONS,MAX_SEGMENT)    ! 0   uniform segment
                                              ! 1   Colonius, Lele and Moin stretching
                                              ! 2   Second order polynomial stretching
                                              ! 3   Third order polynomial stretching
                                              ! 4   Geometric progression
                                              ! 5   Hyperbolic tangent
                                              ! 6   Exponential
    TREAL    vals(MAX_PARAMES,MAX_SEGMENT)
  END TYPE grid_build_dt

  TYPE(grid_build_dt) g_build(3)

CONTAINS
  !########################################################################
  !# The grid spacing dx/ds is given by a hyperbolic tangent function
  !#
  !# dx/ds = f_0 + [(f_1-f_0)/2]*[ 1 + TANH[(s-s_1)/(2*delta_1)] ]
  !#             + [(f_2-f_0)/2]*[ 1 + TANH[(s-s_2)/(2*delta_2)] ]
  !#             + ...
  !########################################################################
  SUBROUTINE BLD_TANH(idir, iseg, x, nmax, work)
    IMPLICIT NONE

    TINTEGER, INTENT(IN   ) :: idir, iseg, nmax
    TREAL,    INTENT(INOUT) :: x(nmax), work(nmax)

    ! -----------------------------------------------------------------------
    TREAL st(3), f(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
    TINTEGER im

    ! #######################################################################
    ! define local variables for readability below; strides of 3
    st    = g_build(idir)%vals(1::3,iseg)     ! transition point in uniform grid
    f     = g_build(idir)%vals(2::3,iseg)     ! ratio dx(i)/dx_0
    delta = g_build(idir)%vals(3::3,iseg)

    ! create mapping from grid s to grid x
    work = C_0_R
    DO im = 1,3                           ! 3 modes at most
      IF ( ABS(delta(im)) > C_0_R ) THEN
        work(:) = work(:) + (f(im)-C_1_R)*delta(im)*LOG(EXP((x(:)-st(im))/delta(im))+C_1_R)
      END IF
    END DO
    work = work - work(1) ! Integration constant

    x = x +work

    RETURN
  END SUBROUTINE BLD_TANH

  !########################################################################
  !# The stretching factor is given by a hyperbolic tangent function, which means (see manual)
  !#
  !# dx/ds = [ 1 + EXP[(s-s_1)/delta_1] ] **(delta_1*f1/h_0)
  !#        *[ 1 + EXP[(s-s_2)/delta_2] ] **(delta_2*f2/h_0)
  !#        *...
  !#
  !# Using compact schemes to integrate this equation
  !########################################################################
  SUBROUTINE BLD_EXP(idir, iseg, x, nmax, w)
    IMPLICIT NONE

    TINTEGER, INTENT(IN   ) :: idir, iseg, nmax
    TREAL,    INTENT(INOUT) :: x(nmax), w(nmax,8)

    ! -----------------------------------------------------------------------
    TREAL st(3), df(3), delta(3)    ! superposition of up to 3 modes, each with 3 parameters
    TINTEGER im, ibc, i1
    TREAL ds

    ! #######################################################################
    ds = ( x(nmax) -x(1) )/ M_REAL(nmax-1)

    ! define local variables for readability below; strides of 3
    st    = g_build(idir)%vals(1::3,iseg)     ! transition point in uniform grid
    df    = g_build(idir)%vals(2::3,iseg) /ds ! normalized maximum stretching factor
    delta = g_build(idir)%vals(3::3,iseg)

    ! create mapping from grid s to grid x
#define jac(j)    w(j,6)
#define rhs(j)    w(j,7)
#define result(j) w(j,8)
    rhs(:) = C_1_R
    DO im = 1,3               ! 3 modes at most
      IF ( ABS(delta(im)) > C_0_R ) THEN
        rhs(:) = rhs(:) *( EXP((x(:)-st(im))/delta(im)) +C_1_R )**(df(im)*delta(im))
      END IF
    END DO
    jac(:) = ds
    ibc = 1                           ! Boundary condition at the first node
    i1 = 1                            ! One equation to be solved
    CALL INT_C1N6_LHS(nmax,     ibc,             w(1,1),w(1,2),w(1,3),w(1,4),w(1,5))
    CALL INT_C1N6_RHS(nmax, i1, ibc,  jac(1), rhs(1), result(1))
    CALL PENTADFS(nmax,     w(2,1),w(2,2),w(2,3),w(2,4),w(2,5))
    CALL PENTADSS(nmax, i1, w(2,1),w(2,2),w(2,3),w(2,4),w(2,5), result(2))
    x(2:nmax) = x(1) +result(2:nmax)  ! Boundary condition
#undef jac
#undef rhs
#undef result

    RETURN
  END SUBROUTINE BLD_EXP

  !########################################################################
  ! Geometric progression and other options
  !########################################################################
  SUBROUTINE BLD_THEREST(idir, iseg, x, nmax)
    IMPLICIT NONE

    TINTEGER, INTENT(IN   ) :: idir, iseg, nmax
    TREAL,    INTENT(INOUT) :: x(nmax)

    ! -----------------------------------------------------------------------
    TREAL a, b, c, d, e
    TREAL eta, deta, dx
    TINTEGER n

    ! #######################################################################
    CALL BLD_CONSTANTS(nmax, x(1),x(nmax),&
        g_build(idir)%opts(1,iseg),g_build(idir)%opts(2,iseg),&
        g_build(idir)%vals(1,iseg),g_build(idir)%vals(2,iseg),&
        g_build(idir)%vals(3,iseg),g_build(idir)%vals(4,iseg), a,b,c,d,e)

    ! Go from computational eta in [0,1] to physical domains
    dx = C_1_R
    IF ( nmax == 1 ) THEN
      deta = C_1_R
    ELSE;
      deta = C_1_R /M_REAL(nmax-1)
    ENDIF

    DO n = 2,nmax
      eta  = M_REAL(n-1)*deta

      SELECT CASE(  g_build(idir)%opts(1,iseg) )
      CASE( 1 )     ! Colonius, Lele and Moin
        x(n) = e + a*eta + d*LOG( EXP(c*(a*eta - b)) + C_1_R - EXP(-b*c) )
      CASE( 2,3 )   ! 2nd- 3rd-Order Polynomial Stretching
        x(n) = a + b*eta + c*eta*eta + d*eta*eta*eta
      CASE( 4 )     ! Geometric prograssion
        dx      = dx*g_build(idir)%vals(1,iseg)
        x(n) = x(n-1) + dx

      END SELECT

    ENDDO

    RETURN
  END SUBROUTINE BLD_THEREST

  ! #######################################################################
  SUBROUTINE BLD_CONSTANTS(imax, vbeg,vend, iopt1_loc,iopt2_loc, val_1,val_2,val_3,val_4, a,b,c,d,e)

    IMPLICIT NONE

    TINTEGER imax, iopt1_loc, iopt2_loc
    TREAL vbeg, vend, val_1, val_2, val_3, val_4
    TREAL a, b, c, d, e

    TREAL x1, x2, x3, x4, valmx
    TREAL z1, z2, z3, z4
    TINTEGER indx1, indx2, indx3, indx4

    ! ######################################
    ! # Colonius, Lele and Moin Stretching #
    ! ######################################

    IF (iopt1_loc .EQ. 1) THEN
      x2 = val_4 - vbeg
      x3 = vend - vbeg

      !
      ! Calculate Constants
      ! -------------------

      a = FLOAT(imax-1) * val_1
      b = (a * (1.0 + val_2/val_1) - x3)/(val_2/val_1)
      c = val_3/val_1 - 1.0
      c = LOG(val_2/(c * val_1))/(b - x2)
      d = val_2/(c * val_1)
      e = vbeg

      !
      ! Force Values to be Exact at Domain Edges (x=0 at Xi=0)
      ! ------------------------------------------------------

      valmx = a + d * LOG(EXP(c*(a - b)) + 1.0 - EXP(-b*c))
      a = a * (x3/valmx)
      b = b * (x3/valmx)
      c = c / (x3/valmx)
      d = d * (x3/valmx)

      ! ######################################
      ! # Second order polynomial stretching #
      ! ######################################

    ELSEIF (iopt1_loc .EQ. 2) THEN

      !
      ! Points to be Used in Calculating Stretching Constants
      ! -----------------------------------------------------

      ! Clustering points at i=1

      IF (iopt2_loc .EQ. 1) THEN
        x1 = vbeg
        indx1 = 1
        x2 = x1 + val_1
        indx2 = 2
        x3 = vend
        indx3 = imax

        ! Clustering points at i=i_max

      ELSEIF (iopt2_loc .EQ. 2) THEN
        x1 = vbeg
        indx1 = 1
        x2 = vend - val_1
        indx2 = imax-1
        x3 = vend
        indx3 = imax

        ! Unknown stretching

      ELSE
        WRITE(*,*) 'Unknown stretching option iopt2_loc =',iopt2_loc
        STOP
      ENDIF

      !
      ! Calculate Stretching Constants
      ! ------------------------------

      z1 = FLOAT(indx1-1)/FLOAT(imax-1)
      z2 = FLOAT(indx2-1)/FLOAT(imax-1)
      z3 = FLOAT(indx3-1)/FLOAT(imax-1)

      a = (-(x3*z1**2*z2) + x3*z1*z2**2 + x2*z1**2*z3 - &
          x1*z2**2*z3 - x2*z1*z3**2 + x1*z2*z3**2)/&
          ((-z1 + z2)*(-z1 + z3)*(-z2 + z3))
      b = (-(x2*z1**2) + x3*z1**2 + x1*z2**2 - x3*z2**2 - x1*z3**2 + &
          x2*z3**2)/((-z1 + z2)*(-z1 + z3)*(-z2 + z3))
      c = (x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3)/&
          ((-z1 + z2)*(-z1 + z3)*(-z2 + z3))
      d = 0.0

      !
      ! Force Values to be Exact at Indx=1
      ! ----------------------------------

      a = a - (a + b*z1 + c*z1*z1 + d*z1*z1*z1 - x1)

      ! #####################################
      ! # Third order polynomial stretching #
      ! #####################################

    ELSEIF (iopt1_loc .EQ. 3) THEN

      !
      ! Points to be Used in Calculating Stretching Constants
      ! -----------------------------------------------------

      ! Clustering Points Both Ends

      IF (iopt2_loc .EQ. 1) THEN
        x1 = vbeg
        indx1 = 1
        x2 = x1 + val_1
        indx2 = 2
        x3 = vend - val_2
        indx3 = imax - 1
        x4 = vend
        indx4 = imax

        ! Clustering Points at an Internal Point

      ELSEIF (iopt2_loc .EQ. 2) THEN
        x1 = vbeg
        indx1 = 1
        x2 = val_2 - val_1/2.0
        indx2 = INT(val_3 * FLOAT(imax))
        x3 = val_2 + val_1/2.0
        indx3 = indx2 + 1
        x4 = vend
        indx4 = imax

        !  Unknown Stretching Option

      ELSE
        WRITE(*,*) 'Unknown stretching option iopt2_loc =',iopt2_loc
        STOP
      ENDIF

      !
      ! Calculate Stretching Constants
      ! ------------------------------

      z1 = FLOAT(indx1-1)/FLOAT(imax-1)
      z2 = FLOAT(indx2-1)/FLOAT(imax-1)
      z3 = FLOAT(indx3-1)/FLOAT(imax-1)
      z4 = FLOAT(indx4-1)/FLOAT(imax-1)

      a = (x4*z1**3*z2**2*z3 - x4*z1**2*z2**3*z3 -&
          x4*z1**3*z2*z3**2 + x4*z1*z2**3*z3**2 + x4*z1**2*z2*z3**3 -&
          x4*z1*z2**2*z3**3 - x3*z1**3*z2**2*z4 + x3*z1**2*z2**3*z4 +&
          x2*z1**3*z3**2*z4 - x1*z2**3*z3**2*z4 - x2*z1**2*z3**3*z4 +&
          x1*z2**2*z3**3*z4 + x3*z1**3*z2*z4**2 - x3*z1*z2**3*z4**2 -&
          x2*z1**3*z3*z4**2 + x1*z2**3*z3*z4**2 + x2*z1*z3**3*z4**2 -&
          x1*z2*z3**3*z4**2 - x3*z1**2*z2*z4**3 + x3*z1*z2**2*z4**3 +&
          x2*z1**2*z3*z4**3 - x1*z2**2*z3*z4**3 - x2*z1*z3**2*z4**3 +&
          x1*z2*z3**2*z4**3)/&
          ((-z1 + z2)*(-z1 + z3)*(-z2 + z3)*(-z1 + z4)*(-z2 + z4)*&
          (-z3 + z4))
      b = (x3*z1**3*z2**2 - x4*z1**3*z2**2 - x3*z1**2*z2**3 +&
          x4*z1**2*z2**3 - x2*z1**3*z3**2 + x4*z1**3*z3**2 +&
          x1*z2**3*z3**2 - x4*z2**3*z3**2 + x2*z1**2*z3**3 -&
          x4*z1**2*z3**3 - x1*z2**2*z3**3 + x4*z2**2*z3**3 +&
          x2*z1**3*z4**2 - x3*z1**3*z4**2 - x1*z2**3*z4**2 +&
          x3*z2**3*z4**2 + x1*z3**3*z4**2 - x2*z3**3*z4**2 -&
          x2*z1**2*z4**3 + x3*z1**2*z4**3 + x1*z2**2*z4**3 -&
          x3*z2**2*z4**3 - x1*z3**2*z4**3 + x2*z3**2*z4**3)/&
          ((-z1 + z2)*(-z1 + z3)*(-z2 + z3)*(-z1 + z4)*(-z2 + z4)*&
          (-z3 + z4))
      c = (-x3*z1**3*z2 + x4*z1**3*z2 + x3*z1*z2**3 - x4*z1*z2**3 +&
          x2*z1**3*z3 - x4*z1**3*z3 - x1*z2**3*z3 + x4*z2**3*z3 -&
          x2*z1*z3**3 + x4*z1*z3**3 + x1*z2*z3**3 - x4*z2*z3**3 -&
          x2*z1**3*z4 + x3*z1**3*z4 + x1*z2**3*z4 - x3*z2**3*z4 -&
          x1*z3**3*z4 + x2*z3**3*z4 + x2*z1*z4**3 - x3*z1*z4**3 -&
          x1*z2*z4**3 + x3*z2*z4**3 + x1*z3*z4**3 - x2*z3*z4**3)/&
          ((-z1 + z2)*(-z1 + z3)*(-z2 + z3)*(-z1 + z4)*(-z2 + z4)*&
          (-z3 + z4))
      d = -(((-z1 + z2)*(x2*z1 - x3*z1 - x1*z2 + x3*z2 &
          + x1*z3 - x2*z3)*&
          (-z1 + z4)*(-z2 + z4) +&
          (-z1 + z2)*(z1 - z3)*(z2 - z3)*&
          (-(x2*z1) + x4*z1 + x1*z2 - x4*z2 - x1*z4 + x2*z4))/&
          ((-z1 + z2)**2*(z1 - z3)*(z2 - z3)*(z1 - z4)*(z2 - z4)*&
          (-z3 + z4)))

      !
      ! Force Values to be Exact at Indx=1
      ! ----------------------------------

      a = a - (a + b*z1 + c*z1*z1 + d*z1*z1*z1 - x1)

    ENDIF

    RETURN
  END SUBROUTINE BLD_CONSTANTS

END MODULE GRID_LOCAL
