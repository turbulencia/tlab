#include "types.h"

SUBROUTINE BLD_GEN(idir, x, imax, scalex)
  USE GRID_LOCAL
  IMPLICIT NONE

  TINTEGER, INTENT(IN   ) :: idir, imax
  TREAL,    INTENT(INOUT) :: x(imax), scalex

  ! -----------------------------------------------------------------------
  TREAL a, b, c, d, e
  TREAL eta, deta, xleft, dx
  TINTEGER iloc, m, iseg

  ! #######################################################################
  iloc = 1
  IF ( g_build(idir)%mirrored ) iloc = imax/2 ! mirrored case; first point in array is imax/2

  ! ##############################
  ! # Build Grid in Each Segment #
  ! ##############################
  x(iloc) = C_0_R ! set origin at zero
  xleft   = x(iloc)
  dx      = C_1_R
  DO iseg = 1,g_build(idir)%nseg
    ! Calculate constants of stretching function, when needed
    CALL BLD_CONSTANTS(g_build(idir)%SIZE(iseg), xleft,g_build(idir)%END(iseg),&
        g_build(idir)%opts(1,iseg),g_build(idir)%opts(2,iseg),&
        g_build(idir)%vals(1,iseg),g_build(idir)%vals(2,iseg),&
        g_build(idir)%vals(3,iseg),g_build(idir)%vals(4,iseg), a,b,c,d,e)

    ! Go from computational eta in [0,1] to physical domains
    IF ( g_build(idir)%SIZE(iseg) .EQ. 1 ) THEN
      deta = C_1_R
    ELSE;
      deta = C_1_R /M_REAL(g_build(idir)%SIZE(iseg)-1)
    ENDIF
    DO m = 2,g_build(idir)%SIZE(iseg)
      eta  = M_REAL(m-1)*deta
      iloc = iloc + 1

      SELECT CASE(  g_build(idir)%opts(1,iseg) )
      CASE( 0 )     ! uniform
        x(iloc) = xleft + eta*(g_build(idir)%END(iseg)-xleft)
      CASE( 1 )     ! Colonius, Lele and Moin
        x(iloc) = e + a*eta + d*LOG( EXP(c*(a*eta - b)) + C_1_R - EXP(-b*c) )
      CASE( 2,3 )   ! 2nd- 3rd-Order Polynomial Stretching
        x(iloc) = a + b*eta + c*eta*eta + d*eta*eta*eta
      CASE( 4 )     ! Geometric prograssion
        dx      = dx*g_build(idir)%vals(1,iseg)
        x(iloc) = x(iloc-1) + dx

      END SELECT

    ENDDO

    xleft = g_build(idir)%END(iseg)

  ENDDO

  ! ################################
  ! # Check Dimensions of new grid #
  ! ################################
  IF ( iloc .NE. imax ) THEN
    WRITE(*,*) '- Error after packing in direction', idir
    WRITE(*,*) '  Processing up to ',iloc,' and dimension is ',imax
    STOP
  ENDIF

  ! ##########################################################
  ! # Set Scale and Force Endpoints to be Exactly as Entered #
  ! ##########################################################
  IF ( g_build(idir)%mirrored ) THEN; iloc = imax/2 ! mirrored case
  ELSE;                                   iloc = 1;     ENDIF

    IF ( g_build(idir)%opts(1,1) .EQ. 4 ) THEN ! Geometric prograssion
      DO m = iloc,imax
        IF ( x(imax) .NE. C_0_R ) THEN
          x(m) = x(m)/x(imax)*scalex
        ENDIF
      ENDDO
    ELSE
      IF ( imax .GT. 1 ) x(imax) = g_build(idir)%END(g_build(idir)%nseg)
    ENDIF

    RETURN
  END SUBROUTINE BLD_GEN

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
