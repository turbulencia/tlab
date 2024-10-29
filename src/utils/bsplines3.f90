!########################################################################
!#
!# 3rd-order B-Spline interpolation for periodic sequence f of imax points,
!# x, between 0 and 1, gives the relative position of the desired point
!# in the interval (i,i+1)
!#
!########################################################################

function BSPLINES3P(f, imax, i, x)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, i
    real(wp) f(imax), x
    real(wp) BSPLINES3P

! -----------------------------------------------------------------------
    real(wp) b1, b2, b3, b4
    integer(wi) i1, i3, i4, iw

! #######################################################################
    b1 = (1.0_wp - x)*(1.0_wp - x)*(1.0_wp - x)
    b2 = 4.0_wp - 6.0_wp*x*x + 3.0_wp*x*x*x
    b3 = 1.0_wp + 3.0_wp*x + 3.0_wp*x*x - 3.0_wp*x*x*x
    b4 = x*x*x

    iw = imax + 1 - i
    i1 = imax - mod(iw, imax)
    i3 = mod(i, imax) + 1
    iw = i + 1
    i4 = mod(iw, imax) + 1

    BSPLINES3P = (b1*f(i1) + b2*f(i) + b3*f(i3) + b4*f(i4))/6.0_wp

end function BSPLINES3P

!########################################################################
!# Non-periodic BCs
!# BCs set by ghost cell with function value equal the first/last point
!########################################################################
function BSPLINES3(f, imax, i, x)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, i
    real(wp) f(imax), x
    real(wp) BSPLINES3

! -----------------------------------------------------------------------
    real(wp) b1, b2, b3, b4
    integer(wi) i1, i3, i4

! #######################################################################
    b1 = (1.0_wp - x)*(1.0_wp - x)*(1.0_wp - x)
    b2 = 4.0_wp - 6.0_wp*x*x + 3.0_wp*x*x*x
    b3 = 1.0_wp + 3.0_wp*x + 3.0_wp*x*x - 3.0_wp*x*x*x
    b4 = x*x*x

    i1 = max(i - 1, 1)
    i3 = i + 1
    i4 = min(i3 + 1, imax)

    BSPLINES3 = (b1*f(i1) + b2*f(i) + b3*f(i3) + b4*f(i4))/6.0_wp

end function BSPLINES3

!########################################################################
!#
!# 3rd-order B-Spline interpolation for periodic sequence f of imax points
!# in non-uniform case.
!# Cox-de Boor algorithm (Pozrikidis, p409)
!########################################################################
function BSPLINES3_NU(f, t, imax, i, tint)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, i
    real(wp), dimension(imax) :: f, t
    real(wp) tint
    real(wp) BSPLINES3_NU

! -----------------------------------------------------------------------
    real(wp) b1, b2, b3, b4, aux1, aux2, aux3
    real(wp) t1, t2, t3, t4, t5, dt, dt1, dt2, dt3
    integer(wi) i1, i3, i4, iw

! #######################################################################
! Coefficients

    dt = (t(imax) - t(1))/real(imax - 1, wp)
    dt1 = t(2) - t(1)
    dt2 = t(3) - t(1)
    dt3 = t(imax) - t(imax - 1)

    if (i == 1) then
        t1 = t(1) - dt - dt3
        t2 = t(1) - dt
        t3 = t(2)
        t4 = t(3)
        t5 = t(4)
    else if (i == 2) then
        t1 = t(1) - dt
        t2 = t(1)
        t3 = t(3)
        t4 = t(4)
        t5 = t(5)
    else if (i == imax - 2) then
        t1 = t(imax - 4)
        t2 = t(imax - 3)
        t3 = t(imax - 1)
        t4 = t(imax)
        t5 = t(imax) + dt
    else if (i == imax - 1) then
        t1 = t(imax - 3)
        t2 = t(imax - 2)
        t3 = t(imax)
        t4 = t(imax) + dt
        t5 = t(imax) + dt + dt1
    else if (i == imax) then
        t1 = t(imax - 2)
        t2 = t(imax - 1)
        t3 = t(imax) + dt
        t4 = t(imax) + dt + dt1
        t5 = t(imax) + dt + dt2
    else
        t1 = t(i - 2)
        t2 = t(i - 1)
        t3 = t(i + 1)
        t4 = t(i + 2)
        t5 = t(i + 3)
    end if

    aux1 = ((tint - t2)*(t3 - tint)/(t3 - t2) + (t4 - tint)*(tint - t(i))/(t4 - t(i)))/ &
           ((t4 - t2)*(t3 - t(i)))
    aux2 = ((t3 - t1)*(t3 - t2)*(t3 - t(i)))
    aux3 = ((t5 - t(i))*(t4 - t(i))*(t3 - t(i)))

    b1 = (t3 - tint)**3.0_wp/aux2
    b2 = (t3 - tint)**2.0_wp*(tint - t1)/aux2 + aux1*(t4 - tint)
    b3 = (tint - t(i))**2.0_wp*(t5 - tint)/aux3 + aux1*(tint - t2)
    b4 = (tint - t(i))**3.0_wp/aux3

! Interpolation

    iw = imax + 1 - i
    i1 = imax - mod(iw, imax)
    i3 = mod(i, imax) + 1
    iw = i + 1
    i4 = mod(iw, imax) + 1

    BSPLINES3_NU = b1*f(i1) + b2*f(i) + b3*f(i3) + b4*f(i4)

end function BSPLINES3_NU
