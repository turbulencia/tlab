subroutine SL_LOWER_BOUNDARY(imax, jmax, kmax, jmin_loc, amin, y, a, at, surface, wrk2d)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, jmin_loc
    real(wp) y(jmax), amin
    real(wp) a(*), at(jmax, imax*kmax)
    real(wp) surface(*), wrk2d(imax*kmax)

! -------------------------------------------------------------------
    integer(wi) i, j, jkmax, ikmax

! ###################################################################
    jkmax = jmax*kmax
    ikmax = imax*kmax

! -------------------------------------------------------------------
! Make y direction the first one; x becomes the last one
! -------------------------------------------------------------------
    call TLab_Transpose(a, imax, jkmax, imax, at, jkmax)

    do i = 1, ikmax
        do j = jmin_loc + 1, jmax
            if (at(j, i) > amin) exit
        end do
        if (j == jmax + 1) then
            wrk2d(i) = y(jmax)
        else
            wrk2d(i) = y(j - 1) + (y(j) - y(j - 1))/(at(j, i) - at(j - 1, i))*(amin - at(j - 1, i))
        end if
    end do

! -------------------------------------------------------------------
! Put array in right order
! -------------------------------------------------------------------
    call TLab_Transpose(wrk2d, kmax, imax, kmax, surface, imax)

    return
end subroutine SL_LOWER_BOUNDARY

!########################################################################
!########################################################################
subroutine SL_UPPER_BOUNDARY(imax, jmax, kmax, jmax_loc, amin, y, a, at, surface, wrk2d)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax, jmax_loc
    real(wp) y(jmax), amin
    real(wp) a(*), at(jmax, imax*kmax)
    real(wp) surface(*), wrk2d(imax*kmax)

! -------------------------------------------------------------------
    integer(wi) i, j, jkmax, ikmax

! ###################################################################
    jkmax = jmax*kmax
    ikmax = imax*kmax

! -------------------------------------------------------------------
! Make y direction the first one; x becomes the last one
! -------------------------------------------------------------------
    call TLab_Transpose(a, imax, jkmax, imax, at, jkmax)

    do i = 1, ikmax
        do j = jmax_loc - 1, 1, -1
            if (at(j, i) > amin) exit
        end do
        if (j == 0) then
            wrk2d(i) = y(1)
        else
            wrk2d(i) = y(j + 1) + (y(j) - y(j + 1))/(at(j, i) - at(j + 1, i))*(amin - at(j + 1, i))
        end if
    end do

! -------------------------------------------------------------------
! Put array in right order
! -------------------------------------------------------------------
    call TLab_Transpose(wrk2d, kmax, imax, kmax, surface, imax)

    return
end subroutine SL_UPPER_BOUNDARY

!########################################################################
!########################################################################
subroutine BOUNDARY_LOWER_INT1(imax, jmax, kmax, avalue, y, a, at, surface, wrk2d)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax
    real(wp), dimension(jmax), intent(IN) :: y
    integer(1), intent(IN) :: avalue
    integer(1), dimension(jmax, imax*kmax), intent(IN) :: a ! real shape is imax*jmax*kmax
    integer(1), dimension(jmax, imax*kmax), intent(INOUT) :: at(jmax, imax*kmax)
    real(wp), dimension(imax*kmax), intent(INOUT) :: wrk2d
    real(wp), dimension(imax*kmax), intent(OUT) :: surface

! -------------------------------------------------------------------
    integer(wi) i, j, jkmax, ikmax

! ###################################################################
    jkmax = jmax*kmax
    ikmax = imax*kmax

! Make y direction the first one; x becomes the last one
    if (imax > 1) then
        call TLab_Transpose_INT1(a, imax, jkmax, imax, at, jkmax)
    else
        at = a
    end if

    do i = 1, ikmax
        do j = 1, jmax
            if (at(j, i) == avalue) exit
        end do
        if (j == jmax + 1) then
            wrk2d(i) = y(jmax) ! in case the avalue is never attained
        else
            wrk2d(i) = y(j)
        end if
    end do

! Put array in right order
    call TLab_Transpose(wrk2d, kmax, imax, kmax, surface, imax)

    return
end subroutine BOUNDARY_LOWER_INT1

! ###################################################################
! ###################################################################

subroutine BOUNDARY_UPPER_INT1(imax, jmax, kmax, avalue, y, a, at, surface, wrk2d)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi) imax, jmax, kmax
    real(wp), dimension(jmax), intent(IN) :: y
    integer(1), intent(IN) :: avalue
    integer(1), dimension(jmax, imax*kmax), intent(IN) :: a ! real shape is imax*jmax*kmax
    integer(1), dimension(jmax, imax*kmax), intent(INOUT) :: at(jmax, imax*kmax)
    real(wp), dimension(imax*kmax), intent(INOUT) :: wrk2d
    real(wp), dimension(imax*kmax), intent(OUT) :: surface

! -------------------------------------------------------------------
    integer(wi) i, j, jkmax, ikmax

! ###################################################################
    jkmax = jmax*kmax
    ikmax = imax*kmax

! Make y direction the first one; x becomes the last one
    if (imax > 1) then
        call TLab_Transpose_INT1(a, imax, jkmax, imax, at, jkmax)
    else
        at = a
    end if

    do i = 1, ikmax
        do j = jmax, 1, -1
            if (at(j, i) == avalue) exit
        end do
        if (j == 0) then; wrk2d(i) = y(1) ! in case the avalue is never attained
        else; wrk2d(i) = y(min(j + 1, jmax)); end if ! +1 so that it coincides with lower boundary
    end do

! Put array in right order
    call TLab_Transpose(wrk2d, kmax, imax, kmax, surface, imax)

    return
end subroutine BOUNDARY_UPPER_INT1

! ###################################################################
! ###################################################################

function UPPER_THRESHOLD(jmax, uc, u, y)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: jmax
    real(wp), dimension(jmax), intent(IN) :: u, y
    real(wp), intent(IN) :: uc
    real(wp) UPPER_THRESHOLD

! -------------------------------------------------------------------
    integer(wi) j

! ###################################################################
    UPPER_THRESHOLD = 0.0_wp
    do j = jmax - 1, 1, -1
        if ((u(j) - uc)*(u(j + 1) - uc) < 0.0_wp) then
            UPPER_THRESHOLD = y(j) + (uc - u(j))*(y(j + 1) - y(j))/(u(j + 1) - u(j))
            exit
        end if
    end do

    return
end function UPPER_THRESHOLD

! ###################################################################
! ###################################################################

function LOWER_THRESHOLD(jmax, uc, u, y)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: jmax
    real(wp), dimension(jmax), intent(IN) :: u, y
    real(wp), intent(IN) :: uc
    real(wp) LOWER_THRESHOLD

! -------------------------------------------------------------------
    integer(wi) j

! ###################################################################
    LOWER_THRESHOLD = 0.0_wp
    do j = 1, jmax - 1
        if ((u(j) - uc)*(u(j + 1) - uc) < 0.0_wp) then
            LOWER_THRESHOLD = y(j) + (uc - u(j))*(y(j + 1) - y(j))/(u(j + 1) - u(j))
            exit
        end if
    end do

    return
end function LOWER_THRESHOLD

! ###################################################################
! ###################################################################

subroutine DELTA_X(imax, jmax, y, a, delta, delta_d, delta_u, A2, eta)
    use TLab_Constants, only: wp, wi
    implicit none

    integer(wi), intent(IN) :: imax, jmax
    real(wp), dimension(jmax), intent(IN) :: y
    real(wp), dimension(imax, jmax), intent(IN) :: a
    real(wp), intent(IN) :: A2, eta
    real(wp), dimension(imax), intent(OUT) :: delta_d, delta_u, delta

! -------------------------------------------------------------------
    integer(wi) i, j
    real(wp) DA, A_05, y_center

! ###################################################################
    y_center = 0.5_wp*(y(jmax/2) + y(jmax/2 + 1))

    do i = 1, imax
        DA = 0.5_wp*(a(i, jmax/2) + a(i, jmax/2 + 1)) - A2
        A_05 = A2 + eta*DA

        do j = 1, jmax/2
            if (a(i, j) <= A_05 .and. a(i, j + 1) > A_05) then
                delta_d(i) = y(j) + (A_05 - a(i, j))*(y(j + 1) - y(j))/(a(i, j + 1) - a(i, j))
            end if
        end do
        delta_d(i) = y_center - delta_d(i)

        do j = jmax/2 + 1, jmax
            if (a(i, j) > A_05 .and. a(i, j + 1) <= A_05) then
                delta_u(i) = y(j) + (A_05 - a(i, j))*(y(j + 1) - y(j))/(a(i, j + 1) - a(i, j))
            end if
        end do
        delta_u(i) = delta_u(i) - y_center

        delta(i) = 0.5_wp*(delta_d(i) + delta_u(i))
    end do

    return
end subroutine DELTA_X
