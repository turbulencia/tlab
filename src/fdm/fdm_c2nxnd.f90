!########################################################################
!#
!# Compact FDMs for non-uniform grids from Shukla and Zhong (2005), JCP, 204, 404–429
!#
!# We set up mmax linear system of size nmax of the form
!#                         A up = B u
!# in this routine. The matrix A is tridiagonal, and B is pentadiagonal,
!# except for the boundary points.
!#
!# We normalize system to get a diagonal 1 in B.
!#
!# 6th-order approximation to 2nd-order derivative:
!#
!# Equations (15) for the interior points (beware that there is a typo in paper).
!# Equations (16) for the first/last points.
!# Table B.2 for the second/second-to-last points.
!#
!########################################################################

!########################################################################
! Create diagonals
!########################################################################
subroutine FDM_C2N6ND_INITIALIZE(nmax, x, lhs, rhs)
    use TLAB_CONSTANTS
    implicit none

    integer(wi), intent(IN) :: nmax
    real(wp), dimension(nmax), intent(IN) :: x
    real(wp), dimension(nmax, 3), intent(OUT) :: lhs ! LHS diagonals (#=3)
    real(wp), dimension(nmax, 4), intent(OUT) :: rhs ! RHS diagonals (#=5-1 because of normalization)

! -------------------------------------------------------------------
    integer(wi) n, j
    real(wp) am1, a, ap1                            ! Left-hand side
    real(wp) bm2, bm1, b, bp1, bp2, bp3             ! Right-hand side
    real(wp) tmp1, tmp2, tmpa, tmpb, tmpd, tmpp, tmpq   ! Intermediate ops
    real(wp) dummy1, dummy2

! #######################################################################
! n = 1
! #######################################################################
    n = 1

! left-hand side
    am1 = 0.0_wp
    a = 1.0_wp
    ap1 = (x(n + 1) - x(n))*(x(n) - x(n + 2)) + (x(n + 1) - x(n))*(x(n) - x(n + 3)) - (x(n) - x(n + 2))*(x(n) - x(n + 3))
    ap1 = ap1/ &
   ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) + (x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

! right-hand side
    bm2 = 0.0_wp
    bm1 = 0.0_wp
    b = (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + 2.0_wp*(x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
    b = b/ &
        ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
    b = b*((x(n) - x(n + 2)) + (x(n) - x(n + 3)))/ &
        (x(n) - x(n + 1)) + 1.0_wp
    b = b/ &
        ((x(n) - x(n + 2))*(x(n) - x(n + 3)))
    tmp1 = ((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))/ &
           ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
    b = (b + tmp1/(x(n + 1) - x(n)))*2.0_wp

    bp1 = (x(n) - x(n + 2)) + (x(n) - x(n + 3)) + ap1*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
    bp1 = bp1*2.0_wp/ &
          ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

    bp2 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 3))) + &
          (x(n) - x(n + 3))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 3)))
    bp2 = bp2/ &
          ((x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) - (x(n + 2) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 3))))
    bp2 = bp2*2.0_wp*(x(n + 1) - x(n))/ &
          ((x(n + 2) - x(n + 1))*(x(n + 2) - x(n))*(x(n + 2) - x(n + 3)))

    bp3 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 2))) + &
          (x(n) - x(n + 2))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 2)))
    bp3 = bp3/ &
          ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) - (x(n + 3) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2))))
    bp3 = bp3*2.0_wp*(x(n + 1) - x(n))/ &
          ((x(n + 3) - x(n + 1))*(x(n + 3) - x(n))*(x(n + 3) - x(n + 2)))

! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
    tmp1 = 1.0_wp/b

    lhs(n, 1) = 0.0_wp
    lhs(n, 2) = a*tmp1
    lhs(n, 3) = ap1*tmp1

    rhs(n, 1) = bp3*tmp1 ! saving here the last term that goes into the 3. superdiagonal
    rhs(n, 2) = 0.0_wp
    rhs(n, 3) = bp1*tmp1
    rhs(n, 4) = bp2*tmp1

! #######################################################################
! n = 2
! #######################################################################
    n = 2

    bm2 = 0.0_wp; bp2 = 0.0_wp

    tmp1 = 1.0_wp/ &
           ((x(n + 1) - x(n))*(x(n) - x(n - 1)) + (x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1)))

    bm1 = (x(n + 1) - x(n))/(x(n + 1) - x(n - 1))
    b = -1.0_wp
    bp1 = (x(n) - x(n - 1))/(x(n + 1) - x(n - 1))

    am1 = (x(n) - x(n - 1))*(x(n) - x(n - 1)) + (x(n + 1) - x(n))*(x(n) - x(n - 1)) - (x(n + 1) - x(n))*(x(n + 1) - x(n))
    am1 = am1*bm1*tmp1
    a = 1.0_wp
    ap1 = (x(n + 1) - x(n))*(x(n + 1) - x(n)) + (x(n + 1) - x(n))*(x(n) - x(n - 1)) - (x(n) - x(n - 1))*(x(n) - x(n - 1))
    ap1 = ap1*bp1*tmp1

    bm1 = bm1*12.0_wp*tmp1
    b = b*12.0_wp*tmp1
    bp1 = bp1*12.0_wp*tmp1

! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
    tmp1 = 1.0_wp/b

    lhs(n, 1) = am1*tmp1
    lhs(n, 2) = a*tmp1
    lhs(n, 3) = ap1*tmp1

    rhs(n, 1) = 0.0_wp
    rhs(n, 2) = bm1*tmp1
    rhs(n, 3) = bp1*tmp1
    rhs(n, 4) = 0.0_wp

! #######################################################################
! 2 < n < nmax-1
! #######################################################################
    do n = 3, nmax - 2

        tmpa = (x(n - 1) - x(n + 1))* &
               (1.0_wp/(x(n - 1) - x(n - 2)) + 1.0_wp/(x(n - 1) - x(n + 2)) + 1.0_wp/(x(n - 1) - x(n)))
        tmpb = (x(n + 1) - x(n - 1))* &
               (1.0_wp/(x(n + 1) - x(n - 2)) + 1.0_wp/(x(n + 1) - x(n + 2)) + 1.0_wp/(x(n + 1) - x(n)))
        tmpd = 1.0_wp/ &                                    ! constant 2/D in Table 2 in Shukla&Zhong(2005)
               ((tmpb + 2.0_wp)*(tmpa + 2.0_wp) - 1.0_wp)
        tmpa = tmpa + 1.0_wp                                 ! constant a in ...
        tmpb = tmpb + 1.0_wp                                 ! constant b in ...

        tmpp = (x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1))* &         ! constant p
               (1.0_wp/((x(n - 1) - x(n))*(x(n - 1) - x(n + 2))) &
                + 1.0_wp/((x(n - 1) - x(n - 2))*(x(n - 1) - x(n + 2))) &
                + 1.0_wp/((x(n - 1) - x(n - 2))*(x(n - 1) - x(n))))
        tmpq = (x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1))* &         ! constant q
               (1.0_wp/((x(n + 1) - x(n))*(x(n + 1) - x(n + 2))) &
                + 1.0_wp/((x(n + 1) - x(n - 2))*(x(n + 1) - x(n + 2))) &
                + 1.0_wp/((x(n + 1) - x(n - 2))*(x(n + 1) - x(n))))

! -------------------------------------------------------------------
        a = 1.0_wp

! -------------------------------------------------------------------
        dummy1 = 1.0_wp + tmpb
        dummy2 = -tmpb

        tmp1 = (x(n) - x(n - 2))*(x(n) - x(n + 2))/((x(n - 1) - x(n - 2))*(x(n - 1) - x(n + 2)))
        tmp2 = dummy1*(1.0_wp + (x(n) - x(n + 1))/(x(n) - x(n - 1))) &
               + dummy2*(2.0_wp*(x(n) - x(n + 1)) + (x(n) - x(n - 1)))/(x(n + 1) - x(n - 1))
        am1 = tmp1*tmp2

        tmp1 = ((x(n) - x(n - 2)) + (x(n) - x(n + 2)))/((x(n - 1) - x(n - 2))*(x(n - 1) - x(n + 2)))*(x(n) - x(n + 1))
        tmp2 = dummy1 &
               + dummy2*(x(n) - x(n - 1))/(x(n + 1) - x(n - 1))
        tmp2 = tmp1*tmp2

        am1 = (am1 + tmp2)*tmpd

! -------------------------------------------------------------------
        dummy1 = 4.0_wp + 2.0_wp*(1.0_wp + tmpb)*(tmpa - 2.0_wp) + 2.0_wp*tmpp*(1.0_wp + tmpb)
        dummy2 = 1.0_wp + (1.0_wp - 2.0_wp*tmpa)*(2.0_wp*tmpb - 1.0_wp) - 2.0_wp*tmpp*tmpb

        tmp1 = (x(n) - x(n - 2))*(x(n) - x(n + 2))/((x(n - 1) - x(n - 2))*(x(n - 1) - x(n + 2)))
        tmp2 = dummy1*(1.0_wp + (x(n) - x(n + 1))/(x(n) - x(n - 1))) &
               + dummy2*(2.0_wp*(x(n) - x(n + 1)) + (x(n) - x(n - 1)))/(x(n + 1) - x(n - 1)) &
               + 2.0_wp/tmpd*(x(n + 1) - x(n - 1))/(x(n) - x(n - 1))
        bm1 = tmp1*tmp2

        tmp1 = ((x(n) - x(n - 2)) + (x(n) - x(n + 2)))/((x(n - 1) - x(n - 2))*(x(n - 1) - x(n + 2)))*(x(n) - x(n + 1))
        tmp2 = dummy1 &
               + dummy2*(x(n) - x(n - 1))/(x(n + 1) - x(n - 1)) &
               + 2.0_wp/tmpd*(x(n + 1) - x(n - 1))/(x(n) - x(n - 1))
        tmp2 = tmp1*tmp2

        bm1 = (bm1 + tmp2)*tmpd/((x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1)))

! -------------------------------------------------------------------
        dummy1 = 1.0_wp + tmpa
        dummy2 = -tmpa

        tmp1 = (x(n) - x(n + 2))*(x(n) - x(n - 2))/((x(n + 1) - x(n + 2))*(x(n + 1) - x(n - 2)))
        tmp2 = dummy1*(1.0_wp + (x(n) - x(n - 1))/(x(n) - x(n + 1))) &
               + dummy2*(2.0_wp*(x(n) - x(n - 1)) + (x(n) - x(n + 1)))/(x(n - 1) - x(n + 1))
        ap1 = tmp1*tmp2

        tmp1 = ((x(n) - x(n - 2)) + (x(n) - x(n + 2)))/((x(n + 1) - x(n - 2))*(x(n + 1) - x(n + 2)))*(x(n) - x(n - 1))
        tmp2 = dummy1 &
               + dummy2*(x(n) - x(n + 1))/(x(n - 1) - x(n + 1))
        tmp2 = tmp1*tmp2

        ap1 = (ap1 + tmp2)*tmpd

! -------------------------------------------------------------------
        dummy1 = 4.0_wp + 2.0_wp*(1.0_wp + tmpa)*(tmpb - 2.0_wp) + 2.0_wp*tmpq*(1.0_wp + tmpa)
        dummy2 = 1.0_wp + (1.0_wp - 2.0_wp*tmpa)*(2.0_wp*tmpb - 1.0_wp) - 2.0_wp*tmpq*tmpa

        tmp1 = (x(n) - x(n - 2))*(x(n) - x(n + 2))/((x(n + 1) - x(n - 2))*(x(n + 1) - x(n + 2)))
        tmp2 = dummy1*(1.0_wp + (x(n) - x(n - 1))/(x(n) - x(n + 1))) &
               + dummy2*(2.0_wp*(x(n) - x(n - 1)) + (x(n) - x(n + 1)))/(x(n - 1) - x(n + 1)) &
               + 2.0_wp/tmpd*(x(n - 1) - x(n + 1))/(x(n) - x(n + 1))
        bp1 = tmp1*tmp2

        tmp1 = ((x(n) - x(n - 2)) + (x(n) - x(n + 2)))/((x(n + 1) - x(n - 2))*(x(n + 1) - x(n + 2)))*(x(n) - x(n - 1))
        tmp2 = dummy1 &
               + dummy2*(x(n) - x(n + 1))/(x(n - 1) - x(n + 1)) &
               + 2.0_wp/tmpd*(x(n - 1) - x(n + 1))/(x(n) - x(n + 1))
        tmp2 = tmp1*tmp2

        bp1 = (bp1 + tmp2)*tmpd/((x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1)))

! -------------------------------------------------------------------
        j = n

        dummy1 = (x(n + 1) - x(j))/(x(n + 1) - x(n - 1))
        dummy2 = (x(n - 1) - x(j))/(x(n + 1) - x(n - 1))

        tmp1 = (dummy1/dummy2 - dummy2/dummy1)*(1.0_wp - tmpa)*(tmpb - 1.0_wp) &
               + (1.0_wp - tmpa)*(2.0_wp/dummy1 + 2.0_wp/dummy2 - 1.0_wp/(dummy1*dummy1)) &
               - (tmpb - 1.0_wp)*(2.0_wp/dummy1 + 2.0_wp/dummy2 + 1.0_wp/(dummy2*dummy2)) &
               - (1.0_wp/dummy1 + 1.0_wp/dummy2)*(5.0_wp + 2.0_wp/(dummy1*dummy2))
        tmp2 = 1.0_wp/(dummy1*dummy1) + 1.0_wp/(dummy2*dummy2) + 1.0_wp/(dummy1*dummy2) &
               - ((tmpb - 1.0_wp)/dummy2*(1.0_wp - tmpa)/dummy1) &
               + ((tmpb - 1.0_wp)/dummy2 - (1.0_wp - tmpa)/dummy1)*(1.0_wp/dummy1 + 1.0_wp/dummy2)

        b = (x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1))/((x(n) - x(n - 2))*(x(n) - x(n + 2))) &
            - (x(n + 1) - x(n - 1))*(1.0_wp/(x(n) - x(n - 2)) + 1.0_wp/(x(n) - x(n + 2)))*(1.0_wp/dummy1 + 1.0_wp/dummy2) &
            + 1.0_wp/(dummy1*dummy2)
        b = b/tmpd &
            + tmp1*((x(n + 1) - x(n - 1))*(1.0_wp/(x(n) - x(n - 2)) + 1.0_wp/(x(n) - x(n + 2))) - (1.0_wp/dummy1 + 1.0_wp/dummy2))
        b = (b + tmp2)*2.0_wp*tmpd/((x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1)))

! -------------------------------------------------------------------
        j = n - 2

        dummy1 = (x(n + 1) - x(j))/(x(n + 1) - x(n - 1))
        dummy2 = (x(n - 1) - x(j))/(x(n + 1) - x(n - 1))

        tmp1 = (dummy1/dummy2 - dummy2/dummy1)*(1.0_wp - tmpa)*(tmpb - 1.0_wp) &
               + (1.0_wp - tmpa)*(2.0_wp/dummy1 + 2.0_wp/dummy2 - 1.0_wp/(dummy1*dummy1)) &
               - (tmpb - 1.0_wp)*(2.0_wp/dummy1 + 2.0_wp/dummy2 + 1.0_wp/(dummy2*dummy2)) &
               - (1.0_wp/dummy1 + 1.0_wp/dummy2)*(5.0_wp + 2.0_wp/(dummy1*dummy2))
        tmp2 = 1.0_wp/(dummy1*dummy1) + 1.0_wp/(dummy2*dummy2) + 1.0_wp/(dummy1*dummy2) &
               - ((tmpb - 1.0_wp)/dummy2*(1.0_wp - tmpa)/dummy1) &
               + ((tmpb - 1.0_wp)/dummy2 - (1.0_wp - tmpa)/dummy1)*(1.0_wp/dummy1 + 1.0_wp/dummy2)

        dummy1 = 1.0_wp/(x(n) - x(n + 1)) + 1.0_wp/(x(n) - x(n - 1))
        dummy1 = tmp1*(1.0_wp + (x(n) - x(j))*dummy1) &
                 + tmp2*(2.0_wp + (x(n) - x(j))*dummy1)*(x(n) - x(j))/(x(n + 1) - x(n - 1)) &
                 + 1.0_wp/tmpd*((x(n + 1) - x(n - 1))*dummy1)
        dummy1 = dummy1*(x(n + 1) - x(n - 1))/(x(n + 2) - x(n - 2))*(x(n) - x(n + 2))/(x(n) - x(n - 2)) ! for n-2

        dummy2 = (x(n) - x(j))/(x(n + 1) - x(n - 1))
        dummy2 = tmp1*dummy2 &
                 + tmp2*dummy2*dummy2 &
                 + 1.0_wp/tmpd
        dummy2 = dummy2*(x(n + 1) - x(n - 1))/(x(n + 2) - x(n - 2))*(x(n + 1) - x(n - 1))/(x(n) - x(n - 2)) ! for n-2

        bm2 = (dummy1 + dummy2)*tmpd/((x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1))) &
              *2.0_wp*(x(n) - x(n + 1))/(x(j) - x(n + 1))*(x(n) - x(n - 1))/(x(j) - x(n - 1))

! -------------------------------------------------------------------
        j = n + 2

        dummy1 = (x(n + 1) - x(j))/(x(n + 1) - x(n - 1))
        dummy2 = (x(n - 1) - x(j))/(x(n + 1) - x(n - 1))

        tmp1 = (dummy1/dummy2 - dummy2/dummy1)*(1.0_wp - tmpa)*(tmpb - 1.0_wp) &
               + (1.0_wp - tmpa)*(2.0_wp/dummy1 + 2.0_wp/dummy2 - 1.0_wp/(dummy1*dummy1)) &
               - (tmpb - 1.0_wp)*(2.0_wp/dummy1 + 2.0_wp/dummy2 + 1.0_wp/(dummy2*dummy2)) &
               - (1.0_wp/dummy1 + 1.0_wp/dummy2)*(5.0_wp + 2.0_wp/(dummy1*dummy2))
        tmp2 = 1.0_wp/(dummy1*dummy1) + 1.0_wp/(dummy2*dummy2) + 1.0_wp/(dummy1*dummy2) &
               - ((tmpb - 1.0_wp)/dummy2*(1.0_wp - tmpa)/dummy1) &
               + ((tmpb - 1.0_wp)/dummy2 - (1.0_wp - tmpa)/dummy1)*(1.0_wp/dummy1 + 1.0_wp/dummy2)

        dummy1 = 1.0_wp/(x(n) - x(n + 1)) + 1.0_wp/(x(n) - x(n - 1))
        dummy1 = tmp1*(1.0_wp + (x(n) - x(j))*dummy1) &
                 + tmp2*(2.0_wp + (x(n) - x(j))*dummy1)*(x(n) - x(j))/(x(n + 1) - x(n - 1)) &
                 + 1.0_wp/tmpd*((x(n + 1) - x(n - 1))*dummy1)
        dummy1 = dummy1*(x(n + 1) - x(n - 1))/(x(n + 2) - x(n - 2))*(x(n) - x(n - 2))/(x(n + 2) - x(n)) ! for n+2

        dummy2 = (x(n) - x(j))/(x(n + 1) - x(n - 1))
        dummy2 = tmp1*dummy2 &
                 + tmp2*dummy2*dummy2 &
                 + 1.0_wp/tmpd
        dummy2 = dummy2*(x(n + 1) - x(n - 1))/(x(n + 2) - x(n - 2))*(x(n + 1) - x(n - 1))/(x(n + 2) - x(n)) ! for n+2

        bp2 = (dummy1 + dummy2)*tmpd/((x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1))) &
              *2.0_wp*(x(n) - x(n + 1))/(x(j) - x(n + 1))*(x(n) - x(n - 1))/(x(j) - x(n - 1))

! -------------------------------------------------------------------
! if uniform, we should have ( 2/11 1 2/11 ) and ( 3/44 12/11 -51/22 12/11 3/44 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
        tmp1 = 1.0_wp/b

        lhs(n, 1) = am1*tmp1
        lhs(n, 2) = a*tmp1
        lhs(n, 3) = ap1*tmp1

        rhs(n, 1) = bm2*tmp1
        rhs(n, 2) = bm1*tmp1
        rhs(n, 3) = bp1*tmp1
        rhs(n, 4) = bp2*tmp1

    end do

! #######################################################################
! n = nmax-1; just a copy of case n = 2
! #######################################################################
    n = nmax - 1

    bm2 = 0.0_wp; bp2 = 0.0_wp

    tmp1 = 1.0_wp/ &
           ((x(n + 1) - x(n))*(x(n) - x(n - 1)) + (x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1)))

    bm1 = (x(n + 1) - x(n))/(x(n + 1) - x(n - 1))
    b = -1.0_wp
    bp1 = (x(n) - x(n - 1))/(x(n + 1) - x(n - 1))

    am1 = (x(n) - x(n - 1))*(x(n) - x(n - 1)) + (x(n + 1) - x(n))*(x(n) - x(n - 1)) - (x(n + 1) - x(n))*(x(n + 1) - x(n))
    am1 = am1*bm1*tmp1
    a = 1.0_wp
    ap1 = (x(n + 1) - x(n))*(x(n + 1) - x(n)) + (x(n + 1) - x(n))*(x(n) - x(n - 1)) - (x(n) - x(n - 1))*(x(n) - x(n - 1))
    ap1 = ap1*bp1*tmp1

    bm1 = bm1*12.0_wp*tmp1
    b = b*12.0_wp*tmp1
    bp1 = bp1*12.0_wp*tmp1

! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
    tmp1 = 1.0_wp/b

    lhs(n, 1) = am1*tmp1
    lhs(n, 2) = a*tmp1
    lhs(n, 3) = ap1*tmp1

    rhs(n, 1) = 0.0_wp
    rhs(n, 2) = bm1*tmp1
    rhs(n, 3) = bp1*tmp1
    rhs(n, 4) = 0.0_wp

! #######################################################################
! n = nmax; same as n = 1, but changing the signs of the increments w.r.t. n
! To understand it, e.g., define a new variable i = -j, where i is the
! discrete variable moving around n
! #######################################################################
    n = nmax

! left-hand side
    am1 = 0.0_wp
    a = 1.0_wp
    ap1 = (x(n - 1) - x(n))*(x(n) - x(n - 2)) + (x(n - 1) - x(n))*(x(n) - x(n - 3)) - (x(n) - x(n - 2))*(x(n) - x(n - 3))
    ap1 = ap1/ &
   ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) + (x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

! right-hand side
    bm2 = 0.0_wp
    bm1 = 0.0_wp
    b = (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + 2.0_wp*(x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
    b = b/ &
        ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
    b = b*((x(n) - x(n - 2)) + (x(n) - x(n - 3)))/ &
        (x(n) - x(n - 1)) + 1.0_wp
    b = b/ &
        ((x(n) - x(n - 2))*(x(n) - x(n - 3)))
    tmp1 = ((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))/ &
           ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
    b = (b + tmp1/(x(n - 1) - x(n)))*2.0_wp

    bp1 = (x(n) - x(n - 2)) + (x(n) - x(n - 3)) + ap1*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
    bp1 = bp1*2.0_wp/ &
          ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

    bp2 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 3))) + &
          (x(n) - x(n - 3))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 3)))
    bp2 = bp2/ &
          ((x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) - (x(n - 2) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 3))))
    bp2 = bp2*2.0_wp*(x(n - 1) - x(n))/ &
          ((x(n - 2) - x(n - 1))*(x(n - 2) - x(n))*(x(n - 2) - x(n - 3)))

    bp3 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 2))) + &
          (x(n) - x(n - 2))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 2)))
    bp3 = bp3/ &
          ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) - (x(n - 3) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2))))
    bp3 = bp3*2.0_wp*(x(n - 1) - x(n))/ &
          ((x(n - 3) - x(n - 1))*(x(n - 3) - x(n))*(x(n - 3) - x(n - 2)))

! if uniform, we should have ( 11 1 0 ) and ( -1 15 -27 13 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
    tmp1 = 1.0_wp/b

    lhs(n, 3) = 0.0_wp
    lhs(n, 2) = a*tmp1
    lhs(n, 1) = ap1*tmp1

    rhs(n, 4) = bp3*tmp1 ! saving here the last term that goes into the 3. subdiagonal
    rhs(n, 3) = 0.0_wp
    rhs(n, 2) = bp1*tmp1
    rhs(n, 1) = bp2*tmp1

    return
end subroutine FDM_C2N6ND_INITIALIZE

! #######################################################################
! Constructing forcing term
! #######################################################################
subroutine FDM_C2NXND_RHS(nmax, mmax, rhs, u, d)
    use TLAB_CONSTANTS
    implicit none

    integer(wi), intent(in) :: nmax, mmax ! m linear systems or size n
    real(wp), intent(in) :: rhs(nmax, 4)       ! RHS diagonals (#=5-1 because of normalization)
    real(wp), intent(in) :: u(mmax, nmax)         ! function
    real(wp), intent(out) :: d(mmax, nmax)         ! RHS

! -------------------------------------------------------------------
    integer(wi) n

! #######################################################################
! Boundary
    n = 1 ! rhs(1,1) contains 3. superdiagonal
    d(:, n) = &
        +u(:, n) &
        + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4) + u(:, n + 3)*rhs(n, 1)

    n = 2
    d(:, n) = u(:, n - 1)*rhs(n, 2) &
              + u(:, n) &
              + u(:, n + 1)*rhs(n, 3)

! Interior points
    do n = 3, nmax - 2
        d(:, n) = u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
                  + u(:, n) &
                  + u(:, n + 1)*rhs(n, 3) + u(:, n + 2)*rhs(n, 4)
    end do

! Boundary
    n = nmax - 1
    d(:, n) = u(:, n - 1)*rhs(n, 2) &
              + u(:, n) &
              + u(:, n + 1)*rhs(n, 3)

    n = nmax ! rhs(1,4) contains 3. subdiagonal
    d(:, n) = u(:, n - 3)*rhs(n, 4) + u(:, n - 2)*rhs(n, 1) + u(:, n - 1)*rhs(n, 2) &
              + u(:, n)

    return
end subroutine FDM_C2NXND_RHS

!########################################################################
!# 4th-order approximation to 2nd-order derivative:
!########################################################################
subroutine FDM_C2N4ND_INITIALIZE(nmax, x, lhs, rhs)
    use TLAB_CONSTANTS
    implicit none

    integer(wi), intent(in) :: nmax
    real(wp), intent(in) :: x(nmax)
    real(wp), intent(out) :: lhs(nmax, 3) ! LHS diagonals (#=3)
    real(wp), intent(out) :: rhs(nmax, 4) ! RHS diagonals (#=5-1 because of normalization)

! -------------------------------------------------------------------
    integer(wi) n
    real(wp) am1, a, ap1                            ! Left-hand side
    real(wp) bm2, bm1, b, bp1, bp2, bp3             ! Right-hand side
    real(wp) tmp1                                   ! Intermediate ops

! #######################################################################
! n = 1
! #######################################################################
    n = 1

! left-hand side
    am1 = 0.0_wp
    a = 1.0_wp
    ap1 = (x(n + 1) - x(n))*(x(n) - x(n + 2)) + (x(n + 1) - x(n))*(x(n) - x(n + 3)) - (x(n) - x(n + 2))*(x(n) - x(n + 3))
    ap1 = ap1/ &
   ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) + (x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

! right-hand side
    bm2 = 0.0_wp
    bm1 = 0.0_wp
    b = (x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + 2.0_wp*(x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
    b = b/ &
        ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
    b = b*((x(n) - x(n + 2)) + (x(n) - x(n + 3)))/ &
        (x(n) - x(n + 1)) + 1.0_wp
    b = b/ &
        ((x(n) - x(n + 2))*(x(n) - x(n + 3)))
    tmp1 = ((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))/ &
           ((x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)) + (x(n + 1) - x(n))*((x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3))))
    b = (b + tmp1/(x(n + 1) - x(n)))*2.0_wp

    bp1 = (x(n) - x(n + 2)) + (x(n) - x(n + 3)) + ap1*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2)) + (x(n + 1) - x(n + 3)))
    bp1 = bp1*2.0_wp/ &
          ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2))*(x(n + 1) - x(n + 3)))

    bp2 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 3))) + &
          (x(n) - x(n + 3))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 3)))
    bp2 = bp2/ &
          ((x(n + 1) - x(n))*(x(n + 1) - x(n + 3)) - (x(n + 2) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 3))))
    bp2 = bp2*2.0_wp*(x(n + 1) - x(n))/ &
          ((x(n + 2) - x(n + 1))*(x(n + 2) - x(n))*(x(n + 2) - x(n + 3)))

    bp3 = (x(n) - x(n + 1))*((x(n + 1) - x(n)) + 2.0_wp*(x(n + 1) - x(n + 2))) + &
          (x(n) - x(n + 2))*(2.0_wp*(x(n + 1) - x(n)) + 3.0_wp*(x(n + 1) - x(n + 2)))
    bp3 = bp3/ &
          ((x(n + 1) - x(n))*(x(n + 1) - x(n + 2)) - (x(n + 3) - x(n + 1))*((x(n + 1) - x(n)) + (x(n + 1) - x(n + 2))))
    bp3 = bp3*2.0_wp*(x(n + 1) - x(n))/ &
          ((x(n + 3) - x(n + 1))*(x(n + 3) - x(n))*(x(n + 3) - x(n + 2)))

! if uniform, we should have ( 0 1 11 ) and ( 13 -27 15 -1 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
    tmp1 = 1.0_wp/b

    lhs(n, 1) = 0.0_wp
    lhs(n, 2) = a*tmp1
    lhs(n, 3) = ap1*tmp1

    rhs(n, 1) = bp3*tmp1 ! saving here the last term that goes into the 3. superdiagonal
    rhs(n, 2) = 0.0_wp
    rhs(n, 3) = bp1*tmp1
    rhs(n, 4) = bp2*tmp1

! #######################################################################
! n = 2
! #######################################################################
    do n = 2, nmax - 1
        bm2 = 0.0_wp; bp2 = 0.0_wp

        tmp1 = 1.0_wp/ &
               ((x(n + 1) - x(n))*(x(n) - x(n - 1)) + (x(n + 1) - x(n - 1))*(x(n + 1) - x(n - 1)))

        bm1 = (x(n + 1) - x(n))/(x(n + 1) - x(n - 1))
        b = -1.0_wp
        bp1 = (x(n) - x(n - 1))/(x(n + 1) - x(n - 1))

        am1 = (x(n) - x(n - 1))*(x(n) - x(n - 1)) + (x(n + 1) - x(n))*(x(n) - x(n - 1)) - (x(n + 1) - x(n))*(x(n + 1) - x(n))
        am1 = am1*bm1*tmp1
        a = 1.0_wp
        ap1 = (x(n + 1) - x(n))*(x(n + 1) - x(n)) + (x(n + 1) - x(n))*(x(n) - x(n - 1)) - (x(n) - x(n - 1))*(x(n) - x(n - 1))
        ap1 = ap1*bp1*tmp1

        bm1 = bm1*12.0_wp*tmp1
        b = b*12.0_wp*tmp1
        bp1 = bp1*12.0_wp*tmp1

! if uniform, we should have ( 1/10 1 1/10 ) and ( 0 12/10 -12/5 12/10 0 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
        tmp1 = 1.0_wp/b

        lhs(n, 1) = am1*tmp1
        lhs(n, 2) = a*tmp1
        lhs(n, 3) = ap1*tmp1

        rhs(n, 1) = 0.0_wp
        rhs(n, 2) = bm1*tmp1
        rhs(n, 3) = bp1*tmp1
        rhs(n, 4) = 0.0_wp

    end do
! #######################################################################
! n = nmax; same as n = 1, but changing the signs of the increments w.r.t. n
! To understand it, e.g., define a new variable i = -j, where i is the
! discrete variable moving around n
! #######################################################################
    n = nmax

! left-hand side
    am1 = 0.0_wp
    a = 1.0_wp
    ap1 = (x(n - 1) - x(n))*(x(n) - x(n - 2)) + (x(n - 1) - x(n))*(x(n) - x(n - 3)) - (x(n) - x(n - 2))*(x(n) - x(n - 3))
    ap1 = ap1/ &
   ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) + (x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

! right-hand side
    bm2 = 0.0_wp
    bm1 = 0.0_wp
    b = (x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + 2.0_wp*(x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
    b = b/ &
        ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
    b = b*((x(n) - x(n - 2)) + (x(n) - x(n - 3)))/ &
        (x(n) - x(n - 1)) + 1.0_wp
    b = b/ &
        ((x(n) - x(n - 2))*(x(n) - x(n - 3)))
    tmp1 = ((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))/ &
           ((x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)) + (x(n - 1) - x(n))*((x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3))))
    b = (b + tmp1/(x(n - 1) - x(n)))*2.0_wp

    bp1 = (x(n) - x(n - 2)) + (x(n) - x(n - 3)) + ap1*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2)) + (x(n - 1) - x(n - 3)))
    bp1 = bp1*2.0_wp/ &
          ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2))*(x(n - 1) - x(n - 3)))

    bp2 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 3))) + &
          (x(n) - x(n - 3))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 3)))
    bp2 = bp2/ &
          ((x(n - 1) - x(n))*(x(n - 1) - x(n - 3)) - (x(n - 2) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 3))))
    bp2 = bp2*2.0_wp*(x(n - 1) - x(n))/ &
          ((x(n - 2) - x(n - 1))*(x(n - 2) - x(n))*(x(n - 2) - x(n - 3)))

    bp3 = (x(n) - x(n - 1))*((x(n - 1) - x(n)) + 2.0_wp*(x(n - 1) - x(n - 2))) + &
          (x(n) - x(n - 2))*(2.0_wp*(x(n - 1) - x(n)) + 3.0_wp*(x(n - 1) - x(n - 2)))
    bp3 = bp3/ &
          ((x(n - 1) - x(n))*(x(n - 1) - x(n - 2)) - (x(n - 3) - x(n - 1))*((x(n - 1) - x(n)) + (x(n - 1) - x(n - 2))))
    bp3 = bp3*2.0_wp*(x(n - 1) - x(n))/ &
          ((x(n - 3) - x(n - 1))*(x(n - 3) - x(n))*(x(n - 3) - x(n - 2)))

! if uniform, we should have ( 11 1 0 ) and ( -1 15 -27 13 )/h^2
! normalize s.t. b = 1 and we only need to store 4 RHS diagonals
! the  RHS diagonals are then O(1), the LHS diagonal O(h^2)
    tmp1 = 1.0_wp/b

    lhs(n, 3) = 0.0_wp
    lhs(n, 2) = a*tmp1
    lhs(n, 1) = ap1*tmp1

    rhs(n, 4) = bp3*tmp1 ! saving here the last term that goes into the 3. subdiagonal
    rhs(n, 3) = 0.0_wp
    rhs(n, 2) = bp1*tmp1
    rhs(n, 1) = bp2*tmp1

    return
end subroutine FDM_C2N4ND_INITIALIZE
