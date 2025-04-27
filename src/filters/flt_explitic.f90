module Filters_Explicit
    use TLab_Constants, only: wp, wi
    use FLT_Base
    implicit none
    private

    integer(wi) i

    public :: FLT_E4, FLT_E4_COEFFS
    public :: FLT_E6
    public :: FLT_ADM

contains
! Explicit 4th order, nonuniform case.
! Coincides with filter as described in Stolz's thesis, imported from  Holger Foysi.
! Coefficient alpha_{-2} determined by setting minimum imaginary part, which implies alpha_m2 = alpha_p2
    subroutine FLT_E4(imax, jkmax, periodic, a, u, uf)
        integer(wi), intent(in) :: imax, jkmax
        logical, intent(in) :: periodic
        real(wp), dimension(jkmax, imax) :: u, uf
        real(wp), dimension(imax, 5) :: a

! boundary points
        if (periodic) then
            i = 1
            uf(:, i) = a(i, 1)*u(:, imax - 1) + a(i, 2)*u(:, imax) + a(i, 3)*u(:, i) &
                       + a(i, 4)*u(:, i + 1) + a(i, 5)*u(:, i + 2)

            i = 2
            uf(:, i) = a(i, 1)*u(:, imax) + a(i, 2)*u(:, i - 1) + a(i, 3)*u(:, i) &
                       + a(i, 4)*u(:, i + 1) + a(i, 5)*u(:, i + 2)

            i = imax - 1
            uf(:, i) = a(i, 1)*u(:, i - 2) + a(i, 2)*u(:, i - 1) + a(i, 3)*u(:, i) &
                       + a(i, 4)*u(:, i + 1) + a(i, 5)*u(:, 1)

            i = imax
            uf(:, i) = a(i, 1)*u(:, i - 2) + a(i, 2)*u(:, i - 1) + a(i, 3)*u(:, i) &
                       + a(i, 4)*u(:, 1) + a(i, 5)*u(:, 2)

        else
            i = 2
            uf(:, 1) = u(:, 1)
            uf(:, i) = a(i, 2)*u(:, i - 1) + a(i, 3)*u(:, i) &
                       + a(i, 4)*u(:, i + 1) + a(i, 5)*u(:, i + 2) + a(i, 1)*u(:, i + 3)

            i = imax - 1
            uf(:, i) = a(i, 1)*u(:, i - 2) + a(i, 2)*u(:, i - 1) + a(i, 3)*u(:, i) &
                       + a(i, 4)*u(:, i + 1) + a(i, 5)*u(:, i - 3)
            uf(:, imax) = u(:, imax)

        end if

! Interior points
        do i = 3, imax - 2
            uf(:, i) = a(i, 1)*u(:, i - 2) + a(i, 2)*u(:, i - 1) + a(i, 3)*u(:, i) &
                       + a(i, 4)*u(:, i + 1) + a(i, 5)*u(:, i + 2)
        end do

        return
    end subroutine FLT_E4

!########################################################################
    subroutine FLT_E4_COEFFS(kmax, periodic, scalez, z, cxi)
        integer(wi), intent(in) :: kmax
        logical, intent(in) :: periodic
        real(wp), intent(IN) :: scalez
        real(wp), dimension(*), intent(IN) :: z
        real(wp), dimension(kmax, 5), intent(out) :: cxi

        ! -----------------------------------------------------------------------
        integer(wi) k, km2, km1, kp1, kp2, kmin_loc, kmax_loc
        real(wp) D0, D1, D2, Dm1, zm2, zm1, zp1, zp2, zp3

#define alpha_m2(k)  cxi(k,1)
#define alpha_m1(k)  cxi(k,2)
#define alpha(k)     cxi(k,3)
#define alpha_p1(k)  cxi(k,4)
#define alpha_p2(k)  cxi(k,5)

        if (kmax <= 1) return

        ! #######################################################################
        ! Coefficients near wall
        ! Vanishing moment of third order
        ! #######################################################################
        if (.not. periodic) then

            ! -----------------------------------------------------------------------
            ! Point 2
            ! -----------------------------------------------------------------------
            k = 2

            zm1 = abs(z(k) - z(k - 1))
            zp1 = abs(z(k + 1) - z(k))
            zp2 = abs(z(k + 2) - z(k))
            zp3 = abs(z(k + 3) - z(k))

            D2 = zp2*(-zp1 + zp2 + zm1)
            D1 = -zp1**2 + zm1**2 + zp2*zp1 + zp2*zm1
            D0 = D2
            Dm1 = D1

            alpha_m2(k) = (zp2**3*zp1*zm1/(2.0_wp*D2) - zp1**3*(zm1**2 + zp2*zm1)/(2.0_wp*D1) + &
                           zm1**3*(zp1**2 - zp2*zp1)/(2.0_wp*Dm1))/ &
                          (zp3**3 - zp2**3*(-zp1*zp3 - zp1*zm1 + zp3**2 + zm1*zp3)/D2 + &
                           zp1**3*(-zm1**2 - zp2*zp3 - zp2*zm1 + 2.0_wp*zp3**2)/D1 - &
                           zm1**3*(-zp2*zp3 - zp1**2 + zp3**2 + zp2*zp1)/Dm1)

            alpha_m1(k) = -0.5_wp*(zp1**2 - zp2*zp1 + 2.0_wp*alpha_m2(k)*(-zp2*zp3 + zp3**2 - zp1**2 + zp2*zp1))/Dm1
            alpha(k) = 0.5_wp*(-zp2*zp1 + zp2**2 + zp2*zm1 + zp1*zm1 - 2.0_wp*alpha_m2(k)*(zp1*zp3 + zp1*zm1 - zp3**2 - zm1*zp3))/D0
            alpha_p1(k) = 0.5_wp*(zm1**2 + zp2*zm1 + 2.0_wp*alpha_m2(k)*(-zm1**2 - zp2*zp3 - zp2*zm1 + zp3**2))/D1
            alpha_p2(k) = -0.5_wp*(2.0_wp*alpha_m2(k)*(-zp1*zp3 - zp1*zm1 + zp3**2 + zm1*zp3) + zp1*zm1)/D2

            ! -----------------------------------------------------------------------
            ! Point N-1
            ! -----------------------------------------------------------------------
            k = kmax - 1

            zm1 = abs(z(k + 1) - z(k))
            zp1 = abs(z(k) - z(k - 1))
            zp2 = abs(z(k) - z(k - 2))
            zp3 = abs(z(k) - z(k - 3))

            D2 = zp2*(-zp1 + zp2 + zm1)
            D1 = -zp1**2 + zm1**2 + zp2*zp1 + zp2*zm1
            D0 = D2
            Dm1 = D1

            alpha_p2(k) = (zp2**3*zp1*zm1/(2.0_wp*D2) - zp1**3*(zm1**2 + zp2*zm1)/(2.0_wp*D1) + &
                           zm1**3*(zp1**2 - zp2*zp1)/(2.0_wp*Dm1))/ &
                          (zp3**3 - zp2**3*(-zp1*zp3 - zp1*zm1 + zp3**2 + zm1*zp3)/D2 + &
                           zp1**3*(-zm1**2 - zp2*zp3 - zp2*zm1 + 2.0_wp*zp3**2)/D1 - &
                           zm1**3*(-zp2*zp3 - zp1**2 + zp3**2 + zp2*zp1)/Dm1)

            alpha_p1(k) = -0.5_wp*(zp1**2 - zp2*zp1 + 2.0_wp*alpha_p2(k)*(-zp2*zp3 + zp3**2 - zp1**2 + zp2*zp1))/Dm1
            alpha(k) = 0.5_wp*(-zp2*zp1 + zp2**2 + zp2*zm1 + zp1*zm1 - 2.0_wp*alpha_p2(k)*(zp1*zp3 + zp1*zm1 - zp3**2 - zm1*zp3))/D0
            alpha_m1(k) = 0.5_wp*(zm1**2 + zp2*zm1 + 2.0_wp*alpha_p2(k)*(-zm1**2 - zp2*zp3 - zp2*zm1 + zp3**2))/D1
            alpha_m2(k) = -0.5_wp*(2.0_wp*alpha_p2(k)*(-zp1*zp3 - zp1*zm1 + zp3**2 + zm1*zp3) + zp1*zm1)/D2

            kmin_loc = 3
            kmax_loc = kmax - 2

        else
            kmin_loc = 1
            kmax_loc = kmax

        end if

        ! #######################################################################
        ! inner points
        ! #######################################################################
        do k = kmin_loc, kmax_loc
            km2 = k - 2; km2 = mod(km2 + kmax - 1, kmax) + 1
            km1 = k - 1; km1 = mod(km1 + kmax - 1, kmax) + 1
            kp1 = k + 1; kp1 = mod(kp1 + kmax - 1, kmax) + 1
            kp2 = k + 2; kp2 = mod(kp2 + kmax - 1, kmax) + 1

            zm2 = z(k) - z(km2); if (zm2 < 0.0_wp) zm2 = zm2 + scalez
            zm1 = z(k) - z(km1); if (zm1 < 0.0_wp) zm1 = zm1 + scalez
            zp1 = z(kp1) - z(k); if (zp1 < 0.0_wp) zp1 = zp1 + scalez
            zp2 = z(kp2) - z(k); if (zp2 < 0.0_wp) zp2 = zp2 + scalez

            D2 = zp2*(zp1 - zp2 - zm1) - (zp1*zm2 + zm2**2 - zm2*zm1)
            D1 = (zp1 - zp2 - zm1)*(zp1 + zm1)
            D0 = zp2*(zp1 - zp2 - zm1)

            alpha_p2(k) = 0.5_wp*(zp1*zm1)/D2
            alpha_m2(k) = alpha_p2(k)
            alpha_p1(k) = -0.5_wp*(2.0_wp*alpha_m2(k)*zm2*(zm2 + zp2) + zm1*(zp2 + zm1))/D1
            alpha_m1(k) = 0.5_wp*(zp1**2 - zp2*zp1 + 2.0_wp*alpha_m2(k)*(zm2**2 + zm2*zp2))/D1
            alpha(k) = 0.5_wp*(1.0_wp - (zp1*zm1 + 2.0_wp*alpha_m2(k)*(zp1*zm2 + zm2**2 - zm1*zm2))/D0) - alpha_m2(k)
        end do

        return
    end subroutine FLT_E4_COEFFS

    !########################################################################
    subroutine FLT_E6(kmax, ijmax, periodic, kflt1bc, kfltmxbc, u, uf)
        integer(wi) kflt1bc, kfltmxbc
        logical periodic
        integer(wi) ijmax, kmax
        real(wp), dimension(ijmax, kmax) :: u, uf

        ! -------------------------------------------------------------------
        integer(wi) k
        integer(wi) kstart, kstop
        real(wp) b0, b1, b2, b3
        real(wp) b_a(7), b_b(7), b_c(7)

        integer(wi) i1, i2, i3, i4, im3, ip3, im2, ip2, im1, ip1

        ! #######################################################################
#ifdef SINGLE_PREC
        b0 = 11.0e0/16.0e0
        b1 = 15.0e0/64.0e0
        b2 = -3.0e0/32.0e0
        b3 = 1.0e0/64.0e0

        b_a(1) = 63.0e0/64.0e0
        b_a(2) = 3.0e0/32.0e0
        b_a(3) = -15.0e0/64.0e0
        b_a(4) = 5.0e0/16.0e0
        b_a(5) = -15.0e0/64.0e0
        b_a(6) = 3.0e0/32.0e0
        b_a(7) = -1.0e0/64.0e0

        b_b(1) = 1.0e0/16.0e0
        b_b(2) = 3.0e0/4.0e0
        b_b(3) = 3.0e0/8.0e0
        b_b(4) = -1.0e0/4.0e0
        b_b(5) = 1.0e0/16.0e0
        b_b(6) = 0.0e0
        b_b(7) = 0.0e0

        b_c(1) = -1.0e0/32.0e0
        b_c(2) = 5.0e0/32.0e0
        b_c(3) = 11.0e0/16.0e0
        b_c(4) = 5.0e0/16.0e0
        b_c(5) = -5.0e0/32.0e0
        b_c(6) = 1.0e0/32.0e0
        b_c(7) = 0.0e0

#else
        b0 = 11.0d0/16.0d0
        b1 = 15.0d0/64.0d0
        b2 = -3.0d0/32.0d0
        b3 = 1.0d0/64.0d0

        b_a(1) = 63.0d0/64.0d0
        b_a(2) = 3.0d0/32.0d0
        b_a(3) = -15.0d0/64.0d0
        b_a(4) = 5.0d0/16.0d0
        b_a(5) = -15.0d0/64.0d0
        b_a(6) = 3.0d0/32.0d0
        b_a(7) = -1.0d0/64.0d0

        b_b(1) = 1.0d0/16.0d0
        b_b(2) = 3.0d0/4.0d0
        b_b(3) = 3.0d0/8.0d0
        b_b(4) = -1.0d0/4.0d0
        b_b(5) = 1.0d0/16.0d0
        b_b(6) = 0.0d0
        b_b(7) = 0.0d0

        b_c(1) = -1.0d0/32.0d0
        b_c(2) = 5.0d0/32.0d0
        b_c(3) = 11.0d0/16.0d0
        b_c(4) = 5.0d0/16.0d0
        b_c(5) = -5.0d0/32.0d0
        b_c(6) = 1.0d0/32.0d0
        b_c(7) = 0.0d0
#endif

        kstart = 1
        kstop = kmax

        if (.not. periodic) then

            k = 1
            uf(:, k) = u(:, k)

            if (kflt1bc == 1) then
                !       If the second plane is included do it separately.

                k = 2
                uf(:, k) = &
                    b_b(1)*u(:, k - 1) &
                    + b_b(2)*u(:, k) &
                    + b_b(3)*u(:, k + 1) &
                    + b_b(4)*u(:, k + 2) &
                    + b_b(5)*u(:, k + 3) &
                    + b_b(6)*u(:, k + 4) &
                    + b_b(7)*u(:, k + 5)

                !       If the third plane is included do it separately.

                k = 3
                uf(:, k) = &
                    b_c(1)*u(:, k - 2) &
                    + b_c(2)*u(:, k - 1) &
                    + b_c(3)*u(:, k) &
                    + b_c(4)*u(:, k + 1) &
                    + b_c(5)*u(:, k + 2) &
                    + b_c(6)*u(:, k + 3) &
                    + b_c(7)*u(:, k + 4)

            else

                do k = 2, 3
                    uf(:, k) = u(:, k)
                end do

            end if

            kstart = 4

            k = kmax
            uf(:, k) = u(:, k)

            if (kfltmxbc == 1) then
                !       If the second plane is included do it separately.

                k = kmax - 1
                uf(:, k) = &
                    b_b(1)*u(:, k + 1) &
                    + b_b(2)*u(:, k) &
                    + b_b(3)*u(:, k - 1) &
                    + b_b(4)*u(:, k - 2) &
                    + b_b(5)*u(:, k - 3) &
                    + b_b(6)*u(:, k - 4) &
                    + b_b(7)*u(:, k - 5)

                !       If the third plane is included do it separately.

                k = kmax - 2
                uf(:, k) = &
                    b_c(1)*u(:, k + 2) &
                    + b_c(2)*u(:, k + 1) &
                    + b_c(3)*u(:, k) &
                    + b_c(4)*u(:, k - 1) &
                    + b_c(5)*u(:, k - 2) &
                    + b_c(6)*u(:, k - 3) &
                    + b_c(7)*u(:, k - 4)
            else

                do k = kmax - 2, kmax - 1
                    uf(:, k) = u(:, k)
                end do

            end if

            kstop = kmax - 3

        end if

        i1 = 1
        i2 = 2
        i3 = 3
        i4 = 4

        do k = kstart, kstop
            im3 = mod(k - i4 + kmax, kmax) + 1
            im2 = mod(k - i3 + kmax, kmax) + 1
            im1 = mod(k - i2 + kmax, kmax) + 1
            ip1 = mod(k + kmax, kmax) + 1
            ip2 = mod(k + i1 + kmax, kmax) + 1
            ip3 = mod(k + i2 + kmax, kmax) + 1
            uf(:, k) = b3*(u(:, im3) + u(:, ip3)) &
                       + b2*(u(:, im2) + u(:, ip2)) &
                       + b1*(u(:, im1) + u(:, ip1)) &
                       + b0*u(:, k)
        end do

        return
    end subroutine FLT_E6

!########################################################################
!# Explicit filter after 3 levels of deconvolution
!#
!# uf = [ I + (I-G) + (I-G)*(I-G) ]*G*u = (3G - 3G^2 + G^3)*u
!#    = G*( 3u - 3G*u + G*G*u )
!#
    subroutine FLT_ADM(imax, jkmax, periodic, a, u, uf, tmp)
        integer(wi) imax, jkmax
        logical periodic
        real(wp), dimension(jkmax, imax) :: u, uf, tmp
        real(wp), dimension(imax, 5) :: a

        ! #######################################################################
        call FLT_E4(imax, jkmax, periodic, u, uf, a)
        call FLT_E4(imax, jkmax, periodic, uf, tmp, a)

        do i = 1, imax
            tmp(:, i) = tmp(:, i) + 3.0_wp*(u(:, i) - uf(:, i))
        end do

        call FLT_E4(imax, jkmax, periodic, tmp, uf, a)

        return
    end subroutine FLT_ADM

end module Filters_Explicit
