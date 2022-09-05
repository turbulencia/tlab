#include "dns_error.h"

module RAND_LOCAL
    use TLAB_TYPES, only: cp, ci
    use TLAB_CONSTANTS, only: efile
    use TLAB_VARS, only: imax, jmax, kmax, isize_field, isize_txc_field
    use TLAB_VARS, only: g
    use TLAB_PROCS
    implicit none
    save

    ! -------------------------------------------------------------------
    integer(ci) :: ispectrum
    real(cp) :: spc_param(5)  ! Fundamental frequency, fmin, fmax, sigma

    integer(ci) :: ipdf
    real(cp) :: ucov(6)

    integer(ci) :: seed          ! Random number generator

    ! -------------------------------------------------------------------
    integer(ci) i

contains

    ! ###################################################################
    subroutine RAND_FIELD(variance, a, tmp1, tmp2, tmp3, wrk2d, wrk3d)
        real(cp), intent(IN) :: variance
        real(cp), dimension(isize_field), intent(OUT) :: a
        real(cp), dimension(isize_txc_field), intent(INOUT) :: tmp1, tmp2, tmp3
        real(cp), dimension(*), intent(INOUT) :: wrk2d, wrk3d

        integer(ci) idim
        real(cp) RAN0, RANG
        external RAN0, RANG

        ! -------------------------------------------------------------------
        select case (ipdf)
        case (1)     ! Uniform distribution
            do i = 1, isize_field
                tmp2(i) = RAN0(seed) - 0.5_cp
            end do

        case (2)     ! Gaussian distribution
            do i = 1, isize_field
                tmp2(i) = RANG(0.0_cp, 1.0_cp, seed)
            end do

        end select

        if (ispectrum > 0) then
            if (g(2)%size == 1) then; idim = 2; ! 2D Fourier transform
            else; idim = 3; end if     ! 3D Fourier transform

            if (ipdf > 0) call OPR_FOURIER_F(idim, imax, jmax, kmax, tmp2, tmp1, tmp3, wrk2d, wrk3d)
            call RAND_PSD(imax, jmax, kmax, tmp1)
            call OPR_FOURIER_B(idim, imax, jmax, kmax, tmp1, tmp2, wrk3d)

        end if

        call RAND_NORMALIZE(variance, tmp2)
        a(1:isize_field) = tmp2(1:isize_field)

        return
    end subroutine RAND_FIELD

    !########################################################################
    subroutine RAND_COVARIANCE(cov, u, v, w)
        real(cp) cov(6)
        real(cp), dimension(isize_field), intent(OUT) :: u, v, w

        ! -------------------------------------------------------------------
        real(cp) trace, lambda1, lambda2, alpha, calpha, salpha
        real(cp) rdummy

#define Rxx cov(1)
#define Ryy cov(2)
#define Rzz cov(3)
#define Rxy cov(4)
#define Rxz cov(5)
#define Ryz cov(6)

        ! ###################################################################
        if (g(3)%size > 1) then
            if (Rxz /= 0.0_cp .or. Ryz /= 0.0_cp) then ! only 2D case developed
                call TLAB_WRITE_ASCII(efile, 'Terms Rxz and Ryz not developed yet.')
                call TLAB_STOP(DNS_ERROR_UNDEVELOP)
            end if

            call RAND_NORMALIZE(Rzz, w)

        end if

        if (Rxy == 0.0_cp) then  ! Diagonal case
            call RAND_NORMALIZE(Rxx, u)
            call RAND_NORMALIZE(Ryy, v)

        else                        ! Nondiagonal case
            ! get eigenvalues
            trace = Rxx + Ryy
            lambda1 = 0.5_cp*(trace + SQRT(trace*trace - 4.0_cp*(Rxx*Ryy - Rxy*Rxy)))
            lambda2 = trace - lambda1

            ! define fields in rotated uncorrelated frame
            call RAND_NORMALIZE(lambda1, u)
            call RAND_NORMALIZE(lambda2, v)

            ! rotate to XY correlated frame
            alpha = ATAN((lambda1 - Rxx)/Rxy)
            calpha = COS(alpha)
            salpha = SIN(alpha)

            do i = 1, isize_field
                rdummy = calpha*u(i) - salpha*v(i)
                v(i) = salpha*u(i) + calpha*v(i)
                u(i) = rdummy
            end do

        end if

        return

#undef Rxx
#undef Ryy
#undef Rzz
#undef Rxy
#undef Rxz
#undef Ryz

    end subroutine RAND_COVARIANCE

    ! ###################################################################
    subroutine RAND_NORMALIZE(variance, a)
        real(cp), intent(IN) :: variance
        real(cp), dimension(imax, jmax, kmax), intent(INOUT) :: a

        ! -------------------------------------------------------------------
        real(cp) AVG1V2D, dummy
        external AVG1V2D

        ! ###################################################################
        dummy = AVG1V2D(imax*jmax, 1, kmax, 1, 1, a) ! 3D average
        a = a - dummy

        dummy = AVG1V2D(imax*jmax, 1, kmax, 1, 2, a) ! 3D average
        if (dummy > 0.0_cp) then
            dummy = SQRT(variance/dummy)
            a = a*dummy
        end if

        return
    end subroutine RAND_NORMALIZE

end module RAND_LOCAL
