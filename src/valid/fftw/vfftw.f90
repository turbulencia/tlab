program VFFTW
    use TLab_Constants, only: wp, wi
    use TLAB_VARS
    use IO_Fields
    use OPR_Partial

    implicit none

#ifdef USE_FFTW
#include "fftw3.f"
#endif

    real(wp), dimension(:, :), allocatable, save, target :: x, y, z
    real(wp), dimension(:, :, :), pointer :: a, b, c
    TCOMPLEX, dimension(:, :, :), pointer :: a1, a2, a3
    real(wp), dimension(:), pointer :: wrk1d, wrk2d, wrk3d

    TCOMPLEX :: Img

    integer*4 i, j, k, ij
    integer(8) fft_plan_fx, fft_plan_fz
    integer(8) fft_plan_bx, fft_plan_bz

!  real(wp) fft_data_x, fft_data_z
    real(wp) dummy, error, params(0)

! ###################################################################
    call DNS_START

    call TLab_Initialize_Parameters(ifile)
    call NavierStokes_Initialize_Parameters(ifile)

! -------------------------------------------------------------------
! allocation of memory space
! -------------------------------------------------------------------
    allocate (wrk1d(isize_wrk1d*10))
    allocate (wrk2d(isize_wrk2d*5))
    allocate (wrk3d(imax*jmax*kmax), a(imax, jmax, kmax), b(imax, jmax, kmax), c(imax, jmax, kmax))
    allocate (a1(imax/2 + 1, jmax, kmax), a2(kmax, imax/2 + 1, jmax), a3(kmax, imax/2 + 1, jmax))

    Img = (-1.0, 0.0)
    Img = sqrt(Img)

    call TLab_Grid_Read(gfile, wrk1d(:, 1), wrk1d(:, 2), wrk1d(:, 3), [g(1)%size, g(2)%size, g(3)%size] )
    call FDM_CreatePlan(wrk1d(:, 1), g(1))
    call FDM_CreatePlan(wrk1d(:, 2), g(2))
    call FDM_CreatePlan(wrk1d(:, 3), g(3))

! ###################################################################
!  Define forcing term
! ###################################################################
!  DO k = 1,kmax
!     DO j = 1,jmax
!        DO i = 1,imax
!           a(i,j,k) = sin(C_2_R*C_PI_R/scalex*x(i)*C_2_R) &
!                    * cos(C_2_R*C_PI_R/scalez*z(k)*C_5_R) &
!                    * M_REAL(j-1)/M_REAL(jmax-1)*C_01_R
!!           c(i,j,k) = sin(C_2_R*C_PI_R/scalex*x(i)*C_2_R) &
!!                    * sin(C_2_R*C_PI_R/scalez*z(k)*C_5_R) * (-C_2_R*C_PI_R/scalez*C_5_R)&
!!                    * M_REAL(j-1)/M_REAL(jmax-1)*C_01_R
!
!        ENDDO
!     ENDDO
!  ENDDO
    call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, a, params)

!  CALL OPR_Partial_X(OPR_P1, imax,jmax,kmax, bcs, g(1), a, c)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), a, c)
    call IO_Write_Fields('field.ref', imax, jmax, kmax, itime, 1, c)

! ###################################################################
    call dfftw_plan_dft_r2c_1d &
        (fft_plan_fx, imax, wrk1d, wrk1d, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d &
        (fft_plan_fz, kmax_total, wrk1d, wrk1d, FFTW_FORWARD, FFTW_ESTIMATE)

    call dfftw_plan_dft_c2r_1d &
        (fft_plan_bx, imax, wrk1d, wrk1d, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d &
        (fft_plan_bz, kmax_total, wrk1d, wrk1d, FFTW_BACKWARD, FFTW_ESTIMATE)

! ###################################################################
! Fourier
! ###################################################################
! Forward Real FFT in x
    do k = 1, kmax
        do j = 1, jmax
            call dfftw_execute_dft_r2c(fft_plan_fx, a(1, j, k), a1(1, j, k))
        end do
    end do

! make j the last index, k the first
    do k = 1, kmax_total
        do ij = 1, (imax/2 + 1)*jmax
            a2(k, ij, 1) = a1(ij, 1, k)
        end do
    end do

! Forward complex FFT in z
    do j = 1, jmax
        do i = 1, imax/2 + 1
            call dfftw_execute_dft(fft_plan_fz, a2(1, i, j), a3(1, i, j))
        end do
    end do

! ###################################################################
! Manipulation
! ###################################################################
    do j = 1, jmax
        do i = 1, imax/2 + 1
            do k = 1, kmax_total
!           dummy = C_2_R*C_PI_R*M_REAL(i-1)/M_REAL(imax)
!           dummy = ( C_14_R/C_9_R*sin(dummy) + C_1_R/C_18_R*sin(2*dummy) ) &
!                 / ( C_1_R + C_2_R/C_3_R*cos(dummy) )
!           a3(k,i,j) = dummy*Img * a3(k,i,j) * M_REAL(imax)/scalex

                if (k <= kmax_total/2) then
                    dummy = C_2_R*C_PI_R*M_REAL(k - 1)/M_REAL(kmax_total)
                else
                    dummy = C_2_R*C_PI_R*M_REAL(k - 1 - kmax_total)/M_REAL(kmax_total)
                end if
                dummy = (C_14_R/C_9_R*sin(dummy) + C_1_R/C_18_R*sin(2*dummy)) &
                        /(C_1_R + C_2_R/C_3_R*cos(dummy))
                a3(k, i, j) = dummy*Img*a3(k, i, j)*M_REAL(kmax_total)/scalez

            end do
        end do
    end do

! ###################################################################
! Back
! ###################################################################
! backwards complex FFT in z
    do j = 1, jmax
        do i = 1, imax/2 + 1
            call dfftw_execute_dft(fft_plan_bz, a3(1, i, j), a2(1, i, j))
        end do
    end do

! make k the last index, kx the first
    do k = 1, kmax_total
        do ij = 1, (imax/2 + 1)*jmax
            a1(ij, 1, k) = a2(k, ij, 1)/M_REAL(imax*kmax_total) ! normalize
        end do
    end do

! backwards real FFT in z
    do k = 1, kmax
        do j = 1, jmax
            call dfftw_execute_dft_c2r(fft_plan_bx, a1(1, j, k), b(1, j, k))
        end do
    end do

    call IO_Write_Fields('field.out', imax, jmax, kmax, itime, 1, b)

! ###################################################################
! Error
! ###################################################################
    error = C_0_R
    dummy = C_0_R
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                error = error + (c(i, j, k) - b(i, j, k))*(c(i, j, k) - b(i, j, k))
                dummy = dummy + c(i, j, k)*c(i, j, k)
            end do
        end do
    end do
    write (*, *) 'Relative error ....: ', sqrt(error)/sqrt(dummy)

    stop
end program VFFTW
