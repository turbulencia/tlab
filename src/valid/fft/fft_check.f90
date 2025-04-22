#include "dns_error.h"
#ifdef USE_MPI

#endif

subroutine FFT_CHECK(check_mode, err_count, case_count, &
                     trans, &
                     trans_ref, &
                     tmp1, &
                     tmp2, &
                     tmp3, &
                     tmp4, &
                     wrk3d, &
                     wrk2d)

    use FDM, only: g
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_Constants, only: lfile, wp, wi
    use OPR_Fourier
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS
#endif

    implicit none

    integer(wi), intent(IN) :: check_mode
    integer(wi), intent(INOUT) :: err_count, case_count
    real(wp), dimension((imax + 2)*(jmax + 2)*kmax) :: tmp1, tmp2, tmp3, tmp4, wrk2d, wrk3d
    real(wp), dimension(imax*jmax*kmax) :: trans, trans_ref

    real(wp) :: scalex, scaley, dummy, field_variance
    real(wp), dimension(jmax) :: spec_variance, dummy_variance

    integer(wi) :: i, j, k, iv, ip, ip_ref, bad_count, good_count, bad, nxz
    integer(wi) :: good_planes, bad_planes
    integer(wi) :: isize_fft3d, isize_wrk3d, isize_trn3d
    integer(wi) :: isize_page_ij, isize_line_i
    real(wp) :: norm, residual, re, im, arg
    character :: fname*32, line*256, check_name*32, label*32

    isize_page_ij = imax*jmax
    isize_line_i = imax

    isize_wrk3d = imax*jmax*kmax
    isize_fft3d = (imax + 2)*(jmax + 2)*kmax
    isize_trn3d = (imax/2)*jmax*(2*kmax)

    call OPR_Fourier_Initialize()

    norm = C_1_R/M_REAL(g(1)%size*g(3)%size)

    tmp1(:) = C_0_R
    tmp3(:) = C_0_R
    tmp4(:) = C_0_R
    wrk2d(:) = C_0_R
    wrk3d(:) = C_0_R; tmp2(:) = C_0_R; trans(:) = C_0_R

    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                ip = (k - 1)*imax*jmax + (j - 1)*imax + i
#ifdef USE_MPI
                tmp1(ip) = setup_check(check_mode, i + ims_offset_i, j, k + ims_offset_k)
#else
                tmp1(ip) = setup_check(check_mode, i, j, k)
#endif
            end do
        end do
    end do

    select case (check_mode)
    case (1)
        label = 'fft_check:    delta'
        nxz = g(1)%size*g(3)%size
        dummy = M_REAL(nxz)/(M_REAL(nxz) - C_1_R)
        field_variance = C_2_R*(M_REAL(nxz) - C_1_R)/nxz/nxz
        field_variance = dummy*field_variance*M_REAL(g(2)%size)/C_2_R
    case (2)
        label = 'fft_check:   cosine'
        field_variance = C_05_R*g(2)%size
    case (3)
        label = 'fft_check:   random'
    case DEFAULT
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
    end select

#ifdef USE_MPI
    do i = 0, ims_npro - 1
        if (ims_pro == i) then
            do k = 1, kmax
                write (*, 1011) ims_pro, 'input', tmp1((k - 1)*imax*jmax + 1:k*imax*jmax)
1011            format(i3, a, 16(g10.3, 1x))
            end do
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    end do
#endif

    tmp2(:) = C_0_R
    trans(:) = C_0_R
    call OPR_Fourier_F(i2, imax, jmax, kmax, tmp1, tmp2, trans)

    tmp4(:) = tmp2(:)
    tmp3(:) = C_0_R
    call OPR_Fourier_B(i2, imax, jmax, kmax, tmp4, tmp3)

    ip = imax*jmax*kmax
    residual = maxval(abs(norm*tmp3(1:ip) - tmp1(1:ip)))
#ifdef USE_MPI
    dummy = residual
    call MPI_Reduce(dummy, residual, 1, MPI_REAL8, MPI_MAX, 0, MPI_COMM_WORLD, ims_err)

    do i = 0, ims_npro - 1
        if (i == ims_pro) then
            write (*, *) 'pro', ims_pro
            do k = 1, kmax
                write (*, 1018) norm*tmp3((k - 1)*imax*jmax + 1:k*imax*jmax)
1018            format(16(G10.3, 1x))
            end do
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    end do

    if (ims_pro == 0) then
#endif

        if (residual > 1.e-14) then
            write (line, 1000) trim(adjustl(label)), 'FAILED', residual
            err_count = err_count + 1
        else
            write (line, 1000) trim(adjustl(label)), 'PASSED', residual
        end if
1000    format(a, 1x, a6, ' Transform-check.    Max. Residual: ', G13.6)
        case_count = case_count + 1
        call TLab_Write_ASCII(lfile, line)

#ifdef USE_MPI
    end if
#endif

    if (check_mode == i3) goto 999

    trans(:) = C_0_R
    spec_variance(:) = C_0_R

!DO i=0,ims_npro-1
!   IF ( i .EQ. ims_pro ) THEN
!      WRITE(*,*) ims_pro, 'after poisson'
!      DO k=1,kmax
!         WRITE(*,*) tmp2((k-1)*2*(imax/2+1)*(jmax+2)+1 : k*2*(imax/2+1)*(jmax+2))
!      ENDDO
!   ENDIF
!   CALL MPI_Barrier(MPI_COMM_WORLD,ims_err)
!ENDDO

    call KXZ_PSD(imax, jmax, kmax, i1, norm, tmp2, trans, spec_variance)

#ifdef USE_MPI
    dummy_variance(:) = spec_variance(:)
    call MPI_Reduce(dummy_variance, spec_variance, jmax, MPI_REAL8, MPI_SUM, i0, MPI_COMM_WORLD, ims_err)
#endif

    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax/2
                ip = (k - 1)*imax/2*jmax + (j - 1)*imax/2 + i
#ifdef USE_MPI
                trans_ref(ip) = power_check(check_mode, i + ims_offset_i/2, j, k + ims_offset_k)
#else
                trans_ref(ip) = power_check(check_mode, i, j, k)
#endif
            end do
        end do
    end do

    ip = imax/2*jmax*kmax
    residual = maxval(abs(trans_ref(1:ip) - trans(1:ip)))

#ifdef USE_MPI
    dummy = residual
    call MPI_Reduce(dummy, residual, 1, MPI_REAL8, 0, MPI_MAX, MPI_COMM_WORLD, ims_err)

    if (ims_pro == 0) then
#endif
        if (maxval(abs(trans_ref(1:ip) - trans(1:ip))) > 1.e-16) then
            write (line, 1017) trim(adjustl(label)), 'FAILED', residual
            err_count = err_count + 1
        else
            write (line, 1017) trim(adjustl(label)), 'PASSED', residual
        end if
        case_count = case_count + 1

1017    format(a, 1x, a6, 1x, 'PSD check.    Max. Residual:', G13.8)
        call TLab_Write_ASCII(lfile, line)
#ifdef USE_MPI
    end if
#endif

#ifdef USE_MPI

    do i = 1, ims_npro
        if (ims_pro == i - 1) then
            do k = 1, kmax
!         WRITE(*,*) 'k=',k
!         WRITE(*,1010) ims_pro, trans(    (k-1)*imax/2*jmax+1 : k*imax/2*jmax)
!         WRITE(*,1010) ims_pro, trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax)
!         WRITE(*,1010) ims_pro, trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax) - trans((k-1)*imax/2*jmax+1 : k*imax/2*jmax)
            end do
        end if
        call MPI_Barrier(MPI_COMM_WORLD, ims_err)
    end do
#else
    do k = 1, kmax
!   WRITE(*,*) 'k=',k
!   WRITE(*,1010) trans(    (k-1)*imax/2*jmax+1 : k*imax/2*jmax)
!   WRITE(*,1010) trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax)
!   WRITE(*,1010) trans_ref((k-1)*imax/2*jmax+1 : k*imax/2*jmax) - trans((k-1)*imax/2*jmax+1 : k*imax/2*jmax)
    end do
#endif
1010 format(16(g12.3, 1x))

#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

        do j = 2, jmax
            spec_variance(1) = spec_variance(1) + spec_variance(j)
        end do

        if (abs(field_variance - spec_variance(1)) < 1.e-07) then
            write (line, 1001) trim(adjustl(label)), 'PASSED', &
                abs(field_variance - spec_variance(1)), field_variance, spec_variance(1)
        else
            write (line, 1001) trim(adjustl(label)), 'FAILED', &
                abs(field_variance - spec_variance(1)), field_variance, spec_variance(1)
            err_count = err_count + 1
        end if
        case_count = case_count + 1

        call TLab_Write_ASCII(lfile, line)
1001    format(a, 1x, a6, 1x, 'PARSEVAL-Identity check. Residual:', &
               1x, G13.8, ' Field:', G13.8, ' Spectrum: ', G13.8)

#ifdef USE_MPI
    end if
#endif

999 continue

contains

    real(wp) function SETUP_CHECK(check_mode, i, j, k)

        implicit none

        integer(wi), intent(IN) :: check_mode, i, j, k

        select case (check_mode)
        case (1) ! delta check
            if (i == i1) then
                if (k <= i2) then
                    SETUP_CHECK = C_1_R*mod(j, i2)
                else
                    SETUP_CHECK = C_0_R
                end if
            else
                SETUP_CHECK = C_0_R
            end if
        case (2) ! cosine check
            SETUP_CHECK = cos(C_2_R*C_PI_R*M_REAL(i - 1)/M_REAL(g(1)%size))
        case (3) ! random check
            call random_number(SETUP_CHECK)
        case default
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end select

    end function SETUP_CHECK

    real(wp) function POWER_CHECK(check_mode, i, j, k)

        integer(wi), intent(IN) :: check_mode, i, j, k

        real(wp) arg, re, im

        select case (check_mode)
        case (1)  ! delta-function check
            if (mod(j, i2) == 1) then
                arg = C_PI_R*C_2_R*M_REAL(k - 1)/M_REAL(g(3)%size)
                re = cos(arg) + cos(C_2_R*arg)
                im = sin(-arg) + sin(-C_2_R*arg)
                power_check = (re*re + im*im)/M_REAL(g(1)%size*g(1)%size*g(3)%size*g(3)%size)
            else
                power_check = C_0_R
            end if
        case (2)  ! cosine check
            if (i == i2 .and. k == i1) then
                power_check = C_025_R
            else
                power_check = C_0_R
            end if
        case (3)
            power_check = -C_1_R
        end select

    end function POWER_CHECK

end subroutine FFT_CHECK
