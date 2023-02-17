subroutine RAND_PSD(nx, ny, nz, u)
    use TLAB_CONSTANTS, only: wp, wi, pi_wp
    use TLAB_VARS, only: isize_txc_dimz
    use TLAB_VARS, only: g
    use RAND_LOCAL
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_offset_i, ims_offset_k
#endif

    implicit none

    integer(wi) nx, ny, nz
    real(wp), dimension(isize_txc_dimz, nz), intent(INOUT) :: u

    ! -----------------------------------------------------------------------
    integer(wi) j, k, iglobal, jglobal, kglobal, ip
    real(wp) pow_dst, pow_org, phase
    real(wp) f, f0, f1, fi, fj, fk
    real(wp) RAN0

    ! #######################################################################
    f0 = spc_param(1); f1 = spc_param(4)

    do k = 1, nz
#ifdef USE_MPI
        kglobal = k + ims_offset_k
#else
        kglobal = k
#endif
        if (kglobal <= g(3)%size/2 + 1) then
            fk = real(kglobal - 1, wp)/g(3)%scale
        else
            fk = -real(g(3)%size + 1 - kglobal, wp)/g(3)%scale
        end if

        do j = 1, ny
            jglobal = j ! No MPI decomposition along Oy
            if (jglobal <= g(2)%size/2 + 1) then
                fj = real(jglobal - 1, wp)/g(2)%scale
            else
                fj = -real(g(2)%size + 1 - jglobal, wp)/g(2)%scale
            end if

            do i = 1, nx/2 + 1
#ifdef USE_MPI
                iglobal = i + ims_offset_i/2
#else
                iglobal = i
#endif
                if (iglobal <= g(1)%size/2 + 1) then
                    fi = real(iglobal - 1, wp)/g(1)%scale
                else
                    fi = -real(g(1)%size + 1 - iglobal, wp)/g(1)%scale
                end if

                f = SQRT(fi**2 + fj**2 + fk**2)

                ! -----------------------------------------------------------------------
                ! target psd
                ! -----------------------------------------------------------------------
                select case (ispectrum)
                case(1)
                    pow_dst = 1.0_wp

                case(2)
                    pow_dst = (f/f0)**4/(1.0_wp+12.0_wp/5.0_wp*(f/f0)**2)**(17.0_wp/6.0_wp)

                case(3)
                    pow_dst = f**4*EXP(-2.0_wp*(f/f0)**2)

                case(4)
                    pow_dst = f**2*EXP(-2.0_wp*f/f0)

                case(5)
                    pow_dst = f**4/16.0_wp*EXP(-2.0_wp*(f/f0)**2)

                case(6)
                    pow_dst = EXP(-0.5_wp*((f - f0)/f1)**2)/(f1*SQRT(2.0_wp*pi_wp))

                case default
                
                end select
                
                if ((f - spc_param(2))*(spc_param(3) - f) > 0.0_wp) pow_dst = 0.0_wp ! Clip

                if (f == 0.0_wp) then
                    pow_dst = 0.0_wp
                else
                    if (g(2)%size == 1 .or. g(3)%size == 1) then ! 2D spectrum
                        pow_dst = pow_dst/(pi_wp*f)
                    else
                        pow_dst = pow_dst/(2*pi_wp*f**2)
                    end if
                end if

                ! -----------------------------------------------------------------------
                ! phase and scaling of complex data
                ! -----------------------------------------------------------------------
                ip = (nx + 2)*(j - 1) + 2*i

                if (ipdf == 0) then
                    if (iglobal == 1 .or. iglobal == g(1)%size/2 + 1) then; phase = 0.0_wp
                    else; phase = (RAN0(seed) - 0.5_wp)*2.0_wp*pi_wp; end if
                    u(ip - 1, k) = SQRT(pow_dst)*COS(phase)
                    u(ip, k) = SQRT(pow_dst)*SIN(phase)

                else
                    pow_org = u(ip - 1, k)**2 + u(ip, k)**2

                    if (pow_org > 0.0_wp) pow_dst = SQRT(pow_dst/pow_org)

                    u(ip - 1, k) = u(ip - 1, k)*pow_dst
                    u(ip, k) = u(ip, k)*pow_dst

                end if

            end do
        end do
    end do

    return
end subroutine RAND_PSD
