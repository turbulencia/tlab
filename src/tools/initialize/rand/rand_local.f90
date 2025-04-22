#include "dns_error.h"

module RAND_LOCAL
    use TLab_Constants, only: wp, wi, efile, lfile
    use TLab_Memory, only: imax, jmax, kmax, isize_field, isize_txc_field, inb_txc
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_pro
#endif
    use FDM, only: g
    use Averages, only: AVG1V2D
    use Distributions
    use OPR_Fourier
    implicit none
    save

    type(distributions_dt) :: psd
    type(distributions_dt) :: pdf
    real(wp) :: ucov(6)
    integer(wi) :: seed          ! Random number generator

    ! -------------------------------------------------------------------
    integer(wi) i

contains

    ! ###################################################################
    subroutine Inirand_Initialize_Parameters(inifile)
        character*(*) inifile

        ! -------------------------------------------------------------------
        character*512 sRes
        character*32 bakfile

        integer(wi) :: idummy
        real(wp) :: rdummy(6)

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        call TLab_Write_ASCII(lfile, 'Reading local input data')
        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#[Broadband]')
        call TLab_Write_ASCII(bakfile, '#Spectrum=<none/uniform/quartic/quadratic/gaussian>')
        call TLab_Write_ASCII(bakfile, '#f0=<frequencies>')
        call TLab_Write_ASCII(bakfile, '#Distribution=<none/uniform/gaussian>')
        call TLab_Write_ASCII(bakfile, '#Covariance=<Rxx,Ryy,Rzz,Rxy,Rxz,Ryz>')
        call TLab_Write_ASCII(bakfile, '#Seed=<random seed>')

        call ScanFile_Int(bakfile, inifile, 'Broadband', 'Seed', '7', seed)
#ifdef USE_MPI
        seed = seed + ims_pro         ! seed for random generator
#endif
        seed = -abs(seed)

        call ScanFile_Char(bakfile, inifile, 'Broadband', 'Spectrum', 'quartic', sRes)
        if (trim(adjustl(sRes)) == 'none') then; psd%type = TYPE_DF_NONE
        else if (trim(adjustl(sRes)) == 'uniform') then; psd%type = TYPE_DF_UNIFORM
        else if (trim(adjustl(sRes)) == 'quartic') then; psd%type = TYPE_DF_QUARTIC
        else if (trim(adjustl(sRes)) == 'quadratic') then; psd%type = TYPE_DF_QUADRATIC
        else if (trim(adjustl(sRes)) == 'gaussian') then; psd%type = TYPE_DF_GAUSSIAN
        end if

        psd%parameters(:) = 0.0_wp
        call ScanFile_Char(bakfile, inifile, 'Broadband', 'f0', '1.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, psd%parameters)
        psd%mean = psd%parameters(1)
        psd%parameters(1) = psd%parameters(2)
        psd%parameters(2) = psd%parameters(3)

        call ScanFile_Real(bakfile, inifile, 'Broadband', 'Sigma', '-1.0', psd%sigma)
        if (psd%sigma < 0.0_wp) psd%sigma = psd%mean/6.0_wp    ! Default value

        call ScanFile_Char(bakfile, inifile, 'Broadband', 'Distribution', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; pdf%type = TYPE_DF_NONE
        else if (trim(adjustl(sRes)) == 'uniform') then; pdf%type = TYPE_DF_UNIFORM
        else if (trim(adjustl(sRes)) == 'gaussian') then; pdf%type = TYPE_DF_GAUSSIAN
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Broadband: Distribution type unknown.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ucov(1:3) = 1.0_wp ! diagonal terms
        ucov(4:6) = 0.0_wp ! off-diagonal terms
        call ScanFile_Char(bakfile, inifile, 'Broadband', 'Covariance', '-1', sRes)
        if (trim(adjustl(sRes)) /= '-1') then
            idummy = 6
            call LIST_REAL(sRes, idummy, rdummy)
            if (idummy == 6) then; ucov(1:6) = rdummy(1:6)
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Broadband: Incorrect number of variances.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if
! CALL TLab_Write_ASCII(bakfile,'Covariance matrix:')
! WRITE(sRes,'(3E11.4)') ucov(1), ucov(4), ucov(5); CALL TLab_Write_ASCII(bakfile,sRes)
! WRITE(sRes,'(3E11.4)') ucov(4), ucov(2), ucov(6); CALL TLab_Write_ASCII(bakfile,sRes)
! WRITE(sRes,'(3E11.4)') ucov(5), ucov(6), ucov(3); CALL TLab_Write_ASCII(bakfile,sRes)

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        inb_txc = 3

        return
    end subroutine Inirand_Initialize_Parameters

    ! ###################################################################
    subroutine RAND_FIELD(variance, a, tmp1, tmp2, tmp3)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc

        real(wp), intent(IN) :: variance
        real(wp), dimension(isize_field), intent(OUT) :: a
        real(wp), dimension(isize_txc_field), intent(INOUT) :: tmp1, tmp2, tmp3

        ! -------------------------------------------------------------------
        integer(wi) idim
        real(wp) RAN0, RANG
        external RAN0, RANG
        complex(wp), pointer :: c_tmp1(:) => null()
        target tmp1

    ! ###################################################################
        select case (pdf%type)
        case (TYPE_DF_UNIFORM)
            do i = 1, isize_field
                tmp2(i) = RAN0(seed) - 0.5_wp
            end do

        case (TYPE_DF_GAUSSIAN)
            do i = 1, isize_field
                tmp2(i) = RANG(0.0_wp, 1.0_wp, seed)
            end do

        end select

        if (psd%type > 0) then
            if (g(2)%size == 1) then        ! 2D Fourier transform
                idim = 2
            else                            ! 3D Fourier transform
                idim = 3
            end if

            call c_f_pointer(c_loc(tmp1), c_tmp1, shape=[isize_txc_field/2])
            if (pdf%type > 0) then
                call OPR_Fourier_F(idim, imax, jmax, kmax, tmp2, tmp1, tmp3)
                call OPR_Fourier_SetPSD(imax, jmax, kmax, c_tmp1, psd)
            else
                do i = 1, isize_txc_field
                    tmp3(i) = RAN0(seed)
                end do
                call OPR_Fourier_SetPSD(imax, jmax, kmax, c_tmp1, psd, locPhase=tmp3)
            end if
            call OPR_Fourier_B(idim, imax, jmax, kmax, tmp1, tmp2)

        end if

        call RAND_NORMALIZE(variance, tmp2)
        a(1:isize_field) = tmp2(1:isize_field)

        return
    end subroutine RAND_FIELD

!########################################################################
    subroutine RAND_COVARIANCE(cov, u, v, w)
        real(wp) cov(6)
        real(wp), dimension(isize_field), intent(OUT) :: u, v, w

        ! -------------------------------------------------------------------
        real(wp) trace, lambda1, lambda2, alpha, calpha, salpha
        real(wp) rdummy

#define Rxx cov(1)
#define Ryy cov(2)
#define Rzz cov(3)
#define Rxy cov(4)
#define Rxz cov(5)
#define Ryz cov(6)

        ! ###################################################################
        if (g(3)%size > 1) then
            if (Rxz /= 0.0_wp .or. Ryz /= 0.0_wp) then ! only 2D case developed
                call TLab_Write_ASCII(efile, 'Terms Rxz and Ryz not developed yet.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            call RAND_NORMALIZE(Rzz, w)

        end if

        if (Rxy == 0.0_wp) then  ! Diagonal case
            call RAND_NORMALIZE(Rxx, u)
            call RAND_NORMALIZE(Ryy, v)

        else                        ! Nondiagonal case
            ! get eigenvalues
            trace = Rxx + Ryy
            lambda1 = 0.5_wp*(trace + sqrt(trace*trace - 4.0_wp*(Rxx*Ryy - Rxy*Rxy)))
            lambda2 = trace - lambda1

            ! define fields in rotated uncorrelated frame
            call RAND_NORMALIZE(lambda1, u)
            call RAND_NORMALIZE(lambda2, v)

            ! rotate to XY correlated frame
            alpha = atan((lambda1 - Rxx)/Rxy)
            calpha = cos(alpha)
            salpha = sin(alpha)

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
        real(wp), intent(IN) :: variance
        real(wp), dimension(imax, jmax, kmax), intent(INOUT) :: a

        ! -------------------------------------------------------------------
        real(wp) dummy

        ! ###################################################################
        dummy = AVG1V2D(imax*jmax, 1, kmax, 1, 1, a) ! 3D average
        a = a - dummy

        dummy = AVG1V2D(imax*jmax, 1, kmax, 1, 2, a) ! 3D average
        if (dummy > 0.0_wp) then
            dummy = sqrt(variance/dummy)
            a = a*dummy
        end if

        return
    end subroutine RAND_NORMALIZE

end module RAND_LOCAL
