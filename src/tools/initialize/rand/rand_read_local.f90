#include "dns_error.h"

subroutine RAND_READ_LOCAL(inifile)
    use TLab_Constants, only: efile, lfile
    use TLab_WorkFlow
    use RAND_LOCAL

    implicit none

    character*(*) inifile

! -------------------------------------------------------------------
    character*512 sRes
    character*32 bakfile

    integer(wi) :: idummy
    real(wp) :: rdummy(6)

! ###################################################################
    bakfile = TRIM(ADJUSTL(inifile))//'.bak'

    call TLab_Write_ASCII(lfile, 'Reading local input data')
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Broadband]')
    call TLab_Write_ASCII(bakfile, '#Spectrum=<none/uniform/quartic/quadratic/gaussian>')
    call TLab_Write_ASCII(bakfile, '#f0=<frequencies>')
    call TLab_Write_ASCII(bakfile, '#Distribution=<none/uniform/gaussian>')
    call TLab_Write_ASCII(bakfile, '#Covariance=<Rxx,Ryy,Rzz,Rxy,Rxz,Ryz>')
    call TLab_Write_ASCII(bakfile, '#Seed=<random seed>')

    call SCANINIINT(bakfile, inifile, 'Broadband', 'Seed', '7', seed)

    call SCANINICHAR(bakfile, inifile, 'Broadband', 'Spectrum', 'quartic', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; ispectrum = 0
    else if (TRIM(ADJUSTL(sRes)) == 'uniform') then; ispectrum = 1
    else if (TRIM(ADJUSTL(sRes)) == 'quartic') then; ispectrum = 3
    else if (TRIM(ADJUSTL(sRes)) == 'quadratic') then; ispectrum = 4
    else if (TRIM(ADJUSTL(sRes)) == 'gaussian') then; ispectrum = 6; end if

    spc_param = 0.0_wp
    call SCANINICHAR(bakfile, inifile, 'Broadband', 'f0', '1.0', sRes)
    idummy = 3
    call LIST_REAL(sRes, idummy, spc_param)

    call SCANINIREAL(bakfile, inifile, 'Broadband', 'Sigma', '-1.0', spc_param(4))
    if (spc_param(4) < 0.0_wp) spc_param(4) = spc_param(1)/6.0_wp    ! Default value

    call SCANINICHAR(bakfile, inifile, 'Broadband', 'Distribution', 'none', sRes)
    if (TRIM(ADJUSTL(sRes)) == 'none') then; ipdf = 0
    else if (TRIM(ADJUSTL(sRes)) == 'uniform') then; ipdf = 1
    else if (TRIM(ADJUSTL(sRes)) == 'gaussian') then; ipdf = 2
    else
        call TLab_Write_ASCII(efile, 'RAND_READ_LOCAL. Broadband: Distribution type unknown.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    ucov(1:3) = 1.0_wp ! diagonal terms
    ucov(4:6) = 0.0_wp ! off-diagonal terms
    call SCANINICHAR(bakfile, inifile, 'Broadband', 'Covariance', '-1', sRes)
    if (TRIM(ADJUSTL(sRes)) /= '-1') then
        idummy = 6
        call LIST_REAL(sRes, idummy, rdummy)
        if (idummy == 6) then; ucov(1:6) = rdummy(1:6)
        else
            call TLab_Write_ASCII(efile, 'RAND_READ_LOCAL. Broadband: Incorrect number of variances.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if
    end if
    ! CALL TLab_Write_ASCII(bakfile,'Covariance matrix:')
    ! WRITE(sRes,'(3E11.4)') ucov(1), ucov(4), ucov(5); CALL TLab_Write_ASCII(bakfile,sRes)
    ! WRITE(sRes,'(3E11.4)') ucov(4), ucov(2), ucov(6); CALL TLab_Write_ASCII(bakfile,sRes)
    ! WRITE(sRes,'(3E11.4)') ucov(5), ucov(6), ucov(3); CALL TLab_Write_ASCII(bakfile,sRes)

    return
end subroutine RAND_READ_LOCAL
