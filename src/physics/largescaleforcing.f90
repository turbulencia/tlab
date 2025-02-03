#include "dns_const.h"
#include "dns_error.h"

module LargeScaleForcing
    use TLab_Constants, only: wp, wi, efile, MAX_PROF, MAX_VARS, MAX_PARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    type term_dt
        sequence
        integer type
        integer scalar(MAX_VARS)                ! fields defining this term
        logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) auxiliar(MAX_PARS)
        real(wp) vector(3)
    end type term_dt
    type(term_dt), public, protected :: subsidenceProps

    public :: LargeScaleForcing_Initialize
    public :: LargeScaleForcing_Subsidence

    integer, parameter :: TYPE_NONE = 0
    integer, parameter, public :: TYPE_SUB_CONSTANT_LOCAL = 1
    integer, parameter :: TYPE_SUB_CONSTANT_GLOBAL = 2

contains
    !########################################################################
    !########################################################################
    subroutine LargeScaleForcing_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes

        integer idummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Subsidence'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/ConstantDivergenceLocal/ConstantDivergenceGlobal>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'TermSubsidence', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; subsidenceProps%type = TYPE_NONE
        else if (trim(adjustl(sRes)) == 'constantdivergencelocal') then; subsidenceProps%type = TYPE_SUB_CONSTANT_LOCAL
        else if (trim(adjustl(sRes)) == 'constantdivergenceglobal') then; subsidenceProps%type = TYPE_SUB_CONSTANT_GLOBAL
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermSubsidence option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        subsidenceProps%active = .false.
        if (subsidenceProps%type /= TYPE_NONE) then
            subsidenceProps%active = .true.

            subsidenceProps%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', '0.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, subsidenceProps%parameters)
        end if

        ! This subsidence type is implemented in OPR_Burgers_y only
        ! to speed up calculation
        if (subsidenceProps%type == TYPE_SUB_CONSTANT_LOCAL) subsidenceProps%active = .false.

        return
    end subroutine LargeScaleForcing_Initialize

    !########################################################################
    !########################################################################
    subroutine LargeScaleForcing_Subsidence(locProps, nx, ny, nz, a, source)
        use TLab_Arrays, only: wrk1d
        use OPR_PARTIAL, only: OPR_PARTIAL_Y
        use Averages, only: AVG1V2D_V
        use FDM, only: g

        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: a(nx, ny, nz)
        real(wp), intent(out) :: source(nx, ny, nz)

        ! -----------------------------------------------------------------------
        integer(wi) bcs(2, 2), j

        !########################################################################
        bcs = 0

        select case (locProps%type)

        case (TYPE_SUB_CONSTANT_LOCAL)
            wrk1d(1:ny, 1) = g(2)%nodes(1:ny)*locProps%parameters(1)      ! subsidence velocity

            call OPR_PARTIAL_Y(OPR_P1, nx, ny, nz, bcs, g(2), a, source)

            do j = 1, ny
                source(:, j, :) = source(:, j, :)*wrk1d(j, 1)
            end do

        case (TYPE_SUB_CONSTANT_GLOBAL)
            wrk1d(1:ny, 1) = g(2)%nodes(1:ny)*locProps%parameters(1)      ! subsidence velocity

            call AVG1V2D_V(nx, ny, nz, 1, a, wrk1d(:, 2), wrk1d(1, 3))     ! Calculate averaged scalar into wrk1d(:,2)
            call OPR_PARTIAL_Y(OPR_P1, 1, ny, 1, bcs, g(2), wrk1d(1, 2), wrk1d(1, 3))

            do j = 1, ny
                source(:, j, :) = wrk1d(j, 3)*wrk1d(j, 1)
            end do

        end select

        return
    end subroutine LargeScaleForcing_Subsidence

end module LargeScaleForcing
