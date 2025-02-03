#include "dns_const.h"
#include "dns_error.h"

module Rotation
    use TLab_Constants, only: wp, wi, efile, lfile, MAX_VARS, MAX_PARS, MAX_PROF
    use NavierStokes, only: rossby
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
    type(term_dt), public :: coriolis

    public :: Rotation_Initialize
    public :: Rotation_Coriolis

contains
    !########################################################################
    !########################################################################
    subroutine Rotation_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=512) sRes
        integer(wi) idummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Rotation'
        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/explicit/normalized>')
        call TLab_Write_ASCII(bakfile, '#Vector=<Fx,Fy,Fz>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'void', sRes)
        if (trim(adjustl(sRes)) == 'void') &
            call ScanFile_Char(bakfile, inifile, 'Main', 'TermCoriolis', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; coriolis%type = EQNS_NONE
        else if (trim(adjustl(sRes)) == 'explicit') then; coriolis%type = EQNS_EXPLICIT
        else if (trim(adjustl(sRes)) == 'normalized') then; coriolis%type = EQNS_COR_NORMALIZED
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Wrong TermCoriolis option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        coriolis%vector(:) = 0.0_wp; coriolis%active = .false.
        if (coriolis%type /= EQNS_NONE) then
            call ScanFile_Char(bakfile, inifile, 'Rotation', 'Vector', '0.0,1.0,0.0', sRes)
            idummy = 3
            call LIST_REAL(sRes, idummy, coriolis%vector)

            if (abs(coriolis%vector(1)) > 0.0_wp) then; coriolis%active(2) = .true.; coriolis%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Ox.'); end if
            if (abs(coriolis%vector(2)) > 0.0_wp) then; coriolis%active(3) = .true.; coriolis%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Oy.'); end if
            if (abs(coriolis%vector(3)) > 0.0_wp) then; coriolis%active(1) = .true.; coriolis%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Oz.'); end if

            if (rossby > 0.0_wp) then
                coriolis%vector(:) = coriolis%vector(:)/rossby ! adding the rossby number into the vector
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Rossby number must be nonzero if coriolis is retained.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            coriolis%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, 'Rotation', 'Parameters', '0.0,1.0', sRes)
            idummy = MAX_PROF
            call LIST_REAL(sRes, idummy, coriolis%parameters)

            if (coriolis%parameters(2) == 0.0_wp) then
                call TLab_Write_ASCII(lfile, __FILE__//'. Default normalized geostrophic velocity set to one.')
                coriolis%parameters(2) = 1.0_wp
            end if

        end if

        ! Consistency check
        if (coriolis%type == EQNS_COR_NORMALIZED) then
            if (coriolis%active(2)) then
                call TLab_Write_ASCII(efile, __FILE__//'. TermCoriolis option only allows for angular velocity along Oy.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if
        end if

        return
    end subroutine Rotation_Initialize

    !########################################################################
    !########################################################################
    subroutine Rotation_Coriolis(locProps, nx, ny, nz, u, r)
        use TLab_OpenMP
        type(term_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz, *)
        real(wp), intent(out) :: r(nx*ny*nz, *)

        ! -----------------------------------------------------------------------
        integer(wi) siz, srt, end    !  Variables for OpenMP Partitioning
        integer(wi) ii, field_sz
        real(wp) dummy, dtr3, dtr1, geo_u, geo_w
        field_sz = nx*ny*nz
        ! -----------------------------------------------------------------------
        ! Coriolis. Remember that coriolis%vector already contains the Rossby #.
        ! -----------------------------------------------------------------------
        select case (locProps%type)
        case (EQNS_EXPLICIT)
            do ii = 1, field_sz
                r(ii, 1) = r(ii, 1) + locProps%vector(3)*u(ii, 2) - locProps%vector(2)*u(ii, 3)
                r(ii, 2) = r(ii, 2) + locProps%vector(1)*u(ii, 3) - locProps%vector(3)*u(ii, 1)
                r(ii, 3) = r(ii, 3) + locProps%vector(2)*u(ii, 1) - locProps%vector(1)*u(ii, 2)
            end do

        case (EQNS_COR_NORMALIZED)
            geo_u = cos(locProps%parameters(1))*locProps%parameters(2)
            geo_w = -sin(locProps%parameters(1))*locProps%parameters(2)

!$omp parallel default( shared ) &
!$omp private( ij, dummy,srt,end,siz )
            call TLab_OMP_PARTITION(field_sz, srt, end, siz)

            dummy = locProps%vector(2)
            dtr3 = 0.0_wp; dtr1 = 0.0_wp
            do ii = srt, end
                r(ii, 1) = r(ii, 1) + dummy*(geo_w - u(ii, 3))
                r(ii, 3) = r(ii, 3) + dummy*(u(ii, 1) - geo_u)
            end do
!$omp end parallel
        end select

    end subroutine Rotation_Coriolis

end module Rotation
