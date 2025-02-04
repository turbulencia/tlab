#include "dns_error.h"
#include "dns_const.h"
!########################################################################
!# HISTORY / AUTHORS
!#
!# 2022/04/01 - J. Kostelecky
!#              Created
!# 2023/10/27 - S. Deshpande
!#              Modified
!########################################################################
!# DESCRIPTION OF SUBROUTINES
!#
!#   called once in dns_main.f90 before time integration starts,
!#   all relevant geometry informations needed for the IBM are
!#   available after calling this routine
!#
!#   options for geometry generation:
!#      0. use existing eps field (ibm_restart==.true.)
!#      1. generate eps field outside tlab and use as existing eps field with
!#         (ibm_restart==.true.)
!#      2. write own routine to generate geometry
!#         (cf. IBM_GENERATE_GEOMETRY_XBARS, ibm_restart==.false.)
!#
!########################################################################
!# ARGUMENTS
!#
!#
!########################################################################
!# REQUIREMENTS
!#
!#
!########################################################################

subroutine IBM_INITIALIZE_GEOMETRY(txc, wrk3d)
    use TLab_Constants, only: efile, wp
    use IBM_VARS
    use TLab_Memory, only: isize_field, inb_txc
    use TLab_WorkFlow, only: stagger_on
    use FDM, only: g
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop

    implicit none

    real(wp), dimension(isize_field, inb_txc), intent(inout) :: txc
    real(wp), dimension(isize_field), intent(inout) :: wrk3d

    target :: txc

    real(wp), dimension(:), pointer :: epsi, epsj, epsk
    real(wp), dimension(:), pointer :: tmp1, tmp2
#ifdef IBM_DEBUG
    real(wp), dimension(:), pointer :: tmp3
#endif
    logical :: flag_epsp

    ! ================================================================== !
    ! assigning pointer to scratch
    txc = 0.0_wp; epsi => txc(:, 1); epsj => txc(:, 2); epsk => txc(:, 3)
    tmp1 => txc(:, 4); tmp2 => txc(:, 5)

    ! eps field (read/create)
    if (ibm_restart) then
        flag_epsp = .false.
        call IBM_IO_READ(wrk3d, flag_epsp)
    else
        if (ibm_geo%name == 'xbars') then
            call IBM_GENERATE_GEOMETRY_XBARS(wrk3d)
        else if (ibm_geo%name == 'hill') then
            call IBM_GENERATE_GEOMETRY_HILL(wrk3d)
        else if (ibm_geo%name == 'valley') then
            call IBM_GENERATE_GEOMETRY_VALLEY(wrk3d)
        else if (ibm_geo%name == 'box') then
            call IBM_GENERATE_GEOMETRY_BOX(wrk3d)
        else
            call TLab_Write_ASCII(efile, 'IBM_GEOMETRY no objects in flow.')
            call TLab_Stop(DNS_ERROR_IBM_MISS_GEO)
        end if
    end if

    ! transpose eps (epsi, epsj, epsk)
    call IBM_GEOMETRY_TRANSPOSE(epsi, epsj, epsk, wrk3d)

    ! generate relevant geometry fields for IBM routines (nobi, nobj, nobk)
    call IBM_GENERATE_GEOMETRY(epsi, epsj, epsk)

    ! verify geometry
    call IBM_VERIFY_GEOMETRY()

    ! epsp field (read/create)
    if (stagger_on) then
        if (ibm_restart) then
            flag_epsp = .true.
            call IBM_IO_READ(wrk3d, flag_epsp)
        else
            if (ibm_geo%name == 'xbars') then; continue
            else if (ibm_geo%name == 'hill') then; continue
            else if (ibm_geo%name == 'valley') then; continue
            else if (ibm_geo%name == 'box') then; continue
            else
                call TLab_Write_ASCII(efile, 'IBM_GEOMETRY epsp field is missing.')
                call TLab_Stop(DNS_ERROR_IBM_MISS_GEO)
            end if
        end if
    end if

    ! genereate ibm_case_xyz - fields
    call IBM_INITIALIZE_CASES(g(1), isize_nobi, isize_nobi_be, nobi, nobi_b, nobi_e, ibm_case_x)
    call IBM_INITIALIZE_CASES(g(2), isize_nobj, isize_nobj_be, nobj, nobj_b, nobj_e, ibm_case_y)
    call IBM_INITIALIZE_CASES(g(3), isize_nobk, isize_nobk_be, nobk, nobk_b, nobk_e, ibm_case_z)

    ! compute gamma_0/1 based on eps-field (volume approach for conditional averages!)
    call IBM_AVG_GAMMA(gamma_0, gamma_1, eps, tmp1)

    ! check idle procs
#ifdef USE_MPI
    call IBM_CHECK_PROCS(epsi, epsj, epsk)
#else
    ! in case of serial mode: one task with full domain, no idle procs
    ims_pro_ibm_x = .true.; ims_pro_ibm_y = .true.; ims_pro_ibm_z = .true.
#endif

#ifdef IBM_DEBUG
    ! io of all geometry fields in debugging mode
    tmp3 => txc(:, 6); tmp3(:) = 0.0_wp
    call IBM_GEOMETRY_DEBUG_IO(epsi, epsj, epsk, tmp1, tmp2, tmp3)
    nullify (tmp3)
#endif

    ! switch to true in routines where the IBM is needed
    ibm_burgers = .false.; ibm_partial = .false.

    ! disassociate pointers
    nullify (epsi, epsj, epsk)
    nullify (tmp1, tmp2)
    wrk3d(:) = 0.0_wp

    return
end subroutine IBM_INITIALIZE_GEOMETRY

!########################################################################

subroutine IBM_IO_READ(wrk3d, flag_epsp)

    use IBM_VARS
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Time, only: itime
    use IO_FIELDS, only: io_fileformat, IO_READ_FIELDS
    use TLab_Constants, only: wp, wi

    implicit none

    real(wp), dimension(isize_field), intent(inout) :: wrk3d
    logical, intent(inout) :: flag_epsp

    character(len=32) :: name
    real(wp) params(0)

    ! ================================================================== !
    wrk3d(:) = 0.0_wp
    select case (io_fileformat)
    case (IO_NOFILE) ! no IO
        if (flag_epsp) then
            epsp = 0.0_wp
        else
            eps = 0.0_wp
        end if
    case (IO_NETCDF) ! not implemented
    case default       ! mpiio - file with header
        select case (ibm_io)
        case (IBM_IO_REAL)
            if (flag_epsp) then
                name = epsp_name_real
                call IO_READ_FIELDS(name, imax, jmax, kmax, itime, 1, 0, epsp, params)
            else
                name = eps_name_real
                call IO_READ_FIELDS(name, imax, jmax, kmax, itime, 1, 0, eps, params)
            end if
        case (IBM_IO_INT)
            call IBM_IO_READ_INT_GEOMETRY(wrk3d, flag_epsp)
        case (IBM_IO_BIT)
            call IBM_IO_READ_BIT_GEOMETRY(wrk3d, flag_epsp)
        end select
    end select

    return
end subroutine IBM_IO_READ

!########################################################################

subroutine IBM_IO_WRITE(wrk3d, flag_epsp)

    use IBM_VARS
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Time, only: itime
    use IO_FIELDS
    use TLab_Constants, only: wp, wi

    implicit none

    real(wp), dimension(isize_field), intent(inout) :: wrk3d
    logical, intent(inout) :: flag_epsp

    character(len=32) :: name

    ! ================================================================== !
    wrk3d(:) = 0.0_wp
    select case (io_fileformat)
    case (IO_NOFILE) ! no IO
    case (IO_NETCDF) ! not implemented
    case default       ! mpiio - file with header
        select case (ibm_io)
        case (IBM_IO_REAL)
            if (flag_epsp) then
                name = epsp_name_real
                call IO_WRITE_FIELDS(name, imax, jmax, kmax, itime, 1, epsp, io_header_q)
            else
                name = eps_name_real
                call IO_WRITE_FIELDS(name, imax, jmax, kmax, itime, 1, eps, io_header_q)
            end if
        case (IBM_IO_INT)
            call IBM_IO_WRITE_INT_GEOMETRY(wrk3d, flag_epsp)
        case (IBM_IO_BIT)
            call IBM_IO_WRITE_BIT_GEOMETRY(wrk3d, flag_epsp)
        end select
    end select

    return
end subroutine IBM_IO_WRITE
