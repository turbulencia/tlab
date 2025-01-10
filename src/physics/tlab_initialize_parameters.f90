#include "dns_error.h"
#include "dns_const.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "TLab_Initialize_Parameters"

!########################################################################
!# Reading general data from file tlab.ini, setting up general parameters
!########################################################################
subroutine TLab_Initialize_Parameters(inifile)
    use TLab_Constants, only: wp, wi, lfile, efile, wfile, MajorVersion, MinorVersion
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, imode_verbosity
    use TLAB_VARS, only: imode_sim
    use TLAB_VARS, only: flow_on, scal_on, fourier_on, stagger_on
    use TLAB_VARS, only: imax, jmax, kmax, isize_field
    use TLAB_VARS, only: isize_wrk1d, isize_wrk2d, isize_wrk3d
    use TLAB_VARS, only: isize_txc_field, isize_txc_dimx, isize_txc_dimz
    use FDM, only: g
    use IO_FIELDS, only: imode_files, imode_precision_files
    ! use Avg_Spatial
#ifdef USE_MPI
    use TLabMPI_VARS
#endif
    use OPR_Filters, only: FilterDomain, FilterDomainActive, FilterDomainBcsFlow, FilterDomainBcsScal, PressureFilter
    use OPR_Filters, only: FILTER_READBLOCK
    implicit none

    character(len=*), intent(in) :: inifile

! -------------------------------------------------------------------
    character*512 sRes
    character*32 bakfile
    integer(wi) ig, idummy
#ifdef USE_MPI
    character*64 lstr
#endif

! ###################################################################
    bakfile = trim(adjustl(inifile))//'.bak'

    call TLab_Write_ASCII(lfile, 'Reading global input data.')

! ###################################################################
! Version Checking
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#[Version]')
    call TLab_Write_ASCII(bakfile, '#Major=<mayor version number>')
    call TLab_Write_ASCII(bakfile, '#Minor=<minor version number>')

    call ScanFile_Int(bakfile, inifile, 'Version', 'Major', '0', idummy)
    if (MajorVersion /= idummy) then
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Major version error.')
        call TLab_Stop(DNS_ERROR_VERSION)
    end if
    call ScanFile_Int(bakfile, inifile, 'Version', 'Minor', '0', idummy)
    if (MinorVersion /= idummy) then
        write (sRes, '(I5)') MinorVersion
        call TLab_Write_ASCII(wfile, 'DNS_REAL_GLOBAL. Minor version warning. Expected : '//sRes)
    end if

! ###################################################################
! Global information
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Main]')
    call TLab_Write_ASCII(bakfile, '#FileFormat=<mpiio/NetCDF/None>')
    call TLab_Write_ASCII(bakfile, '#FileType=<Double/Single>')
    call TLab_Write_ASCII(bakfile, '#VerbosityLevel=<0/1/2>')
    call TLab_Write_ASCII(bakfile, '#Type=<temporal/spatial>')
    call TLab_Write_ASCII(bakfile, '#CalculateFlow=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#CalculateScalar=<yes/no>')
    !
    call TLab_Write_ASCII(bakfile, '#SpaceOrder=<CompactJacobian4/CompactJacobian6/CompactJacobian6Penta/CompactDirect6>')
    !
    call TLab_Write_ASCII(bakfile, '#ComModeITranspose=<none,asynchronous,sendrecv>')
    call TLab_Write_ASCII(bakfile, '#ComModeKTranspose=<none,asynchronous,sendrecv>')

! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, inifile, 'Main', 'FileFormat', 'MpiIO', sRes)
    if (trim(adjustl(sRes)) == 'mpiio') then; imode_files = IO_MPIIO
    elseif (trim(adjustl(sRes)) == 'netcdf') then; imode_files = IO_NETCDF
    elseif (trim(adjustl(sRes)) == 'none') then; imode_files = IO_NOFILE
    else
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Wrong Main.FileFormat.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'FileType', 'Double', sRes)
    if (trim(adjustl(sRes)) == 'double') then; imode_precision_files = IO_TYPE_DOUBLE
    elseif (trim(adjustl(sRes)) == 'single') then; imode_precision_files = IO_TYPE_SINGLE
    else
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Wrong Main.FileType.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
    end if

    call ScanFile_Int(bakfile, inifile, 'Main', 'VerbosityLevel', '1', imode_verbosity)

    call ScanFile_Char(bakfile, inifile, 'Main', 'Type', 'temporal', sRes)
    if (trim(adjustl(sRes)) == 'temporal') then; imode_sim = DNS_MODE_TEMPORAL
    elseif (trim(adjustl(sRes)) == 'spatial') then; imode_sim = DNS_MODE_SPATIAL
    else
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Entry Main.Type must be temporal or spatial')
        call TLab_Stop(DNS_ERROR_SIMTYPE)
    end if

    if (imode_sim == DNS_MODE_TEMPORAL) fourier_on = .true.

    call ScanFile_Char(bakfile, inifile, 'Main', 'CalculateFlow', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; flow_on = .true.
    elseif (trim(adjustl(sRes)) == 'no') then; flow_on = .false.
    else
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Entry Main.CalculateFlow must be yes or no')
        call TLab_Stop(DNS_ERROR_CALCFLOW)
    end if

    call ScanFile_Char(bakfile, inifile, 'Main', 'CalculateScalar', 'yes', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; scal_on = .true.
    elseif (trim(adjustl(sRes)) == 'no') then; scal_on = .false.
    else
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Entry Main.CalculateScalar must be yes or no')
        call TLab_Stop(DNS_ERROR_CALCSCALAR)
    end if

! -------------------------------------------------------------------
    call ScanFile_Char(bakfile, inifile, 'Main', 'SpaceOrder', 'void', sRes)
    if (trim(adjustl(sRes)) == 'compactjacobian4') then; g(1:3)%mode_fdm1 = FDM_COM4_JACOBIAN; 
    elseif (trim(adjustl(sRes)) == 'compactjacobian6') then; g(1:3)%mode_fdm1 = FDM_COM6_JACOBIAN; 
    elseif (trim(adjustl(sRes)) == 'compactjacobian6hyper') then; g(1:3)%mode_fdm1 = FDM_COM6_JACOBIAN_HYPER; 
    elseif (trim(adjustl(sRes)) == 'compactjacobian6penta') then; g(1:3)%mode_fdm1 = FDM_COM6_JACOBIAN_PENTA; 
    elseif (trim(adjustl(sRes)) == 'compactdirect4') then; g(1:3)%mode_fdm1 = FDM_COM4_DIRECT; 
    elseif (trim(adjustl(sRes)) == 'compactdirect6') then; g(1:3)%mode_fdm1 = FDM_COM6_DIRECT; 
    else
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Wrong SpaceOrder option.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if
    g(1:3)%mode_fdm2 = g(1:3)%mode_fdm1

! ###################################################################
! Pressure staggering
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Staggering]')
    call TLab_Write_ASCII(bakfile, '#StaggerHorizontalPressure=<yes/no>')

    call ScanFile_Char(bakfile, inifile, 'Staggering', 'StaggerHorizontalPressure', 'no', sRes)
    if (trim(adjustl(sRes)) == 'yes') then; stagger_on = .true.; call TLab_Write_ASCII(lfile, 'Horizontal staggering of the pressure along Ox and Oz.')
    elseif (trim(adjustl(sRes)) == 'no') then; stagger_on = .false.
    else
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Entry Main. StaggerHorizontalPressure must be yes or no')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

! ###################################################################
! Grid Parameters
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Grid]')
    call TLab_Write_ASCII(bakfile, '#Imax=<imax>')
    call TLab_Write_ASCII(bakfile, '#Imax(*)=<imax_proc>')
    call TLab_Write_ASCII(bakfile, '#Jmax=<jmax>')
    call TLab_Write_ASCII(bakfile, '#Kmax=<kmax>')
    call TLab_Write_ASCII(bakfile, '#Kmax(*)=<kmax_proc>')
    call TLab_Write_ASCII(bakfile, '#XUniform=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#YUniform=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#ZUniform=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#XPeriodic=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#YPeriodic=<yes/no>')
    call TLab_Write_ASCII(bakfile, '#ZPeriodic=<yes/no>')

    call ScanFile_Int(bakfile, inifile, 'Grid', 'Imax', '0', g(1)%size)
    call ScanFile_Int(bakfile, inifile, 'Grid', 'Jmax', '0', g(2)%size)
    call ScanFile_Int(bakfile, inifile, 'Grid', 'Kmax', '0', g(3)%size)

! default
    imax = g(1)%size
    jmax = g(2)%size
    kmax = g(3)%size

    g(1)%name = 'x'
    g(2)%name = 'y'
    g(3)%name = 'z'

    do ig = 1, 3
        call ScanFile_Char(bakfile, inifile, 'Grid', g(ig)%name(1:1)//'Uniform', 'void', sRes)
        if (trim(adjustl(sRes)) == 'yes') then; g(ig)%uniform = .true.
        else if (trim(adjustl(sRes)) == 'no') then; g(ig)%uniform = .false.
        else
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Error in Uniform '//g(ig)%name(1:1)//' grid')
            call TLab_Stop(DNS_ERROR_UNIFORMX)
        end if

        call ScanFile_Char(bakfile, inifile, 'Grid', g(ig)%name(1:1)//'Periodic', 'void', sRes)
        if (trim(adjustl(sRes)) == 'yes') then; g(ig)%periodic = .true.
        else if (trim(adjustl(sRes)) == 'no') then; g(ig)%periodic = .false.
        else
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Error in Periodic '//g(ig)%name(1:1)//' grid')
            call TLab_Stop(DNS_ERROR_IBC)
        end if

        ! consistency check
        if (g(ig)%periodic .and. (.not. g(ig)%uniform)) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Grid must be uniform in periodic direction.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

    end do

    ! consistency check
    select case (imode_sim)
    case (DNS_MODE_TEMPORAL)
        if (.not. g(1)%periodic) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Grid must be uniform and periodic in direction X for temporal simulation')
            call TLab_Stop(DNS_ERROR_CHECKUNIFX)
        end if
    case (DNS_MODE_SPATIAL)
    end select

! -------------------------------------------------------------------
! Domain decomposition in parallel mode
! -------------------------------------------------------------------
#ifdef USE_MPI
    if (ims_npro > 1) then
        call ScanFile_Int(bakfile, inifile, 'Grid', 'Kmax(*)', '-1', kmax)
        if (kmax > 0 .and. mod(g(3)%size, kmax) == 0) then
            ims_npro_k = g(3)%size/kmax
        else
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Input kmax incorrect')
            call TLab_Stop(DNS_ERROR_KMAXTOTAL)
        end if

        call ScanFile_Int(bakfile, inifile, 'Grid', 'Imax(*)', '-1', imax)
        if (imax > 0 .and. mod(g(1)%size, imax) == 0) then
            ims_npro_i = g(1)%size/imax
        else
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Input imax incorrect')
            call TLab_Stop(DNS_ERROR_KMAXTOTAL)
        end if

        ! consistency check
        if (ims_npro_i*ims_npro_k == ims_npro) then
            write (lstr, *) ims_npro_i; write (sRes, *) ims_npro_k
            lstr = trim(adjustl(lstr))//'x'//trim(adjustl(sRes))
            call TLab_Write_ASCII(lfile, 'Initializing domain partition '//trim(adjustl(lstr)))
        else
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Inconsistency in total number of PEs')
            call TLab_Stop(DNS_ERROR_KMAXTOTAL)
        end if

    end if

#endif

! ###################################################################
! Filters
! ###################################################################
! Domain
    call FILTER_READBLOCK(bakfile, inifile, 'Filter', FilterDomain)
    FilterDomainActive(:) = .true.                      ! Variable to eventually allow for control field by field
    FilterDomainBcsFlow(:) = FilterDomain(2)%BcsMin     ! To allow a difference between flow and scalar
    FilterDomainBcsScal(:) = FilterDomain(2)%BcsMin     ! Can be modified later depending on BCs of equations

! Pressure
    call FILTER_READBLOCK(bakfile, inifile, 'PressureFilter', PressureFilter)

! ###################################################################
! specific data for spatial mode
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')
    call TLab_Write_ASCII(bakfile, '#[Statistics]')
    call TLab_Write_ASCII(bakfile, '#IAvera=<plane1,plane2,...>')

    call ScanFile_Char(bakfile, inifile, 'Statistics', 'IAvera', 'void', sRes)
    if (trim(adjustl(sRes)) /= 'void') then
        ! nstatavg = MAX_STATS_SPATIAL
        ! call LIST_INTEGER(sRes, nstatavg, statavg)
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Initialization of spatial statistics to be updated.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
    end if

! ###################################################################
! Initialization of array sizes
! ###################################################################
    call TLab_Write_ASCII(bakfile, '#')

    isize_field = imax*jmax*kmax

! auxiliar array txc for intermediate calculations
    isize_txc_field = imax*jmax*kmax
    if (fourier_on) then
        isize_txc_dimz = (imax + 2)*(jmax + 2)
        isize_txc_dimx = kmax*(jmax + 2)
        isize_txc_field = isize_txc_dimz*kmax ! space for FFTW lib
#ifdef USE_MPI
        if (ims_npro_k > 1) then
            if (mod(isize_txc_dimz, (2*ims_npro_k)) /= 0) then ! add space for MPI transposition
                isize_txc_dimz = isize_txc_dimz/(2*ims_npro_k)
                isize_txc_dimz = (isize_txc_dimz + 1)*(2*ims_npro_k)
            end if
            isize_txc_field = max(isize_txc_field, isize_txc_dimz*kmax)
        end if
        if (ims_npro_i > 1) then
            if (mod(isize_txc_dimx, (2*ims_npro_i)) /= 0) then ! add space for MPI transposition
                isize_txc_dimx = isize_txc_dimx/(2*ims_npro_i)
                isize_txc_dimx = (isize_txc_dimx + 1)*(2*ims_npro_i)
            end if
            isize_txc_field = max(isize_txc_field, isize_txc_dimx*(imax + 2))
        end if
#endif
        if (mod(imax, 2) /= 0) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Imax must be a multiple of 2 for the FFT operations.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if
    end if

! scratch arrays
    isize_wrk1d = max(g(1)%size, max(g(2)%size, g(3)%size))
    isize_wrk2d = max(imax*jmax, max(imax*kmax, jmax*kmax))
    isize_wrk3d = max(isize_field, isize_txc_field)

    return
end subroutine TLab_Initialize_Parameters
