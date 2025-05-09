#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI

#endif

#define C_FILE_LOC "TRANSFORM"

program TRANSFIELDS

    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal, tag_part
    use FDM, only: fdm_dt
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_txc, isize_wrk3d, inb_flow, inb_scal, isize_txc_field, isize_wrk1d, inb_wrk1d
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, scal_on, flow_on, fourier_on
    use TLab_Memory, only: TLab_Initialize_Memory, TLab_Allocate_Real
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_k
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use FDM, only: g, FDM_Initialize
    use Thermodynamics
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use TLab_Background, only: TLab_Initialize_Background, qbg, sbg
    use Gravity, only: Gravity_Initialize
    use Rotation, only: Rotation_Initialize
    use Thermo_Anelastic
    use IO_Fields
    use OPR_FILTERS
    use OPR_INTERPOLATORS
    use OPR_Fourier
    use TLab_Grid

    implicit none

    ! Parameter definitions
    integer(wi), parameter :: itime_size_max = 3000
    integer(wi), parameter :: iopt_size_max = 512

    ! -------------------------------------------------------------------
    ! Additional local arrays
    type(grid_dt), target :: g_dst(3)
    type(grid_dt), pointer :: x_dst => g_dst(1), y_dst => g_dst(2), z_dst => g_dst(3)       ! to be cleaned
    real(wp), allocatable, save :: q_dst(:, :), s_dst(:, :)

    real(wp), allocatable, save :: x_aux(:), y_aux(:), z_aux(:)
    real(wp), pointer :: txc_aux(:, :, :) => null()

    ! -------------------------------------------------------------------
    ! Local variables
    ! -------------------------------------------------------------------
    integer(wi) opt_main, opt_function
    integer(wi) iq, is, ig, ip, j, k
    integer(wi) idummy
    logical iread_flow, iread_scal
    character*32 bakfile, flow_file, scal_file
    character*64 str
    character*512 sRes
    integer(wi) subdomain(6)

    integer(wi) imax_dst, jmax_dst, kmax_dst

    logical :: flag_crop = .false., flag_extend = .false.
    integer(wi) jmax_aux, inb_scal_dst
    real(wp) dummy, tolerance

    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi) iopt_size
    real(wp) opt_vec(iopt_size_max)

    real(wp) params(1)!, scales(3)

    ! ###################################################################
    bakfile = trim(adjustl(ifile))//'.bak'

    call TLab_Start

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
    call Rotation_Initialize(ifile)

    call TLab_Consistency_Check()

    ! -------------------------------------------------------------------
    ! File names
    ! -------------------------------------------------------------------
#include "dns_read_times.h"

    ! -------------------------------------------------------------------
    ! Read local options
    ! -------------------------------------------------------------------
    opt_main = -1 ! default values
    opt_function = 0

    call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'ParamTransform', '-1', sRes)
    iopt_size = iopt_size_max
    call LIST_REAL(sRes, iopt_size, opt_vec)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, '(A)') 'Option ?'
        write (*, '(A)') '1. Crop fields'
        write (*, '(A)') '2. Extend fields in Ox and Oy'
        write (*, '(A)') '3. Remesh fields'
        write (*, '(A)') '4. Linear combination of fields'
        write (*, '(A)') '5. Filter fields'
        write (*, '(A)') '6. Mapping of scalar fields'
        write (*, '(A)') '7. Blend fields'
        write (*, '(A)') '8. Add mean profiles'
        write (*, '(A)') '9. Extrude fields in Oz'
        write (*, '(A)') '10. Change to single precision'
        read (*, *) opt_main
#endif
    else
        opt_main = int(opt_vec(1))
    end if

    call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'Subdomain', '-1', sRes)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, *) 'Subdomain limits ?'
        read (*, '(A64)') sRes
#endif
    end if
    idummy = 6
    call LIST_INTEGER(sRes, idummy, subdomain)

    if (idummy < 6) then ! default
        subdomain(1) = 1; subdomain(2) = g(1)%size
        subdomain(3) = 1; subdomain(4) = g(2)%size
        subdomain(5) = 1; subdomain(6) = g(3)%size
    end if

    ! -------------------------------------------------------------------
    select case (opt_main)
    case (:0)
        call TLab_Write_ASCII(efile, C_FILE_LOC//'. Missing input [ParamTransform] in tlab.ini.')
        call TLab_Stop(DNS_ERROR_INVALOPT)

    case (1)
        if (subdomain(1) /= 1 .or. subdomain(2) /= g(1)%size .or. &
            subdomain(5) /= 1 .or. subdomain(6) /= g(3)%size) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Cropping only in Oy.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (subdomain(3) < 1 .or. subdomain(4) > g(2)%size) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Cropping out of bounds in Oy.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

    case (4)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Coefficients ?'
            read (*, '(A512)') sRes
            iopt_size = iopt_size_max - 1
            call LIST_REAL(sRes, iopt_size, opt_vec(2))
#endif
        else
            iopt_size = iopt_size - 1
        end if
        if (iopt_size == 0) then
            call TLab_Write_ASCII(lfile, C_FILE_LOC//'. Performing arithmetic mean of fields.')
            iopt_size = itime_size
            opt_vec(2:) = 1.0_wp/real(itime_size, wp)
        end if
        if (iopt_size /= itime_size) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Number of coefficient incorrect.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if
        if (iopt_size > iopt_size_max) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Array opt_vec too small.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

    case (5)
        if (FilterDomain(1)%type == DNS_FILTER_NONE) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Filter information needs to be provided in block [Filter].')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

    case (6)
        flow_on = .false. ! Force not to process the flow fields

        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, '(A)') 'Function type ?'
            write (*, '(A)') '1. Vapor Saturation Humidity'
            write (*, '(A)') '2. Linear'
            read (*, *) opt_function
#endif
        else
            opt_function = int(opt_vec(2))
        end if

        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Coefficients ?'
            read (*, '(A512)') sRes
            iopt_size = iopt_size_max - 2
            call LIST_REAL(sRes, iopt_size, opt_vec(3))
#endif
        else
            iopt_size = iopt_size - 2
        end if

    case (7) ! 2nd and 3rd entries in opt_vec contain coeffs.
        if (sRes == '-1') then
            write (*, *) 'Coefficients ?'
            read (*, '(A512)') sRes
            iopt_size = 2
            call LIST_REAL(sRes, iopt_size, opt_vec(2))
        else
            iopt_size = iopt_size - 1
        end if
        if (iopt_size /= 2) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Number of blend coefficient incorrect.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

    end select

    ! -------------------------------------------------------------------
    select case (opt_main)
    case (1, 3, 9) ! Crop, remesh, extrude
        g_dst(1)%size = subdomain(2) - subdomain(1) + 1
        g_dst(2)%size = subdomain(4) - subdomain(3) + 1
        g_dst(3)%size = subdomain(6) - subdomain(5) + 1

    case (2)   ! Extend
        g_dst(1)%size = g(1)%size + subdomain(2) + subdomain(1)
        g_dst(2)%size = g(2)%size + subdomain(4) + subdomain(3)
        g_dst(3)%size = g(3)%size

    case DEFAULT
        g_dst(1)%size = g(1)%size
        g_dst(2)%size = g(2)%size
        g_dst(3)%size = g(3)%size

    end select

#ifdef USE_MPI
    imax_dst = g_dst(1)%size/ims_npro_i
    jmax_dst = g_dst(2)%size
    kmax_dst = g_dst(3)%size/ims_npro_k
#else
    imax_dst = g_dst(1)%size
    jmax_dst = g_dst(2)%size
    kmax_dst = g_dst(3)%size
#endif

    iread_flow = flow_on
    iread_scal = scal_on

    ! #######################################################################
    ! Initialize memory space and grid data
    ! #######################################################################
    inb_txc = 0
    inb_scal_dst = inb_scal
    select case (opt_main)
    case (3)            ! Remesh
        isize_txc_field = max(isize_txc_field, imax_dst*jmax_dst*kmax_dst)
        inb_txc = 5
    case (5)            ! Filter
        inb_txc = 4
    case (6)            ! Scalar field mapping
        inb_txc = 5
        inb_scal_dst = 1
    end select
    isize_wrk3d = max(isize_wrk3d, isize_txc_field)
    isize_wrk3d = max(isize_wrk3d, imax_dst*jmax_dst*kmax_dst)
    if (fourier_on) inb_txc = max(inb_txc, 1)

    call TLab_Initialize_Memory(C_FILE_LOC)

    call TLab_Initialize_Background(ifile)

    ! Further allocation
    if (flow_on) call TLab_Allocate_Real(__FILE__, q_dst, [imax_dst*jmax_dst*kmax_dst, inb_flow], 'flow-dst')
    if (scal_on) call TLab_Allocate_Real(__FILE__, s_dst, [imax_dst*jmax_dst*kmax_dst, inb_scal_dst], 'scal-dst')

    ! ###################################################################
    ! Initialize operators and reference data
    ! ###################################################################
    if (opt_main == 5) then
        call OPR_Filter_Initialize_Parameters(ifile)
        do ig = 1, 3
            call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig))
        end do
    end if

    if (fourier_on) call OPR_Fourier_Initialize()

    call OPR_CHECK()

    ! -------------------------------------------------------------------
    ! Initialize cumulative field
    ! -------------------------------------------------------------------
    if (opt_main == 4 .or. opt_main == 7) then
        if (flow_on) q_dst = 0.0_wp
        if (scal_on) s_dst = 0.0_wp
    end if

    ! -------------------------------------------------------------------
    ! Initialize remeshing
    ! -------------------------------------------------------------------
    if (opt_main == 3) then
        call TLab_Grid_Read('grid.trn', x_dst, y_dst, z_dst, [g_dst(1)%size, g_dst(2)%size, g_dst(3)%size])
        g_dst(1:3)%periodic = g(1:3)%periodic
        
        tolerance = 0.001_wp    ! percentage of grid spacing

#define xn(i) x%nodes(i)
#define yn(i) y%nodes(i)
#define zn(i) z%nodes(i)

        ! Check grids; Ox and Oz directions are assumed to be periodic
        dummy = (g_dst(1)%scale - g(1)%scale)/(xn(g(1)%size) - xn(g(1)%size - 1))
        if (abs(dummy) > tolerance) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Ox scales are not equal at the end.')
            call TLab_Stop(DNS_ERROR_GRID_SCALE)
        end if

        dummy = (g_dst(3)%scale - g(3)%scale)/(zn(g(3)%size) - zn(g(3)%size - 1))
        if (abs(dummy) > tolerance) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Oz scales are not equal')
            call TLab_Stop(DNS_ERROR_GRID_SCALE)
        end if

        ! In the Oy direction, we allow to have a different box
        jmax_aux = g(2)%size; subdomain = 0

        dummy = (y_dst%nodes(g_dst(2)%size) - yn(g(2)%size))/(yn(g(2)%size) - yn(g(2)%size - 1))
        if (dummy > tolerance) then                 ! Extend
            flag_extend = .true.
            subdomain(4) = int(dummy) + 1           ! # planes to add at the top
            jmax_aux = jmax_aux + subdomain(4)
        else if (dummy < -tolerance) then           ! Crop
            flag_crop = .true.
            do j = jmax - 1, 1, -1
                if ((g(2)%nodes(j) - y_dst%nodes(g_dst(2)%size))*(g(2)%nodes(j + 1) - y_dst%nodes(g_dst(2)%size)) < 0.0_wp) exit
            end do
            subdomain(4) = j + 1                    ! top plane of cropped region
            jmax_aux = subdomain(4)
            subdomain(3) = 1
        end if

        dummy = (y_dst%nodes(1) - yn(1))/(yn(2) - yn(1))
        if (dummy < -tolerance) then                ! Extend
            flag_extend = .true.
            subdomain(3) = int(abs(dummy)) + 1      ! # planes to add at the bottom
            jmax_aux = jmax_aux + subdomain(3)
        else if (dummy > tolerance) then            ! Crop
            flag_crop = .true.
            do j = 1, jmax - 1, 1
                if ((g(2)%nodes(j) - y_dst%nodes(1))*(g(2)%nodes(j + 1) - y_dst%nodes(1)) < 0.0_wp) exit
            end do
            subdomain(3) = j                        ! bottom plane of cropped region
            jmax_aux = jmax_aux - subdomain(3) + 1
        end if

        if (flag_extend .and. flag_crop) then
            call TLab_Write_ASCII(efile, C_FILE_LOC//'. Simultaneous extend and crop is undeveloped.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! Reallocating memory space because jmax_aux can be larger than jmax, jmax_dst
        isize_wrk1d = max(isize_wrk1d, jmax_aux)

        idummy = max(imax, imax_dst)*max(jmax_aux, max(jmax, jmax_dst))*max(kmax, kmax_dst)
        isize_txc_field = max(isize_txc_field, idummy)
#ifdef USE_MPI
        idummy = kmax*jmax_aux
        if (mod(idummy, ims_npro_i) /= 0) then  ! add space for MPI transposition
            idummy = idummy/ims_npro_i
            idummy = (idummy + 1)*ims_npro_i
        end if
        idummy = idummy*max(imax, imax_dst)
        isize_txc_field = max(isize_txc_field, idummy)
#endif
        isize_wrk3d = max(isize_wrk3d, isize_txc_field)

        deallocate (txc, wrk1d, wrk3d)
        call TLab_Allocate_Real(C_FILE_LOC, txc, [isize_txc_field, inb_txc], 'txc')
        call TLab_Allocate_Real(C_FILE_LOC, wrk1d, [isize_wrk1d, inb_wrk1d], 'wrk1d')
        call TLab_Allocate_Real(C_FILE_LOC, wrk3d, [isize_wrk3d], 'wrk3d')
        txc_aux(1:imax, 1:jmax_aux, 1:kmax) => txc(1:imax*jmax_aux*kmax, 1)

        allocate (x_aux(g(1)%size + 1))         ! need extra space in cubic splines
        allocate (z_aux(g(3)%size + 1))
        allocate (y_aux(jmax_aux + 1))

        x_aux(1:g(1)%size) = xn(1:g(1)%size)  ! need extra space in cubic splines
        z_aux(1:g(3)%size) = zn(1:g(3)%size)  ! need extra space in cubic splines

        ! Creating grid
        if (flag_crop) then
            write (str, '(I3)') subdomain(4)
            call TLab_Write_ASCII(lfile, 'Croping above '//trim(adjustl(str))//' for remeshing...')
            write (str, '(I3)') subdomain(3)
            call TLab_Write_ASCII(lfile, 'Croping below '//trim(adjustl(str))//' for remeshing...')
            call TRANS_CROP(1, jmax, 1, subdomain, g(2)%nodes, y_aux)

            y_aux(1) = y_dst%nodes(1)                 ! Using min and max of new grid
            y_aux(jmax_aux) = y_dst%nodes(g_dst(2)%size)

        else
            y_aux(1 + subdomain(3):g(2)%size + subdomain(3)) = yn(1:g(2)%size) ! we need extra space

            if (subdomain(4) > 0) then
                write (str, '(I3)') subdomain(4)
                call TLab_Write_ASCII(lfile, 'Adding '//trim(adjustl(str))//' planes at the top for remeshing...')
                dummy = (y_dst%nodes(g_dst(2)%size) - yn(g(2)%size))/real(subdomain(4)) ! distributing the points uniformly
                do ip = g(2)%size + subdomain(3) + 1, g(2)%size + subdomain(3) + subdomain(4)
                    y_aux(ip) = y_aux(ip - 1) + dummy
                end do
            end if

            if (subdomain(3) > 0) then
                write (str, '(I3)') subdomain(3)
                call TLab_Write_ASCII(lfile, 'Adding '//trim(adjustl(str))//' planes at the bottom for remeshing...')
                dummy = (y_dst%nodes(1) - yn(1))/real(subdomain(3))
                do ip = subdomain(3), 1, -1
                    y_aux(ip) = y_aux(ip + 1) + dummy ! dummy is negative
                end do
            end if

        end if

        call TLab_Write_ASCII(efile, 'Changing grid variable. Code to be fixed.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
        ! g(2)%scale = g_dst(2)%scale     ! watch out, overwriting grid information
        ! g(2)%size = jmax_aux

    end if

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))
        call TLab_Write_ASCII(lfile, sRes)

        if (iread_flow) then ! Flow variables
            write (flow_file, *) itime; flow_file = trim(adjustl(tag_flow))//trim(adjustl(flow_file))
            call IO_Read_Fields(flow_file, imax, jmax, kmax, itime, inb_flow, 0, q, params)
            rtime = params(1)
        end if

        if (iread_scal) then ! Scalar variables
            write (scal_file, *) itime; scal_file = trim(adjustl(tag_scal))//trim(adjustl(scal_file))
            call IO_Read_Fields(scal_file, imax, jmax, kmax, itime, inb_scal, 0, s, params)
            rtime = params(1)
        end if

        select case (opt_main)
            ! ###################################################################
            ! Cropping
            ! ###################################################################
        case (1)
            if (flow_on) then
                do iq = 1, inb_flow
                    call TLab_Write_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_CROP(imax, jmax, kmax, subdomain, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (scal_on) then
                do is = 1, inb_scal
                    call TLab_Write_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_CROP(imax, jmax, kmax, subdomain, s(1, is), s_dst(1, is))
                end do
            end if

            ! ###################################################################
            ! Extension
            ! ###################################################################
        case (2)
            if (flow_on) then
                do iq = 1, inb_flow
                    call TLab_Write_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_EXTEND(imax, jmax, kmax, subdomain, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (scal_on) then
                do is = 1, inb_scal
                    call TLab_Write_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_EXTEND(imax, jmax, kmax, subdomain, s(1, is), s_dst(1, is))
                end do
            end if

            ! ###################################################################
            ! Change grid
            ! ###################################################################
        case (3)

            if (flow_on) then
                do iq = 1, inb_flow
                    call TLab_Write_ASCII(lfile, 'Transfering data to new array...')
                    if (flag_crop) then
                        call TRANS_CROP(imax, jmax, kmax, subdomain, q(:, iq), txc_aux)
                        do k = 1, kmax
                            txc_aux(:, 1, k) = txc_aux(:, 1, k) &
                                               + (y_aux(1) - yn(subdomain(3)))*(txc_aux(:, 2, k) - txc_aux(:, 1, k))/(yn(subdomain(3) + 1) - yn(subdomain(3)))
                            txc_aux(:, jmax_aux, k) = txc_aux(:, jmax_aux - 1, k) &
                    + (y_aux(jmax_aux) - yn(subdomain(4) - 1))*(txc_aux(:, jmax_aux, k) - txc_aux(:, jmax_aux - 1, k))/(yn(subdomain(4)) - yn(subdomain(4) - 1))
                        end do
                    else
                        call TRANS_EXTEND(imax, jmax, kmax, subdomain, q(:, iq), txc_aux)
                    end if
                    call OPR_INTERPOLATE(imax, jmax_aux, kmax, &
                                         x_aux, y_aux, z_aux, x_dst, y_dst, z_dst, &
                                         txc_aux, q_dst(:, iq), txc(:, 2))
                end do
            end if

            if (scal_on) then
                do is = 1, inb_scal
                    call TLab_Write_ASCII(lfile, 'Transfering data to new array...')
                    if (flag_crop) then
                        call TRANS_CROP(imax, jmax, kmax, subdomain, s(:, is), txc_aux)
                        do k = 1, kmax
                            txc_aux(:, 1, k) = txc_aux(:, 1, k) &
                                               + (y_aux(1) - yn(subdomain(3)))*(txc_aux(:, 2, k) - txc_aux(:, 1, k))/(yn(subdomain(3) + 1) - yn(subdomain(3)))
                            txc_aux(:, jmax_aux, k) = txc_aux(:, jmax_aux - 1, k) &
                    + (y_aux(jmax_aux) - yn(subdomain(4) - 1))*(txc_aux(:, jmax_aux, k) - txc_aux(:, jmax_aux - 1, k))/(yn(subdomain(4)) - yn(subdomain(4) - 1))
                        end do
                    else
                        call TRANS_EXTEND(imax, jmax, kmax, subdomain, s(:, is), txc_aux)
                    end if
                    call OPR_INTERPOLATE(imax, jmax_aux, kmax,  &
                                         x_aux, y_aux, z_aux, x_dst, y_dst, z_dst, &
                                         txc_aux, s_dst(:, is), txc(:, 2))
                end do
            end if

            ! ###################################################################
            ! Linear combination of fields
            ! ###################################################################
        case (4)
            if (flow_on) then
                q_dst = q_dst + q*opt_vec(it + 1)
            end if

            if (scal_on) then
                s_dst = s_dst + s*opt_vec(it + 1)
            end if

            ! ###################################################################
            ! Filter
            ! ###################################################################
        case (5)
            if (flow_on) then
                do iq = 1, inb_flow
                    call TLab_Write_ASCII(lfile, 'Filtering...')
                    q_dst(:, iq) = q(:, iq) ! in-place operation
                    if (FilterDomain(1)%type == DNS_FILTER_HELMHOLTZ) &  ! Bcs depending on field
                        FilterDomain(2)%BcsMin = FilterDomainBcsFlow(iq)
                    call OPR_FILTER(imax, jmax, kmax, FilterDomain, q_dst(1, iq), txc)
                end do
            end if

            if (scal_on) then
                do is = 1, inb_scal
                    call TLab_Write_ASCII(lfile, 'Filtering...')
                    s_dst(:, is) = s(:, is) ! in-place operation
                    if (FilterDomain(1)%type == DNS_FILTER_HELMHOLTZ) & ! Bcs depending on field
                        FilterDomain(2)%BcsMin = FilterDomainBcsScal(is)
                    call OPR_FILTER(imax, jmax, kmax, FilterDomain, s_dst(1, is), txc)
                end do
            end if

            ! ###################################################################
            ! Transformation
            ! ###################################################################
        case (6)
            if (opt_function == 1) then
                call TRANS_FUNCTION(imax, jmax, kmax, s, s_dst, txc)

            else if (opt_function == 2) then
                s_dst(:, 1) = 0.0_wp
                do is = 1, min(inb_scal, iopt_size)
                    s_dst(:, 1) = s_dst(:, 1) + opt_vec(2 + is)*s(:, is)
                end do

            end if

            ! ###################################################################
            ! Blend
            ! ###################################################################
        case (7)
            if (it == 1) opt_vec(2) = yn(1) + opt_vec(2)*g(2)%scale
            write (sRes, *) opt_vec(2), opt_vec(3); sRes = 'Blending with '//trim(adjustl(sRes))
            call TLab_Write_ASCII(lfile, sRes)

            if (scal_on) then
                do is = 1, inb_scal
                    call TRANS_BLEND(imax, jmax, kmax, opt_vec(2), y%nodes, s(1, is), s_dst(1, is))
                end do
            end if

            if (flow_on) then ! Blended fields have rtime from last velocity field
                do iq = 1, inb_flow
                    call TRANS_BLEND(imax, jmax, kmax, opt_vec(2), y%nodes, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (it == 1) opt_vec(3) = -opt_vec(3) ! flipping blending shape

            ! ###################################################################
            ! Adding mean profiles
            ! ###################################################################
        case (8)
            if (flow_on) then
                do iq = 1, inb_flow
                    call TLab_Write_ASCII(lfile, 'Adding mean flow profiles...')
                    call TRANS_ADD_MEAN(0, iq, imax, jmax, kmax, y%nodes, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (scal_on) then
                do is = 1, inb_scal
                    call TLab_Write_ASCII(lfile, 'Adding mean scal profiles...')
                    call TRANS_ADD_MEAN(1, is, imax, jmax, kmax, y%nodes, s(1, is), s_dst(1, is))
                end do
            end if

            ! ###################################################################
            ! Extrude
            ! ###################################################################
        case (9)
            if (flow_on) then
                do iq = 1, inb_flow
                    call TLab_Write_ASCII(lfile, 'Extruding along Oz...')
                    call TRANS_EXTRUDE(imax, jmax, kmax, subdomain, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (scal_on) then
                do is = 1, inb_scal
                    call TLab_Write_ASCII(lfile, 'Extruding along Oz...')
                    call TRANS_EXTRUDE(imax, jmax, kmax, subdomain, s(1, is), s_dst(1, is))
                end do
            end if

            ! ###################################################################
            ! Change to single precision
            ! ###################################################################
        case (10)
            q_dst(:, 1:inb_flow) = q(:, 1:inb_flow)
            s_dst(:, 1:inb_scal) = s(:, 1:inb_scal)

            idummy = io_datatype
            io_datatype = IO_TYPE_SINGLE

        end select

        ! ###################################################################
        ! Writing transform fields
        ! ###################################################################
        if (opt_main /= 4 .and. opt_main /= 7) then
            if (flow_on) then
                flow_file = trim(adjustl(flow_file))//'.trn'
                call IO_Write_Fields(flow_file, imax_dst, jmax_dst, kmax_dst, itime, inb_flow, q_dst, io_header_q(1:1))
            end if
            if (scal_on) then
                scal_file = trim(adjustl(scal_file))//'.trn'
                call IO_Write_Fields(scal_file, imax_dst, jmax_dst, kmax_dst, itime, inb_scal_dst, s_dst, io_header_s(1:inb_scal_dst))
            end if
        end if

        if (opt_main == 10) then
            io_datatype = idummy
        end if

    end do

    ! ###################################################################
    ! Final operations
    ! ###################################################################
    if (opt_main == 4 .or. opt_main == 7) then
        if (flow_on) then
            flow_file = trim(adjustl(flow_file))//'.trn'
            call IO_Write_Fields(flow_file, imax_dst, jmax_dst, kmax_dst, itime, inb_flow, q_dst, io_header_q(1:1))
        end if
        if (scal_on) then
            scal_file = trim(adjustl(scal_file))//'.trn'
            call IO_Write_Fields(scal_file, imax_dst, jmax_dst, kmax_dst, itime, inb_scal_dst, s_dst, io_header_s(1:inb_scal_dst))
        end if
    end if

    call TLab_Stop(0)

contains
    !########################################################################
    !# Crop array a into array b in the first two indices
    !########################################################################
    subroutine TRANS_CROP(nx, ny, nz, subdomain, a, b)
        integer(wi) nx, ny, nz, subdomain(6)
        real(wp), dimension(nx, ny, nz) :: a
        real(wp), dimension(nx, subdomain(4) - subdomain(3) + 1, nz) :: b

        ! #######################################################################
        do k = 1, nz
            do j = subdomain(3), subdomain(4)
                b(:, j - subdomain(3) + 1, k) = a(:, j, k)
            end do
        end do

        return
    end subroutine TRANS_CROP

    !########################################################################
    !# Crop array a into array b in the first two indices
    !########################################################################
    subroutine TRANS_EXTRUDE(nx, ny, nz, subdomain, a, b)
        integer(wi) nx, ny, nz, subdomain(6)
        real(wp), intent(IN) :: a(nx, ny, nz)
        real(wp), intent(OUT) :: b(subdomain(2) - subdomain(1) + 1, subdomain(4) - subdomain(3) + 1, subdomain(6) - subdomain(5) + 1)

        ! #######################################################################
        do k = 1, subdomain(6) - subdomain(5) + 1
            b(:, :, k) = a(subdomain(1):subdomain(2), subdomain(3):subdomain(4), 1)
        end do

        return
    end subroutine TRANS_EXTRUDE

    !########################################################################
    !# Extend array a into array b in the first two indices
    !########################################################################
    subroutine TRANS_EXTEND(nx, ny, nz, planes, a, b)
        integer(wi) nx, ny, nz, planes(6)
        real(wp), dimension(nx, ny, nz) :: a
        real(wp), dimension(planes(1) + nx + planes(2), planes(3) + ny + planes(4), nz) :: b

        ! -----------------------------------------------------------------------
        integer(wi) j, k

        ! #######################################################################
        do k = 1, nz
            b(1 + planes(1):nx + planes(1), 1 + planes(3):ny + planes(3), k) = a(1:nx, 1:ny, k)

            ! extension in i
            do j = 1, ny
                b(1:planes(1), j, k) = b(1 + planes(1), j, k)
                b(nx + planes(1) + 1:nx + planes(1) + planes(2), j, k) = b(nx + planes(1), j, k)
            end do

            ! extension in j; corners are now written
            do j = 1, planes(3)
                b(:, j, k) = b(:, 1 + planes(3), k)
            end do

            do j = ny + planes(3) + 1, ny + planes(3) + planes(4)
                b(:, j, k) = b(:, ny + planes(3), k)
            end do

        end do

        return
    end subroutine TRANS_EXTEND

    !########################################################################
    !########################################################################
    subroutine TRANS_ADD_MEAN(flag_mode, is, nx, ny, nz, y, a, b)
        use Profiles, only: Profiles_Calculate

        integer(wi) flag_mode, is, nx, ny, nz
        real(wp), dimension(*), intent(IN) :: y
        real(wp), dimension(nx, ny, nz), intent(IN) :: a
        real(wp), dimension(nx, ny, nz), intent(OUT) :: b

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy

        ! #######################################################################
        if (flag_mode == 0) then ! Velocity
            if (is == 1) then ! Only the mean velocity
                do j = 1, ny
                    dummy = Profiles_Calculate(qbg(1), y(j))
                    b(:, j, :) = dummy + a(:, j, :)
                end do
            else
                b = a
            end if

        else                         ! Scalars
            do j = 1, ny
                dummy = Profiles_Calculate(sbg(is), y(j))
                b(:, j, :) = dummy + a(:, j, :)
            end do

        end if

        return
    end subroutine TRANS_ADD_MEAN

    !########################################################################
    !# Calculate b = f(a)
    !########################################################################
    subroutine TRANS_FUNCTION(nx, ny, nz, a, b, txc)
        use Thermodynamics, only: imixture
        use Thermodynamics, only: rd_ov_rv, Lvl
        use THERMO_ANELASTIC, only: pbackground

        integer(wi) nx, ny, nz
        real(wp), dimension(nx*ny*nz) :: a, b
        real(wp), dimension(nx*ny*nz, *) :: txc

        ! -----------------------------------------------------------------------
        real(wp) qt_0, qt_1, h_0, h_1, p(1)

        ! #######################################################################
        imixture = MIXT_TYPE_AIRWATER
        call Thermodynamics_Initialize_Parameters(ifile)
        call Gravity_Initialize(ifile)
        call Rotation_Initialize(ifile)
        inb_scal = 1

        qt_0 = 9.0d-3; qt_1 = 1.5d-3
        h_0 = 0.955376d0; h_1 = 0.981965d0
        p = 0.940d0

        txc(:, 1) = h_0 + a(:)*(h_1 - h_0)          ! total enthalpy
        txc(:, 2) = qt_0 + a(:)*(qt_1 - qt_0)       ! total water, space for q_l
        txc(:, 3) = 0.0_wp
        txc(:, 4) = p                               ! pressure

        pbackground = p
        call Thermo_Anelastic_PH(nx, ny, nz, txc(1, 2), txc(1, 1))        ! Calculate q_l
        call Thermo_Anelastic_TEMPERATURE(nx, ny, nz, txc(1, 1), txc(1, 5))

        ! Calculate saturation specific humidity
        call Thermo_Psat_Polynomial(nx*ny*nz, txc(1, 5), txc(1, 1))
        txc(:, 1) = 1.0_wp/(txc(:, 4)/txc(:, 1) - 1.0_wp)*rd_ov_rv
        txc(:, 1) = txc(:, 1)/(1.0_wp + txc(:, 1))

        ! Calculate parameter eta (assuming c_p = c_p,d)
        txc(:, 3) = rd_ov_rv*Lvl*Lvl/(txc(:, 5)*txc(:, 5))

        ! Calculate s
        b(:) = txc(:, 2) - txc(:, 1)*(1.0_wp + txc(:, 3)*txc(:, 2))/(1.0_wp + txc(:, 3)*txc(:, 1))

        return
    end subroutine TRANS_FUNCTION

    !########################################################################
    !# b <- b + f(y)*a
    !########################################################################
    subroutine TRANS_BLEND(nx, ny, nz, params, y, a, b)
        integer(wi) nx, ny, nz
        real(wp), dimension(*) :: params
        real(wp), dimension(ny) :: y
        real(wp), dimension(nx, ny, nz) :: a, b

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) shape, xi

        ! #######################################################################
        do j = 1, ny
            xi = (y(j) - params(1))/params(2)
            shape = 0.5_wp*(1.0_wp + tanh(-0.5_wp*xi))
            b(:, j, :) = b(:, j, :) + shape*a(:, j, :)
        end do

        return

    end subroutine TRANS_BLEND

end program TRANSFIELDS
