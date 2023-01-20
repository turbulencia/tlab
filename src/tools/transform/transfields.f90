#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "TRANSFORM"

program TRANSFIELDS

    use TLAB_TYPES, only: filter_dt, grid_dt
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_VARS, only: ims_npro_i, ims_npro_k
    use TLAB_MPI_PROCS
#endif
    use IO_FIELDS
    use OPR_FILTERS
    use OPR_INTERPOLATORS
    use OPR_FOURIER

    implicit none

#include "integers.h"

    ! Parameter definitions
    TINTEGER, parameter :: itime_size_max = 3000
    TINTEGER, parameter :: iopt_size_max = 512

    ! -------------------------------------------------------------------
    ! Additional local arrays
    TREAL, dimension(:), allocatable, save :: x_dst, y_dst, z_dst
    TREAL, dimension(:, :), allocatable, save :: q_dst, s_dst

    TREAL, dimension(:), allocatable, save :: y_aux
    TREAL, dimension(:, :, :), allocatable, save :: txc_aux

    ! -------------------------------------------------------------------
    ! Local variables
    ! -------------------------------------------------------------------
    TINTEGER opt_main, opt_function
    TINTEGER iq, is, ig, ip, j, k
    TINTEGER idummy, iread_flow, iread_scal
    character*32 bakfile, flow_file, scal_file
    character*64 str
    character*512 sRes
    TINTEGER subdomain(6)

    type(grid_dt), dimension(3) :: g_dst
    TINTEGER imax_dst, jmax_dst, kmax_dst

    logical flag_crop, flag_extend
    TINTEGER jmax_aux, inb_scal_dst
    TREAL dummy

    TINTEGER itime_size, it
    TINTEGER itime_vec(itime_size_max)

    TINTEGER iopt_size
    TREAL opt_vec(iopt_size_max)

    ! ###################################################################
    bakfile = trim(adjustl(ifile))//'.bak'

    call TLAB_START

    call IO_READ_GLOBAL(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    ! -------------------------------------------------------------------
    ! File names
    ! -------------------------------------------------------------------
#include "dns_read_times.h"

    ! -------------------------------------------------------------------
    ! Read local options
    ! -------------------------------------------------------------------
    opt_main = -1 ! default values
    opt_function = 0

    call SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamTransform', '-1', sRes)
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
        write (*, '(A)') '6. Transform scalar fields'
        write (*, '(A)') '7. Blend fields'
        write (*, '(A)') '8. Add mean profiles'
        write (*, '(A)') '9. Extrude fields in Oz'
        read (*, *) opt_main
#endif
    else
        opt_main = int(opt_vec(1))
    end if

    select case (opt_main)
    case (:0)
        call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Missing input [ParamTransform] in dns.ini.')
        call TLAB_STOP(DNS_ERROR_INVALOPT)

    case (1)
        if (subdomain(1) /= 1 .or. subdomain(2) /= g(1)%size .or. &
            subdomain(5) /= 1 .or. subdomain(6) /= g(3)%size) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Cropping only in Oy.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (subdomain(3) < 1 .or. subdomain(4) > g(2)%size) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Cropping out of bounds in Oy.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
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
            call TLAB_WRITE_ASCII(lfile, C_FILE_LOC//'. Performing arithmetic mean of fields.')
            iopt_size = itime_size
            opt_vec(2:) = C_1_R/M_REAL(itime_size)
        end if
        if (iopt_size /= itime_size) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Number of coefficient incorrect.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if
        if (iopt_size > iopt_size_max) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Array opt_vec too small.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

    case (5)
        if (FilterDomain(1)%type == DNS_FILTER_NONE) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Filter information needs to be provided in block [Filter].')
            call TLAB_STOP(DNS_ERROR_OPTION)
        end if

    case (6)
        icalc_flow = 0 ! Force not to process the flow fields

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
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Number of blend coefficient incorrect.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

    end select

    ! -------------------------------------------------------------------
    call SCANINICHAR(bakfile, ifile, 'PostProcessing', 'Subdomain', '-1', sRes)

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

    ! -------------------------------------------------------------------
    ! Further allocation of memory space
    ! -------------------------------------------------------------------
    inb_txc = 0

    inb_scal_dst = inb_scal

    iread_flow = icalc_flow
    iread_scal = icalc_scal

    if (opt_main == 3) then ! Remesh
        isize_txc_field = max(isize_txc_field, imax_dst*jmax_dst*kmax_dst)
        inb_txc = 5
    else if (opt_main == 5) then ! Filter
        inb_txc = 4
    else if (opt_main == 6) then
        inb_txc = 5
        inb_scal_dst = 1
    end if

    if (ifourier == 1) inb_txc = max(inb_txc, 1)

    ! -------------------------------------------------------------------
    isize_wrk3d = max(isize_txc_field, imax_dst*jmax_dst*kmax_dst)

    ! -------------------------------------------------------------------
    if (icalc_flow == 1) allocate (q_dst(imax_dst*jmax_dst*kmax_dst, inb_flow))
    if (icalc_scal == 1) allocate (s_dst(imax_dst*jmax_dst*kmax_dst, inb_scal_dst))

    if (opt_main == 3) then
        allocate (x_dst(g_dst(1)%size))
        allocate (y_dst(g_dst(2)%size))
        allocate (z_dst(g_dst(3)%size))
    end if

    call TLAB_ALLOCATE(C_FILE_LOC)

    ! -------------------------------------------------------------------
    ! Read the grid
    ! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE(wrk1d)

    if (opt_main == 5) then
        do ig = 1, 3
            call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig), wrk1d)
        end do
    end if

    ! -------------------------------------------------------------------
    ! Initialize Poisson solver
    ! -------------------------------------------------------------------
    if (ifourier == 1) call OPR_FOURIER_INITIALIZE()
    
    call OPR_CHECK()

    ! -------------------------------------------------------------------
    ! Initialize cumulative field
    ! -------------------------------------------------------------------
    if (opt_main == 4 .or. opt_main == 7) then
        if (icalc_flow == 1) q_dst = C_0_R
        if (icalc_scal == 1) s_dst = C_0_R
    end if

    ! -------------------------------------------------------------------
    ! Initialize remeshing
    ! -------------------------------------------------------------------
    if (opt_main == 3) then
        call IO_READ_GRID('grid.trn', g_dst(1)%size, g_dst(2)%size, g_dst(3)%size, &
                          g_dst(1)%scale, g_dst(2)%scale, g_dst(3)%scale, x_dst, y_dst, z_dst, dummy)

        ! Check grids; Ox and Oz directions are assumed to be periodic
        dummy = (g_dst(1)%scale - g(1)%scale)/(x(g(1)%size, 1) - x(g(1)%size - 1, 1))
        if (abs(dummy) > C_1EM3_R) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Ox scales are not equal at the end.')
            call TLAB_STOP(DNS_ERROR_GRID_SCALE)
        end if
        wrk1d(1:g(1)%size, 1) = x(1:g(1)%size, 1) ! we need extra space

        dummy = (g_dst(3)%scale - g(3)%scale)/(z(g(3)%size, 1) - z(g(3)%size - 1, 1))
        if (abs(dummy) > C_1EM3_R) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Oz scales are not equal')
            call TLAB_STOP(DNS_ERROR_GRID_SCALE)
        end if
        wrk1d(1:g(3)%size, 3) = z(1:g(3)%size, 1) ! we need extra space

        ! In the Oy direction, we allow to have a different box
        jmax_aux = g(2)%size; subdomain = 0

        dummy = (y_dst(g_dst(2)%size) - y(g(2)%size, 1))/(y(g(2)%size, 1) - y(g(2)%size - 1, 1))
        if (dummy > C_1EM3_R) then ! Extend
            flag_extend = .true.
            subdomain(4) = int(dummy) + 1       ! # planes to add at the top
            jmax_aux = jmax_aux + subdomain(4)
        else if (dummy < -C_1EM3_R) then ! Crop
            flag_crop = .true.
            do j = jmax - 1, 1, -1
                if ((g(2)%nodes(j) - y_dst(g_dst(2)%size))*(g(2)%nodes(j + 1) - y_dst(g_dst(2)%size)) < C_0_R) exit
            end do
            subdomain(4) = j + 1                ! top plane of cropped region
            jmax_aux = subdomain(4)
            subdomain(3) = 1
        end if

        dummy = (y_dst(1) - y(1, 1))/(y(2, 1) - y(1, 1))
        if (dummy < -C_1EM3_R) then ! Extend
            flag_extend = .true.
            subdomain(3) = int(abs(dummy)) + 1       ! # planes to add at the bottom
            jmax_aux = jmax_aux + subdomain(3)
        else if (dummy > C_1EM3_R) then ! Crop
            flag_crop = .true.
            do j = 1, jmax - 1, 1
                if ((g(2)%nodes(j) - y_dst(1))*(g(2)%nodes(j + 1) - y_dst(1)) < C_0_R) exit
            end do
            subdomain(3) = j                   ! bottom plane of cropped region
            jmax_aux = jmax_aux - subdomain(3) + 1
        end if

        if (flag_extend .and. flag_crop) then
            call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Simultaneous extend and crop is undeveloped.')
            call TLAB_STOP(DNS_ERROR_UNDEVELOP)
        end if

        ! Reallocating memory space
        idummy = max(jmax_aux, max(g(1)%size, g(3)%size))
        isize_wrk1d = max(isize_wrk1d, idummy)
        isize_wrk1d = isize_wrk1d + 1

        inb_txc = inb_txc - 1    ! Creating txc_aux
        idummy = max(imax, imax_dst)*max(jmax_aux, max(jmax, jmax_dst))*max(kmax, kmax_dst)
        isize_txc_field = max(isize_txc_field, idummy)
#ifdef USE_MPI
        idummy = kmax*jmax_aux
        if (mod(idummy, ims_npro_i) /= 0) then ! add space for MPI transposition
            idummy = idummy/ims_npro_i
            idummy = (idummy + 1)*ims_npro_i
        end if
        idummy = idummy*max(imax, imax_dst)
        isize_txc_field = max(isize_txc_field, idummy)
#endif
        isize_wrk3d = isize_txc_field

        idummy = isize_wrk1d*7 + (isize_wrk1d + 10)*36
        isize_wrk3d = max(isize_wrk3d, idummy)

        deallocate (txc, wrk3d)

        allocate (txc(isize_txc_field, inb_txc))
        allocate (wrk3d(isize_wrk3d))
        allocate (y_aux(isize_wrk1d))
        allocate (txc_aux(imax, jmax_aux, kmax))

        ! Creating grid
        if (flag_crop) then
            write (str, '(I3)') subdomain(4)
            call TLAB_WRITE_ASCII(lfile, 'Croping above '//trim(adjustl(str))//' for remeshing...')
            write (str, '(I3)') subdomain(3)
            call TLAB_WRITE_ASCII(lfile, 'Croping below '//trim(adjustl(str))//' for remeshing...')
            call TRANS_CROP(i1, jmax, 1, subdomain, g(2)%nodes, y_aux)

            y_aux(1) = y_dst(1)             ! Using min and max of new grid
            y_aux(jmax_aux) = y_dst(g_dst(2)%size)

        else
            y_aux(1 + subdomain(3):g(2)%size + subdomain(3)) = y(1:g(2)%size, 1) ! we need extra space

            if (subdomain(4) > 0) then
                write (str, '(I3)') subdomain(4)
                call TLAB_WRITE_ASCII(lfile, 'Adding '//trim(adjustl(str))//' planes at the top for remeshing...')
                dummy = (y_dst(g_dst(2)%size) - y(g(2)%size, 1))/real(subdomain(4)) ! distributing the points uniformly
                do ip = g(2)%size + subdomain(3) + 1, g(2)%size + subdomain(3) + subdomain(4)
                    y_aux(ip) = y_aux(ip - 1) + dummy
                end do
            end if

            if (subdomain(3) > 0) then
                write (str, '(I3)') subdomain(3)
                call TLAB_WRITE_ASCII(lfile, 'Adding '//trim(adjustl(str))//' planes at the bottom for remeshing...')
                dummy = (y_dst(1) - y(1, 1))/real(subdomain(3))
                do ip = subdomain(3), 1, -1
                    y_aux(ip) = y_aux(ip + 1) + dummy ! dummy is negative
                end do
            end if

        end if

        g(2)%scale = g_dst(2)%scale  ! watch out, overwriting grid information
        g(2)%size = jmax_aux

    end if

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))
        call TLAB_WRITE_ASCII(lfile, sRes)

        if (iread_flow == 1) then ! Flow variables
            write (flow_file, *) itime; flow_file = trim(adjustl(tag_flow))//trim(adjustl(flow_file))
            call IO_READ_FIELDS(flow_file, IO_FLOW, imax, jmax, kmax, inb_flow, i0, q, wrk3d)
        end if

        if (iread_scal == 1) then ! Scalar variables
            write (scal_file, *) itime; scal_file = trim(adjustl(tag_scal))//trim(adjustl(scal_file))
            call IO_READ_FIELDS(scal_file, IO_SCAL, imax, jmax, kmax, inb_scal, i0, s, wrk3d)
        end if

        ! ###################################################################
        ! Cropping
        ! ###################################################################
        if (opt_main == 1) then
            if (icalc_flow > 0) then
                do iq = 1, inb_flow
                    call TLAB_WRITE_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_CROP(imax, jmax, kmax, subdomain, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (icalc_scal > 0) then
                do is = 1, inb_scal
                    call TLAB_WRITE_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_CROP(imax, jmax, kmax, subdomain, s(1, is), s_dst(1, is))
                end do
            end if

            ! ###################################################################
            ! Extension
            ! ###################################################################
        else if (opt_main == 2) then
            if (icalc_flow > 0) then
                do iq = 1, inb_flow
                    call TLAB_WRITE_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_EXTEND(imax, jmax, kmax, subdomain, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (icalc_scal > 0) then
                do is = 1, inb_scal
                    call TLAB_WRITE_ASCII(lfile, 'Transfering data to new array...')
                    call TRANS_EXTEND(imax, jmax, kmax, subdomain, s(1, is), s_dst(1, is))
                end do
            end if

            ! ###################################################################
            ! Change grid
            ! ###################################################################
        else if (opt_main == 3) then

            if (icalc_flow > 0) then
                do iq = 1, inb_flow
                    call TLAB_WRITE_ASCII(lfile, 'Transfering data to new array...')
                    if (flag_crop) then
                        call TRANS_CROP(imax, jmax, kmax, subdomain, q(:, iq), txc_aux)
                        do k = 1, kmax
                            txc_aux(:, 1, k) = txc_aux(:, 1, k) &
               + (y_aux(1) - y(subdomain(3), 1))*(txc_aux(:, 2, k) - txc_aux(:, 1, k))/(y(subdomain(3) + 1, 1) - y(subdomain(3), 1))
                            txc_aux(:, jmax_aux, k) = txc_aux(:, jmax_aux - 1, k) &
   + (y_aux(jmax_aux)-y(subdomain(4)-1,1)) *(txc_aux(:,jmax_aux,k)-txc_aux(:,jmax_aux-1,k)) /(y(subdomain(4),1)-y(subdomain(4)-1,1))
                        end do
                    else
                        call TRANS_EXTEND(imax, jmax, kmax, subdomain, q(:, iq), txc_aux)
                    end if
                    call OPR_INTERPOLATE(imax, jmax_aux, kmax, imax_dst, jmax_dst, kmax_dst, &
                                         g, wrk1d(:, 1), y_aux, wrk1d(:, 3), x_dst, y_dst, z_dst, &
                                         txc_aux, q_dst(:, iq), txc)
                end do
            end if

            if (icalc_scal > 0) then
                do is = 1, inb_scal
                    call TLAB_WRITE_ASCII(lfile, 'Transfering data to new array...')
                    if (flag_crop) then
                        call TRANS_CROP(imax, jmax, kmax, subdomain, s(:, is), txc_aux)
                        do k = 1, kmax
                            txc_aux(:, 1, k) = txc_aux(:, 1, k) &
               + (y_aux(1) - y(subdomain(3), 1))*(txc_aux(:, 2, k) - txc_aux(:, 1, k))/(y(subdomain(3) + 1, 1) - y(subdomain(3), 1))
                            txc_aux(:, jmax_aux, k) = txc_aux(:, jmax_aux - 1, k) &
   + (y_aux(jmax_aux)-y(subdomain(4)-1,1)) *(txc_aux(:,jmax_aux,k)-txc_aux(:,jmax_aux-1,k)) /(y(subdomain(4),1)-y(subdomain(4)-1,1))
                        end do
                    else
                        call TRANS_EXTEND(imax, jmax, kmax, subdomain, s(:, is), txc_aux)
                    end if
                    call OPR_INTERPOLATE(imax, jmax_aux, kmax, imax_dst, jmax_dst, kmax_dst, &
                                         g, wrk1d(:, 1), y_aux, wrk1d(:, 3), x_dst, y_dst, z_dst, &
                                         txc_aux, s_dst(:, is), txc)
                end do
            end if

            ! ###################################################################
            ! Linear combination of fields
            ! ###################################################################
        else if (opt_main == 4) then
            if (icalc_flow > 0) then
                q_dst = q_dst + q*opt_vec(it + 1)
            end if

            if (icalc_scal > 0) then
                s_dst = s_dst + s*opt_vec(it + 1)
            end if

            ! ###################################################################
            ! Filter
            ! ###################################################################
        else if (opt_main == 5) then
            if (icalc_flow > 0) then
                do iq = 1, inb_flow
                    call TLAB_WRITE_ASCII(lfile, 'Filtering...')
                    q_dst(:, iq) = q(:, iq) ! in-place operation
                    if (FilterDomain(1)%type == DNS_FILTER_HELMHOLTZ) &  ! Bcs depending on field
                        FilterDomain(2)%BcsMin = FilterDomainBcsFlow(iq)
                    call OPR_FILTER(imax, jmax, kmax, FilterDomain, q_dst(1, iq), wrk1d, wrk2d, txc)
                end do
            end if

            if (icalc_scal > 0) then
                do is = 1, inb_scal
                    call TLAB_WRITE_ASCII(lfile, 'Filtering...')
                    s_dst(:, is) = s(:, is) ! in-place operation
                    if (FilterDomain(1)%type == DNS_FILTER_HELMHOLTZ) & ! Bcs depending on field
                        FilterDomain(2)%BcsMin = FilterDomainBcsScal(is)
                    call OPR_FILTER(imax, jmax, kmax, FilterDomain, s_dst(1, is), wrk1d, wrk2d, txc)
                end do
            end if

            ! ###################################################################
            ! Transformation
            ! ###################################################################
        else if (opt_main == 6) then
            if (opt_function == 1) then
                call TRANS_FUNCTION(imax, jmax, kmax, s, s_dst, txc)

            else if (opt_function == 2) then
                s_dst(:, 1) = C_0_R
                do is = 1, min(inb_scal, iopt_size)
                    s_dst(:, 1) = s_dst(:, 1) + opt_vec(2 + is)*s(:, is)
                end do

            end if

            ! ###################################################################
            ! Blend
            ! ###################################################################
        else if (opt_main == 7) then
            if (it == 1) opt_vec(2) = y(1, 1) + opt_vec(2)*g(2)%scale
            write (sRes, *) opt_vec(2), opt_vec(3); sRes = 'Blending with '//trim(adjustl(sRes))
            call TLAB_WRITE_ASCII(lfile, sRes)

            if (icalc_scal > 0) then
                do is = 1, inb_scal
                    call TRANS_BLEND(imax, jmax, kmax, opt_vec(2), y, s(1, is), s_dst(1, is))
                end do
            end if

            if (icalc_flow > 0) then ! Blended fields have rtime from last velocity field
                do iq = 1, inb_flow
                    call TRANS_BLEND(imax, jmax, kmax, opt_vec(2), y, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (it == 1) opt_vec(3) = -opt_vec(3) ! flipping blending shape

            ! ###################################################################
            ! Adding mean profiles
            ! ###################################################################
        else if (opt_main == 8) then
            if (icalc_flow > 0) then
                do iq = 1, inb_flow
                    call TLAB_WRITE_ASCII(lfile, 'Adding mean flow profiles...')
                    call TRANS_ADD_MEAN(i0, iq, imax, jmax, kmax, y, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (icalc_scal > 0) then
                do is = 1, inb_scal
                    call TLAB_WRITE_ASCII(lfile, 'Adding mean scal profiles...')
                    call TRANS_ADD_MEAN(i1, is, imax, jmax, kmax, y, s(1, is), s_dst(1, is))
                end do
            end if

            ! ###################################################################
            ! Extrude
            ! ###################################################################
        else if (opt_main == 9) then
            if (icalc_flow > 0) then
                do iq = 1, inb_flow
                    call TLAB_WRITE_ASCII(lfile, 'Extruding along Oz...')
                    call TRANS_EXTRUDE(imax, jmax, kmax, subdomain, q(1, iq), q_dst(1, iq))
                end do
            end if

            if (icalc_scal > 0) then
                do is = 1, inb_scal
                    call TLAB_WRITE_ASCII(lfile, 'Extruding along Oz...')
                    call TRANS_EXTRUDE(imax, jmax, kmax, subdomain, s(1, is), s_dst(1, is))
                end do
            end if

        end if

        ! ###################################################################
        ! Writing transform fields
        ! ###################################################################
        if (opt_main /= 4 .and. opt_main /= 7) then
            if (icalc_flow > 0) then
                flow_file = trim(adjustl(flow_file))//'.trn'
                call IO_WRITE_FIELDS(flow_file, IO_FLOW, imax_dst, jmax_dst, kmax_dst, inb_flow, q_dst, wrk3d)
            end if
            if (icalc_scal > 0) then
                scal_file = trim(adjustl(scal_file))//'.trn'
                call IO_WRITE_FIELDS(scal_file, IO_SCAL, imax_dst, jmax_dst, kmax_dst, inb_scal_dst, s_dst, wrk3d)
            end if
        end if

    end do

    ! ###################################################################
    ! Final operations
    ! ###################################################################
    if (opt_main == 4 .or. opt_main == 7) then
        if (icalc_flow > 0) then
            flow_file = trim(adjustl(flow_file))//'.trn'
            call IO_WRITE_FIELDS(flow_file, IO_FLOW, imax_dst, jmax_dst, kmax_dst, inb_flow, q_dst, wrk3d)
        end if
        if (icalc_scal > 0) then
            scal_file = trim(adjustl(scal_file))//'.trn'
            call IO_WRITE_FIELDS(scal_file, IO_SCAL, imax_dst, jmax_dst, kmax_dst, inb_scal_dst, s_dst, wrk3d)
        end if
    end if

    call TLAB_STOP(0)

contains
    !########################################################################
    !# Crop array a into array b in the first two indices
    !########################################################################
    subroutine TRANS_CROP(nx, ny, nz, subdomain, a, b)
        implicit none

        TINTEGER nx, ny, nz, subdomain(6)
        TREAL, dimension(nx, ny, nz) :: a
        TREAL, dimension(nx, subdomain(4) - subdomain(3) + 1, nz) :: b

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
        implicit none

        TINTEGER nx, ny, nz, subdomain(6)
        TREAL, intent(IN) :: a(nx, ny, nz)
        TREAL, intent(OUT) :: b(subdomain(2) - subdomain(1) + 1, subdomain(4) - subdomain(3) + 1, subdomain(6) - subdomain(5) + 1)

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

        implicit none

        TINTEGER nx, ny, nz, planes(6)
        TREAL, dimension(nx, ny, nz) :: a
        TREAL, dimension(planes(1) + nx + planes(2), planes(3) + ny + planes(4), nz) :: b

        ! -----------------------------------------------------------------------
        TINTEGER j, k

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
    !# DESCRIPTION
    !#
    !########################################################################
    subroutine TRANS_ADD_MEAN(flag_mode, is, nx, ny, nz, y, a, b)

        use TLAB_CONSTANTS, only: efile
        use TLAB_VARS, only: sbg, qbg
        use PROFILES
        implicit none

        TINTEGER flag_mode, is, nx, ny, nz
        TREAL, dimension(*), intent(IN) :: y
        TREAL, dimension(nx, ny, nz), intent(IN) :: a
        TREAL, dimension(nx, ny, nz), intent(OUT) :: b

        ! -----------------------------------------------------------------------
        TINTEGER j
        TREAL dummy

        ! #######################################################################
        if (flag_mode == 0) then ! Velocity
            if (is == 1) then ! Only the mean velocity
                do j = 1, ny
                    dummy = PROFILES_CALCULATE(qbg(1), y(j))
                    b(:, j, :) = dummy + a(:, j, :)
                end do
            else
                b = a
            end if

        else                         ! Scalars
            do j = 1, ny
                dummy = PROFILES_CALCULATE(sbg(is), y(j))
                b(:, j, :) = dummy + a(:, j, :)
            end do

        end if

        return
    end subroutine TRANS_ADD_MEAN

    !########################################################################
    !# DESCRIPTION
    !#
    !# Calculate b = f(a)
    !#
    !########################################################################
    subroutine TRANS_FUNCTION(nx, ny, nz, a, b, txc)

        use TLAB_VARS, only: inb_scal, epbackground
        use THERMO_VARS, only: imixture, MRATIO, GRATIO, dsmooth
        use THERMO_VARS, only: THERMO_AI, WGHT_INV

        implicit none

        TINTEGER nx, ny, nz
        TREAL, dimension(nx*ny*nz) :: a, b
        TREAL, dimension(nx*ny*nz, *) :: txc

        ! -----------------------------------------------------------------------
        TREAL qt_0, qt_1, h_0, h_1, p
        TREAL LATENT_HEAT

        ! #######################################################################
        imixture = MIXT_TYPE_AIRWATER
        call THERMO_INITIALIZE
        MRATIO = C_1_R
        dsmooth = C_0_R
        inb_scal = 1

        LATENT_HEAT = THERMO_AI(6, 1, 1) - THERMO_AI(6, 1, 3)

        qt_0 = 9.0d-3; qt_1 = 1.5d-3
        h_0 = 0.955376d0; h_1 = 0.981965d0
        p = 0.940d0

        txc(:, 1) = h_0 + a(:)*(h_1 - h_0) ! total enthalpy
        txc(:, 2) = qt_0 + a(:)*(qt_1 - qt_0) ! total water, space for q_l
        txc(:, 3) = C_0_R
        txc(:, 4) = p                       ! pressure

        call THERMO_AIRWATER_PH(nx, ny, nz, txc(1, 2), txc(1, 1), epbackground, p)        ! Calculate q_l
        call THERMO_ANELASTIC_TEMPERATURE(nx, ny, nz, txc(1, 1), epbackground, txc(1, 5))

        ! Calculate saturation specific humidity
        call THERMO_POLYNOMIAL_PSAT(nx, ny, nz, txc(1, 5), txc(1, 1))
        txc(:, 1) = C_1_R/(MRATIO*txc(:, 4)/txc(:, 1) - C_1_R)*WGHT_INV(2)/WGHT_INV(1)
        txc(:, 1) = txc(:, 1)/(C_1_R + txc(:, 1))

        ! Calculate parameter \beta (assuming c_p = c_p,d)
        txc(:, 3) = WGHT_INV(2)/WGHT_INV(1)/GRATIO*LATENT_HEAT*LATENT_HEAT/(txc(:, 5)*txc(:, 5))

        ! Calculate s
        b(:) = txc(:, 2) - txc(:, 1)*(C_1_R + txc(:, 3)*txc(:, 2))/(C_1_R + txc(:, 3)*txc(:, 1))

        return
    end subroutine TRANS_FUNCTION

    !########################################################################
    !# DESCRIPTION
    !#
    !# b <- b + f(y)*a
    !#
    !########################################################################
    subroutine TRANS_BLEND(nx, ny, nz, params, y, a, b)

        implicit none

        TINTEGER nx, ny, nz
        TREAL, dimension(*) :: params
        TREAL, dimension(ny) :: y
        TREAL, dimension(nx, ny, nz) :: a, b

        ! -----------------------------------------------------------------------
        TINTEGER j
        TREAL shape, xi

        ! #######################################################################
        do j = 1, ny
            xi = (y(j) - params(1))/params(2)
            shape = C_05_R*(C_1_R + tanh(-C_05_R*xi))
            b(:, j, :) = b(:, j, :) + shape*a(:, j, :)
        end do

        return

    end subroutine TRANS_BLEND

end program TRANSFIELDS
