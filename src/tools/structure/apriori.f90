#include "types.h"
#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI
#include "dns_const_mpi.h"
#endif

#define C_FILE_LOC "APRIORI"

program APRIORI

    use TLAB_TYPES, only: pointers_dt
    use TLAB_TYPES, only: filter_dt
    use TLAB_CONSTANTS
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_PROCS
#endif
    use Thermodynamics
    use IO_FIELDS
    use OPR_FILTERS
    use OPR_FOURIER
    use OPR_PARTIAL

    implicit none

! Parameter definitions
    TINTEGER, parameter :: itime_size_max = 512
    TINTEGER, parameter :: iopt_size_max = 512

    ! -------------------------------------------------------------------
    ! Additional local arrays
    TREAL, dimension(:, :), allocatable, save, target :: qf, sf
    TREAL, dimension(:), allocatable, save :: mean, y_aux
    type(pointers_dt), dimension(16) :: vars

    integer, parameter :: i1 = 1
    
! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
    TINTEGER opt_main, opt_block, opt_order, opt_format
    TINTEGER iq, is, ig, ij, bcs(2, 2)
    TINTEGER nfield, idummy, jmax_aux, MaskSize
    logical iread_flow, iread_scal
    character*32 fname, bakfile, flow_file, scal_file, plot_file, time_str
    TINTEGER subdomain(6)

    integer(1) opt_gate
    integer(1), dimension(1) :: gate

! Reading variables
    character*512 sRes
    TINTEGER itime_size, it
    TINTEGER itime_vec(itime_size_max)

    TINTEGER iopt_size
    TREAL opt_vec(iopt_size_max)

! ###################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    bakfile = trim(adjustl(ifile))//'.bak'

    call TLAB_START()

    call IO_READ_GLOBAL(ifile)
    call Thermodynamics_Initialize(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

! -------------------------------------------------------------------
! Allocating memory space
! -------------------------------------------------------------------
    allocate (y_aux(g(2)%size)) ! Reduced vertical grid

! -------------------------------------------------------------------
! File names
! -------------------------------------------------------------------
#include "dns_read_times.h"

! -------------------------------------------------------------------
! Read local options
! -------------------------------------------------------------------
    opt_main = -1 ! default values
    opt_format = 2

    opt_block = 1 ! not used yet
    opt_gate = 0
    opt_order = 1

    call SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamStructure', '-1', sRes)
    iopt_size = iopt_size_max
    call LIST_REAL(sRes, iopt_size, opt_vec)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, '(A)') 'Option ?'
        write (*, '(A)') '1. Subgrid stress'
        write (*, '(A)') '2. Velocity derivatives'
        read (*, *) opt_main
#endif
    else
        opt_main = int(opt_vec(1))
    end if

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

    MaskSize = 6

! -------------------------------------------------------------------
! Further allocation of memory space
! -------------------------------------------------------------------
    nfield = 1

    iread_flow = flow_on
    iread_scal = scal_on

    inb_txc = 4
    if (fourier_on) inb_txc = max(inb_txc, 1)

    select case (opt_main)

    case (1)
        inb_txc = max(inb_txc, 8)
        nfield = 6

    case (2)
        inb_txc = max(inb_txc, 11)
        nfield = 9

    end select

! -------------------------------------------------------------------
    isize_wrk3d = isize_txc_field
#ifdef USE_MPI
    isize_wrk3d = isize_wrk3d + isize_field ! more space in wrk3d array needed in IO_WRITE_VISUALS
#endif

    jmax_aux = g(2)%size/opt_block

! -------------------------------------------------------------------
    if (flow_on) allocate (qf(imax*jmax*kmax, inb_flow))
    if (scal_on) allocate (sf(imax*jmax*kmax, inb_scal))

    allocate (mean(2*opt_order*nfield))

    call TLAB_ALLOCATE(C_FILE_LOC)

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

! ------------------------------------------------------------------------
! Define size of blocks
! ------------------------------------------------------------------------
    y_aux(:) = 0
    do ij = 1, jmax
        is = (ij - 1)/opt_block + 1
        y_aux(is) = y_aux(is) + y(ij, 1)/M_REAL(opt_block)
    end do

! -------------------------------------------------------------------
! Initialize filters
! -------------------------------------------------------------------
    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig))
    end do

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
    if (fourier_on) call OPR_FOURIER_INITIALIZE()

    call OPR_CHECK()

! ###################################################################
! Postprocess given list of files
! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))
        call TLAB_WRITE_ASCII(lfile, sRes)

        if (iread_flow) then ! Flow variables
            write (flow_file, *) itime; flow_file = trim(adjustl(tag_flow))//trim(adjustl(flow_file))
            call IO_READ_FIELDS(flow_file, IO_FLOW, imax, jmax, kmax, inb_flow, 0, q)
        end if

        if (iread_scal) then ! Scalar variables
            write (scal_file, *) itime; scal_file = trim(adjustl(tag_scal))//trim(adjustl(scal_file))
            call IO_READ_FIELDS(scal_file, IO_SCAL, imax, jmax, kmax, inb_scal, 0, s)
        end if

! -------------------------------------------------------------------
! define time string
! -------------------------------------------------------------------
        do ij = MaskSize, 1, -1
            time_str(ij:ij) = '0'
        end do
        write (plot_file, '(I10)') itime
        time_str(MaskSize - len_trim(adjustl(plot_file)) + 1:Masksize) = trim(adjustl(plot_file))

        select case (opt_main)

! ###################################################################
! Subgrid stress
! ###################################################################
        case (1)
            do iq = 1, 3
                qf(1:isize_field, iq) = q(1:isize_field, iq)
                call OPR_FILTER(imax, jmax, kmax, FilterDomain, qf(1, iq), txc)
            end do

            nfield = 0
            nfield = nfield + 1; vars(nfield)%field => txc(:, 1); vars(nfield)%tag = 'Tauxx'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 2); vars(nfield)%tag = 'Tauyy'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 3); vars(nfield)%tag = 'Tauzz'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 4); vars(nfield)%tag = 'Tauxy'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 5); vars(nfield)%tag = 'Tauxz'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 6); vars(nfield)%tag = 'Tauyz'

            txc(1:isize_field, 1) = q(1:isize_field, 1)*q(1:isize_field, 1)
            txc(1:isize_field, 2) = q(1:isize_field, 2)*q(1:isize_field, 2)
            txc(1:isize_field, 3) = q(1:isize_field, 3)*q(1:isize_field, 3)
            txc(1:isize_field, 4) = q(1:isize_field, 1)*q(1:isize_field, 2)
            txc(1:isize_field, 5) = q(1:isize_field, 1)*q(1:isize_field, 3)
            txc(1:isize_field, 6) = q(1:isize_field, 2)*q(1:isize_field, 3)

            do is = 1, nfield
                call OPR_FILTER(imax, jmax, kmax, FilterDomain, txc(1, is), txc(1, 7))
            end do

            txc(1:isize_field, 1) = txc(1:isize_field, 1) - qf(1:isize_field, 1)*qf(1:isize_field, 1)
            txc(1:isize_field, 2) = txc(1:isize_field, 2) - qf(1:isize_field, 2)*qf(1:isize_field, 2)
            txc(1:isize_field, 3) = txc(1:isize_field, 3) - qf(1:isize_field, 3)*qf(1:isize_field, 3)
            txc(1:isize_field, 4) = txc(1:isize_field, 4) - qf(1:isize_field, 1)*qf(1:isize_field, 2)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) - qf(1:isize_field, 1)*qf(1:isize_field, 3)
            txc(1:isize_field, 6) = txc(1:isize_field, 6) - qf(1:isize_field, 2)*qf(1:isize_field, 3)

            if (jmax_aux*opt_block /= g(2)%size) then
                do is = 1, nfield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
            end if

            write (fname, *) itime; fname = 'tau'//trim(adjustl(fname))
            call AVG_N_XZ(fname, itime, rtime, imax*opt_block, jmax_aux, kmax, &
                          nfield, opt_order, vars, opt_gate, gate, y_aux, mean)

            do is = 1, nfield
                plot_file = trim(adjustl(vars(is)%tag))//time_str(1:MaskSize)
                ! to be written in terms of subarrays
                ! call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, is), wrk3d)
            end do

! ###################################################################
! Velocity derivatives
! ###################################################################
        case (2)
            nfield = 0
            nfield = nfield + 1; vars(nfield)%field => txc(:, 1); vars(nfield)%tag = 'Ux'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 2); vars(nfield)%tag = 'Uy'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 3); vars(nfield)%tag = 'Uz'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 4); vars(nfield)%tag = 'Vx'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 5); vars(nfield)%tag = 'Vy'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 6); vars(nfield)%tag = 'Vz'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 7); vars(nfield)%tag = 'Wx'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 8); vars(nfield)%tag = 'Wy'
            nfield = nfield + 1; vars(nfield)%field => txc(:, 9); vars(nfield)%tag = 'Wz'

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(:, 1), txc(1, 1))
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), q(:, 1), txc(1, 2))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(:, 1), txc(1, 3))

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(:, 2), txc(1, 4))
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), q(:, 2), txc(1, 5))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(:, 2), txc(1, 6))

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(:, 3), txc(1, 7))
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), q(:, 3), txc(1, 8))
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(:, 3), txc(1, 9))

            do is = 1, nfield
                call OPR_FILTER(imax, jmax, kmax, FilterDomain, txc(1, is), txc(1, 10))
            end do

            if (jmax_aux*opt_block /= g(2)%size) then
                do is = 1, nfield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
            end if

            write (fname, *) itime; fname = 'gradU'//trim(adjustl(fname))
            call AVG_N_XZ(fname, itime, rtime, imax*opt_block, jmax_aux, kmax, &
                          nfield, opt_order, vars, opt_gate, gate, y_aux, mean)

            do is = 1, nfield
                plot_file = trim(adjustl(vars(is)%tag))//time_str(1:MaskSize)
                ! to be written in terms of subarrays
                ! call IO_WRITE_VISUALS(plot_file, opt_format, imax, jmax, kmax, i1, subdomain, txc(1, is), wrk3d)
            end do

        end select

    end do

    call TLAB_STOP(0)
end program APRIORI
