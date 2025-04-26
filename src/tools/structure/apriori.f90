#include "dns_const.h"
#include "dns_error.h"
#ifdef USE_MPI

#endif

#define C_FILE_LOC "APRIORI"

program APRIORI
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal, tag_part
    use TLab_Pointers, only: pointers_dt
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal, inb_txc, isize_field
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, flow_on, scal_on, fourier_on
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
    use TLab_Memory, only: isize_wrk3d
#endif
    use FDM, only: g, FDM_Initialize
    use Thermodynamics
    use NavierStokes, only: NavierStokes_Initialize_Parameters
    use Gravity, only: Gravity_Initialize
    use Rotation, only: Rotation_Initialize
    use TLab_Grid
    use IO_Fields
    use OPR_FILTERS
    use OPR_Fourier
    use OPR_Partial

    implicit none

! Parameter definitions
    integer(wi), parameter :: itime_size_max = 512
    integer(wi), parameter :: iopt_size_max = 512

    ! -------------------------------------------------------------------
    ! Additional local arrays
    real(wp), dimension(:, :), allocatable, save, target :: qf, sf
    real(wp), dimension(:), allocatable, save :: mean, y_aux
    type(pointers_dt), dimension(16) :: vars

    integer, parameter :: i1 = 1

! -------------------------------------------------------------------
! Local variables
! -------------------------------------------------------------------
    integer(wi) opt_main, opt_block, opt_order, opt_format
    integer(wi) iq, is, ig, ij, bcs(2, 2)
    integer(wi) nfield, idummy, jmax_aux, MaskSize
    logical iread_flow, iread_scal
    character*32 fname, bakfile, flow_file, scal_file, plot_file, time_str
    integer(wi) subdomain(6)

    integer(1) opt_gate
    integer(1), dimension(1) :: gate

! Reading variables
    character*512 sRes
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi) iopt_size
    real(wp) opt_vec(iopt_size_max)

    real(wp) params(1)

! ###################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    bakfile = trim(adjustl(ifile))//'.bak'

    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif
    call NavierStokes_Initialize_Parameters(ifile)
    call Thermodynamics_Initialize_Parameters(ifile)
    call Gravity_Initialize(ifile)
    call Rotation_Initialize(ifile)

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

    call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'ParamStructure', '-1', sRes)
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
#ifdef USE_MPI
    isize_wrk3d = isize_wrk3d + isize_field ! more space in wrk3d array needed in IO_WRITE_VISUALS
#endif

    jmax_aux = g(2)%size/opt_block

! -------------------------------------------------------------------
    if (flow_on) allocate (qf(imax*jmax*kmax, inb_flow))
    if (scal_on) allocate (sf(imax*jmax*kmax, inb_scal))

    allocate (mean(2*opt_order*nfield))

    call TLab_Initialize_Memory(C_FILE_LOC)

! -------------------------------------------------------------------
! Read the grid
! -------------------------------------------------------------------
    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

! ------------------------------------------------------------------------
! Define size of blocks
! ------------------------------------------------------------------------
    y_aux(:) = 0
    do ij = 1, jmax
        is = (ij - 1)/opt_block + 1
        y_aux(is) = y_aux(is) + g(2)%nodes(ij)/real(opt_block, wp)
    end do

! -------------------------------------------------------------------
! Initialize filters
! -------------------------------------------------------------------
    call OPR_Filter_Initialize_Parameters(ifile)
    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), FilterDomain(ig))
    end do

! -------------------------------------------------------------------
! Initialize Poisson solver
! -------------------------------------------------------------------
    if (fourier_on) call OPR_Fourier_Initialize()

    call OPR_CHECK()

! ###################################################################
! Postprocess given list of files
! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))
        call TLab_Write_ASCII(lfile, sRes)

        if (iread_flow) then ! Flow variables
            write (flow_file, *) itime; flow_file = trim(adjustl(tag_flow))//trim(adjustl(flow_file))
            call IO_Read_Fields(flow_file, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:1))
            rtime = params(1)
        end if

        if (iread_scal) then ! Scalar variables
            write (scal_file, *) itime; scal_file = trim(adjustl(tag_scal))//trim(adjustl(scal_file))
            call IO_Read_Fields(scal_file, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
            rtime = params(1)
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

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(:, 1), txc(1, 1))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), q(:, 1), txc(1, 2))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(:, 1), txc(1, 3))

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(:, 2), txc(1, 4))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), q(:, 2), txc(1, 5))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(:, 2), txc(1, 6))

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(:, 3), txc(1, 7))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), q(:, 3), txc(1, 8))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(:, 3), txc(1, 9))

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

    call TLab_Stop(0)
end program APRIORI
