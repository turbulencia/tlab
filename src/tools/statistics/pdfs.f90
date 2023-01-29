#include "types.h"
#include "dns_error.h"
#include "dns_const.h"

#define C_FILE_LOC "PDFS"

program PDFS

    use TLAB_TYPES, only: pointers_dt
    use TLAB_CONSTANTS, only: ifile, efile, lfile, gfile, tag_flow, tag_scal, wp
    use TLAB_VARS
    use TLAB_ARRAYS
    use TLAB_PROCS
#ifdef USE_MPI
    use TLAB_MPI_PROCS
#endif
    use FI_SOURCES, only: FI_BUOYANCY
    use THERMO_VARS, only: imixture
    use IO_FIELDS
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use OPR_FILTERS
    use OPR_FOURIER
    use OPR_PARTIAL

    implicit none

    TINTEGER, parameter :: itime_size_max = 512
    TINTEGER, parameter :: iopt_size_max = 20
    TINTEGER, parameter :: igate_size_max = 8
    TINTEGER, parameter :: params_size_max = 2

    ! -------------------------------------------------------------------
    ! Additional local arrays
    TREAL, allocatable :: pdf(:), y_aux(:)
    integer(1), allocatable :: gate(:)
    type(pointers_dt) :: vars(16)
    integer, parameter :: i1 = 1
    
    ! -------------------------------------------------------------------
    ! Local variables
    ! -------------------------------------------------------------------
    character*512 sRes
    character*32 fname, bakfile
    character*64 str

    TINTEGER opt_main, opt_block, opt_bins(2)
    TINTEGER opt_cond, opt_cond_scal, opt_cond_relative
    TINTEGER nfield, ifield, ij, is, bcs(2, 2), isize_pdf, ig
    TREAL dummy, eloc1, eloc2, eloc3, cos1, cos2, cos3
    TINTEGER jmax_aux, iread_flow, iread_scal, idummy
    TINTEGER ibc(16)
    TREAL vmin(16), vmax(16)
    logical reduce_data

    ! Gates for the definition of the intermittency function (partition of the fields)
    TINTEGER igate_size
    TREAL gate_threshold(igate_size_max)
    integer(1) gate_level

    TINTEGER itime_size, it
    TINTEGER itime_vec(itime_size_max)

    TINTEGER iopt_size
    TINTEGER opt_vec(iopt_size_max)
    TREAL opt_vec2(iopt_size_max)

    TINTEGER params_size
    TREAL params(params_size_max)

    !########################################################################
    !########################################################################
    bcs = 0 ! Boundary conditions for derivative operator set to biased, non-zero

    bakfile = trim(adjustl(ifile))//'.bak'

    call TLAB_START()

    call IO_READ_GLOBAL(ifile)

#ifdef USE_MPI
    call TLAB_MPI_INITIALIZE
#endif

    ! -------------------------------------------------------------------
    ! File names
    ! -------------------------------------------------------------------
#include "dns_read_times.h"

    ! -------------------------------------------------------------------
    ! Additional options
    ! -------------------------------------------------------------------
    opt_main = -1 ! default values
    opt_block = 1
    gate_level = 0
    opt_bins = 16

    call SCANINICHAR(bakfile, ifile, 'PostProcessing', 'ParamPdfs', '-1', sRes)
    iopt_size = iopt_size_max
    call LIST_INTEGER(sRes, iopt_size, opt_vec)

    if (sRes == '-1') then
#ifdef USE_MPI
#else
        write (*, *) 'Option ?'
        write (*, *) ' 1. Main variables'
        write (*, *) ' 2. Scalar gradient G_iG_i/2 equation'
        write (*, *) ' 3. Enstrophy W_iW_i/2 equation'
        write (*, *) ' 4. Strain 2S_ijS_ij/2 equation'
        write (*, *) ' 5. Velocity gradient invariants'
        write (*, *) ' 6. Flamelet equation'
        write (*, *) ' 7. Joint enstrophy and strain'
        write (*, *) ' 9. Joint scalar and scalar gradient'
        write (*, *) '10. Scalar gradient components'
        write (*, *) '11. Eigenvalues of rate-of-strain tensor'
        write (*, *) '12. Eigenframe of rate-of-strain tensor'
        write (*, *) '13. Longitudinal velocity derivatives'
        write (*, *) '14. Potential vorticity'
        write (*, *) '15. Analysis of B and V'
        read (*, *) opt_main

        write (*, *) 'Planes block size ?'
        read (*, *) opt_block

        write (*, *) 'Gate level to be used ?'
        read (*, *) gate_level

        write (*, *) 'Number of PDF bins ?'
        read (*, '(A)') sRes
        idummy = 2
        call LIST_INTEGER(sRes, idummy, opt_bins)

#endif
    else
        opt_main = opt_vec(1)
        if (iopt_size >= 2) opt_block = opt_vec(2)
        if (iopt_size >= 3) gate_level = int(opt_vec(3), KIND=1)
        if (iopt_size >= 4) opt_bins = opt_vec(4:5)

    end if

    if (opt_main < 0) then ! Check
        call TLAB_WRITE_ASCII(efile, 'PDFS. Missing input [ParamPdfs] in dns.ini.')
        call TLAB_STOP(DNS_ERROR_INVALOPT)
    end if

    if (opt_block < 1) then
        call TLAB_WRITE_ASCII(efile, 'PDFS. Invalid value of opt_block.')
        call TLAB_STOP(DNS_ERROR_INVALOPT)
    end if

    ! -------------------------------------------------------------------
    iread_flow = 0
    iread_scal = 0
    if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then; inb_txc = 6; 
    else; inb_txc = 1
    end if
    if (ifourier == 1) inb_txc = max(inb_txc, 1)
    nfield = 2

    select case (opt_main)
    case (1)
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 6)
        nfield = 4 + inb_scal_array; isize_pdf = opt_bins(1) + 2
        if (imode_eqns == DNS_EQNS_INTERNAL .or. imode_eqns == DNS_EQNS_TOTAL) nfield = nfield + 2
    case (2) ! Scalar gradient equation
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 6)
        nfield = 5; isize_pdf = opt_bins(1) + 2
    case (3) ! Enstrophy equation
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 8)
        nfield = 7; isize_pdf = opt_bins(1) + 2
    case (4) ! Strain equation
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 8)
        nfield = 5; isize_pdf = opt_bins(1) + 2
    case (5) ! Invariants
        iread_flow = 1; inb_txc = max(inb_txc, 6)
        nfield = 3; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
    case (6) ! Chi-flamelet
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 6)
        nfield = 2; isize_pdf = opt_bins(1) + 2
    case (7)
        iread_flow = 1; inb_txc = max(inb_txc, 4)
        nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
    case (9)
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 3)
        nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
    case (10)
        iread_scal = 1; inb_txc = max(inb_txc, 4)
        nfield = 5; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
    case (11) ! eigenvalues
        iread_flow = 1; inb_txc = max(inb_txc, 9)
        nfield = 3; isize_pdf = opt_bins(1) + 2
    case (12) ! eigenframe
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 9)
        nfield = 6; isize_pdf = opt_bins(1) + 2
    case (13) ! longitudinal velocity derivatives
        iread_flow = 1; inb_txc = max(inb_txc, 3)
        nfield = 3; isize_pdf = opt_bins(1) + 2
    case (14) ! potential vorticity
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 6)
        nfield = 2; isize_pdf = opt_bins(1) + 2
    case (15) ! joint s and v
        iread_scal = 1; iread_flow = 1; inb_txc = max(inb_txc, 8)
        nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
    end select

    ! -------------------------------------------------------------------
    ! Defining gate levels for conditioning
    ! -------------------------------------------------------------------
    opt_cond = 0 ! default values
    opt_cond_relative = 0
    igate_size = 0

    if (gate_level /= 0) then
#include "dns_read_partition.h"
        if (opt_cond > 1) inb_txc = max(inb_txc, 5)
    end if

    ! -------------------------------------------------------------------
    ! Allocating memory space
    ! -------------------------------------------------------------------
    allocate (gate(isize_field))

    ! in case g(2)%size is not divisible by opt_block, drop the upper most planes
    jmax_aux = g(2)%size/opt_block
    allocate (y_aux(jmax_aux)) ! Reduced vertical grid

    ! Space for the min and max of sampling variable at opt_bins+1,opt_bins+2
    ! Space for the 3D pdf at jmax_aux+1
    allocate (pdf(isize_pdf*(jmax_aux + 1)*nfield))

    inb_wrk2d = max(inb_wrk2d, 4)
    isize_wrk3d = max(isize_field, isize_txc_field)
    call TLAB_ALLOCATE(C_FILE_LOC)

    ! -------------------------------------------------------------------
    ! Initialize
    ! -------------------------------------------------------------------
    call IO_READ_GRID(gfile, g(1)%size, g(2)%size, g(3)%size, g(1)%scale, g(2)%scale, g(3)%scale, x, y, z, area)
    call FDM_INITIALIZE(x, g(1), wrk1d)
    call FDM_INITIALIZE(y, g(2), wrk1d)
    call FDM_INITIALIZE(z, g(3), wrk1d)

    call FI_BACKGROUND_INITIALIZE()  ! Initialize thermodynamic quantities

    do ig = 1, 3
        call OPR_FILTER_INITIALIZE(g(ig), Dealiasing(ig))
    end do

    if (ifourier == 1) then         ! For Poisson solver
        call OPR_FOURIER_INITIALIZE()
    end if

    call OPR_CHECK()

    y_aux(:) = 0                        ! Reduced vertical grid
    do ij = 1, jmax_aux*opt_block
        is = (ij - 1)/opt_block + 1
        y_aux(is) = y_aux(is) + y(ij, 1)/M_REAL(opt_block)
    end do

    ibc(1:nfield) = 1

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (sRes, *) itime; sRes = 'Processing iteration It'//trim(adjustl(sRes))
        call TLAB_WRITE_ASCII(lfile, sRes)

        if (iread_scal == 1) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            call IO_READ_FIELDS(fname, IO_SCAL, imax, jmax, kmax, inb_scal, 0, s)
        end if

        if (iread_flow == 1) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            call IO_READ_FIELDS(fname, IO_FLOW, imax, jmax, kmax, inb_flow, 0, q)
        end if

        call FI_DIAGNOSTIC(imax, jmax, kmax, q, s)

        ! -------------------------------------------------------------------
        ! Calculate intermittency
        ! -------------------------------------------------------------------
        if (opt_cond == 1) then ! External file
            write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname)); params_size = 2
            call IO_READ_FIELD_INT1(fname, 1, imax, jmax, kmax, itime, params_size, params, gate)
            igate_size = int(params(2))

        else if (opt_cond > 1) then
            opt_cond_scal = 1 ! Scalar field to use for the conditioning
            if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
                opt_cond_scal = inb_scal_array
            end if

            call TLAB_WRITE_ASCII(lfile, 'Calculating partition...')
            call FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
                         imax, jmax, kmax, igate_size, gate_threshold, q, s, txc, gate)

            if (jmax_aux*opt_block /= g(2)%size) then
                call REDUCE_BLOCK_INPLACE_INT1(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, gate, wrk1d)
            end if

        end if

        ! -------------------------------------------------------------------
        ! Type of PDFs
        ! -------------------------------------------------------------------
        ifield = 0
        reduce_data = .true.

        select case (opt_main)

            ! ###################################################################
            ! Main variable 2D-PDF
            ! ###################################################################
        case (1)
            ifield = ifield + 1; vars(ifield)%field => q(:, 1); vars(ifield)%tag = 'u'
            ifield = ifield + 1; vars(ifield)%field => q(:, 2); vars(ifield)%tag = 'v'
            ifield = ifield + 1; vars(ifield)%field => q(:, 3); vars(ifield)%tag = 'w'
            if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
                call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'p'
            else
                ifield = ifield + 1; vars(ifield)%field => q(:, 6); vars(ifield)%tag = 'p'
                ifield = ifield + 1; vars(ifield)%field => q(:, 5); vars(ifield)%tag = 'r'
                ifield = ifield + 1; vars(ifield)%field => q(:, 7); vars(ifield)%tag = 't'
            end if

            do is = 1, inb_scal_array
                ifield = ifield + 1; vars(ifield)%field => s(:, is); vars(ifield)%tag = 's'
                write (str, *) is; vars(ifield)%tag = trim(adjustl(vars(ifield)%tag))//trim(adjustl(str))
            end do

            do is = 1, ifield ! In case we want same interval for all heights
                if (ibc(is) == 0) call MINMAX(imax, jmax, kmax, vars(is)%field, vmin(is), vmax(is))
            end do

            ! ###################################################################
            ! Scalar gradient equation
            ! ###################################################################
        case (2)
            call TLAB_WRITE_ASCII(lfile, 'Computing scalar gradient equation...')

            call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s, q(1, 1), q(1, 2), q(1, 3), &
                                        txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
            call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s, & ! array q used as auxiliar
                                       txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), q(1, 1), wrk2d, wrk3d)
            txc(1:isize_field, 2) = txc(1:isize_field, 2)*visc/schmidt(inb_scal)
            call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 3), txc(1, 4), wrk2d, wrk3d)
            txc(1:isize_field, 5) = txc(1:isize_field, 1)/txc(1:isize_field, 3)
            txc(1:isize_field, 4) = log(txc(1:isize_field, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'GiGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'LnGiGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionMsGiGjSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuGiLapGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'StrainAMsNiNjSij'; ibc(ifield) = 2

            ! ###################################################################
            ! Enstrophy equation
            ! ###################################################################
        case (3)
            call TLAB_WRITE_ASCII(lfile, 'Computing enstrophy equation...')

            if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
                if (buoyancy%type == EQNS_NONE) then
                    txc(:, 4) = C_0_R; txc(:, 5) = C_0_R; txc(:, 6) = C_0_R
                else
                    if (buoyancy%type == EQNS_EXPLICIT) then
                        call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, epbackground, pbackground, rbackground, wrk3d)
                    else
                        wrk1d(1:jmax, 1) = C_0_R
                        call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, wrk3d, wrk1d)
                    end if
                    s(1:isize_field, 1) = wrk3d(1:isize_field)*buoyancy%vector(2)

                    call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s, txc(1, 4), wrk3d, wrk2d, wrk3d)
                    txc(:, 4) = -txc(:, 4)
                    txc(:, 5) = C_0_R
                    call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s, txc(1, 6), wrk3d, wrk2d, wrk3d)
                end if

            else
                call FI_VORTICITY_BAROCLINIC(imax, jmax, kmax, q(1, 5), q(1, 6), txc(1, 4), txc(1, 3), txc(1, 7), wrk2d, wrk3d)
            end if
            ! result vector in txc1, txc2, txc3
            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 7))
            ! scalar product, store in txc8
      txc(1:isize_field,8) = txc(1:isize_field,1)*txc(1:isize_field,4) + txc(1:isize_field,2)*txc(1:isize_field,5) + txc(1:isize_field,3)*txc(1:isize_field,6)

            call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                         txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), wrk2d, wrk3d)

            call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 2), &
                                        txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), wrk2d, wrk3d)
            txc(1:isize_field, 2) = visc*txc(1:isize_field, 2)

            call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 3), txc(1, 4), txc(1, 5), wrk2d, wrk3d)

            call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 4), txc(1, 5))

            txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 3) ! -w^2 div(u)
            txc(1:isize_field, 4) = txc(1:isize_field, 1)/txc(1:isize_field, 3) ! production rate
            txc(1:isize_field, 6) = log(txc(1:isize_field, 3))                  ! ln(w^2)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'WiWi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 6); vars(ifield)%tag = 'LnWiWi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionWiWjSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuWiLapWi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'DilatationMsWiWiDivU'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Baroclinic'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'RateANiNjSij'; ibc(ifield) = 2

            ! ###################################################################
            ! Strain equation
            ! ###################################################################
        case (4)
            call TLAB_WRITE_ASCII(lfile, 'Computing strain equation...')

            if (imode_eqns == DNS_EQNS_INCOMPRESSIBLE .or. imode_eqns == DNS_EQNS_ANELASTIC) then
                call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            else
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), q(1, 6), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            end if
            txc(1:isize_field, 1) = C_2_R*txc(1:isize_field, 2)

            call FI_STRAIN_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                      txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
            txc(1:isize_field, 2) = C_2_R*txc(1:isize_field, 2)

            call FI_STRAIN_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                     txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8))
            txc(1:isize_field, 3) = C_2_R*visc*txc(1:isize_field, 3)

            call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            txc(1:isize_field, 4) = C_2_R*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = log(txc(1:isize_field, 4))

            ifield = ifield + 1; vars(1)%field => txc(:, 4); vars(ifield)%tag = '2SijSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(2)%field => txc(:, 5); vars(ifield)%tag = 'Ln2SijSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(3)%field => txc(:, 2); vars(ifield)%tag = 'ProductionMs2SijSjkS_ki'; ibc(ifield) = 2
            ifield = ifield + 1; vars(4)%field => txc(:, 3); vars(ifield)%tag = 'DiffusionNuSijLapSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(5)%field => txc(:, 1); vars(ifield)%tag = 'Pressure2SijPij'; ibc(ifield) = 2

            ! ###################################################################
            ! Velocity gradient invariants
            ! ###################################################################
        case (5)
            call TLAB_WRITE_ASCII(lfile, 'Computing velocity gradient invariants...')

  call FI_INVARIANT_R(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call FI_INVARIANT_Q(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
            call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 3), txc(1, 4))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'InvP'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'InvQ'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'InvR'; ibc(ifield) = 2

            if (jmax_aux*opt_block /= g(2)%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.RQ'
            call PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, y_aux, txc(1, 1), txc(1, 2), pdf)

            ! ###################################################################
            ! Chi flamelet equation PDF
            ! ###################################################################
        case (6)
            call TLAB_WRITE_ASCII(lfile, 'Computing flamelet equation...')

            call FI_STRAIN_A(imax, jmax, kmax, s, q(1, 1), q(1, 2), q(1, 3), &
                             txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'StrainAGiGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'StrainA'; ibc(ifield) = 2

            ! ###################################################################
            ! Joint PDF W^2 and 2S^2
            ! ###################################################################
        case (7)
            call TLAB_WRITE_ASCII(lfile, 'Computing enstrophy-strain pdf...')

            call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), wrk2d, wrk3d)
            call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 2), txc(1, 3), txc(1, 4))
            txc(1:isize_field, 2) = C_2_R*txc(1:isize_field, 2)

            if (jmax_aux*opt_block /= g(2)%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.WS'
            call PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, y_aux, txc(1, 1), txc(1, 2), pdf)

            ! ###################################################################
            ! Joint PDF Scalar and Scalar Gradient
            ! ###################################################################
        case (9)
            call TLAB_WRITE_ASCII(lfile, 'Computing scalar-scalar--gradient pdf...')

            call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 1), txc(1, 2), wrk2d, wrk3d)
            txc(1:isize_field, 2) = log(txc(1:isize_field, 1))

            ifield = ifield + 1; vars(1)%field => s(:, 1); vars(ifield)%tag = 's'; ibc(ifield) = 1
            ifield = ifield + 1; vars(2)%field => txc(:, 1); vars(ifield)%tag = 'GiGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(3)%field => txc(:, 2); vars(ifield)%tag = 'LnGiGi'; ibc(ifield) = 3

            if (jmax_aux*opt_block /= g(2)%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.SLnG'
            call PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, s(1, 1), txc(1, 2), y_aux, pdf)

            write (fname, *) itime; fname = 'cavgGiGi'//trim(adjustl(fname))
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          1, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 1), y_aux, pdf)

            write (fname, *) itime; fname = 'cavgLnGiGi'//trim(adjustl(fname))
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          1, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 2), y_aux, pdf)

            ! ###################################################################
            ! Scalar gradient components
            ! ###################################################################
        case (10)
            call TLAB_WRITE_ASCII(lfile, 'Computing scalar gradient components...')

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s, txc(1, 1), wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s, txc(1, 2), wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s, txc(1, 3), wrk3d, wrk2d, wrk3d)
            ! Angles; s array is overwritten to save space
            do ij = 1, isize_field
                dummy = txc(ij, 2)/sqrt(txc(ij, 1)*txc(ij, 1) + txc(ij, 2)*txc(ij, 2) + txc(ij, 3)*txc(ij, 3))
                txc(ij, 4) = asin(dummy)                 ! with Oy
                s(ij, 1) = atan2(txc(ij, 3), txc(ij, 1))  ! with Ox in plane xOz
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Gx'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Gy'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Gz'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => s(:, 1); vars(ifield)%tag = 'Gtheta'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Gphi'; ibc(ifield) = 2

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.GphiS'
            call PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, s(1, 1), txc(1, 4), y_aux, pdf)

            ! ###################################################################
            ! eigenvalues of rate-of-strain tensor
            ! ###################################################################
        case (11)
            call TLAB_WRITE_ASCII(lfile, 'Computing eigenvalues of Sij...')

     CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7)) ! txc7-txc9

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'Lambda1'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Lambda2'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'Lambda3'; ibc(ifield) = 2

            ! ###################################################################
            ! eigenframe of rate-of-strain tensor
            ! ###################################################################
        case (12)
            call TLAB_WRITE_ASCII(lfile, 'Computing eigenframe of Sij...')

     CALL FI_STRAIN_TENSOR(imax,jmax,kmax, q(1,1),q(1,2),q(1,3), txc(1,1),txc(1,2),txc(1,3),txc(1,4),txc(1,5),txc(1,6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7)) ! txc7-txc9
            call TENSOR_EIGENFRAME(imax, jmax, kmax, txc(1, 1), txc(1, 7)) ! txc1-txc6
            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 7), txc(1, 8), txc(1, 9), txc(1, 10))

            do ij = 1, isize_field ! local direction cosines of vorticity vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                q(ij, 1) = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                q(ij, 2) = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                q(ij, 3) = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
            end do

            ifield = ifield + 1; vars(ifield)%field => q(:, 1); vars(ifield)%tag = 'cos(w,lambda1)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => q(:, 2); vars(ifield)%tag = 'cos(w,lambda2)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => q(:, 3); vars(ifield)%tag = 'cos(w,lambda3)'; ibc(ifield) = 2

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s, txc(1, 7), wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s, txc(1, 8), wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s, txc(1, 9), wrk3d, wrk2d, wrk3d)

            do ij = 1, isize_field ! local direction cosines of scalar gradient vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                cos1 = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                cos2 = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                cos3 = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
                txc(ij, 7) = cos1; txc(ij, 8) = cos2; txc(ij, 9) = cos3
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'cos(G,lambda1)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'cos(G,lambda2)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'cos(G,lambda3)'; ibc(ifield) = 2

            ! ###################################################################
            ! Longitudinal velocity derivatives
            ! ###################################################################
        case (13)
            call TLAB_WRITE_ASCII(lfile, 'Computing longitudinal velocity derivatives...')

            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), q(1, 1), txc(1, 1), wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), q(1, 2), txc(1, 2), wrk3d, wrk2d, wrk3d)
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), q(1, 3), txc(1, 3), wrk3d, wrk2d, wrk3d)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Sxx'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Syy'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Szz'; ibc(ifield) = 2

            ! ###################################################################
            ! Potential vorticity
            ! ###################################################################
        case (14)
            call TLAB_WRITE_ASCII(lfile, 'Computing potntial vorticity...')

            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
            txc(1:isize_field, 6) = txc(1:isize_field, 1)**2 + txc(1:isize_field, 2)**2 + txc(1:isize_field, 3)**2 ! Enstrophy
            call OPR_PARTIAL_X(OPR_P1, imax, jmax, kmax, bcs, g(1), s(1, 1), txc(1, 4), wrk3d, wrk2d, wrk3d)
            txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), s(1, 1), txc(1, 4), wrk3d, wrk2d, wrk3d)
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_PARTIAL_Z(OPR_P1, imax, jmax, kmax, bcs, g(3), s(1, 1), txc(1, 4), wrk3d, wrk2d, wrk3d)
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 3)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b

            txc(1:isize_field, 5) = sqrt(txc(1:isize_field, 5) + C_SMALL_R)
            txc(1:isize_field, 6) = sqrt(txc(1:isize_field, 6) + C_SMALL_R)
            txc(1:isize_field, 2) = txc(1:isize_field, 1)/(txc(1:isize_field, 5)*txc(1:isize_field, 6)) ! Cosine of angle between 2 vectors

            txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 1) ! Squared of the potential voticity
            txc(1:isize_field, 1) = log(txc(1:isize_field, 1) + C_SMALL_R)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'LnPotentialEnstrophy'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'CosPotentialEnstrophy'; ibc(ifield) = 2

            ! ###################################################################
            ! Analysis of B and V
            ! ###################################################################
        case (15)
            call TLAB_WRITE_ASCII(lfile, 'Computing analysis of B and V...')

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'b'; ibc(ifield) = 1
            if (buoyancy%type == EQNS_EXPLICIT) then
                call THERMO_ANELASTIC_BUOYANCY(imax, jmax, kmax, s, epbackground, pbackground, rbackground, txc(1, 1))
            else
                wrk1d(1:jmax, 1) = C_0_R
                call FI_BUOYANCY(buoyancy, imax, jmax, kmax, s, txc(1, 1), wrk1d)
            end if
            dummy = C_1_R/froude
            txc(1:isize_field, 1) = txc(1:isize_field, 1)*dummy

            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'v'; ibc(ifield) = 1
            txc(1:isize_field, 2) = q(1:isize_field, 2)

            ! I need tmp1 w/o reduction to calculate derivatives
            call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), txc(1, 1), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
            call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), txc(1, 1), txc(1, 4), txc(1, 6), wrk2d, wrk3d)
            call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), txc(1, 1), txc(1, 3), txc(1, 6), wrk2d, wrk3d)
            txc(1:isize_field, 3) = txc(1:isize_field, 3) + txc(1:isize_field, 4) + txc(1:isize_field, 5)

            ! -------------------------------------------------------------------
            if (jmax_aux*opt_block /= g(2)%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.bv'
            call PDF2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1, 1), txc(1, 2), y_aux, pdf)

            ! -------------------------------------------------------------------
            write (fname, *) itime; fname = 'cavgB'//trim(adjustl(fname))
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 1), y_aux, pdf)

            write (fname, *) itime; fname = 'cavgBii'//trim(adjustl(fname))
            if (jmax_aux*opt_block /= g(2)%size) then
                call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, txc(1, 3), wrk1d)
            end if
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 3), y_aux, pdf)
            fname = trim(adjustl(fname))//'.bv'
            call CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1, 1), txc(1, 2), txc(1, 3), y_aux, pdf)

            ! -------------------------------------------------------------------
            write (fname, *) itime; fname = 'cavgU'//trim(adjustl(fname))
            txc(1:isize_field, 3) = q(1:isize_field, 1)
            if (jmax_aux*opt_block /= g(2)%size) then
                call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, txc(1, 3), wrk1d)
            end if
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 3), y_aux, pdf)
            fname = trim(adjustl(fname))//'.bv'
            call CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1, 1), txc(1, 2), txc(1, 3), y_aux, pdf)

            write (fname, *) itime; fname = 'cavgW'//trim(adjustl(fname))
            txc(1:isize_field, 3) = q(1:isize_field, 3)
            if (jmax_aux*opt_block /= g(2)%size) then
                call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, txc(1, 3), wrk1d)
            end if
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 3), y_aux, pdf)
            fname = trim(adjustl(fname))//'.bv'
            call CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1, 1), txc(1, 2), txc(1, 3), y_aux, pdf)

            write (fname, *) itime; fname = 'cavgVii'//trim(adjustl(fname))
            call OPR_PARTIAL_Z(OPR_P2, imax, jmax, kmax, bcs, g(3), q(1, 2), txc(1, 5), txc(1, 6), wrk2d, wrk3d)
            call OPR_PARTIAL_Y(OPR_P2, imax, jmax, kmax, bcs, g(2), q(1, 2), txc(1, 4), txc(1, 6), wrk2d, wrk3d)
            call OPR_PARTIAL_X(OPR_P2, imax, jmax, kmax, bcs, g(1), q(1, 2), txc(1, 3), txc(1, 6), wrk2d, wrk3d)
            txc(1:isize_field, 3) = txc(1:isize_field, 3) + txc(1:isize_field, 4) + txc(1:isize_field, 5)
            if (jmax_aux*opt_block /= g(2)%size) then
                call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, txc(1, 3), wrk1d)
            end if
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 3), y_aux, pdf)
            fname = trim(adjustl(fname))//'.bv'
            call CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1, 1), txc(1, 2), txc(1, 3), y_aux, pdf)

            ! -------------------------------------------------------------------
            bbackground = C_0_R
            call FI_PRESSURE_BOUSSINESQ(q, s, txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call OPR_PARTIAL_Y(OPR_P1, imax, jmax, kmax, bcs, g(2), txc(1, 3), txc(1, 4), wrk3d, wrk2d, wrk3d)
            if (jmax_aux*opt_block /= g(2)%size) then
                call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, txc(1, 3), wrk1d)
                call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, txc(1, 4), wrk1d)
            end if

            write (fname, *) itime; fname = 'cavgP'//trim(adjustl(fname))
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 3), y_aux, pdf)
            fname = trim(adjustl(fname))//'.bv'
            call CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1, 1), txc(1, 2), txc(1, 3), y_aux, pdf)

            write (fname, *) itime; fname = 'cavgPy'//trim(adjustl(fname))
            call CAVG1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                          ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 4), y_aux, pdf)
            fname = trim(adjustl(fname))//'.bv'
            call CAVG2V(fname, rtime, imax*opt_block, jmax_aux, kmax, opt_bins, txc(1, 1), txc(1, 2), txc(1, 4), y_aux, pdf)

        end select

        ! ###################################################################
        if (ifield > 0) then
            if (nfield < ifield) then
                call TLAB_WRITE_ASCII(efile, C_FILE_LOC//'. Array space nfield incorrect.')
                call TLAB_STOP(DNS_ERROR_WRKSIZE)
            end if

            if (jmax_aux*opt_block /= g(2)%size .and. reduce_data) then
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, i1, i1, i1, imax, jmax_aux*opt_block, kmax, vars(is)%field, wrk1d)
                end do
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))
            call PDF1V_N(fname, rtime, imax*opt_block, jmax_aux, kmax, &
                         ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, y_aux, pdf)

        end if

    end do

    call TLAB_STOP(0)
end program PDFS
